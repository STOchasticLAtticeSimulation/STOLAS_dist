#ifndef INCLUDED_simd_rk4_hpp_
#define INCLUDED_simd_rk4_hpp_

// Explicit SIMD RK4 drift-step kernel for the MODEL==2 (USR) piecewise-linear
// potential, i.e. what used to be the ndiv-substep loop in evolution()
// (STOLAS.hpp, the "for (int dn = 0; dn < ndiv; ++dn) stepper_noise.do_step(...)"
// block). Measured to be >99% of the per-step lattice cost, and the compiler
// does not auto-vectorize it: a branchless rewrite of VV/Vphi alone made the
// scalar loop *slower* (extra arithmetic, no vectorization), while hand
// intrinsics gave a real ~1.8x on NEON (verified in
// STOLAS/noIO scratch benchmark). Lattice points are therefore batched by
// hand: NEON (2 doubles/vector, tested on Apple Silicon) and SVE (runtime
// vector length, for Fugaku's A64FX -- untested here, no SVE hardware
// available; verify on Fugaku before trusting it) share the same structure.
//
// dphidN for MODEL==2 does not depend on N (the "t" argument is unused), so
// this is an autonomous ODE and the classic 4th-order RK4 update below is
// exactly what boost::numeric::odeint::runge_kutta4::do_step computes.
//
// phievol itself is left untouched -- calPphi/calPpi/RecalPphipi need the
// pre-drift state -- results are written to driftedPhievol instead.

#if defined(__ARM_FEATURE_SVE)
#include <arm_sve.h>
#elif defined(__ARM_NEON)
#include <arm_neon.h>
#endif

static std::array<state_type,NLnoiseAll> driftedPhievol{};

inline void scalar_rk4_drift_tail(const std::array<state_type,NLnoiseAll>& phievol,
                                   std::array<state_type,NLnoiseAll>& out,
                                   int from, int to, int ndivLocal, double dNsubLocal) {
  for (int i = from; i < to; i++) {
    double ph = phievol[i][0], pp = phievol[i][1];
    for (int dn = 0; dn < ndivLocal; dn++) {
      auto rhs = [](double ph_, double pp_, double &dph, double &dpp) {
        double vv    = (ph_>phi1) ? (V0+B1*(ph_-phi1)-B2*phi2)
                     : (ph_>phi2) ? (V0+B2*(ph_-phi1-phi2))
                                  : (V0+B3*(ph_-phi2)-B2*phi1);
        double vphi_ = (ph_>phi1) ? B1 : (ph_>phi2) ? B2 : B3;
        double H = std::sqrt((pp_*pp_*0.5 + vv)/3.0);
        dph = pp_/H;
        dpp = -3.0*pp_ - vphi_/H;
      };
      double k1p,k1i,k2p,k2i,k3p,k3i,k4p,k4i;
      rhs(ph, pp, k1p, k1i);
      rhs(ph+0.5*dNsubLocal*k1p, pp+0.5*dNsubLocal*k1i, k2p, k2i);
      rhs(ph+0.5*dNsubLocal*k2p, pp+0.5*dNsubLocal*k2i, k3p, k3i);
      rhs(ph+dNsubLocal*k3p,     pp+dNsubLocal*k3i,     k4p, k4i);
      ph += dNsubLocal/6.0*(k1p+2*k2p+2*k3p+k4p);
      pp += dNsubLocal/6.0*(k1i+2*k2i+2*k3i+k4i);
    }
    out[i][0] = ph; out[i][1] = pp;
  }
}

#if defined(__ARM_FEATURE_SVE)
// ---------------------------------------------------------------------------
// SVE path (Fugaku / A64FX). NOT tested on real hardware -- there is no SVE
// available on the Mac this was written on. Verify numerically against the
// scalar path (see scalar_rk4_drift_tail above) on Fugaku before trusting
// production runs on this path. Written as plain functions rather than
// lambdas because sizeless SVE vector types cannot be captured by value into
// a closure.
// ---------------------------------------------------------------------------

inline void simd_rhs_sve(svbool_t pg, svfloat64_t vphi1, svfloat64_t vphi2,
                          svfloat64_t vB1, svfloat64_t vB2, svfloat64_t vB3,
                          svfloat64_t vV0mB2phi1, svfloat64_t vV0mB2phi1phi2, svfloat64_t vV0mB1phi1mB2phi2,
                          svfloat64_t vthird, svfloat64_t vhalf, svfloat64_t vnegthree,
                          svfloat64_t ph, svfloat64_t pp, svfloat64_t &dph, svfloat64_t &dpp) {
  svbool_t gt2 = svcmpgt_f64(pg, ph, vphi2);
  svbool_t gt1 = svcmpgt_f64(pg, ph, vphi1);

  svfloat64_t base = svmla_f64_x(pg, vV0mB2phi1, svsub_f64_x(pg, ph, vphi2), vB3); // V0-B2*phi1 + (ph-phi2)*B3
  svfloat64_t r2   = svmla_f64_x(pg, vV0mB2phi1phi2, ph, vB2);                     // V0-B2*(phi1+phi2) + ph*B2
  svfloat64_t r1   = svmla_f64_x(pg, vV0mB1phi1mB2phi2, ph, vB1);                  // V0-B1*phi1-B2*phi2 + ph*B1
  svfloat64_t vv = svsel_f64(gt2, r2, base);
  vv = svsel_f64(gt1, r1, vv);

  svfloat64_t vphi_ = svsel_f64(gt2, vB2, vB3);
  vphi_ = svsel_f64(gt1, vB1, vphi_);

  svfloat64_t Hsq = svmul_f64_x(pg, svmla_f64_x(pg, vv, svmul_f64_x(pg, pp, pp), vhalf), vthird); // ((pp*pp)*0.5+vv)/3
  svfloat64_t H = svsqrt_f64_x(pg, Hsq);
  svfloat64_t invH = svdiv_f64_x(pg, svdup_f64(1.0), H);

  dph = svmul_f64_x(pg, pp, invH);
  dpp = svsub_f64_x(pg, svmul_f64_x(pg, vnegthree, pp), svmul_f64_x(pg, vphi_, invH));
}

inline int simd_width() { return (int)svcntd(); }

inline void simd_rk4_drift_batch(const std::array<state_type,NLnoiseAll>& phievol,
                                  std::array<state_type,NLnoiseAll>& out,
                                  int ndivLocal, double dNsubLocal) {
  const svbool_t pg = svptrue_b64();
  const svfloat64_t vphi1 = svdup_f64(phi1), vphi2 = svdup_f64(phi2);
  const svfloat64_t vB1 = svdup_f64(B1), vB2 = svdup_f64(B2), vB3 = svdup_f64(B3);
  const svfloat64_t vV0mB2phi1        = svdup_f64(V0 - B2*phi1);
  const svfloat64_t vV0mB2phi1phi2    = svdup_f64(V0 - B2*(phi1+phi2));
  const svfloat64_t vV0mB1phi1mB2phi2 = svdup_f64(V0 - B1*phi1 - B2*phi2);
  const svfloat64_t vthird = svdup_f64(1.0/3.0), vhalf = svdup_f64(0.5), vnegthree = svdup_f64(-3.0);
  const svfloat64_t vh = svdup_f64(dNsubLocal), vh2 = svdup_f64(0.5*dNsubLocal), vh6 = svdup_f64(dNsubLocal/6.0);
  const svfloat64_t vtwo = svdup_f64(2.0);

  const int VW = (int)svcntd();
  const int nVec = (NLnoiseAll/VW)*VW;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < nVec; i += VW) {
    svfloat64x2_t pair = svld2_f64(pg, &phievol[i][0]);
    svfloat64_t phi0 = svget2_f64(pair, 0);
    svfloat64_t pi0  = svget2_f64(pair, 1);

    for (int dn = 0; dn < ndivLocal; dn++) {
      svfloat64_t k1p,k1i,k2p,k2i,k3p,k3i,k4p,k4i;
      simd_rhs_sve(pg,vphi1,vphi2,vB1,vB2,vB3,vV0mB2phi1,vV0mB2phi1phi2,vV0mB1phi1mB2phi2,vthird,vhalf,vnegthree,
                   phi0, pi0, k1p, k1i);
      simd_rhs_sve(pg,vphi1,vphi2,vB1,vB2,vB3,vV0mB2phi1,vV0mB2phi1phi2,vV0mB1phi1mB2phi2,vthird,vhalf,vnegthree,
                   svmla_f64_x(pg, phi0, vh2, k1p), svmla_f64_x(pg, pi0, vh2, k1i), k2p, k2i);
      simd_rhs_sve(pg,vphi1,vphi2,vB1,vB2,vB3,vV0mB2phi1,vV0mB2phi1phi2,vV0mB1phi1mB2phi2,vthird,vhalf,vnegthree,
                   svmla_f64_x(pg, phi0, vh2, k2p), svmla_f64_x(pg, pi0, vh2, k2i), k3p, k3i);
      simd_rhs_sve(pg,vphi1,vphi2,vB1,vB2,vB3,vV0mB2phi1,vV0mB2phi1phi2,vV0mB1phi1mB2phi2,vthird,vhalf,vnegthree,
                   svmla_f64_x(pg, phi0, vh,  k3p), svmla_f64_x(pg, pi0, vh,  k3i), k4p, k4i);

      svfloat64_t sump = svadd_f64_x(pg, svadd_f64_x(pg, k1p, k4p), svmul_f64_x(pg, vtwo, svadd_f64_x(pg, k2p, k3p)));
      svfloat64_t sumi = svadd_f64_x(pg, svadd_f64_x(pg, k1i, k4i), svmul_f64_x(pg, vtwo, svadd_f64_x(pg, k2i, k3i)));
      phi0 = svmla_f64_x(pg, phi0, vh6, sump);
      pi0  = svmla_f64_x(pg, pi0,  vh6, sumi);
    }

    svst2_f64(pg, &out[i][0], svcreate2_f64(phi0, pi0));
  }

  scalar_rk4_drift_tail(phievol, out, nVec, NLnoiseAll, ndivLocal, dNsubLocal);
}

#elif defined(__ARM_NEON)
// ---------------------------------------------------------------------------
// NEON path (dev machine, Apple Silicon). 2 doubles/vector. Verified: bit
// consistent with the scalar path to ~1e-18, ~1.8x faster in isolation.
// ---------------------------------------------------------------------------

inline float64x2_t simd_VV_neon(float64x2_t ph, float64x2_t vphi1, float64x2_t vphi2,
                                 float64x2_t vB1, float64x2_t vB2, float64x2_t vB3) {
  uint64x2_t gt2 = vcgtq_f64(ph, vphi2);
  uint64x2_t gt1 = vcgtq_f64(ph, vphi1);
  float64x2_t base = vaddq_f64(vdupq_n_f64(V0 - B2*phi1), vmulq_f64(vsubq_f64(ph, vphi2), vB3));
  float64x2_t r2   = vaddq_f64(vdupq_n_f64(V0 - B2*(phi1+phi2)), vmulq_f64(ph, vB2));
  float64x2_t r1   = vaddq_f64(vdupq_n_f64(V0 - B1*phi1 - B2*phi2), vmulq_f64(ph, vB1));
  float64x2_t v = vbslq_f64(gt2, r2, base);
  v = vbslq_f64(gt1, r1, v);
  return v;
}

inline float64x2_t simd_Vphi_neon(float64x2_t ph, float64x2_t vphi1, float64x2_t vphi2,
                                   float64x2_t vB1, float64x2_t vB2, float64x2_t vB3) {
  uint64x2_t gt2 = vcgtq_f64(ph, vphi2);
  uint64x2_t gt1 = vcgtq_f64(ph, vphi1);
  float64x2_t v = vbslq_f64(gt2, vB2, vB3);
  v = vbslq_f64(gt1, vB1, v);
  return v;
}

inline void simd_rhs_neon(float64x2_t vphi1, float64x2_t vphi2, float64x2_t vB1, float64x2_t vB2, float64x2_t vB3,
                           float64x2_t ph, float64x2_t pp, float64x2_t &dph, float64x2_t &dpp) {
  float64x2_t Hsq = vmulq_f64(vaddq_f64(vmulq_n_f64(vmulq_f64(pp,pp), 0.5),
                                         simd_VV_neon(ph,vphi1,vphi2,vB1,vB2,vB3)),
                               vdupq_n_f64(1.0/3.0));
  float64x2_t H = vsqrtq_f64(Hsq);
  float64x2_t invH = vdivq_f64(vdupq_n_f64(1.0), H);
  dph = vmulq_f64(pp, invH);
  dpp = vsubq_f64(vmulq_n_f64(pp, -3.0), vmulq_f64(simd_Vphi_neon(ph,vphi1,vphi2,vB1,vB2,vB3), invH));
}

inline int simd_width() { return 2; }

inline void simd_rk4_drift_batch(const std::array<state_type,NLnoiseAll>& phievol,
                                  std::array<state_type,NLnoiseAll>& out,
                                  int ndivLocal, double dNsubLocal) {
  constexpr int VW = 2;
  constexpr int nVec = (NLnoiseAll/VW)*VW;

  const float64x2_t vphi1 = vdupq_n_f64(phi1), vphi2 = vdupq_n_f64(phi2);
  const float64x2_t vB1 = vdupq_n_f64(B1), vB2 = vdupq_n_f64(B2), vB3 = vdupq_n_f64(B3);
  const float64x2_t vh = vdupq_n_f64(dNsubLocal), vh2 = vdupq_n_f64(0.5*dNsubLocal), vh6 = vdupq_n_f64(dNsubLocal/6.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < nVec; i += VW) {
    float64x2x2_t pair = vld2q_f64(&phievol[i][0]);
    float64x2_t phi0 = pair.val[0];
    float64x2_t pi0  = pair.val[1];

    for (int dn = 0; dn < ndivLocal; dn++) {
      float64x2_t k1p,k1i,k2p,k2i,k3p,k3i,k4p,k4i;
      simd_rhs_neon(vphi1,vphi2,vB1,vB2,vB3, phi0, pi0, k1p, k1i);
      simd_rhs_neon(vphi1,vphi2,vB1,vB2,vB3, vaddq_f64(phi0, vmulq_f64(vh2,k1p)), vaddq_f64(pi0, vmulq_f64(vh2,k1i)), k2p, k2i);
      simd_rhs_neon(vphi1,vphi2,vB1,vB2,vB3, vaddq_f64(phi0, vmulq_f64(vh2,k2p)), vaddq_f64(pi0, vmulq_f64(vh2,k2i)), k3p, k3i);
      simd_rhs_neon(vphi1,vphi2,vB1,vB2,vB3, vaddq_f64(phi0, vmulq_f64(vh, k3p)), vaddq_f64(pi0, vmulq_f64(vh, k3i)), k4p, k4i);

      float64x2_t sump = vaddq_f64(vaddq_f64(k1p,k4p), vmulq_n_f64(vaddq_f64(k2p,k3p), 2.0));
      float64x2_t sumi = vaddq_f64(vaddq_f64(k1i,k4i), vmulq_n_f64(vaddq_f64(k2i,k3i), 2.0));
      phi0 = vaddq_f64(phi0, vmulq_f64(vh6, sump));
      pi0  = vaddq_f64(pi0,  vmulq_f64(vh6, sumi));
    }

    float64x2x2_t outpair;
    outpair.val[0] = phi0; outpair.val[1] = pi0;
    vst2q_f64(&out[i][0], outpair);
  }

  scalar_rk4_drift_tail(phievol, out, nVec, NLnoiseAll, ndivLocal, dNsubLocal);
}

#else
// ---------------------------------------------------------------------------
// Portable fallback: no NEON/SVE detected, just run the scalar kernel over
// everything (still OpenMP-parallel across lattice points).
// ---------------------------------------------------------------------------

inline int simd_width() { return 1; }

inline void simd_rk4_drift_batch(const std::array<state_type,NLnoiseAll>& phievol,
                                  std::array<state_type,NLnoiseAll>& out,
                                  int ndivLocal, double dNsubLocal) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    scalar_rk4_drift_tail(phievol, out, i, i+1, ndivLocal, dNsubLocal);
  }
}

#endif

#endif
