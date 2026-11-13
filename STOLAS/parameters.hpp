#ifndef INCLUDED_parameters_hpp_
#define INCLUDED_parameters_hpp_

// double sigma = std::stod(std::getenv("SIGMA"));

// Parameters of STOLAS
const double sigma = 1./16; //0.5;// //ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
constexpr double dn = 1; // Thickness of nsigma sphere shell
constexpr int NLpower = 8; // Box size L = 2^NLpower lattice points per side
constexpr int NLnoise = 1 << NLpower; // (derived) number of lattice points per side
constexpr double dN = 0.001;//*M_LN2; // e-folds step
constexpr double Nprec = 1e-7; // Precision of e-foldings
constexpr double dlogn = 0.1; // Width of bin in power spectrum
// constexpr double LL = 1.; // Box size L

double nsigmareset = 8.;
constexpr int aninum = 20; // dividing number for animation
constexpr bool EoI_noise = false; // Adding the noise until EoN
constexpr int MeanNumber = 0; // Number of averageing in each
constexpr int internumber = 0; // Number of zoom-in

// FFTW setup
constexpr bool splan = false;
constexpr bool FFTWwisdom = false;

// Outputs
constexpr bool sfield = false; // Output field
constexpr bool strajectory = false; // Output the trajectory
constexpr bool spower = true; // Output the power spectrum (ignored when sbias = true)
constexpr bool sanimation = false; // Output the animation
constexpr bool snoisemap = false; // Output the noisemap
constexpr bool susrlength = false; // Output the noisemap

// Importance sampling
constexpr bool sbias = false; // Use the importance sampling (false: no bias, and no weight/compaction)
constexpr bool sweight = true; // Output the weight (needs sbias)
constexpr bool scompaction = true; // Output the compaction (needs sbias)
const double nbias = 16.; // wavenumber of the bias
const double Nbias = log(nbias/sigma); // Time of the bias
constexpr double dNbias = 0.1; // Variance of the bias
const double bias = 25.*sqrt(dNbias); // Amplitude of the bias

// Directory name of saved data, you can change after "make clean" in your terminal.
const std::string sdatadir = "data";

#endif
