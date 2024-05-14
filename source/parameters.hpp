#ifndef INCLUDED_parameters_hpp_
#define INCLUDED_parameters_hpp_

#define BrokenPoint 1 // 0 -> No broken point, 1 -> Exist broken point

// Parameters of STOLAS
const double sigma = 0.1; // ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
const double dn = 1; // Thickness of nsigma sphere shell
const int NL = pow(2,6); // Box size L
const double dN = 0.01; // e-folds step
const double Nprec = 1e-7; // Precision of e-foldings
const double dlnk = 0.1; // Width of bin in power spectrum

// Importance sampling
const double Nbias = 3.8; // Time of the bias
const double dNbias = 0.1; // Variance of the bias
const double bias = 0*sqrt(dNbias); // Amplitude of the bias

// Outputs
const bool szeta = true; // Output zeta
const bool spower = true; // Output power spectrum
const bool scompaction = true; // Output compaction
const bool sanimation = false; // Output animation

// Directory name of saved data, you can change after "make clean" in your terminal.
const std::string sdatadir = "data";
const std::string sourcedir = "../source";
const std::string noisefilename = "noisedata/noisemap_";
const std::string biasfilename = "biasdata/biasmap.dat";

#endif