<p align="center">
  <img src="assets/STOLAS_LOGO.png" width="400">
</p>

# How to use STOLAS

STOLAS is a stochastic lattice simulation code of inflation: it solves the stochastic
evolution of the inflaton on a 3D lattice and computes the curvature perturbation
by the delta-N formalism. The noise and the importance-sampling bias are generated
on the fly, so no large noise/bias files are needed.

## Download
Download from the GitHub page [STOLAS](https://github.com/STOchasticLAtticeSimulation/STOLAS_dist) to the local folder.

## Requirements
- C++20 compiler with OpenMP (GCC, or Apple clang + `libomp`)
- [FFTW3](https://www.fftw.org/) (with the threads library)
- [Boost](https://www.boost.org/) (header only; `odeint` is used)

On macOS with Homebrew: `brew install libomp fftw boost`

## Build
```
cd STOLAS
make
```
This builds two executables: `STOLAS` (the simulation) and `functions` (post-processing).
If FFTW or Boost is installed in a non-standard place, give the paths:
`make FFTW_ROOT=/path/to/fftw BOOST_ROOT=/path/to/boost`.
The compiler can be changed with `make CXX=g++-14`.

## Run
```
SIGMA=100. ./STOLAS 0
```
- The argument is the **noise seed** (here `0`). Different seeds give different realizations.
- The environment variable **SIGMA** is required for the USR model (`MODEL 2`), where it sets
  the amplitude of the power spectrum `calPzeta = SIGMA * 1e-5`.
  For the other models it is optional (default 0) and is only used as a label in the output file names.

The progress of each stage (lattice evolution, noise evolution, delta-N map, ...) is printed
to the standard output with the elapsed time and the ETA.
In a terminal the progress line is updated in place; when the output is redirected to a file
(e.g. batch jobs) a new line is printed every 10%.
Output directories under `data/` are created automatically.

## Choose the model
Give `MODEL` to `make`:
```
make MODEL=2
```
The build is redone automatically when the model is changed. Without `MODEL`, the default
written in *STOLAS/model.hpp* (`#define MODEL`) is used. For the hybrid model the number of
waterfall fields is given by `NFIELDS`, e.g. `make MODEL=3 NFIELDS=2`.

| `MODEL` | model |
|---|---|
| 0 | chaotic inflation |
| 1 | Starobinsky (piecewise linear potential) |
| 2 | ultra-slow-roll (USR) |
| 3 | hybrid inflation (`NFIELDS` waterfall fields) |

The parameters of each model (potential, initial conditions, ...) are defined below the
corresponding `#if MODEL==...` in the same file.

## Set the parameters
Parameters of STOLAS are in *STOLAS/parameters.hpp*.
- **sigma**: ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
- **dn**: Thickness of nsigma sphere shell
- **NLpower**: Box size L = 2^NLpower lattice points per side (`NLnoise = 2^NLpower` is derived automatically)
- **dN**: e-folds step
- **Nprec**: Precision of e-foldings
- **dlogn**: Width of bin in power spectrum
- **nsigmareset**, **internumber**: Settings of the zoom-in (interpolation) procedure; **internumber** is the number of zoom-in
- **EoI_noise**, **MeanNumber**: Add the noise until the end of noise (EoN), and the number of averaging
- **splan**, **FFTWwisdom**: FFTW plan/wisdom setup

You can choose the output of the result.
- **sfield**: Output the field
- **strajectory**: Output the trajectory
- **spower**: Output the power spectrum (only without importance sampling; ignored when `sbias = true`)
- **sanimation**: Output the animation
- **snoisemap**: Output the noise map
- **susrlength**: Output the USR length (`MODEL 2` only)
- **sdatadir**: Name of the output directory

Text outputs (power spectrum `power_*.dat` with columns `n` and `P`, `compaction_*.dat`, `probabilities.dat`, `trajectory_*.dat`)
are whitespace-separated tables whose first line is a `#` header naming the columns
(the weight file `logw_*.dat` has no header).
Output files are written to `data/<model>/`, e.g. the delta-N map is
`Nmap_<NLnoise>[_<NFIELDS>]_<seed>[_<zoom-in number>][_<SIGMA>].bin`
(binary, `NLnoise^3` doubles of N), where the bracketed parts are added only when they matter:
`NFIELDS` for the hybrid model, the zoom-in number when `internumber != 0`, and `SIGMA` (as an integer) for the USR model.
For example `Nmap_64_0.bin` (seed 0, no zoom-in), or `Nmap_64_0_100.bin` for USR with `SIGMA=100.`.

## Importance sampling
In *parameters.hpp*, you can switch the importance sampling and change the detail of the bias.
- **sbias**: Use the importance sampling. If `false`, the bias is not generated (saving one FFT per step) and the weight and the compaction function are neither computed nor output
- **sweight**: Output the weight (ln W of the importance sampling)
- **scompaction**: Output the compaction function
- **nbias**: Wavenumber of the bias (**Nbias** = log(nbias/sigma) is the time of the bias)
- **dNbias**: Variance of the bias
- **bias**: Amplitude of the bias

## Post-processing
`functions` reads the delta-N map (`data/<model>/Nfile/`) and `ln W` (`data/<model>/logwfile/`)
saved by `STOLAS` (so it needs `sbias = true`) and writes the compaction function recentered on the peak
(`data/<model>/cmpfile_recenter/`):
```
SIGMA=100. ./functions 0
```

## Tests
The noise generation is optimized to touch only the lattice points on the shell of each step.
The original (slow) implementation is kept in *STOLAS/tests/noise_bias_ref.hpp* for regression tests:
- `make test` compares the generated noise and bias fields of the two implementations
  (several seeds, steps and fields) for the current `MODEL` and `parameters.hpp`.
  Try several `NLpower` (e.g. 4 to 7).
- `tests/compare_nmap.sh [MODEL] [seed] [NLpower] [dN]` runs STOLAS with both implementations on a
  small lattice (copies of the sources are used, defaults: `1 0 5 0.01`) and compares all output files.
  `tests/compare.cpp` is the comparison tool (max difference, relative to the largest value).

## Batch jobs
*STOLAS/bstolas.sh* and *STOLAS/bfunc.sh* are examples of job scripts which run many seeds
(they are written for an SGE-like queue; edit the `#$` lines for your system).
