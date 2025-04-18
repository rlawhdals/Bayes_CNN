# About

This GitHub repository contains simulation studies and data analyses presented in the paper:

Jong Hyeon Lee, Jongmin Kim, Heesang Lee, and Jaewoo Park (2025). **A Delayed Acceptance Auxiliary Variable MCMC for Spatial Models with Intractable Likelihood Functions.**

Instructions for running the emulator and DA-AVM algorithms are provided in their respective directories.

# 0. Package Versions

(Write here if you’re listing versions of R, Rcpp, etc.)

# 1. First-Stage Kernel Construction

- `kernel.R` : Constructs the first-stage kernel using subsampling, function emulation, and a frequentist estimator.
- `kernel.cpp` : C++ functions used to construct the kernels.

# 2. DA-AVM Algorithms

- `potts.R` : Simulates the Potts model and implements the DA-AVM algorithm.
- `pp.R` : Implements DA-AVM in the point process example for the RSV dataset.
- `ergm.R` : Implements DA-AVM in the ergm example using the Faux Mesa high school network data.
