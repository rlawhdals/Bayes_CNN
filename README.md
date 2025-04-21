# About

This GitHub repository contains simulation studies and data analyses presented in the paper:  
Jong Hyeon Lee, Jongmin Kim, Heesang Lee, and Jaewoo Park (2025).  
**A Delayed Acceptance Auxiliary Variable MCMC for Spatial Models with Intractable Likelihood Functions.**

Instructions for running the emulator and DA-AVM algorithms are provided in their respective directories.

# 0. Software Environment

This analysis was conducted using the following environment:

- **R version**: 4.4.0 (2024-04-24)
- **RStudio version**: 2023.12.1+402
- **Operating System**: macOS Sonoma 15.3.2 (Apple Silicon)
- **Key packages**:
  - `Rcpp` 1.0.13
  - `RcppArmadillo` 14.0.0.1
  - `ergm` 4.6.0
  - `potts` 0.5.11

# 1. First-Stage Kernel Construction

- `kernel.cpp` : C++ functions used to construct the kernels and our algorithms.

# 2. DA-AVM Algorithms

- `potts.R` : Simulates the Potts model and implements the DA-AVM algorithm.
- `pp.R` : Implements DA-AVM in the point process example for the RSV dataset.
- `ergm.R` : Implements DA-AVM in the ergm example using the Faux Mesa high school network data.
