# About

This GitHub repository contains simulation studies and data analyses presented in the paper:  
Jong Hyeon Lee, Jongmin Kim, Heesang Lee, and Jaewoo Park (2025).  
**A Delayed Acceptance Auxiliary Variable MCMC for Spatial Models with Intractable Likelihood Functions.**

# 0. Software Environment

This analysis was conducted using the following environment:

- **R version**: 4.4.0 (2024-04-24)
- **RStudio version**: 2023.12.1+402
- **Operating System**: macOS Sonoma 15.3.2 (Apple Silicon)
- **Key packages**:
  - `Rcpp` 1.0.13
  - `RcppArmadillo` 14.0.0.1
  - `potts` 0.5.11
  - `spatstat` 3.1.1
  - `ergm` 4.6.0

# 1. Our Algorithms

- `da_avm.cpp` : C++ functions for implementing our algorithms, including the construction of first-stage kernels.
- `packages.R` : R script for loading the required packages.

# 2. Applications with DA-AVM Algorithms

- `potts.R`: Simulates the Potts model and applies the DA-AVM algorithm.
- `pp.R`: Applies the DA-AVM algorithm to the point process example using the RSV dataset.
- `ergm.R`: Applies the DA-AVM algorithm to the ERGM example using the Faux Mesa high school network data.
