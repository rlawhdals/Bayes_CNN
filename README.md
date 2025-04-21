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
  - `DiceKriging` 1.6.0
  - `potts` 0.5.11
  - `spatstat` 3.1.1
  - `ergm` 4.6.0

# 1. Our Algorithms

- `potts.cpp` : C++ functions for implementing our algorithms in the potts example, including the construction of first-stage kernels.
- `pointprocess.cpp` : C++ functions for implementing our algorithms in the attraction-repulsion point process example, including the construction of first-stage kernels.
- `ergm.cpp` : C++ functions for implementing our algorithms in the ergm example, including the construction of first-stage kernels.
- `packages.R` : R script for loading the required packages.

# 2. Applications with DA-AVM Algorithms

- `rsvb2.txt` : RSV dataset
- `potts.R`: Simulates the Potts model and applies the DA-AVM algorithm.
- `pointprocess.R`: Applies the DA-AVM algorithm to the point process example using the RSV dataset.
- `ergm.R`: Applies the DA-AVM algorithm to the ERGM example using the Faux Mesa high school network data.
