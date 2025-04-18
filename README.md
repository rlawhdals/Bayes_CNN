# About

This GitHub repository includes some of the simulation studies and data analysis included in the paper: Jong Hyeon Lee, Jongmin Kim, Heesang Lee and Jaewoo Park. (2025). A Delayed Acceptance Auxiliary Variable MCMC for Spatial Models with Intractable Likelihood Function.

Instructions for running the emulator and DA-AVM algorithms can be found in their respective directories.

# 0. Package version

# 1. Emulator
- `data_sim.R` : simulate Ising, potts, and the autologistic model
- `data_sim.cpp` : C++ function to generate Ising, and autologistic model
# 2. DA-AVM
- `potts.R` : simulate the Potts model and implement DA-AVM
- `pp.R`: implement DA-AVM using the RSV dataset
- `ergm.R` : implement DA-AVM using the Faux Mesa high school data
