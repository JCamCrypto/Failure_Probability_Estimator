# Failure Probability Estimator

Estimate the failure probability in a lattice-based scheme.

## Organization

This repository contains code for CPU and for GPU.  
It is organized as follows:
* The [distributions](src/distributions.hpp) module
  defines optimized operations on symmetric distributions.
* For each supported scheme,
  a module under [src/schemes](src/schemes) computes the failure probability,
  using the `distributions` module.  
  We currently support:
  * [CKKS19](https://eprint.iacr.org/2019/1468)
  * [DLP14](https://eprint.iacr.org/2014/794)
  * [ID-ML-KEM_MNTRU](https://eprint.iacr.org/2025/2143)
  * [ML-KEM](https://csrc.nist.gov/pubs/fips/203/final)
* The [main](src/main.cpp) file displays the failure probability of some scheme,
  computed by the corresponding module.
* The compilation results (if any) are located under the [build](build) directory.
* Some results can be stored under the [saved](saved) directory.

## Extension

To add support for a new scheme:
1. Create a module under [src/schemes](src/schemes) with the name of your scheme.
1. Include the [distributions](src/distributions.hpp) module.
1. Define a function which computes the failure probability in your scheme.

## Usage

### CPU version

To estimate the failure probability of a supported scheme:
1. Modify [main](src/main.cpp) such that it includes the module corresponding to that scheme.
1. In [src/distributions-cpu.cpp](src/distributions-cpu.cpp), set `C_LEN_THREADS` to the appropriate number of threads.
1. Use `make DEVICE=CPU SCHEME=<name> run` to create the executable `build/main` and run it.  
   `<name>` must be such that `src/schemes/<name>.(cpp|hpp)` are the corresponding files.

### GPU version

To estimate the failure probability of a supported scheme:
1. Modify [main](src/main.cpp) such that it includes and uses the appropriate function.
1. In [src/distributions-gpu.cu](src/distributions-gpu.cu), set `C_NB_DEVICES` to the appropriate number of GPUs.
1. Use `make DEVICE=GPU SCHEME=<name> run` to create the executable `build/main` and run it.  
   `<name>` must be such that `src/schemes/<name>.(cpp|hpp)` are the corresponding files.
