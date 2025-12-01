# Failure Probability Estimator

Estimate the failure probability in a lattice-based scheme.

## Organization

This repository is organized as follows:
* The [distributions](src/distributions.h) module
  defines optimized operations on symmetric distributions.
* For each supported scheme,
  a module under [src/schemes](src/schemes) computes the failure probability,
  using the `distributions` module.  
  We currently support:
  * [CKKS19](https://eprint.iacr.org/2019/1468)
  * [DLP14](https://eprint.iacr.org/2014/794)
  * [ID-ML-KEM_MNTRU](https://eprint.iacr.org/2025/2143)
  * [ML-KEM](https://csrc.nist.gov/pubs/fips/203/final)
* The [main](src/main.c) file displays the failure probability of some scheme,
  computed by the corresponding module.
* The compilation results (if any) are located under the [build](build) directory.
* Some intermediate results can be stored under the [saved](saved) directory.

## Usage

To add support for a new scheme:
1. Create a module under [src/schemes](src/schemes) with the name of your scheme.
1. Include the [distributions](src/distributions.h) module.
1. Define a function which computes the failure probability in your scheme.

To estimate the failure probability of a supported scheme:
1. Modify [main](src/main.c) such that it includes and uses the appropriate function.
1. Modify the `SCHEME_NAME` variable in [Makefile](Makefile) with the name of the scheme.
1. In [src/distributions.c](src/distributions.c), set `C_LEN_THREADS` to the appropriate number of threads.
1. Use `make main` to create the executable `build/main`, or `make run` to additionally run it.

**Intermediate results are stored under a `saved` directory.**  
Therefore, there must exist a directory `${PWD}/saved` when executing the program.  
This directory is automatically created by `make run`, but may have to be created manually otherwise.