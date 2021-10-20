# UFS CMake API

This directory contains a simplified implementation of CMake
functions and macros required to build GOCART within the UFS.

The included APIs eliminate GOCART's dependency on ESMA_cmake
and ecBuild packages required to build GEOS-5. The code has
been adapted from the following package versions:

 - ESMA_cmake v3.4.2 (May 17, 2021)
   Retrieved from: https://github.com/GEOS-ESM/ESMA_cmake

 - ecBuild 3.3.2.jcsda3 (Aug 31, 2020)
   Retrieved from: https://github.com/JCSDA-internal/ecbuild
