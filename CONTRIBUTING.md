# Contributing to GOCART

Contributing code to GOCART should be done via a Pull Request (PR) to this repository. 

Prior to submitting a PR, testing must be performed for start-stop and layout regression. Regression testing can be accomplished by running gcm_regress.j in the regress subdirectory of an experiment setup. Regression testing will take about thirty minutes to complete, and text files will be generated that indicate either PASS or FAIL for each regression test.

For zero-diff PRs, it must be documented that the change is zero-diff by running three sets of one day simulations using the executable before and after the changes were made: 1) an AMIP, 2) a replay, and 3) a zero-increment replay. The zero-increment replay should also be identical to the AMIP simulation. More information on how to run these types of simulations can be found at https://github.com/GEOS-ESM/GEOSgcm_GridComp/wiki/How-to-Run-GEOS-Zero-Diff-Tests.

The impacts of non zero-diff PRs must be clearly documented for consideration.


## Contributor License Agreement (CLA)

All external developers contributing to GEOS-ESM projects must complete a [Contributor License
Agreement](https://github.com/GEOS-ESM/cla).

NOTE: Internal NASA contributors associated with GEOS, including contractors,
are covered by other agreements and do not need to sign a CLA.

## License

By contributing to GEOS-ESM projects, you agree your contributions will be
licensed under the [Apache 2.0 License](LICENSE)
