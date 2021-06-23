# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Updated `CODEOWNERS` to reflect changes in staffing
- Updated `components.yaml`
  - Added fixture block
  - Added ESMA_env 3.2.1
  - Updated ESMA_cmake to 3.4.0
  - Updated GMAO_Shared to 1.4.0
  - Updated MAPL to v2.6.7
- Updated CircleCI to use `large` resource
- Rename BUILD_UFS CMake flag as UFS_GOCART
- Rename UFS target as UFS_Aerosols
- Add CMake macros replacing ecBuild and ESMA_cmake for the UFS
- Relax PFLOGGER dependency requirement outside Baselibs
- Refactored UFS Aerosols: introduced dynamic tracer mapping and limited unit conversion

### Fixed

- Fixes to CMake to build with UFS and with no Baselibs
- Fixed return code handling
- Fixed uninitialized rc in Cubic in process library
- Fixed OH unit conversion

### Removed

### Added

- Compute and export PM2.5 and PM10 diagnostic tracers in UFS interface
- Add NOAA/ARL FENGSHA dust scheme

## [1.0.1] - 2021-03-22

### Fixed

- Corrected CircleCI configuration
- Added CHANGELOG enforcer

## [1.0.0] - 2021-03-17

This release of GOCART is aimed at moving GOCART from GEOSchem_GridComp to its own repository. At the moment the only code here that is verified as correct is the "GOCART Legacy" code which is code that came from GEOSchem_GridComp v1.4.3.

To use this release, you need to use GEOSchem_GridComp v1.4.4 which has the changes necessary to move GOCART Legacy to this separate repository.

Identical in code to v1.0.0-beta.1
 
## [1.0.0-beta.1] - 2021-03-01

### Changed

- Pull in updates from `develop`. In GEOSgcm setup, this is zero-diff for both GOCART actual and climatological (e.g., not GOCART2G).

## [1.0.0-beta] - 2021-02-23

### Added

- Added GOCART legacy source code which is been pulled out of the GEOSchem_GridComp v1.4.2.

## [0.9.1] - 2020-09-29

### Changed

- Switched StrTemplate from using the one provided by `GMAO_MPEU` to one provided by `MAPL`

### Fixed

- Fixed bugs in the process library
- Bug with how some optional `rc` arguments where handled in the Process Library.
 
### Added

- CircleCI configuration added
 
## [0.9.0] - 2020-09-03

### Added

- Initial release of GOCART repository
