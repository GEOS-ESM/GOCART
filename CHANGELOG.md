# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Add `.editorconfig` file
  - This matches the styles currently used in MAPL (2 space indents in CMake and yaml, 4 spaces for Python)

### Fixed
- Added protection guard for pointer DU_SRC. fixed issue #148
- Initialized pointers by allocation instead of assignment. fixed issue #127
- Removed ExtData2G yaml files from all AMIP.20C directories as these are not needed anymore 

### Changed

- Updated ExtData2G yaml files to handle AMIP magic date
- Fix bug in getAerosolSum
- more hard-coded name changes for Issue #93
- Fixed bug in MieQuery.H, shape of not present variable is used
- remove logic dinosaur "goto" 
- Removed some print statements that have been commented out.
- Update `CODEOWNERS` file to make approvals less restrictive
- Updated the CircleCI to use circleci-tools 0.13.0 orb
  - Moves CI to use Baselibs 6.2.13 needed by MAPL development
- Update `components.yaml` to be in line with GEOSgcm v10.22.1
- Updates to support Spack
- Changed the handling of state variable names in multiple instances of component (see Issue #93)
- Major refactoring of Mie table class. (see Issue #96)
   - Renamed Chem_MieTableMod.F90 --> GOCART2G_Mie2GMod.F90
   - Renamed module Chem_MieTableMod2G --> GOCART2G_Mie2GMod
   - Introduced object oriented design with type-bound methods
   - renamed some components/arguments for clarity
   - eliminated extraneous container data type that is not needed under new GOCART design.
- Cleaned up optional keyword arguments in call to Mie calculator for aerosol
  radiative forcing calculation; zero diff change
- Simplified loading of radiation MieTables.


## [2.0.7] - 2021-04-29

### Fixed

- Added a restart column to all statespec.rc files in GOCART2G. Variables being provided by ExtData must be MAPL_RestartSkip

## [2.0.6] - 2021-04-28

### Fixed

- Initialize allocatable variables in Process Library. Fixes #130

## [2.0.5] - 2021-03-14

### Added

- Added AMIP.20C ExtData configs to allow AMIP GOCART runs to work before Y2000 (during the transition from HFED to QFED)
  - Note 1: This is not a *new* scenario but rather a stopgap until Extdata is updated to allow time ranges to be specified.
  - Note 2: Temporarily, this will allow runs before Y2000 using magic-data scripting in `gcm_run.j` (a la MERRA2 in GOCART1).

### Fixed

- Fix for bad entries in `DU2G_GridComp_ExtData.rc` (see Issue #109)

## [2.0.4] - 2021-02-18

### Fixed

- Fix for layout regression failure (see Issue #103)

## [2.0.3] - 2021-02-14

### Fixed

- Fixes to allow GNU to run GOCART2G
  - Changes to `if (present(foo)) .and. associated(foo))` construct which is ambiguous Fortran

### Changed

- Updated FENGSHA dust flux according to Webb et al., Aeolian Res. 42 (2020) 100560

## [2.0.2] - 2021-01-07

### Changed

- Update CI to Baselibs 6.2.8 and GCC 11.2
- Updated dust emissions tuning to be consistent with current Ginoux default.

## [2.0.1] - 2021-10-20

### Changed

- Several updates in preparation for release for GEOS-FP, x-experiments,
etc.
- Updated rc files to reflect new ExtData structure

## Removed

- Removed aerosols from legacy GOCART

### Fixed

- Added `CONFIGURE_DEPENDS` to the `GLOB` calls in GOCART2G

## [2.0.0-rc.1] - 2021-08-06

### Changed

- Added callbacks needed for GAAS2G. Removed diagnostic prints statements. Fixed bug with
  the creation of Diagnostic Mie tables.

- Updated `CODEOWNERS` to reflect changes in staffing
- Updated `components.yaml`
  - Added fixture block
  - Added ESMA_env 3.3.0
  - Updated ESMA_cmake to 3.5.0
  - Updated GMAO_Shared to 1.4.3
  - Updated MAPL to v2.7.0
- Updated CircleCI to use `large` resource and v6.2.4 Baselibs
- Rename BUILD_UFS CMake flag as UFS_GOCART
- Rename UFS target as UFS_Aerosols
- Add CMake macros replacing ecBuild and ESMA_cmake for the UFS
- Relax PFLOGGER dependency requirement outside Baselibs
- Refactored UFS Aerosols: introduced dynamic tracer mapping and limited unit conversion

### Fixed

- Fixes to CMake to build with UFS and with no Baselibs
- Fixed return code handling
- Fixed uninitialized rc in Cubic in process library
- Fixed build issue with GNU 9.2.0

### Removed

- Removed Soil Erosion Potential Distribution from FENGSHA dust scheme

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
