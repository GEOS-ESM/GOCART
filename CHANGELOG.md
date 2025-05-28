# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Removed

### Changed

- The pressure lid change associated with the introduction of run0 to set 0 above the lid
- Fwet value in dust modified from 0.8 to 1.0
- Dust and Sea salt Emission scale factors updated for L181
- Several changes in the DT alarm logic
  - GOCART reference time removed in `GOCART2G_GridCompMod.F90`
  - Heartbeat time step removed and `timeToWork` logical added

### Fixed

### Added

## [v2.4.1] - 2025-05-28

### Changed

- Update `components.yaml` to match that of GEOSgcm v11.7.1
  - ESMA_env v4.35.0
  - GMAO_Shared v2.0.0
  - MAPL v2.54.1
  - HEMCO geos/v2.3.0
- Updates to couple GOCART to the NOAA/UFS system after the ESMF 8.8.0 and MAPL 2.53.0 update
- Modified filepaths for the optics files to no longer link to a personal nobackup directory

### Fixed

- fixed a bug that avoids dividing nSubsteps = 0

### Added

- Added GitHub Action CI tests
- Added new drag partition options for FENGSHA to use modeled GVF or LAI for drag partition following Darmenova (2011) and Martecorena (2005) (see https://github.com/GEOS-ESM/GOCART/pull/305)

## [v2.4.0] - 2025-03-26

### Removed

- Removed all ExtData.rc files

### Changed

- Modified the file paths in carbon, sulfate, and nitrate ExtData.yaml files to used the revised version of the CEDS anthropogenic emissions. Note the previous version has an incorrect seasonal cycle.
- Sulfate surface area density calculation in SU_Compute_Diags was incorrectly being passed the effective radius used for settling along with the sigma width of the number distribution.
  Properly it should be passed the number median radius, also present in the RC file. Added a hook to read that field from the RC file ("particle_radius_number"), store in SU grid comp, and pass to SU_Compute_Diags.
  This change is zero-diff to the SU internal state. It changes value of export SO4AREA.
- Changed DMS concentration data holder from ExtData provided (SU_DMSO) to local copy (dmso_conc).
  This is relevant since if we run source tagged instances where we don't want DMS emissions we would zero out dmso_conc and that is what should be passed to DMSemission subroutine. This is zero diff except in that case.
- Changed SU2G_instance_SU.rc to now have separate filename inputs for explosive and degassing volcanoes
- It changes the formulation of the hydrophobic to hydrophilic conversion for carbon species, now defined by a time scale specified in the instance RC file.
  This is now specified by providing an e-folding time in days. This moves the time constant from outside the fortran to the run-time configurable RC file.
  This is not quite zero-diff with original code because of the precision of the specification, but testing shows nearly zero-diff result.
- Also now present in the carbon instance RC files is a run-time configurable optional parameterized loss rate (e-folding time in days) per species and per mode.
  Default value for all is set to "-1" which means no use of this function.

### Fixed

- Use 'CA' component name to identify carbonaceous contributions to PM in UFS diagnostic calculations. These contributions were missing due to changes in field names.
- Add replay import patch for ozone.
- corrected reading variable 'rhod' from files ( it was mispelled as 'rhop') - it is reverted now
- Silenced unwarranted error messages from wavelength/channel retrieval functions occurring when 470nm and/or 870nm channels are not included in GOCART resource file.
- Add explicit `find_package()` calls for missing dependencies for MAPL for builds with spack-stack. Will eventually be fixed in MAPL in later versions
- Corrected the units of the gravimetric soil moisture to percent instead of fractional in the FENGSHA dust scheme.
- Fix issue of GOCART/GEOSgcm circular CMake dependencies when used as external project
- Fix UFS/Standalone CMake issue
- Fix type of `k` in `SUvolcanicEmissions`

### Added

- Additional tuning parameters for the soil moisture and drylimit calculations for application specific tuning.
- Required attributes for the 2D GOCART export fields in AERO_DP bundle have been set in subroutine append_to_bundle in Chem_AeroGeneric.F90.
  These export fields are imported by OBIO via Surface GC, and the missing of the attributes was causing the writing of surface import checkpoint to fail.
  The issue has been explained in detail on https://github.com/GEOS-ESM/GOCART/issues/258
- Added export line to GOCART2G_GridCompMod to couple allow use of GOCART
  SU sulfate production tendency elsewhere in Chemistry, specifically for
  CARMA
- Update ESMF CMake target to `ESMF::ESMF`
- Moved present volcanic emission inventories to one or the other line for these new entries; set other
  line /dev/null; this is stop gap until next time we update volcanic emission inventories, at which
  point will provide (for AMIP and AMIP.20C) separate explosive and degassing emissions
- Made accommodating changes for above in SU2G_GridCompMod.F90 and in the Process Library
- Verified zero diff in current configuration (this is true of tracers and restarts, but not diagnostics:
  until an actual split is made in the input emissions then the volcanic emissions are being assigned to
  one or the other emission diagnostics (explosive or degassing).
- Changed Chem_SettlingSimple in the process library to call Mie Query for radius and rhop inputs to the settling velocity calculation.
  The calls to Chem_SettlingSimple were then changed accordingly in each of the species' grid comps.
  Since the RH flag is no longer needed, it was removed from GA_EnvironmentMod.F90 and each of the instance RC files.
- State Spec RC files for GOCART2G, CA, DU, NI, SU, and SS were updated such that the long names for AOD are more intuitive
- Modified ExtData.yaml files to persist as climatological anthropogenic emissions after the end of the CEDS dataset in 2019.
  Analogous rc files removed as this capability is only available with ExtData2G
- Update `components.yaml` to match that of GEOSgcm v11.6.1
  - ESMA_env v4.29.0 (Baselibs 7.24.0, Updates for SLES15 at NCCS, various fixes)
  - ESMA_cmake v3.48.0 (Fixes for NAS, debug flags, Updates for SLES15 at NCCS, MPI detection, ESMF and MPI CMake fixes for Spack)
  - GMAO_Shared v1.9.8 (Bug fix for MITgcm, CI fixes, SLES15 Updates)
  - MAPL 2.47.1 (Various fixes and features, support for Spack)
- Update CI to use Baselibs by default from CircleCI Orb
- Correct soil moisture parameterization in FENGSHA
- Add `soil_moisture_factor` to the DU2G_instance_DU.rc (same name used in the K14 scheme) and DU2G_GridCompMod.F90 files for FENGSHA
- Add `soil_drylimit_factor` to the DU2G_instance_DU.rc and DU2G_GridCompMod.F90 files for FENGSHA
- Moved process library macros to header file.

## [v2.3.0] - 2025-01-16

### Changed

- Filepath for CEDS has been updated in the ExtData yaml and rc files. Note the old version had an incorrect seasonal cycle.
- State Spec RC files for GOCART2G, CA, DU, NI, SU, and SS were updated such that the long names for AOD are more intuitive
- Updated CI to use v4 orb
- Update components.yaml to match GEOSgcm `main` as of 2025-01-16

## [v2.2.1] - 2023-05-30

### Fixed

- In dust and sea-salt, changed dimensions back to `globalCellCountPerDim` since these are needed to determine emission tuning parameters, not to allocate arrays.

## [v2.2.0] - 2023-05-18

### Fixed

- Made needed code changes in `SS2G_GridCompMod.F90` and `CA2G_GridCompMod.F90` to permit data instances of of GOCART aerosols to run
- Added missing brown carbon (BR) climatology hooks to yaml and rc files for data driven instances
- Changed pointers to climatological deposition inputs in yaml and rc files to `/dev/null` since the files pointed to didn't provide them anyway, and in any case they are being used presently in the model
- Changed pointers to climatological nitrate inputs in yaml and rc files to `/dev/null` since pointing to FP files was inconsistent with MERRA-2 files used for other species
- Ensured zero-diff in performance of yaml vs. rc files for ExtData2G vs. ExtData1g for data driven aerosols
- To do: remove hooks to old (legacy) GOCART.data instances in CHEM and setup scripts
- Fixed rc file in legacy O3 component.
- Fixed issue #223 where Global dimension was being used for allocating a local array
- This fixes a long standing issue that one can not start and stop the model in anything less than 3 hour increments to test start/stop regression because of GOCART.
- Fix issue with scattering coefficient calculation with oc
- Fix a long standing issue that one can not start and stop the model in anything less than 3 hour increments to test start/stop regression because of GOCART.

### Changed

- Comment out ASSERT to allow `GOCART_DT` to not match the `HEARTBEAT_DT`
- Single-moment moist changes from Donifan
- Change names of microphysics schemes to match refactored physics
- Set `SS_SCALE` default to 0.0
- Updates in CA2G for OpenMP
- Updates for CI
  - Update BCs version
  - Update components to match GEOSgcm v11.0.0

## [2.1.4] - 2023-05-12

### Fixed

- Fix in GOCART2G parent so that it can run with nitrates turned off. This patch of general utility was contributed by NOAA.

## [2.1.3] - 2023-02-27

### Added

- Added `*` to CA State specs file to allow for ACG to substitute in the long name
- Changes were made so GOCART2G and its children can be run with component level
OpenMP threading. The key change is to create the data structure ThreadWorkspace
to hold variables that should be private to each thread to avoid race conditions.
Additionally spatially distributed arrays that are not in any of the ESMF states
were added to the ESMF internal state so they could be properly handled when
the 'mini' ESMF sates are created. Those arrays are xhno3 for NI2G, h202_init
for SU2G, and deep_lakes_mask for SS2G. All of these arrays have MAPL_RestartSkip
option so they are not written to restart.
- Aerosol single scattering backscatter coefficient for each instances and total at wavelengths_profile
- Total (molecular + aerosols) attenuated backscatter coefficient from TOA and sfc at 532nm

### Changed

- Moved to use GitHub Action for label enforcement
- For OPS configuration: removal of links, change of QFED paths from vNRT/ to v2.5r1-nrt/ (note after November 2021, files are v2.6r1)
- For AMIP configuration: update of QFED from v2.5r1 to v2.6r1
- Update of climatological paths from MERRAero to MERRA-2
- Updated CircleCI image to use Baselibs 7.7.0
- Update `components.yaml` to reflect GEOSgcm
- CA restarts will have a change in longname for `philic` and `phobic` variables due to addition of `*` in the CA State specs file
  for the Internal state variables

## [2.1.2] - 2022-10-07

## Added

- Extinction/Scattering profile exports at model RH at wavelengths_profile
- Extinction/Scattering profile exports with RH=20% and RH=80% at wavelengths_profile

## [2.1.1] - 2022-09-16

### Fixed

- Remove GOCART requirement for gFTL, pFlogger and yaFyaml. These are requirements of MAPL. (#184)
- Update BCs version in CI for GEOSgcm run

### Added

- Added CI test for building GOCART like UFS does

## [2.1.0] - 2022-08-24

### Added

- Add `.editorconfig` file
  - This matches the styles currently used in MAPL (2 space indents in CMake and yaml, 4 spaces for Python)
- Add YAML validator GitHub Action
  - This action makes sure all YAML files are valid (to a relaxed standard)
- Initial implementation of offline simulator for 3D profiles of extinction based on the new GOCART2G-Mie interface. Work for the 2D simulator in progress.

### Fixed

- Added protection guard for pointer DU_SRC. fixed issue #148
- Initialized pointers by allocation instead of assignment. fixed issue #127
- Removed ExtData2G yaml files from all AMIP.20C directories as these are not needed anymore
- Removed declaration of Disable_Convection, no longer needed
- Added extra `esmf` to CMake files for UFS

### Changed
- Fixed typo in PM2.5 calculation (check if nitrate is active for not double counting ammonium)
- Removed nbinsDU and nbinsSS arguments from subroutine NIheterogenousChem
- Updated ExtData2G yaml files to handle AMIP magic date
- Fix bug in getAerosolSum
- more hard-coded name changes for Issue #93
- Fixed bug in MieQuery.H, shape of not present variable is used
- remove logic dinosaur "goto"
- Removed some print statements that have been commented out.
- Update `CODEOWNERS` file to make approvals less restrictive
- Updated the CircleCI to use circleci-tools v1 orb
  - Moves CI to use Baselibs 6.2.13 needed by MAPL development
- Update `components.yaml` to be in line with GEOSgcm v10.22.4
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
- Added SetVM in UFS Aerosol Cap for ESMF managed threading

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
