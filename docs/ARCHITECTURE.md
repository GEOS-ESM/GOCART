# GOCART: architecture

**Audience:** coding agents and developers working on the GOCART repository (Goddard Chemistry, Aerosol, Radiation, and Transport model).
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GOCART` @ `0d19ada` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`; repository `GOCART`, mounted at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GOCART`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths below are relative to the GOCART repository root unless stated otherwise.
Component detail: [GOCART2G](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md). Agent rules: [AGENTS.md](../AGENTS.md). Decision log: [decisions.md](decisions.md).
Up:
- Parent repository `GEOSchem_GridComp`: `../../docs/ARCHITECTURE.md` (written separately).
- Physics: `../../../ARCHITECTURE.md`.
- Moist, a consumer of the aerosol activation method: `../../../GEOSmoist_GridComp/ARCHITECTURE.md`.
- Repository `GEOSgcm_GridComp`: `../../../../../docs/ARCHITECTURE.md`.
- Fixture: `../../../../../../../../docs/ARCHITECTURE.md`.

## 1. Purpose
GOCART holds the aerosol components of GEOS. The code that runs in a standard GEOSgcm configuration is `GOCART2G`, a MAPL2 gridded component with five aerosol-species children (dust, sea salt, sulfate, carbonaceous, nitrate) plus two shared libraries. The repository also holds an older component (`GOCART_GridComp`, trace gases and legacy aerosol) and cap code for UFS (the NOAA Unified Forecast System). Chemistry adds GOCART2G as a child and passes its aerosol optical and microphysical states up to Radiation, Moist and Surface (`../GEOS_ChemGridComp.F90`, `SetServices`; `../GEOS_ChemGridComp.rc`, `ENABLE_GOCART2G`, outside this repository).

## 2. Place in the hierarchy (GEOSgcm fixture)
```
GCM (GEOSgcm_GridComp)
  AGCM -> Physics (GEOSphysics_GridComp)
            CHEM (GEOSchem_GridComp, GEOS_ChemGridComp.F90)
              GOCART2G   (this repository, ESMF/GOCART2G_GridComp)
                DU, SS, CA.*, SU, NI   (one child per instance name)
              GOCART, GOCART.data  (this repository, ESMF/GOCART_GridComp; legacy)
              other chemistry children (ACHEM, GAAS, GMICHEM, CARMA, HEMCO, ...)
```
- Chemistry adds `GOCART2G` with `MAPL_AddChild(GC, NAME='GOCART2G', SS=GOCART2G_SetServices)` when `ENABLE_GOCART2G` is true (`../GEOS_ChemGridComp.F90`, `SetServices`). The code default is `.FALSE.`; `../GEOS_ChemGridComp.rc` sets `ENABLE_GOCART2G: .TRUE.`.
- Chemistry adds the legacy `GOCART` and `GOCART.data` children under `ENABLE_GOCART` and `ENABLE_GOCART_DATA`, using `GOCART_GridCompMod.F90` (`../GEOS_ChemGridComp.F90`, `SetServices`; `../GEOS_ChemGridComp.rc`). The rc enables all three; `gcm_setup` in the fixture (`Applications/@GEOSgcm_App/gcm_setup`) sets `ENABLE_GOCART_DATA` false when `DATA_DRIVEN` is true (outside this repository).
- Which provider supplies `AERO` and `AERO_DP` is chosen by `AERO_PROVIDER` (`../GEOS_ChemGridComp.F90`, `GetProvider_`; `Applications/@GEOSgcm_App/gcm_setup`). See section 5.
- The parent component repository `GEOSchem_GridComp` has its own architecture document, written separately (link above).

## 3. Module map
| Directory | Content | Status in the GEOSgcm fixture |
|---|---|---|
| `ESMF/GOCART2G_GridComp` | Component `GOCART2G`: parent (`GOCART2G_GridCompMod.F90`) and children `DU2G_GridComp`, `SS2G_GridComp`, `SU2G_GridComp`, `CA2G_GridComp`, `NI2G_GridComp`; `GA_Environment` (base type for children). [Detail](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md) | Built (`ESMF/CMakeLists.txt`); live |
| `Process_Library` | Library `Process_Library`: physical process routines (`GOCART2G_Process.F90`) and Mie-table code (`GOCART2G_MieMod.F90`) (`Process_Library/CMakeLists.txt`) | Built; used by GOCART2G and by `@GMI` |
| `ESMF/Shared` | Library `Chem_Shared2G`: `Chem_AeroGeneric.F90`, `ReplenishAlarm.F90`, `GOCART2G_AerosolActivationMod.F90`; `g2g_rst_generator.py` (`ESMF/Shared/CMakeLists.txt`) | Built; live |
| `ESMF/Aerosol_GridComp` | Module `Aerosol_GridComp_mod` (`Aerosol_GridComp.F90`): a wrapper that adds `GOCART2G` as its only child (`MAPL_AddChild`) and runs all phases of all children; ACG file `Aerosol_StateSpecs.rc` | Built in both modes (`ESMF/CMakeLists.txt`). Not added by Chemistry in this fixture: *Inferred:* no `AddChild` of it appears in `../GEOS_ChemGridComp.F90`, and no CMake file in the fixture other than GOCART's references it |
| `ESMF/UFS` | NUOPC cap (`Aerosol_Cap.F90`, which uses `Aerosol_GridComp_mod`) and `Aerosol_*` support modules; library `UFS_Aerosols` | Built only when option `UFS_GOCART` is on (`CMakeLists.txt`, `ESMF/CMakeLists.txt`) |
| `ESMF/GOCART_GridComp` | Legacy `GOCART`: `GOCART_GridCompMod.F90`, `Aero_GridCompMod.F90`, `GOCART_StubCompMod.F90`; trace-gas children `O3_GridComp`, `CO_GridComp`, `CO2_GridComp`, `CFC_GridComp`, `Rn_GridComp`, `CH4_GridComp`; `GOCARTdata_AerRegistry.rc` | Built in non-UFS mode (`ESMF/CMakeLists.txt`); added by Chemistry as `GOCART` and `GOCART.data` |
| `ESMF/Apps` | Offline tool `aop_calculator.x` (aerosol optical properties) using `GOCART2G_AopMod`, `GOCART2G_SimpleBundleMod`, `GOCART2G_MieMod` (`ESMF/Apps/CMakeLists.txt`) | Built in non-UFS mode |
| `ESMF/HEMCO_GridComp` | Not in this repository tree at the verified commit; `components.yaml` fetches HEMCO there (`components.yaml`, `HEMCO`). In the fixture, `HEMCO_GridComp` is a directory of the parent repository (`../HEMCO_GridComp/CMakeLists.txt`, relative to the GOCART root) | Not documented here |
| `CCPP` | `GOCART_CCPP.F90`: empty file (0 bytes, tracked) | Not built; no CMake reference |
| `ESMF/NUOPC` | `NUOPC_Cap.F90`: empty file (0 bytes, tracked) | Not built; no CMake reference |
| `config` | `GMAO_Shared.sparse`: sparse checkout list for the `GMAO_Shared` dependency | Used by `components.yaml` |
| `.github` | `CODEOWNERS`, workflows (changelog, labels, YAML lint, CI build) | Repository tooling |

`README.md` states the framework interfaces as "MAPL, NUOPC, and CCPP" (`README.md`). *Unverified:* the CCPP and NUOPC directories are placeholders for planned work; their files are empty.

## 4. Build structure
- Top level (`CMakeLists.txt`): in the fixture the repository is a subproject. A standalone build is possible only when an ESMA cmake directory exists; option `UFS_GOCART` (default OFF) forces the standalone path and a 32-bit build. It sets `DOING_GEOS5`, then `add_subdirectory(ESMF)` and `add_subdirectory(Process_Library)`. `CHANGELOG.md` records that standalone GOCART builds were removed. *Unverified:* the extent to which the standalone path still works.
- `ESMF/CMakeLists.txt`: UFS mode builds `Aerosol_GridComp`, `GOCART2G_GridComp`, `Shared`, `UFS`. Otherwise it builds `Shared`, `Apps`, `Aerosol_GridComp`, `GOCART_GridComp`, `GOCART2G_GridComp`, `HEMCO_GridComp`. It sets `HEMCO_EXTERNAL_CONFIG`, `MAPL_ESMF` and `BUILD_GEOS_INTERFACE`. (`ESMF/CMakeLists.txt`)
- Libraries and dependencies (CMake files in each directory):
  - `Process_Library`: depends on ESMF and NetCDF Fortran.
  - `Chem_Shared2G` (`esma_set_this(OVERRIDE Chem_Shared2G)`): built from three files; used by GOCART2G.
  - `GOCART2G_GridComp`: `esma_set_this` with SUBCOMPONENTS `GA_Environment DU2G SS2G CA2G NI2G SU2G`; depends on MAPL, `Chem_Shared2G`, `Process_Library`, ESMF; links OpenMP unless the compiler is NAG; installs `GOCART2G_GridComp.rc` to `etc`.
  - Each child's CMake runs `mapl_acg` (the MAPL Automatic Code Generator) on its `*_StateSpecs.rc` and installs every `*.rc` and `*.yaml` in its directory, including the `AMIP/` and `AMIP.20C/` subdirectories, to `etc/<directory>`.
- Direct dependency declarations for the repository (`components.yaml`): `env`, `cmake`, `ecbuild`, `HEMCO` (`ESMF/HEMCO_GridComp@`), `GMAO_Shared` (sparse, `config/GMAO_Shared.sparse`), `MAPL` (`ESMF/Shared/MAPL@`, tag `v2.71.0`). *Inferred:* in the GEOSgcm fixture the fixture's own `components.yaml` supplies MAPL and GMAO_Shared and these entries are used only when GOCART is cloned alone.
- The parent `@GEOSchem_GridComp/CMakeLists.txt` depends on `GOCART_GridComp GOCART2G_GridComp` and calls `esma_add_subdirectories(Shared GOCART TR ...)` (outside this repository).
- Tests: no test is registered in any CMake file of this repository; `ESMF/GOCART_GridComp/ut_GOCART.F90` is not referenced.
- Continuous integration (`.github/workflows/workflow.yml`, `.circleci/config.yml`) builds the GEOSgcm fixture with several compilers.

## 5. Interfaces to the rest of GEOS
Everything below is on the Chemistry side unless noted. Chemistry file: `../GEOS_ChemGridComp.F90`.

### 5.1 Aerosol states exported upward
- GOCART2G exports `AERO` (an ESMF state) and `AERO_DP` (a field bundle, used for dry-deposition-related bins) (`ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `SetServices`).
- Chemistry promotes both from the provider child: `MAPL_AddExportSpec(..., CHILD_ID=AERO_PROVIDER)` (`../GEOS_ChemGridComp.F90`, `SetServices`).
- Physics connects `AERO` from CHEM to Radiation and to Moist, and `AERO_DP` to Surface (`../../GEOS_PhysicsGridComp.F90`, `SetServices`; see the Physics document (`../../../ARCHITECTURE.md` in the fixture)).
- `AERO` carries ESMF methods that consumers call, not only data: `run_aerosol_optics` (Radiation, GMI photolysis), `get_mixRatioSum`, `serialize_bundle`, `get_monochromatic_aop` (GAAS), `aerosol_activation_properties` (Moist) (`GOCART2G_GridCompMod.F90`, `Initialize`; `ESMF/Shared/GOCART2G_AerosolActivationMod.F90`, `setup_aerosol_activation`). Details in the [GOCART2G document](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md).
- Consumers of the methods (outside this repository): `GEOSsolar_GridComp/GEOS_SolarGridComp.F90` and `GEOSirrad_GridComp/GEOS_IrradGridComp.F90` (optics); `GEOSmoist_GridComp/aer_actv_single_moment.F90` (activation; see the Moist document (`../../../GEOSmoist_GridComp/ARCHITECTURE.md` in the fixture)); `@GAAS/GAAS_GridCompMod.F90`; `@GMI/GmiPhotolysis_GridCompClassMod.F90`.

### 5.2 Provider selection
- `AERO_PROVIDER` selects the child that supplies `AERO`: `GetProvider_` maps `GOCART2G` to the child index (`../GEOS_ChemGridComp.F90`, `GetProvider_`). The code default in Chemistry is `GOCART2G`.
- If the provider is `none`, Chemistry creates its own `AERO` and `AERO_DP` and marks the optics and activation methods as not implemented; a source comment calls this "an architecture violation" (`../GEOS_ChemGridComp.F90`, `SetServices`).
- The fixture's `gcm_setup` writes `AERO_PROVIDER = GOCART2G` into the run configuration (`Applications/@GEOSgcm_App/gcm_setup`, `AGCM.rc.tmpl`; comment lists options `GOCART2G, MAM, none`).

### 5.3 Connections between GOCART2G and other Chemistry children
Made in `../GEOS_ChemGridComp.F90`, `SetServices`, with `MAPL_AddConnectivity`:
| From | To GOCART2G import / state | Fields |
|---|---|---|
| `CHEMENV` | `GOCART2G` | `DELP`, `AIRDENS`, `NCN_PRCP` |
| `GMICHEM` | `GOCART2G` | `OH`, `H2O2`, `NO3`, `HNO3` (renamed `GMI_*`) |
| `ACHEM` | `GOCART2G` | `pSO2_OCS` (only if `doing_OCS`), `pSOA_ANTHRO_VOC`, `pSOA_BIOB_VOC` |
| `HEMCO` | `GOCART2G` | `OC_ISOPRENE`, `OC_MTPA`, `OC_MTPO`, `OC_LIMO` |
| `GOCART2G` | `GAAS` | `AERO` (GAAS requires GOCART2G; `__raise__` otherwise) |
| `GOCART2G` | `GMICHEM` | `AERO` |
| `GOCART2G` | `CARMA` | `PSO4` as `CARMA_PSO4TOT` |
The connections follow MAPL2 semantics: a connection between two components with the same spec shares the field (no copy). Whether a lag arises depends on the run order of the two components in Physics, which is set by the Physics and Chemistry Run routines (see the Physics document (`../../../ARCHITECTURE.md` in the fixture)); this document does not assert any lag.

### 5.4 Run phases and timing
- Parent: two run phases (`Run1`, `Run2`); each non-data child: three (`Run`, `Run2`, `Run0`); each data child: one (`Run`). (`GOCART2G_GridCompMod.F90`, `SetServices`; each child's `SetServices`). Chemistry calls phase 1 only for children with more than one run phase, and in its own `Run2` calls phase 2 for those and phase 1 for the others (`../GEOS_ChemGridComp.F90`, `Run1`, `Run2`). Physics runs Chemistry phase 1 after Surface phase 1 and phase 2 before Radiation (see the Physics document (`../../../ARCHITECTURE.md` in the fixture)).
- Time step: children read `RUN_DT` and `GOCART2G_DT` (default `RUN_DT`) (`ESMF/GOCART2G_GridComp/*/*2G_GridCompMod.F90`, `Initialize`). The fixture template sets `GOCART2G_DT: @CHEM_DT` (`Applications/@GEOSgcm_App/AGCM.rc.tmpl`); `CHEM_DT` is set per resolution in `set_resolution_params.csh`.

### 5.5 Tracer friendliness
Chemistry collects "friendly" tracers from its children with `MAPL_GridCompGetFriendlies(..., AddGCPrefix=.true.)` for the lists ANALYSIS, TURBULENCE, DYNAMICS, MOIST, STRATCHEM, GMICHEM, GEOSCHEMCHEM, TR (`../GEOS_ChemGridComp.F90`, Initialize). GOCART2G internal species are declared friendly to `DYNAMICS:TURBULENCE:MOIST` in the ACG files (for example `ESMF/GOCART2G_GridComp/DU2G_GridComp/DU2G_StateSpecs.rc`, `friendlyto`). Recursion into grandchildren happens only with `AddGCPrefix` (`base/MAPL_Generic.F90`, comment).

### 5.6 Surface coupling by bin count
Surface expects a fixed number of bins per aerosol species in `AERO_DP` and stops if the count differs (`GEOSsurface_GridComp/GEOS_SurfaceGridComp.F90`, `Run2`; counts `NUM_*` in `GEOSsurface_GridComp/Shared/StieglitzSnow.F90`, relative to the GEOSgcm_GridComp repository root: DU 5, BC 2, OC 2, SU 1, SS 5). The GOCART2G children set `nbins` in their instance rc files (`ESMF/GOCART2G_GridComp/*/*_instance_*.rc`, `nbins`). See Surface (`../../../GEOSsurface_GridComp/ARCHITECTURE.md` in the fixture).

### 5.7 Other users of GOCART code
- `@GMI` calls `ReadPointEmissions` from `GOCART2G_Process` and links `Process_Library` (`@GMI/GmiEmiss_GridCompClassMod.F90`, `@GMI/CMakeLists.txt`; outside this repository).
- `@CARMA` and `@QuickChem` mention `GOCART2G_GridCompMod.F90` or `Chem_AeroGeneric.F90` only in comments.

## 6. Configuration at repository level
Only structure-changing keys; per-child keys are in the [GOCART2G document](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md).
| Key | Read in / set by | Effect |
|---|---|---|
| `ENABLE_GOCART2G`, `ENABLE_GOCART`, `ENABLE_GOCART_DATA` | `../GEOS_ChemGridComp.rc` (read in `../GEOS_ChemGridComp.F90`) | Which GOCART-family children Chemistry creates |
| `AERO_PROVIDER` | `AGCM.rc` (from `AGCM.rc.tmpl`) | Which child supplies `AERO` and `AERO_DP` |
| `ACTIVE_INSTANCES_{DU,SS,SU,CA,NI}`, `PASSIVE_INSTANCES_*` | `ESMF/GOCART2G_GridComp/GOCART2G_GridComp.rc` | Which aerosol child instances exist; names containing `data` select data-driven instances |
| `use_threads` | `ESMF/GOCART2G_GridComp/GOCART2G_GridComp.rc` | OpenMP threading of children |
| `EMISSIONS` (`AMIP`, `OPS`) | `gcm_setup`, `gcm_run.j` in the fixture | Copies the child rc/yaml variants from `AMIP/` to the run directory |
| `UFS_GOCART` | CMake option, `CMakeLists.txt` | UFS build mode |

## 7. Restarts and data
- Internal species fields of DU, SS, SU, CA, NI are `MAPL_RestartOptional` (`ESMF/GOCART2G_GridComp/*/*_StateSpecs.rc`, `restart`). Restart file names come from the fixture's `AGCM.rc.tmpl` (`du_internal_rst`, `ss_internal_rst`, `caoc_/cabc_/cabr_internal_rst`, `ni_internal_rst`, `su_internal_rst`, `gocart_internal_rst`, `gocartdata_internal_rst`).
- `ESMF/Shared/g2g_rst_generator.py` converts a legacy `gocart_internal_rst` into GOCART2G restarts (script header; *Unverified:* current use).
- Emission and climatology inputs are declared per child in `*_GridComp_ExtData.yaml`; the fixture lists them in `Applications/@GEOSgcm_App/construct_extdata_yaml_list.py`.

## 8. Known traps (repository level)
- Several directories look live but are not built or not used in the fixture: `CCPP`, `ESMF/NUOPC` (empty files), `ESMF/UFS` (only with `UFS_GOCART`), `ESMF/Aerosol_GridComp` (built, not added by Chemistry). Check `CMakeLists.txt` and `../GEOS_ChemGridComp.F90` before editing.
- There are two component families with similar names: `GOCART_GridComp` (legacy, uses `Chem_Registry`) and `GOCART2G_GridComp` (ACG-based). Their registries, restart files and `AERO` behaviour differ.
- Resource files are copied to the run directory by the fixture scripts, so edits to shipped `*.rc` values do not take effect in an existing experiment directory.
- `Process_Library` is also linked by `@GMI`; changing the interface of `GOCART2G_Process` routines affects code outside this repository.
- The `AERO_PROVIDER` default differs between layers (section 9).

## 9. Observations (noted, not fixed)
- `ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90` (`Run2`) and `Process_Library/GOCART2G_Process.F90` (`DustAerosolDistributionMeng2022`, `Chem_CalcVsettle`, `Chem_CalcVsettle2Gorig`) set `pi = 3.141529265`. The standard value is 3.14159265. Question for maintainers: is the digit order intentional?
- Default of `AERO_PROVIDER`: `../GEOS_ChemGridComp.F90` defaults to `GOCART2G`; `GEOS_GcmGridComp.F90` (outside this repository) defaults to `GOCART.data` when the value is `NULL`; `gcm_setup` writes `GOCART2G` both in its `DATA_ATMOS` branch and later unconditionally; the `AGCM.rc.tmpl` comment lists the options as "GOCART2G, MAM, none" and omits `GOCART.data`. Question for maintainers: which value is intended when the key is missing?
- Legacy path under `GOCART.data`: Chemistry connects `AERO` from `GOCARTdata` to `GMICHEM` (`../GEOS_ChemGridComp.F90`, `SetServices`). *Unverified:* how the legacy `GOCART_GridComp` supplies `AERO` to Physics when `AERO_PROVIDER` is `GOCART.data`.
- The revision history of `Process_Library/GOCART2G_Process.F90` (header comment) says the library "only uses intrinsic Fortran functions"; `Process_Library/GOCART2G_MieMod.F90` uses `netcdf`, and `Process_Library/CMakeLists.txt` depends on ESMF and NetCDF Fortran. Question for maintainers: is the header statement still current?
- `CCPP/GOCART_CCPP.F90` and `ESMF/NUOPC/NUOPC_Cap.F90` are empty tracked files. Question for maintainers: are they placeholders?
- `AMIP.20C/` rc variants are installed by each child's CMake, but a search of the fixture's `gcm_setup` found handling only for `AMIP` and `OPS` (*Unverified:* whether another script selects `AMIP.20C`).
- `ESMF/CMakeLists.txt` includes `HEMCO_GridComp`, which is not in this tree at the verified commit; `ESMF/Aerosol_GridComp/CMakeLists.txt` builds only `if EXISTS Aerosol_GridComp.F90`. *Inferred:* both guard against missing optional sources.
- Chemistry comment: `AERO_ACI` is "currently disabled in Physics ... waiting for patches in the ESMF" (`../GEOS_ChemGridComp.F90`, under `#ifdef ENABLE_AERO_ACI`).
- More component-level observations are in the [GOCART2G document](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md).
