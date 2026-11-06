# GOCART: architecture

**Audience:** coding agents and developers working on the GOCART repository (Goddard Chemistry, Aerosol, Radiation, and Transport model).
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GOCART` @ `96b3c9a` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`; repository `GOCART`, mounted at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GOCART`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths below are relative to the GOCART repository root unless stated otherwise.
Component detail: [GOCART2G](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md). Agent rules: [AGENTS.md](../AGENTS.md). Decision log: [decisions.md](decisions.md).
> **Build status on this line:** built, partly ported. `@GEOSchem_GridComp/CMakeLists.txt` (`../CMakeLists.txt`) calls `esma_add_subdirectories(GOCART)` under the comment "Only GOCART is needed on release/MAPL-v3 (for GOCART2G_GridComp)". Within this repository, `ESMF/CMakeLists.txt` builds `Shared` and `GOCART2G_GridComp`; `Apps`, `Aerosol_GridComp`, `GOCART_GridComp` and `UFS` are commented out as "not yet ported to MAPL3 - restore when ported"; `ESMF/GOCART2G_GridComp/CMakeLists.txt` builds `GA_Environment`, `DU2G_GridComp`, `SS2G_GridComp`, `SU2G_GridComp` and comments out `CA2G_GridComp` and `NI2G_GridComp`. The Chemistry top-level library (`GEOS_ChemGridComp.F90`, which adds this component) is present but not compiled on this line: the `esma_add_library` block that builds it is commented out in `../CMakeLists.txt`; it says "Restore subdirs as components are ported to MAPL3". Statements below about Chemistry describe that file as it is (MAPL2 code, `MAPL_AddChild`, ESMF Config), not running code. The Chemistry `Shared` directory (`Chem_Base`, `Chem_Shared`) is also not built (`../CMakeLists.txt`: `Shared` appears only in a commented-out list).

Up (on this line the parent's documents are not written yet; the MAPL2-line versions are on its `develop` branch):
- Parent repository `GEOSchem_GridComp`: `../../docs/ARCHITECTURE.md`.
- Physics: `../../../ARCHITECTURE.md`.
- Moist, a consumer of the aerosol activation method: `../../../GEOSmoist_GridComp/ARCHITECTURE.md`.
- Repository `GEOSgcm_GridComp`: `../../../../../docs/ARCHITECTURE.md`.
- Fixture: `../../../../../../../../docs/ARCHITECTURE.md`.

## 1. Purpose
GOCART holds the aerosol components of GEOS. The code that is built on this line is `GOCART2G`, a MAPL3 gridded component, with three ported aerosol-species children (dust, sea salt, sulfate) plus two shared libraries (`ESMF/CMakeLists.txt`; `ESMF/GOCART2G_GridComp/CMakeLists.txt`). The carbonaceous and nitrate children, the older component (`GOCART_GridComp`, trace gases and legacy aerosol) and the cap code for UFS (the NOAA Unified Forecast System) are in the tree but not built. Where the component is added to a model, it passes its aerosol optical and microphysical states up through the exported `AERO` state (`ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `SetServices`). The code that adds it in the GEOSgcm fixture, Chemistry (`../GEOS_ChemGridComp.F90`), is not compiled on this line (see the build status above).

## 2. Place in the hierarchy (GEOSgcm fixture)
```
GCM (GEOSgcm_GridComp)
  AGCM -> Physics (GEOSphysics_GridComp)
            CHEM (GEOSchem_GridComp, GEOS_ChemGridComp.F90; not compiled on this line)
              GOCART2G   (this repository, ESMF/GOCART2G_GridComp)
                DU, SS, SU   (one child per instance name; CA.* and NI not created on this line)
              GOCART, GOCART.data  (this repository, ESMF/GOCART_GridComp; legacy, not built)
              other chemistry children (not built on this line)
```
- On this line GOCART2G is instantiated by name: the shared-library entry point is `gocart2g_setservices_` (the external wrapper `GOCART2G_SetServices` at the end of `ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`), and the regression cases load it from `libGOCART2G_GridComp.so` with a YAML child list (`ESMF/GOCART2G_GridComp/regression/gocart+seasalt+dust+sulfate/run1/root.yaml`, `children: GOCART2G`).
- The Chemistry file `../GEOS_ChemGridComp.F90` is the MAPL2 version: it `use`s `GOCART2G_SetServices` and calls `MAPL_AddChild(GC, NAME='GOCART2G', SS=GOCART2G_SetServices)` when `ENABLE_GOCART2G` is true (`SetServices`). It is not compiled here, so no running component adds GOCART2G in the fixture on this line.
- The legacy `GOCART` and `GOCART.data` children are added by the same uncompiled Chemistry file under `ENABLE_GOCART` and `ENABLE_GOCART_DATA` (`../GEOS_ChemGridComp.F90`, `SetServices`).
- Provider selection (`AERO_PROVIDER`, `GetProvider_`) is also in that uncompiled file. See section 5.
- The parent component repository `GEOSchem_GridComp` has no architecture document on this line yet.

## 3. Module map
| Directory | Content | Status in the GEOSgcm fixture |
|---|---|---|
| `ESMF/GOCART2G_GridComp` | Component `GOCART2G`: parent (`GOCART2G_GridCompMod.F90`) and children `DU2G_GridComp`, `SS2G_GridComp`, `SU2G_GridComp`, `CA2G_GridComp`, `NI2G_GridComp`; `GA_Environment` (base type for children). [Detail](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md) | Built (`ESMF/CMakeLists.txt`). Ported to MAPL3 for `GA_Environment`, `DU2G`, `SS2G`, `SU2G`; `CA2G_GridComp` and `NI2G_GridComp` directories hold MAPL2 code and are commented out of `GOCART2G_GridComp/CMakeLists.txt` |
| `Process_Library` | Library `Process_Library`: physical process routines (`GOCART2G_Process.F90`) and Mie-table code (`GOCART2G_MieMod.F90`) (`Process_Library/CMakeLists.txt`) | Built (`CMakeLists.txt`: `add_subdirectory (Process_Library)`); used by GOCART2G; `@GMI` also links it but is not built on this line |
| `ESMF/Shared` | Library `Chem_Shared2G`: `Chem_AeroGeneric.F90`, `ReplenishAlarm.F90`, `GOCART2G_AerosolActivationMod.F90`; `g2g_rst_generator.py` (`ESMF/Shared/CMakeLists.txt`) | Built; `Chem_AeroGeneric.F90`, `ReplenishAlarm.F90`, `GOCART2G_AerosolActivationMod.F90` use the MAPL3 API |
| `ESMF/Aerosol_GridComp` | Module `Aerosol_GridComp_mod` (`Aerosol_GridComp.F90`): a wrapper that adds `GOCART2G` as its only child (`MAPL_AddChild`) and runs all phases of all children; ACG file `Aerosol_StateSpecs.rc` | Not built: commented out in `ESMF/CMakeLists.txt` ("not yet ported to MAPL3"). The source is the MAPL2 code |
| `ESMF/UFS` | NUOPC cap (`Aerosol_Cap.F90`, which uses `Aerosol_GridComp_mod`) and `Aerosol_*` support modules; library `UFS_Aerosols` | Not built: in `ESMF/CMakeLists.txt` the `UFS` entry inside the `UFS_GOCART` block is commented out ("not yet ported to MAPL3"), as are the other UFS-mode subdirectories listed there |
| `ESMF/GOCART_GridComp` | Legacy `GOCART`: `GOCART_GridCompMod.F90`, `Aero_GridCompMod.F90`, `GOCART_StubCompMod.F90`; trace-gas children `O3_GridComp`, `CO_GridComp`, `CO2_GridComp`, `CFC_GridComp`, `Rn_GridComp`, `CH4_GridComp`; `GOCARTdata_AerRegistry.rc` | Not built: commented out in `ESMF/CMakeLists.txt`. The source is the MAPL2 code |
| `ESMF/Apps` | Offline tool `aop_calculator.x` (aerosol optical properties) using `GOCART2G_AopMod`, `GOCART2G_SimpleBundleMod`, `GOCART2G_MieMod` (`ESMF/Apps/CMakeLists.txt`) | Not built: commented out in `ESMF/CMakeLists.txt` |
| `ESMF/HEMCO_GridComp` | Not in this repository tree at the verified commit; `components.yaml` fetches HEMCO there (`components.yaml`, `HEMCO`, branch `geos/release/MAPL-v3`). `ESMF/CMakeLists.txt` no longer lists it. In the fixture, `HEMCO_GridComp` is a directory of the parent repository, commented out of its build (`../CMakeLists.txt`) | Not built; not documented here |
| `CCPP` | `GOCART_CCPP.F90`: empty file (0 bytes, tracked) | Not built; no CMake reference |
| `ESMF/NUOPC` | `NUOPC_Cap.F90`: empty file (0 bytes, tracked) | Not built; no CMake reference |
| `config` | `GMAO_Shared.sparse`: sparse checkout list for the `GMAO_Shared` dependency | Used by `components.yaml` |
| `.github` | `CODEOWNERS`, workflows (changelog, labels, YAML lint, CI build) | Repository tooling |

`README.md` states the framework interfaces as "MAPL, NUOPC, and CCPP" (`README.md`). *Unverified:* the CCPP and NUOPC directories are placeholders for planned work; their files are empty.

## 4. Build structure
- Top level (`CMakeLists.txt`): in the fixture the repository is a subproject. A standalone build is possible only when an ESMA cmake directory exists; option `UFS_GOCART` (default OFF) forces the standalone path and a 32-bit build. It sets `DOING_GEOS5`, then `add_subdirectory(ESMF)` and `add_subdirectory(Process_Library)`. `CHANGELOG.md` records that standalone GOCART builds were removed. *Unverified:* the extent to which the standalone path still works.
- `ESMF/CMakeLists.txt`: the only active entries are `Shared` and `GOCART2G_GridComp`; `Apps`, `Aerosol_GridComp`, `GOCART_GridComp` and (in the `UFS_GOCART` block) `UFS` are commented out with "not yet ported to MAPL3 - restore when ported". It sets `HEMCO_EXTERNAL_CONFIG`, `MAPL_ESMF` and `BUILD_GEOS_INTERFACE`. (`ESMF/CMakeLists.txt`)
- Libraries and dependencies (CMake files in each directory):
  - `Process_Library`: depends on ESMF and NetCDF Fortran.
  - `Chem_Shared2G` (`esma_set_this(OVERRIDE Chem_Shared2G)`): built from three files; depends on MAPL and ESMF; used by GOCART2G.
  - `GOCART2G_GridComp`: `esma_set_this` with SUBCOMPONENTS `GA_Environment DU2G_GridComp SS2G_GridComp SU2G_GridComp` (`CA2G_GridComp`, `NI2G_GridComp` commented out); depends on `Chem_Shared2G`, `Process_Library`, ESMF; type `SHARED`; links OpenMP unless the compiler is NAG; installs `GOCART2G_GridComp.yaml` to `etc`. It also adds the `regression` subdirectory.
  - Each ported child's CMake runs `mapl_acg` (the MAPL Automatic Code Generator) on its `*_StateSpecs.rc` with `IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS`, and installs every `*.rc` and `*.yaml` in its directory, including the `AMIP/` and `AMIP.20C/` subdirectories, to `etc/<directory>`.
- Regression tests: `ESMF/GOCART2G_GridComp/regression/CMakeLists.txt` calls `esma_add_regression_tests(NAME gocart2g ...)` over the cases in `cases.txt` (`dust-sa`, `seasalt-sa`, `sulfate-sa`, `gocart+seasalt+dust+sulfate`, each `run1` and `run2`). `run_case.cmake` skips a case when its regression data directory is absent, runs the case on 6 processes, and compares checkpoints with baselines only for the IntelLLVM compiler.
- Direct dependency declarations for the repository (`components.yaml`): `env`, `cmake`, `ecbuild`, `HEMCO` (`ESMF/HEMCO_GridComp@`, branch `geos/release/MAPL-v3`), `GMAO_Shared` (sparse, branch `release/MAPL-v3`), `MAPL` (`ESMF/Shared/MAPL@`, branch `develop`). *Inferred:* in the GEOSgcm fixture the fixture's own `components.yaml` supplies MAPL and GMAO_Shared and these entries are used only when GOCART is cloned alone.
- The parent `@GEOSchem_GridComp/CMakeLists.txt` adds only `GOCART` on this line (`esma_add_subdirectories(GOCART)`); the block that would build the Chemistry library, with dependencies `Chem_Shared Chem_Base GOCART_GridComp GOCART2G_GridComp ...`, is commented out (outside this repository).
- Tests: no unit test is registered; `ESMF/GOCART_GridComp/ut_GOCART.F90` is not referenced.
- Continuous integration (`.github/workflows/workflow.yml`) builds the GEOSgcm fixture at `release/MAPL-v3` with the `ifx` and `gfortran-15` compilers (Debug), and a spack build.
- Porting record: `ESMF/GOCART2G_GridComp/SU2G_GridComp/mapl3-porting-notes.md` lists the SU2G port steps and the items left after it (advection, diffusion and convection services, NI2G inputs).

## 5. Interfaces to the rest of GEOS
Everything below is on the Chemistry side unless noted. Chemistry file: `../GEOS_ChemGridComp.F90`, which on this line is the MAPL2 code (`MAPL_AddChild`, `MAPL_AddConnectivity`, ESMF Config) and is not compiled. The tables describe that file as it is; they are not connections that exist in a built MAPL3 model. The GOCART2G side below is the MAPL3 code.

### 5.1 Aerosol states exported upward
- GOCART2G declares the exports `AERO` (an ESMF state item, `MAPL_STATEITEM_STATE`) and `AERO_DP` (`MAPL_STATEITEM_FIELDBUNDLE`) with `MAPL_GridCompAddSpec`; it also declares per-child `<name>_AERO` and `<name>_AERO_DP` exports and connects each child's items to itself (`ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `SetServices`).
- Chemistry (uncompiled, MAPL2 code) promotes both from the provider child: `MAPL_AddExportSpec(..., CHILD_ID=AERO_PROVIDER)` (`../GEOS_ChemGridComp.F90`, `SetServices`).
- Physics (MAPL2-style `../../GEOS_PhysicsGridComp.F90` in the fixture; see the Physics document, if present on the line) connects `AERO` from CHEM to Radiation and to Moist, and `AERO_DP` to Surface (`../../GEOS_PhysicsGridComp.F90`, `SetServices`; see the Physics document (`../../../ARCHITECTURE.md` in the fixture)).
- `AERO` carries ESMF methods that consumers call, not only data: `run_aerosol_optics` (Radiation, GMI photolysis), `get_mixRatioSum`, `serialize_bundle`, `get_monochromatic_aop` (GAAS), `aerosol_activation_properties` (Moist) (`GOCART2G_GridCompMod.F90`, `Initialize`, `ESMF_MethodAdd`; `ESMF/Shared/GOCART2G_AerosolActivationMod.F90`, `setup_aerosol_activation`). Details in the [GOCART2G document](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md).
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
These are MAPL2 connections in an uncompiled file; this document does not assert how they behave on the MAPL3 line. Inside GOCART2G, the MAPL3 connections between siblings are made with `MAPL_GridCompAddConnection` (section 2 of the GOCART2G document).

### 5.4 Run phases and timing
- Parent: two run phases, named `Run1` and `Run2` (`MAPL_GridCompSetEntryPoint(..., phase_name=...)`); each non-data child: `Run1`, `Run2`, `Run0`. The parent runs the children by phase name with `MAPL_GridCompRunChild`: `Run1` in its `run1`, `Run0` then `Run2` in its `run2` (`GOCART2G_GridCompMod.F90`, `SetServices`, `run1`, `run2`; `DU2G_GridCompMod.F90`, `SetServices`). Data-driven children are not supported on this line: each ported child fails with "data driv... section has not been activated yet" in `SetServices` (`DU2G_GridCompMod.F90`, `SS2G_GridCompMod.F90`, `SU2G_GridCompMod.F90`). In the uncompiled Chemistry file, Chemistry calls phase 1 only for children with more than one run phase, and in its own `Run2` calls phase 2 for those and phase 1 for the others (`../GEOS_ChemGridComp.F90`, `Run1`, `Run2`). Physics runs Chemistry phase 1 after Surface phase 1 and phase 2 before Radiation (see the Physics document (`../../../ARCHITECTURE.md` in the fixture)).
- Time step: ported children read `GOCART2G_DT` with default the heartbeat step (`MAPL_GridCompGetResource(gc, "GOCART2G_DT", CDT, default=real(HDT))`) (`ESMF/GOCART2G_GridComp/DU2G_GridComp/DU2G_GridCompMod.F90`, `Initialize`).

### 5.5 Tracer friendliness
In the uncompiled Chemistry file, Chemistry collects "friendly" tracers from its children with `MAPL_GridCompGetFriendlies(..., AddGCPrefix=.true.)` for the lists ANALYSIS, TURBULENCE, DYNAMICS, MOIST, STRATCHEM, GMICHEM, GEOSCHEMCHEM, TR (`../GEOS_ChemGridComp.F90`, Initialize). GOCART2G internal species still carry `FRIENDLYTO` `DYNAMICS:TURBULENCE:MOIST` in the ACG files (for example `ESMF/GOCART2G_GridComp/DU2G_GridComp/DU2G_StateSpecs.rc`, `FRIENDLYTO`). In addition, each ported child declares an import `TRADV` of item type `MAPL_STATEITEM_SERVICE` listing its advected tracers, with the comment "Replaces MAPL2's FRIENDLYTO=DYNAMICS" (`DU2G_GridCompMod.F90`, `SS2G_GridCompMod.F90`, `SU2G_GridCompMod.F90`, `SetServices`, `TRADV`). `ESMF/GOCART2G_GridComp/SU2G_GridComp/mapl3-porting-notes.md` says the `TRADV` subscription takes effect once the AGCM and Physics are ported. *Unverified:* how the `TURBULENCE` and `MOIST` parts of `FRIENDLYTO` are used on this line (the notes say diffusion and convection need MAPL3 service providers).

### 5.6 Surface coupling by bin count
Surface expects a fixed number of bins per aerosol species in `AERO_DP` and stops if the count differs (`GEOSsurface_GridComp/GEOS_SurfaceGridComp.F90`, `Run2`; counts `NUM_*` in `GEOSsurface_GridComp/Shared/StieglitzSnow.F90`, relative to the GEOSgcm_GridComp repository root: DU 5, BC 2, OC 2, SU 1, SS 5). The ported GOCART2G children set `nbins` in their instance YAML files (`ESMF/GOCART2G_GridComp/*/*_instance_*.yaml`, `nbins`). See Surface (`../../../GEOSsurface_GridComp/ARCHITECTURE.md` in the fixture).

### 5.7 Other users of GOCART code
- `@GMI` calls `ReadPointEmissions` from `GOCART2G_Process` and links `Process_Library` (`@GMI/GmiEmiss_GridCompClassMod.F90`, `@GMI/CMakeLists.txt`; outside this repository; `@GMI` is not built on this line).
- `@CARMA` and `@QuickChem` mention `GOCART2G_GridCompMod.F90` or `Chem_AeroGeneric.F90` only in comments.

## 6. Configuration at repository level
Only structure-changing keys; per-child keys are in the [GOCART2G document](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md). On this line the GOCART2G configuration is YAML, read through the component's `hconfig` and `MAPL_GridCompGetResource`.
| Key | Read in / set by | Effect |
|---|---|---|
| `ENABLE_GOCART2G`, `ENABLE_GOCART`, `ENABLE_GOCART_DATA` | `../GEOS_ChemGridComp.rc` (read in `../GEOS_ChemGridComp.F90`; not compiled on this line) | Which GOCART-family children Chemistry would create |
| `AERO_PROVIDER` | `AGCM.rc` (from `AGCM.rc.tmpl`), read by the uncompiled Chemistry file | Which child supplies `AERO` and `AERO_DP` |
| `ACTIVE_INSTANCES`, `PASSIVE_INSTANCES` (YAML maps keyed `DU`, `SS`, `SU`, `CA`, `NI`, each a list of instance names) | `ESMF/GOCART2G_GridComp/GOCART2G_GridComp.yaml` (read in `GOCART2G_GridCompMod.F90`, `setup_constituents_`) | Which aerosol child instances exist. Names containing `data` select data-driven instances, which fail in `SetServices` on this line; `CA` and `NI` instances are marked inactive and not created |
| `use_threads` | none | The read of `use_threads` is commented out in the parent ("Do we re-implement threading?", `GOCART2G_GridCompMod.F90`, `SetServices`) |
| `EMISSIONS` (`AMIP`, `OPS`) | `gcm_setup`, `gcm_run.j` in the fixture | Copies the child rc/yaml variants from `AMIP/` to the run directory |
| `UFS_GOCART` | CMake option, `CMakeLists.txt` | UFS build mode |

## 7. Restarts and data
- Internal species fields of DU, SS, SU are `MAPL_RestartOptional` (`ESMF/GOCART2G_GridComp/*/*_StateSpecs.rc`, `RESTART`). The parent's configuration sets `mapl: misc: restart: import: true` and `checkpoint: export: true, import: true` (`ESMF/GOCART2G_GridComp/GOCART2G_GridComp.yaml`). The regression cases read restarts and compare checkpoints (`ESMF/GOCART2G_GridComp/regression/run_case.cmake`, `copy_restarts`, `compare_results`). Restart file naming in a full model is set outside this repository.
- `ESMF/Shared/g2g_rst_generator.py` converts a legacy `gocart_internal_rst` into GOCART2G restarts (script header; *Unverified:* current use).
- Emission and climatology inputs are declared per child in `*_GridComp_ExtData.yaml`; the fixture lists them in `Applications/@GEOSgcm_App/construct_extdata_yaml_list.py`.

## 8. Known traps (repository level)
- Several directories look live but are not built on this line: `CCPP`, `ESMF/NUOPC` (empty files), `ESMF/UFS`, `ESMF/Aerosol_GridComp`, `ESMF/GOCART_GridComp`, `ESMF/Apps`, and the `CA2G_GridComp` and `NI2G_GridComp` directories (MAPL2 code). Check `ESMF/CMakeLists.txt` and `ESMF/GOCART2G_GridComp/CMakeLists.txt` before editing.
- Two files with the same role coexist in each ported child directory: `<X>2G_instance_<name>.yaml` (read on this line) and the MAPL2 `<X>2G_instance_<name>.rc` (not read by the ported code; `GOCART2G_GridComp.rc` likewise beside `GOCART2G_GridComp.yaml`).
- There are two component families with similar names: `GOCART_GridComp` (legacy, uses `Chem_Registry`) and `GOCART2G_GridComp` (ACG-based). Their registries, restart files and `AERO` behaviour differ.
- Resource files may be copied to a run directory by scripts outside this repository, so edits to shipped values may not take effect in an existing experiment directory. *Unverified* on this line.
- `Process_Library` is also linked by `@GMI` (not built on this line); changing the interface of `GOCART2G_Process` routines affects code outside this repository.
- The `AERO_PROVIDER` default differs between layers (section 9; files outside this repository, not compiled on this line).

## 9. Observations (noted, not fixed)
- `ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90` (`run2`, next to a comment "TODO - use MAPL_PI instead??") and `Process_Library/GOCART2G_Process.F90` (`DustAerosolDistributionMeng2022`, `Chem_CalcVsettle`, `Chem_CalcVsettle2Gorig`) set `pi = 3.141529265`. The standard value is 3.14159265. Question for maintainers: is the digit order intentional?
- Default of `AERO_PROVIDER`: `../GEOS_ChemGridComp.F90` defaults to `GOCART2G`; `GEOS_GcmGridComp.F90` (outside this repository) defaults to `GOCART.data` when the value is `NULL`; `gcm_setup` writes `GOCART2G` both in its `DATA_ATMOS` branch and later unconditionally; the `AGCM.rc.tmpl` comment lists the options as "GOCART2G, MAM, none" and omits `GOCART.data`. Question for maintainers: which value is intended when the key is missing?
- Legacy path under `GOCART.data`: Chemistry connects `AERO` from `GOCARTdata` to `GMICHEM` (`../GEOS_ChemGridComp.F90`, `SetServices`). *Unverified:* how the legacy `GOCART_GridComp` supplies `AERO` to Physics when `AERO_PROVIDER` is `GOCART.data`.
- The revision history of `Process_Library/GOCART2G_Process.F90` (header comment) says the library "only uses intrinsic Fortran functions"; `Process_Library/GOCART2G_MieMod.F90` uses `netcdf`, and `Process_Library/CMakeLists.txt` depends on ESMF and NetCDF Fortran. Question for maintainers: is the header statement still current?
- `CCPP/GOCART_CCPP.F90` and `ESMF/NUOPC/NUOPC_Cap.F90` are empty tracked files. Question for maintainers: are they placeholders?
- `AMIP.20C/` rc variants are installed by each child's CMake, but a search of the fixture's `gcm_setup` found handling only for `AMIP` and `OPS` (*Unverified:* whether another script selects `AMIP.20C`).
- `ESMF/CMakeLists.txt` still sets `HEMCO_EXTERNAL_CONFIG TRUE` although no HEMCO subdirectory is listed; `components.yaml` still fetches HEMCO into `ESMF/HEMCO_GridComp@`.
- `ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90` (`SetServices`) reads `aerosol_photolysis_wavelengths_in_nm_from_LUT`, while `GOCART2G_GridComp.yaml` defines `aerosol_photolysis_wavelength_in_nm_from_LUT` and the children read the singular form (`DU2G_GridCompMod.F90`, `Initialize`). Question for maintainers: is the parent's plural key intended?
- Chemistry comment (in the uncompiled file): `AERO_ACI` is "currently disabled in Physics ... waiting for patches in the ESMF" (`../GEOS_ChemGridComp.F90`, under `#ifdef ENABLE_AERO_ACI`).
- More component-level observations are in the [GOCART2G document](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md).
