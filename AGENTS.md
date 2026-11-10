# GOCART: guide for agents

**Describes:** the MAPL3 line, branch `release/MAPL-v3`; verified against `GOCART` @ `40316d5` (2026-11-10).
**Placement:** `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GOCART` in the GEOSgcm fixture (`components.yaml`); other fixtures may place it elsewhere. Paths below are relative to this repository.

> **Only GOCART2G is supported on the MAPL3 line.** `ESMF/GOCART2G_GridComp` (the parent, its DU, SS, SU and NI children, `GA_Environment`), `ESMF/Shared` (`Chem_Shared2G`) and `Process_Library` are ported and built. Everything else is MAPL2 code kept for reference and is not built: `ESMF/GOCART_GridComp`, `ESMF/Aerosol_GridComp`, `ESMF/Apps`, `ESMF/UFS` (`ESMF/CMakeLists.txt`, "not yet ported to MAPL3"), `CA2G_GridComp` (commented out in `ESMF/GOCART2G_GridComp/CMakeLists.txt`), `CCPP` and `ESMF/NUOPC` (empty files). Do not port or fix those unless asked. The Chemistry parent that would add GOCART2G (`GEOS_ChemGridComp.F90` in `GEOSchem_GridComp`) is not compiled on this line; GOCART2G runs only in the regression cases, which load it from a YAML child list.

## Rules for agents

- Read [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) (what GOCART2G is and how it is put together: structure, contracts, invariants, the background to the traps below) before changing children, connections, exports, `AERO` or configuration. Update it when the structure, a contract or an invariant changes; it describes only what changes rarely, so keep status and to-dos out of it.
- Note bugs and oddities; do not fix them unless asked. Record a bug or an open question as a GitHub issue, not in a document; a document may link to the issue.
- Do not edit the generated `*___.h` files (`*_Import___.h`, `*_Export___.h`, `*_Internal___.h`, `*_GetPointer___.h`, `*_DeclarePointer___.h`); edit the `*_StateSpecs.rc` table. Only state the table cannot express is added in `SetServices`.
- Keep this file current: update it when status, a rule or a trap changes, and update its "verified against" commit.
- Cite files and routines, never line numbers. Label unchecked statements *Inferred:* or *Unverified:*.
- Error handling in code you write or touch uses the MAPL3 macros (`@MAPL/include/MAPL_ErrLog.h`): `_RC` on calls, `_STAT` on `allocate`, `_ASSERT(cond, msg)` with a message, and `_RETURN(_SUCCESS)` at every successful exit. Convert `RC=STATUS`/`VERIFY_`/`ASSERT_`/`RETURN_` in any routine you change. `_RC` expands to two statements (`rc=status);_VERIFY(status`), so it cannot end a one-line `if`; use a block `if`. Do not add decorative `!----` separator lines.

## Status and to-dos

- CA instances are deactivated with a warning (`deactivate_instances__`), because CA2G is not built.
- `ESMF/GOCART2G_GridComp/SU2G_GridComp/mapl3-porting-notes.md` and `ESMF/GOCART2G_GridComp/NI2G_GridComp/mapl3-porting-notes.md` record the ports and what is left.
- Services: each computational DU, SS, SU and NI instance subscribes its tracers to `TRADV` (`SetServices`).
- **TODO: subscribe to `MTR` and `TRDIFF`.** On `develop` (MAPL2), these same tracers (and CA2G's `*phobic`/`*philic`) are `FRIENDLYTO` `DYNAMICS:TURBULENCE:MOIST` in each `*_StateSpecs.rc`: Moist's GF convectively transports and scavenges them, and Turbulence diffuses them. On this line only the `DYNAMICS` part is carried over (`TRADV`). Each child should also subscribe the same `service_items` to Moist's `MTR` (`tracers_for_moist`; Moist provides it, GEOSgcm_GridComp#1566) and to `TRDIFF` (`turbulently_diffused_quantities`, the name Moist uses for its own subscription). Neither needs to wait for a provider: an unconnected subscription is accepted. A `TODO` comment marks the spot after each `TRADV` subscription in `SetServices`.

## Known traps

- **Phases are run by name.** A new or renamed child phase must also be run by name from the parent; a phase registered without `phase_name` is not found.
- **Adding or reordering instances reorders `AERO` and `AERO_DP`** (`create_instances_`, `add_aero_states_`); check what depends on the order.
- **NI's connections do not check that DU, SS and SU exist.** The parent's `SetServices` reads `DU%instances(1)`, `SS%instances(1)` and `SU%instances(1)` without a size check, so a configuration that lists NI but no DU, SS or SU instance fails there. `run2` also uses the first NI instance for `nifactor` (sulfate particulate-matter sums: 1.0 when it is active, otherwise 132.14/96.06).
- **`AERO_DP` must match what Surface expects**, or Surface stops (`GEOS_SurfaceGridComp.F90`, `Run2`; counts `NUM_*` in `StieglitzSnow.F90`). Changing `nbins`, or what a child appends to its `<name>_AERO_DP`, can break it.
- **Do not remove the data-instance code paths.** Data instances will have to be supported, but on this line a child named with `data` stops in its `SetServices` (`_FAIL("data driven section has not been activated yet")` in SS2G, SU2G, NI2G; DU2G says "data driver section"). Keep `run_data`, the parent's `index(child_name, "data")` checks and the SU2G `AERO_DP` branch for that port.
- **Do not compare MAPL2 inputs with `MAPL_UNDEF`** (`-huge` in MAPL3); use the children's `MAPL2_UNDEF`. *Unverified:* what MAPL3 ExtData writes for missing data; if it is `-huge`, these masks miss it (SU2G porting notes, "Risks to settle early").
- **Wavelength lists.** `TOTANGSTR` needs 550 nm in `wavelengths_for_vertically_integrated_aop_in_nm` (`run2` asserts it); the backscatter exports need 532 nm.
- **Change a `.yaml` and its `.rc` together**, including the `AMIP/` and `AMIP.20C/` variants. Only the form may differ: the YAML adds a `mapl:` block, and `GOCART2G_GridComp.yaml` uses `ACTIVE_INSTANCES`/`PASSIVE_INSTANCES` maps where the `.rc` has one `ACTIVE_INSTANCES_<X>` key per species.

## Build and test

- Builds as part of the GEOSgcm fixture (`@GEOSchem_GridComp/CMakeLists.txt` adds only GOCART on this line) or standalone: `mepo clone` in a GOCART checkout fetches `components.yaml` (ESMA_env, ESMA_cmake, ecbuild, HEMCO, GMAO_Shared, MAPL `develop`), and `CMakeLists.txt` then builds its own project (`GOCART_STANDALONE`). The standalone build was last checked at `40316d5` (2026-11-10, ifx Debug; it also registers the regression tests), but CI builds only the fixture. Libraries: `Chem_Shared2G`, `Process_Library`, `GOCART2G_GridComp` with `GA_Environment`, `DU2G`, `SS2G`, `SU2G`, `NI2G`.
- Regression tests (label `REGRESSION`): `ESMF/GOCART2G_GridComp/regression`, one test per line of `cases.txt` (`dust-sa`, `seasalt-sa`, `sulfate-sa`, `nitrate-sa`, `gocart+seasalt+dust+sulfate`, `gocart+seasalt+dust+sulfate+nitrate`, each `run1` and `run2`). Its `README.md` describes the cases, the `DuSsStub` component used by `nitrate-sa`, how to build a case from MAPL2 captures, and `cmpchk.py`. Each case directory has a `logging.yaml` (keep it in new cases). A case whose data directory is missing passes without running (`run_case.cmake`).
- Pull requests for this line target `release/MAPL-v3`, not `develop` (the MAPL2 line), and need a `CHANGELOG.md` entry and labels (`.github/workflows/changelog-enforcer.yml`, `enforce-labels.yml`); YAML files are linted (`.github/workflows/validate_yaml_files.yml`, `.yamllint.yml`). Code owners: `.github/CODEOWNERS`. CI builds the GEOSgcm fixture at `release/MAPL-v3` (`.github/workflows/workflow.yml`).
