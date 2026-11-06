# GOCART: guide for agents

**Audience:** coding agents and developers working on the GOCART repository (Goddard Chemistry, Aerosol, Radiation, and Transport model).
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GOCART` @ `96b3c9a` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`; repository `GOCART`, mounted at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GOCART`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

## What this repository is
Aerosol (and a few trace-gas) components for GEOS. On this line the aerosol component GOCART2G has been ported to MAPL3: specs are declared by ACG `*_StateSpecs.rc` files plus `MAPL_GridCompAddSpec` calls, configuration is YAML, children and connections use `MAPL_GridCompAddChild` and `MAPL_GridCompAddConnection` (`ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `SetServices`). Only the dust, sea-salt and sulfate children are ported and built. See `docs/ARCHITECTURE.md`.

> **Build status on this line:** built, partly ported. `@GEOSchem_GridComp/CMakeLists.txt` calls `esma_add_subdirectories(GOCART)`, with the comment "Only GOCART is needed on release/MAPL-v3 (for GOCART2G_GridComp)". Inside this repository `ESMF/CMakeLists.txt` builds `Shared` and `GOCART2G_GridComp`; `Apps`, `Aerosol_GridComp`, `GOCART_GridComp` and `UFS` are commented out as "not yet ported to MAPL3", and `GOCART2G_GridComp/CMakeLists.txt` comments out `CA2G_GridComp` and `NI2G_GridComp`. The Chemistry top-level library (`GEOS_ChemGridComp.F90`, which adds GOCART2G as a child) is present but not compiled on this line: the `esma_add_library` block that builds it is commented out in `../CMakeLists.txt` ("Restore subdirs as components are ported to MAPL3"). On the parent's side, `Chem_Base` and `Chem_Shared` (Chemistry `Shared` directory) are not built either. The documents of the parent repository do not exist yet on this line; the MAPL2-line versions are on its `develop` branch.

## Documents
| File | Content |
|---|---|
| `docs/ARCHITECTURE.md` | Placement, module map, build structure, interfaces to the rest of GEOS, links up and down |
| `docs/decisions.md` | Previous architectural decisions (pattern, likely reason, open questions) |
| `ESMF/GOCART2G_GridComp/ARCHITECTURE.md` | GOCART2G parent, the five aerosol children, `Process_Library`, `GA_Environment` |
| `../docs/ARCHITECTURE.md` | Parent repository `GEOSchem_GridComp` (on this line the parent's documents are not written yet; the MAPL2-line versions are on its `develop` branch) |
| `../../ARCHITECTURE.md` | Physics, the grandparent component |
| `../../../../docs/ARCHITECTURE.md` | Repository `GEOSgcm_GridComp` overview |
| `../../../../../../../docs/ARCHITECTURE.md` | GEOSgcm fixture overview |

## Rules for agents
- Read `docs/ARCHITECTURE.md` and the component document before editing.
- Note bugs and oddities; do not fix them unless asked. Record them under "Observations (noted, not fixed)" in the relevant document.
- Keep the documents current: when you change a routine, resource key or file that a document cites, update the citation in the same change.
- Cite as `(path, routine or resource key)`, never line numbers. Label unchecked statements *Inferred:* or *Unverified:*.
- GOCART2G children are `AERO`-state providers with a hand-ordered coupling to each other (MAPL2 line: the nitrate child reads the first dust, sea-salt and sulfate instances; on this line NI is not built and those connections are commented out in the parent's `SetServices`). Read "Known traps" in the GOCART2G document before adding or reordering instances. On this line CA and NI instances are not created (ignored with a warning).
- Resource files (`*.yaml`, and the `*.rc` that remain beside them) are installed to `etc/`; the instance `*.yaml` files are the ones read on this line (`ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `add_children__`). The `*.rc` files beside them are MAPL2 leftovers that the ported code does not read.
- `CCPP/GOCART_CCPP.F90` and `ESMF/NUOPC/NUOPC_Cap.F90` are empty files; the UFS cap is in `ESMF/UFS`, not built on this line.
- `ESMF/GOCART2G_GridComp/SU2G_GridComp/mapl3-porting-notes.md` records the SU2G port steps and what is left after it.

## Build and test
- Not built standalone in the GEOSgcm fixture: built as part of GEOSgcm (`CMakeLists.txt`, `esma_add_subdirectories`; `CHANGELOG.md` records that standalone builds were removed). Libraries built: `Chem_Shared2G`, `Process_Library`, `GOCART2G_GridComp` with `GA_Environment`, `DU2G`, `SS2G`, `SU2G`.
- Regression tests: `ESMF/GOCART2G_GridComp/regression` registers cases from `cases.txt` (dust, seasalt, sulfate stand-alone and a combined case, each with two runs) through `esma_add_regression_tests`; a case is skipped when its data directory is absent (`run_case.cmake`). No unit tests are registered (`ESMF/GOCART_GridComp/ut_GOCART.F90` is not referenced).
- Pull requests need a `CHANGELOG.md` entry and labels (`.github/workflows/changelog-enforcer.yml`, `enforce-labels.yml`); YAML files are linted (`.github/workflows/validate_yaml_files.yml`, `.yamllint.yml`). Code owners: `.github/CODEOWNERS`. CI builds the GEOSgcm fixture at `release/MAPL-v3` (`.github/workflows/workflow.yml`).
