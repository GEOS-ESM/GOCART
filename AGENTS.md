# GOCART: guide for agents

**Audience:** coding agents and developers working on the GOCART repository (Goddard Chemistry, Aerosol, Radiation, and Transport model).
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GOCART` @ `0d19ada` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`; repository `GOCART`, mounted at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GOCART`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

## What this repository is
Aerosol (and a few trace-gas) components for GEOS, written as MAPL2 gridded components. In the GEOSgcm fixture the live code is `ESMF/GOCART2G_GridComp` (parent plus five aerosol children), `Process_Library` and `ESMF/Shared`. Chemistry (`GEOS_ChemGridComp.F90` in the parent repository `GEOSchem_GridComp`) adds GOCART2G as a child and exports its `AERO` and `AERO_DP` states to Physics. See `docs/ARCHITECTURE.md`.

## Documents
| File | Content |
|---|---|
| `docs/ARCHITECTURE.md` | Placement, module map, build structure, interfaces to the rest of GEOS, links up and down |
| `docs/decisions.md` | Previous architectural decisions (pattern, likely reason, open questions) |
| `ESMF/GOCART2G_GridComp/ARCHITECTURE.md` | GOCART2G parent, the five aerosol children, `Process_Library`, `GA_Environment` |
| `../docs/ARCHITECTURE.md` | Parent repository `GEOSchem_GridComp` (written separately; may not exist yet) |
| `../../ARCHITECTURE.md` | Physics, the grandparent component |
| `../../../../docs/ARCHITECTURE.md` | Repository `GEOSgcm_GridComp` overview |
| `../../../../../../../docs/ARCHITECTURE.md` | GEOSgcm fixture overview |

## Rules for agents
- Read `docs/ARCHITECTURE.md` and the component document before editing.
- Note bugs and oddities; do not fix them unless asked. Record them under "Observations (noted, not fixed)" in the relevant document.
- Keep the documents current: when you change a routine, resource key or file that a document cites, update the citation in the same change.
- Cite as `(path, routine or resource key)`, never line numbers. Label unchecked statements *Inferred:* or *Unverified:*.
- GOCART2G children are `AERO`-state providers with a hand-ordered coupling to each other (the nitrate child reads the first dust, sea-salt and sulfate instances). Read "Known traps" in the GOCART2G document before adding or reordering instances.
- Resource files (`*.rc`, `*.yaml`) are installed to `etc/` and may be overridden by the run directory; check `Applications/@GEOSgcm_App` in the fixture before assuming a file's contents are used as shipped.
- `CCPP/GOCART_CCPP.F90` and `ESMF/NUOPC/NUOPC_Cap.F90` are empty files; the UFS cap is in `ESMF/UFS`.

## Build and test
- Not built standalone in the GEOSgcm fixture: the repository is built as part of GEOSgcm (`CMakeLists.txt`, `esma_add_subdirectories`; `CHANGELOG.md` records that standalone builds were removed).
- No unit tests are registered in this repository's CMake files (`ESMF/GOCART_GridComp/ut_GOCART.F90` is not referenced).
- Pull requests need a `CHANGELOG.md` entry and labels (`.github/workflows/changelog-enforcer.yml`, `enforce-labels.yml`); YAML files are linted (`.github/workflows/validate_yaml_files.yml`, `.yamllint.yml`). Code owners: `.github/CODEOWNERS`.
