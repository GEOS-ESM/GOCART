# GOCART2G MAPL3 regression tests

Each test runs a single run phase of a GOCART2G component for one 20-minute
time step, starting from state captured in a MAPL2 GEOSgcm run, and compares
the resulting MAPL3 checkpoint against a stored baseline.

## Running

`CMakeLists.txt` registers one ctest per line of `cases.txt`, all with the
label `REGRESSION`:

```
ctest -L REGRESSION                # all regression tests
ctest -R 'dust-sa/run1'            # a single case
```

Every case depends on the `gocart2g/sync_data` fixture, which syncs
regression data from
`s3://gmao-usw2-si-team/regression-data/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/GOCART2G_GridComp`
into `$LOCAL_REGRESSION_DATA_DIR/<same path>`.

## What a test does

`run_case.cmake` runs the following steps for a case `<case>/runN`:

1. Copies the config directory `<case>/runN/` into a temporary experiment
   directory.
2. Copies the restarts `<data>/<case>/runN/checkpoints/2000-04-14T21:20:00/`
   into `checkpoints/` and symlinks `<data>/ExtData` as `ExtData`.
3. Runs `GEOS.x mapl.yaml` on 6 processes. The clock starts at
   `2000-04-14T21:20:00` (from `cap_restart.yaml`), and `dt` and
   `segment_duration` are both `PT20M`, so exactly one step runs.
4. On IntelLLVM builds only, compares every `*.nc` file in
   `<data>/<case>/runN/checkpoints/last/` against the new
   `checkpoints/last/` with `nccmp` (or `cmp` if `nccmp` is not found).
   NaNs compare equal, and `lons`, `lats`, `corner_lons` and `corner_lats`
   are ignored.
5. Deletes the experiment directory.

### Common setup

All cases share the following settings:

- **Grid:** C12 cubed sphere with 181 model levels, where `PLE` gives the
  pressure (set in `root.yaml`).
- **Root component:** `configurable_gridcomp`, which runs each child's
  `run_phases` from `root.yaml`/`fake-gocart.yaml` in order.
- **`DataFor<XX>` component:** a `configurable_gridcomp` whose *internal*
  state holds every import the component under test needs. Its restart
  file, `DataFor<XX>_internal.nc`, is the MAPL2 import state captured
  before the run phase. `root.yaml` connects each of those internal fields
  to the matching import of the component under test.
- **Checkpoints:** the component under test sets
  `activate_all_exports: true` and checkpoints its import, internal and
  export states. These are the files compared in step 4.
- **Logging:** `mapl.yaml` loads `logging.yaml` (`pflogger: logging.yaml`), which logs to the console only, sets the `MAPL` logger to ERROR and keeps `MAPL.profiler` at INFO on the root rank. The cases have no field dictionary, so MAPL warns once per variable spec per rank, and its default logging writes each warning to `warnings_and_errors.log` under an MPI lock, which serializes the ranks.

### `run1` vs `run2`

Each case has a `run1` and a `run2` directory, each testing one GOCART2G
run phase:

- **Run1:** emissions and sources (in `DU2G_GridCompMod.F90`, "Computes
  emissions/sources for Dust").
- **Run2:** settling, dry deposition, wet removal, chemistry (SU) and
  diagnostics.

The two directories differ only in the selected phase and, because the
phases read different fields, in the imports that `DataFor<XX>` provides.

## Cases

### Standalone (`*-sa`) cases

In these cases one aerosol component runs on its own:

```
Root
├── DataFor<XX>
└── FakeGOCART            (configurable_gridcomp)
    └── <XX>              (lib<XX>2G_GridComp.so, run_phases: [RunN])
```

`FakeGOCART` stands in for the real GOCART2G parent, so that
`fake-gocart.yaml` can select a single run phase of the aerosol component.
Because MAPL connections can only target a component's direct children,
`root.yaml` connects to `FakeGOCART`. MAPL passes the aerosol component's
unsatisfied imports up to `FakeGOCART`, and that is how the connections
reach the aerosol component.

| Case | Component | Instance config | Imports |
|---|---|---|---|
| `dust-sa` | DU2G (`DU`) | `DU_instance_DU.yaml` | `data-for-du.yaml` |
| `seasalt-sa` | SS2G (`SS`) | `SS_instance_SS.yaml` | `data-for-ss.yaml` |
| `sulfate-sa` | SU2G (`SU`) | `SU_instance_SU.yaml` | `data-for-su.yaml` |
| `nitrate-sa` | NI2G (`NI`) | `NI_instance_NI.yaml` | `data-for-ni.yaml` |
| `carbon-bc-sa` | CA2G (`CA_bc`) | `CA_instance_CA_bc.yaml` | `data-for-ca.yaml` |

The restarts are `DataFor<XX>_internal.nc` and `<XX>_internal.nc`. The
baseline consists of `<XX>_{import,internal,export}.nc`.

`carbon-bc-sa` is the one case where the `DataFor` component is not named
after the child: the child is the CA2G *instance* `CA_bc`, so the restarts
are `DataForCA_internal.nc` and `CA_bc_internal.nc`. Only the black-carbon
instance is exercised; `CA_oc` and `CA_br` are ported and runnable but are
not yet part of the suite.

`nitrate-sa` has one more child, `DuSsStub` (`../DuSsStub_GridComp`, built only for this case):

```
Root
├── DataForNI
├── DuSsStub              (libDuSsStub_GridComp.so)
└── FakeGOCART            (configurable_gridcomp)
    └── NI                (libNI2G_GridComp.so, run_phases: [RunN])
```

NI2G's Initialize reads the `radius` and `fnum` ESMF Info attributes from its `DU` and `SS` imports. Only DU2G's and SS2G's own Initialize set them, and a `configurable_gridcomp` cannot. `DuSsStub` imports `DU` and `SS` from `DataForNI`, sets the attributes as DU2G and SS2G do, and exports them; `root.yaml` connects those exports to `FakeGOCART`.

### `gocart+seasalt+dust+sulfate`

This case runs the real GOCART2G parent with its DU, SS and SU children:

```
Root
├── DataForGOCART2G       (run_phases: [run])
└── GOCART2G              (libGOCART2G_GridComp.so, run_phases: [RunN])
    ├── DU
    ├── SS
    └── SU
```

- **Configs:** `GOCART2G_GridComp.yaml`, plus `DU_instance_DU.yaml`,
  `SS_instance_SS.yaml` and `SU_instance_SU.yaml`.
- **Imports:** `data-for-gocart2g.yaml` holds the combined imports of all
  three aerosols.
- **Restarts:** `DataForGOCART2G_internal.nc`, `DU_internal.nc`,
  `SS_internal.nc` and `SU_internal.nc`.
- **Baseline:** `{DU,SS,SU}_{import,internal,export}.nc` and
  `GOCART2G_{import,export}.nc`.

### `gocart+seasalt+dust+sulfate+nitrate`

The same as `gocart+seasalt+dust+sulfate`, with an NI child added:

- **Configs:** `GOCART2G_GridComp.yaml` lists `NI: [NI]` (an empty list in the case above); `NI_instance_NI.yaml` is added.
- **Imports:** `data-for-gocart2g.yaml` and `root.yaml` add NI's imports, such as the ammonia emissions `EMI_NH3_AG`, `EMI_NH3_BB` and `EMI_NH3_EN`.
- **Restarts:** adds `NI_internal.nc`.
- **Baseline:** adds `NI_{import,internal,export}.nc`.

### `gocart+seasalt+dust+sulfate+nitrate+carbon-bc`

The same as `gocart+seasalt+dust+sulfate+nitrate`, with a CA child added:

- **Configs:** `GOCART2G_GridComp.yaml` lists `CA: [CA_bc]` (an empty list
  in the case above); `CA_instance_CA_bc.yaml` is added. Only the
  black-carbon instance is active, so the child is `CA_bc` rather than a
  bare species name.
- **Imports:** `data-for-gocart2g.yaml` and `root.yaml` add CA's imports,
  such as the carbon emissions `BC_BIOMASS`, `BC_BIOFUEL` and
  `BC_ANTEBC1`.
- **Restarts:** adds `CA_bc_internal.nc`.
- **Baseline:** adds `CA_bc_{import,internal,export}.nc`.

This is the only case that exercises the parent's carbonaceous connection
block (`GOCART2G_GridCompMod.F90`), including `CA_bc_AERO` and
`CA_bc_AERO_DP`; `carbon-bc-sa` runs against `FakeGOCART` and so does not.

## Creating a case from MAPL2 captures

The input data comes from a MAPL2 run instrumented to dump each
component's states before and after each run phase. The dumps are then
converted to the new cubed-sphere layout with
`GEOS_Restart_Utilities/convertOldCStoNewCS.py`.

### Baseline GEOSgcm

The captures come from GEOSgcm `v11.11.0` (MAPL `v2.71.0`), built with ifx
rather than the default ifort, with these main components:

| Component | Version |
|---|---|
| GEOSgcm | `v11.11.0` |
| MAPL | `v2.71.0` |
| env | `v5.26.0`, with `g5_modules` edited for ifx (see below) |
| cmake | `v4.51.0` |
| ecbuild | `geos/v3.15.2` |
| GEOSgcm_GridComp | `v3.1.0` |
| GEOSchem_GridComp | `v2.1.0` |
| GOCART | `fd894af`, plus local changes kept in a stash |
| FVdycoreCubed_GridComp / fvdycore | `v3.1.0` / `geos/v3.1.0` |
| GEOSradiation_GridComp | `develop` at `0b402f6`, with a local edit to `GEOS_SolarGridComp.F90` (see below) |
| RRTMGP | `geos/v1.8+1.0.0` |
| GMAO_Shared | `v3.0.3` |
| GEOS_Util | `v3.0.4` |
| NCEP_Shared | `v1.5.0` |
| GEOSgcm_App | `v3.1.0` |

All other components were at the tags pinned by GEOSgcm `v11.11.0`.

The run was on discover (NCCS), where the `g5_modules` edit switches from
ifort 2021.13 / Intel MPI 2021.13 to ifx 2025.3.0 / Intel MPI 2021.17:

```
basedir = /discover/swdev/gmao_SIteam/Baselibs/ESMA-Baselibs-8.33.0/x86_64-pc-linux-gnu/ifx_2025.3.0-intelmpi_2021.17
comp/gcc/12.3.0
comp/intel/2025.3.0
mpi/impi/2021.17
other/jemalloc/5.3.0
other/ISSM/geos/1.0.3/ifx_2025.3.0-intelmpi_2021.17
```

The local edit to `GEOSsolar_GridComp/GEOS_SolarGridComp.F90` applies fix 1
from
[GEOSradiation_GridComp#110](https://github.com/GEOS-ESM/GEOSradiation_GridComp/issues/110):
it deletes the two `have_rrtmgp_wavenums = .true.` lines. Without the fix,
`UPDATE_EXPORT` reads the RRTMGP band limits uninitialized on every call
after the first. With it, the limits are reloaded on every call.

### Mapping captures to case files

For component `<XX>` and phase `N`, the converted files go into a case as
follows:

| MAPL2 capture | Destination |
|---|---|
| `<XX>_import_before_runPhaseN.nc` | `runN/checkpoints/2000-04-14T21:20:00/DataFor<XX>_internal.nc` |
| `<XX>_internal_before_runPhaseN.nc` | `runN/checkpoints/2000-04-14T21:20:00/<XX>_internal.nc` |
| `<XX>_*_after_runPhaseN.nc` | `baseline/` (for `cmpchk.py`) |

The ctest baseline in `checkpoints/last/` is not the MAPL2 capture but a
MAPL3 checkpoint from a known-good IntelLLVM run.

The aerosol optics tables go under
`ExtData/chemistry/AerosolOptics/v2.x.x/x/`, where `aerosol_*_optics_file`
in each instance yaml points.

## Comparing against MAPL2 by hand (`cmpchk.py`)

Each run directory has a `cmpchk.py`, not used by the ctest, that compares
a MAPL3 checkpoint with the MAPL2 states captured *after* the same run
phase. Run it from the experiment directory, which must contain:

- `baseline/<XX>_<state>_after_runPhaseN.nc`: the MAPL2 captures.
- `checkpoints/last/<XX>_<state>.nc`: the MAPL3 checkpoint.

For each of the import, internal and export states, it prints the min, max
and norm of the difference for every variable present in both files.

### Validating a case before updating regression-data

Before a new or changed case goes into regression-data, check it against
MAPL2:

1. Copy the MAPL2 "before" captures, renamed as in the mapping table
   above, into `checkpoints/2000-04-14T21:20:00/` and run the case by hand.
2. Copy the MAPL2 "after" captures into `baseline/` and run `cmpchk.py`.
3. If every difference is zero, add the renamed MAPL2 inputs to
   regression-data as `runN/checkpoints/2000-04-14T21:20:00/`, and the
   MAPL3 outputs as `runN/checkpoints/last/`.
