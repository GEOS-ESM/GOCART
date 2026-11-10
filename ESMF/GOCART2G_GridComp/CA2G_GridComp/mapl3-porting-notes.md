# CA2G_GridComp MAPL3 porting plan

Status: **code-complete, built, and run-tested end to end.** Ported on top
of `release/MAPL-v3`, following the SU2G port
(`../SU2G_GridComp/mapl3-porting-notes.md`) as the primary template for the
single-state-bundle `_AERO`/`_AERO_DP` export pattern, with NI2G
(`../NI2G_GridComp/mapl3-porting-notes.md`) as the most recent precedent for
multi-instance parent wiring and the "GOCART2G already anticipated this
species" starting point.

Scope of this pass, agreed up front:
- Full port: code + GOCART2G parent wiring + a new standalone `carbon-bc-sa`
  regression harness (not code-only, not harness-only), mirroring the NI2G
  pass end to end.
- Attempt an actual build (`cmake --build`) and an actual `GEOS.x` run in
  this environment (compilers and a partial baseline-data tree are
  available here, unlike the NI2G porting session).

## Starting point

- `CA2G_GridCompMod.F90`: fully MAPL2 (`MAPL_GenericSetServices`,
  `ESMF_UserCompSetInternalState`, `MAPL_AddImportSpec`/
  `MAPL_AddInternalSpec`, `MAPL_GetPointer`, `.rc`-based config, a
  `pflogger` logger field, a `deactivate_instances__`-style dead code path).
- `CA2G_StateSpecs.rc`: MAPL2 column layout (`FRIENDLYTO`, literal `2`
  ungridded sizes, `MAPL_RestartOptional`).
- `GOCART2G_GridCompMod.F90` already anticipated CA in several places
  (`ACTIVE_INSTANCES: CA: [CA_oc, CA_bc, CA_br]` already present in
  `GOCART2G_GridComp.yaml`, `get_mixRatioSum`/`serialize_bundle`/
  `run_aerosol_optics`/`getAerosolSum` already branching on `"CA_oc"`/
  `"CA_bc"`/`"CA_br"`), but `create_instances_` was unconditionally calling a
  `deactivate_instances__(gc, self%CA, _RC)` stub that logged a warning and
  never added CA as a real child (CA2G wasn't ported, same pattern NI2G's
  notes describe for NI before that port).

## Reference files

- Template component (closest: single aerosol species, `_AERO`+`_AERO_DP`
  export pattern, phobic/philic-style internal split):
  `../SU2G_GridComp/SU2G_GridCompMod.F90`
- Template state spec: `../SU2G_GridComp/SU2G_StateSpecs.rc`
- Template instance config: `../SU2G_GridComp/SU2G_instance_SU.yaml`
- Parent: `../GOCART2G_GridCompMod.F90`, `../GOCART2G_StateSpecs.rc`
- Shared helpers: `../../Shared/Chem_AeroGeneric.F90` (`add_aero`,
  `append_to_bundle`, `get_mixR`, `findKlid`, `setZeroKlid`,
  `determine_data_driven`)
- Multi-instance parent-wiring precedent: `../NI2G_GridComp/
  mapl3-porting-notes.md` (step 7, "Parent wiring")

## CA's structural difference from DU/SS/SU/NI

CA2G is **one Fortran module shared by three independently-configured
GridComp instances** (`CA_oc`, `CA_bc`, `CA_br` -- organic carbon, black
carbon, brown carbon), each with its own instance YAML/`.rc` (own
`nbins=2`, `fscav`/`fwet`/particle properties, aerosol-optics LUT file,
etc.), unlike DU/SS/SU/NI which each have exactly one active computational
instance. This drives most of the differences from the SU2G template:

- **`internal_variable_name`**: set to the scalar `trim(comp_name)` (e.g.
  `"CA_bc"`), used only as an identifying token for
  `get_mixR`/`get_mixRatioSum`, not an actual field name -- same idiom NI2G
  introduced for its `"NI"` token.
- **`GCsuffix`**: a local `character(len=2)` derived by comparing
  `trim(comp_name)` against the literal strings `'CA_oc'`/`'CA_bc'`/
  `'CA_br'` (`if/else if` chain, three places: `Run1`, the AERO-optics
  routine, and `Run_data`), mapping to `'OC'`/`'BC'`/`'BR'` respectively.
  This drives which emission import fields (`OC_*`/`BC_*`/`BRC_*`), which
  process-library behavior (`GOCART2G_Process`'s settling/wet-removal calls
  take a `GCsuffix` argument), and which climatology import names
  (`climOCphobic`, `climBCphobic`, ... in the still-inactive data-driven
  path) apply to a given instance.
- **`determine_data_driven(comp_name, data_driven, _RC)`**: the same shared
  helper NI2G uses (`Chem_AeroGeneric.F90`), called identically in
  `SetServices`/`Initialize`/`Run` to decide the `_data`-suffix branch.
- **No replenishment alarm, no `TRADV`/`FRIENDLYTO` subscription, no
  `pflogger`**: unlike SU2G/NI2G, CA2G's own internal fields
  (`<comp_name>phobic`/`<comp_name>philic`) were never friendly to
  DYNAMICS/TURBULENCE/MOIST in the MAPL2 source, and CA has no oxidant
  replenishment mechanism, so none of SU2G's/NI2G's `MAPL_STATEITEM_SERVICE
  TRADV` or `createReplenishAlarm` machinery applies. The MAPL2 source's
  `use pflogger, only: logger_t => logger` import and a vestigial logger
  field were dropped (dead code, never used for actual logging calls).
- **`_AERO_DP` bundle IS exported** (unlike NI2G, which has none): CA2G
  keeps the `MAPL_STATEITEM_FIELDBUNDLE` deposition-bundle export, same as
  DU/SS/SU, built via `append_to_bundle` calls in `Initialize`.
- **Dead scaffolding removed**: the MAPL2 source's commented-out pointer
  declarations in `Run2` and the entire `deactivate_instances__`
  subroutine/call in the *parent* (see "Parent wiring" below) were CA-only
  leftovers from before this port and are now obsolete; removed rather than
  carried forward as comments.
- **Module-level parameters**: `n_ca_bins = 2`, `instanceComputational = 1`,
  `instanceData = 2`, `MAPL2_UNDEF = 1.0e15` -- same rationale as SU2G/NI2G
  (MAPL3's `MAPL_UNDEF` is `-huge`, which would break `1.01*x > undef`-style
  masks carried over verbatim from MAPL2).

## Steps

### 1. State specs and build

- `CA2G_StateSpecs.rc` rewritten column-by-column following SU2G: `nbins`
  (`=2`, per-instance from each `CA2G_instance_CA.{oc,bc,br}.yaml`) drives
  `ungrd_nbins = MAPL_UngriddedDim(self%nbins, name="nbins", ...)`, matching
  the `"nbins"` naming convention NI2G's porting notes flagged as required
  for cross-component ungridded-dim connection matching (not an issue here
  since CA has no incoming ungridded connections from other species, but
  kept consistent).
- Wavelength-indexed exports use `ungrd_wavelengths_profile`/
  `ungrd_wavelengths_vertint`, as in SU2G/NI2G.
- `CMakeLists.txt`: `TYPE SHARED`, `mapl_acg` call with
  `IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS`,
  matching the other species' CMakeLists.

### 2. Instance configuration

- Three instance YAMLs: `CA2G_instance_CA_oc.yaml`, `CA2G_instance_CA_bc.yaml`,
  `CA2G_instance_CA_br.yaml` (plus matching `.rc` files kept alongside, as
  DU/SS/SU/NI do), each with its own `aerosol_radBands_optics_file`/
  `aerosol_monochromatic_optics_file` pointing at
  `ExtData/chemistry/AerosolOptics/v2.x.x/x/opticsBands_{BC,OC,BR}...`/
  `optics_{BC,OC,BR}...`, own `fscav`/`fwet`/particle-property lists sized
  to `nbins: 2`, and `point_emissions_srcfilen: /dev/null` (disabling point
  emissions by default, matching SU2G's regression default).
- Each instance yaml's `mapl:` block carries `misc: activate_all_imports:
  true, activate_all_exports: true, restart: {import: true, internal:
  true}, checkpoint: {export/import/internal: true}` -- used as-is for the
  `carbon-bc-sa` regression copy except with `activate_all_imports` dropped
  and `restart: import` dropped, matching SU2G's regression-vs-production
  instance yaml diff (imports are wired explicitly via `root.yaml`
  connections in the regression harness instead).

### 3. SetServices

Rewritten after SU2G:

- `MAPL_GridCompGet(gc, name=comp_name)`; private state via
  `_SET_NAMED_PRIVATE_STATE`/`_GET_NAMED_PRIVATE_STATE` with
  `PRIVATE_STATE = "CA2G_GridComp"`.
- `ungrd_nbins`/`ungrd_wavelengths_profile`/`ungrd_wavelengths_vertint`
  declared before the ACG includes, as in SU2G/NI2G.
- CA-specific resources read via `MAPL_GridCompGetResource`:
  `my_day_of_week`, `hydrophobic_fraction`, `time_days_hydrophobic_to_hydrophilic`,
  `time_days_chemical_destruction` (asserted `size == self%nbins`),
  `pom_ca_ratio`, `monoterpenes_emission_fraction`,
  `isoprene_emission_fraction`, `aircraft_fuel_emission_factor`,
  `aviation_vertical_layers` (asserted `size == 4`),
  `point_emissions_srcfilen` (default `/dev/null`, used to toggle
  `self%doing_point_emissions`).
- Entry points: `Run -> "Run1"` always; `Run2 -> "Run2"` and
  `Run0 -> "Run0"` only when `.not. data_driven` (the data-driven branch
  has no Run0/Run2 phases, matching the still-`_FAIL`-guarded MAPL2
  behavior).
- Data-driven import/internal specs (`PLE`, `RH2`, `<comp_name>phobic`/
  `philic`, `clim<comp_name>phobic`/`philic`, per-bin `climCADP*`/
  `climCAWT*`/`climCASD*`/`climCASV*`) converted to `MAPL_GridCompAddSpec`
  calls behind `_FAIL("data driven section has not been activated yet")`,
  same inactive-by-construction pattern as SU2G/DU2G/NI2G's data-driven
  guards -- never exercised in MAPL2 either.
- Computational-instance specs come from the generated
  `CA2G_Import___.h`/`CA2G_Export___.h`/`CA2G_Internal___.h`.
- `<comp_name>_AERO` (`MAPL_STATEITEM_STATE`) **and** `<comp_name>_AERO_DP`
  (`MAPL_STATEITEM_FIELDBUNDLE`) export specs added unconditionally (both
  data-driven and computational instances), matching SU2G.

### 4. Initialize

Rewritten after SU2G's Initialize:

- `MAPL_GridCompGet(gc, name=, geom=, num_levels=)`; private/internal state
  via `_GET_NAMED_PRIVATE_STATE`/`MAPL_GridCompGetInternalState`.
- `determine_data_driven` branches `provider_state`/`prefix` (`import`/
  `'clim'` vs `export`/`''`), used to fetch phobic/philic fields from the
  right state.
- `<comp_name>_AERO`/`<comp_name>_AERO_DP` fetched via `ESMF_StateGet`;
  phobic/philic internal fields added to the AERO state with
  `ESMF_StateAdd`; deposition fields added to the `_AERO_DP` bundle via
  `append_to_bundle`, following SU2G/DU2G's pattern.
- `add_aero`/`add_aero_named_alias` calls for `effective_radius_in_microns`
  and other diagnostic aliases, as in SU2G.
- `internal_variable_name` set to the scalar `trim(comp_name)` (see
  "structural difference" above).
- Mie-table setup (`ESMF_InfoGetAlloc(info, key='mieTable_pointer', ...)`)
  per the SU2G/NI2G `ESMF_InfoGetAlloc`-for-allocatable-array-Info fix.

### 5. Run routines

Run0, Run (dispatch to `Run_data` or `Run1`), Run1, Run2, following SU2G:

- `MAPL_GridCompGet`/`MAPL_GridCompGetInternalState`,
  `_GET_NAMED_PRIVATE_STATE`; `rc` non-optional (except the internal
  `Run_data`).
- Fields come from the generated `CA2G_GetPointer___.h`/
  `CA2G_DeclarePointer___.h`.
- `GCsuffix` derived once per routine via the `trim(comp_name) ==
  'CA.{oc,bc,br}'` chain (see "structural difference" above), then passed
  into `GOCART2G_Process`'s `Chem_SettlingSimple`/`WetRemovalGOCART2G`/
  `WetRemovalUFS`/emissions calls exactly as the MAPL2 source did.
- Edge-level fields (`PLE`, `ZLE`, `PFL_LSAN`, `PFI_LSAN`) remapped to
  0-based lower bounds before any process-library call, matching the
  SU2G/DU2G/NI2G fix for the MAPL2-vs-MAPL3 edge-indexing convention
  change.

### 6. AERO methods

`aerosol_optics`/`monochromatic_aerosol_optics` rewritten after SU2G's
MAPL3 versions: `ESMF_InfoGet` on the state's Info (replacing
`ESMF_AttributeGet`), 0-based `ple0` remap, `ESMF_InfoGetAlloc` for the
`mieTable_pointer` integer array, `get_mixR` dispatch on
`internal_variable_name` (the existing `"DU"`/`"SS"` ungridded-sum /
else-branch single-field logic in `Chem_AeroGeneric::get_mixR` already
covers CA's case -- phobic+philic summed as a plain 2-bin ungridded field,
no new branch needed, unlike NI which needed a dedicated named-bins
branch).

### 7. Parent wiring (GOCART2G)

- `GOCART2G_GridCompMod.F90`:
  - `use CA2G_GridCompMod`; `add_children__(gc, self%CA, CA2G_SetServices,
    _RC)` enabled in `create_instances_`; the `deactivate_instances__`
    subroutine (and its call site) removed entirely, along with the now
    -unused `use pflogger, only: logger_t => logger` -- this was dead
    scaffolding left over from before CA2G was ported, not something NI2G
    or any other species still needs.
  - SetServices: `CA_AERO`/`CA_AERO_DP` import specs added per active CA
    instance, alongside the existing DU/SS/SU/NI loops (using `child%name`
    as the literal short-name prefix, since CA's three instances don't
    share one constant name the way DU/SS/SU/NI's single instance does).
  - SetServices: **per-instance dynamic connection strings** -- unlike
    DU/SS/SU/NI (one constant `child_items` string built once outside the
    loop), CA's `child_items` must be rebuilt *inside* the loop over
    `self%CA%instances`, using `child%name` (e.g. `"CA_bc"`) as the literal
    prefix for each of the 13 export fields (`EXTTAU`, `STEXTTAU`,
    `SCATAU`, `STSCATAU`, `EXTCOEF`, `EXTCOEFRH20`, `EXTCOEFRH80`,
    `SCACOEF`, `SCACOEFRH20`, `SCACOEFRH80`, `BCKCOEF`, `ANGSTR`, `SMASS`),
    because the three instances have different, non-interchangeable names.
    `_AERO`/`_AERO_DP` connections follow the same per-instance pattern.
  - Run2: the totals block (`totexttau = totexttau + <CA instance>exttau`
    etc.) now live, summing across whichever CA instances are active,
    following the DU/SS/SU/NI style.
- `GOCART2G_StateSpecs.rc`: **ALIAS column added** (last column of the
  single IMPORT header) plus 39 new CA import rows
  (`#-imports-from-carbon-` section, 13 fields x 3 instances:
  `CA_bcEXTTAU`/`CA_ocEXTTAU`/`CA_brEXTTAU`, etc.). The ALIAS column was
  introduced **because of CA**: the ACG (`MAPL_GridCompSpecs_ACG.py`)
  derives a field's Fortran pointer variable name from its `NAME`
  column, and the original dotted spelling (`CA.bcEXTTAU`) contained a
  literal `.`, which is not a legal Fortran identifier character. Adding
  `ALIAS` as the last column
  (e.g. `bcexttau` for `CA_bcEXTTAU`) makes the generated
  `GOCART2G_DeclarePointer___.h`/`GetPointer___.h` use the alias verbatim
  as the Fortran variable name while keeping `CA_bcEXTTAU` as the actual
  ESMF field/state name used for connections. Since the rename to
  underscores, `CA_bcEXTTAU` is itself a legal Fortran identifier, so the
  ALIAS column is **no longer strictly required**; it is retained only
  because the pointer names (`bcexttau` etc.) are what
  `GOCART2G_GridCompMod.F90` refers to. Dropping it is a safe follow-up
  that would also require renaming those pointer uses. Verified the existing
  DU/SS/SU/NI rows are unaffected (the ACG's row parser truncates/ignores
  trailing columns that are absent, so their un-aliased rows parse
  unchanged). A single IMPORT header for the whole file is required for
  this to work: `MAPL_GridCompSpecs_ACG.py` reads one column header per
  `category: IMPORT` block and multiple separate `category: IMPORT` blocks
  in the same file would not safely append/merge an extra column across
  blocks.
- `../CMakeLists.txt`: `CA2G_GridComp` enabled in `alldirs`.
- `GOCART2G_GridComp.yaml`: `ACTIVE_INSTANCES: CA: [CA_oc, CA_bc, CA_br]`
  was already present (pre-dating this port, left disabled only by the
  parent's `deactivate_instances__` stub); no yaml change needed in the
  production config.

### 8. `carbon-bc-sa` regression harness

Added a new standalone case, `regression/carbon-bc-sa/`, mirroring
`sulfate-sa`'s single-instance structure (one `DataForCA` +
`FakeGOCART(CA_bc)` pair) rather than NI2G's combined-case structure,
because -- unlike NI2G, which unconditionally reads `DU`/`SS`/`SO4` import
fields and their custom `radius`/`fnum` `Info` attributes at `Initialize`
-- a `grep` of `CA2G_GridCompMod.F90` confirmed **no cross-species
`ESMF_StateGet` calls on `import` at all**: CA2G only needs generic
met/emissions fields, the same standalone-friendly shape as SU2G. Picked
the `CA_bc` (black carbon) instance arbitrarily as the one to exercise,
since all three instances share one code path and only differ by
`GCsuffix`/instance-config values.

- `run{1,2}/data-for-ca.yaml`: full 58-field import set from
  `CA2G_StateSpecs.rc`'s IMPORT section (met/surface fields, emissions for
  OC/BC/BRC, aircraft/VOC-SOA fields), as a `configurable_gridcomp`
  `internal` state.
- `CA2G_instance_CA_bc.yaml` (top-level, regression copy of
  `CA2G_GridComp/CA2G_instance_CA_bc.yaml`): same resource values, `mapl:`
  block adapted to the SU2G regression pattern (`setServices: sharedObj:
  libCA2G_GridComp, userRoutine: ca2g_setservices_`;
  `activate_all_imports` dropped since imports are wired explicitly by
  `root.yaml`; `restart: internal: true` only).
- `run{1,2}/fake-gocart.yaml`: one child, `CA_bc:` (dso
  `libCA2G_GridComp.so`, `setservices: ca2g_setservices_`, `config_file:
  CA_instance_CA_bc.yaml`, `run_phases: [Run1]`/`[Run2]`).
- `run{1,2}/root.yaml`: `DataForCA` + `FakeGOCART` children, `im_world: 12`/
  `num_levels: 181` geometry (matching SU2G/NI2G), and an explicit
  `connections:` list for all 58 fields (`src_comp: DataForCA, src_intent:
  internal` -> `dst_comp: FakeGOCART`, same name both sides).
- `run{1,2}/cap_gridcomp.yaml`, `cap_driver.yaml`, `cap_restart.yaml`,
  `mapl.yaml`, `cmpchk.py`: near-verbatim copies of `sulfate-sa`'s,
  substituting `CA_bc` for `SU`.
- `run{1,2}/CA_instance_CA_bc.yaml`: symlinks to `../CA2G_instance_CA_bc.yaml`,
  matching the `SU_instance_SU.yaml -> ../SU2G_instance_SU.yaml` convention.
- Registered `carbon-bc-sa/run1`/`carbon-bc-sa/run2` in `regression/cases.txt`.
- Updated the stale `CA: []  # [CA_oc, CA_bc, CA_br] not ported to MAPL3
  yet` comment in both existing combined cases'
  (`gocart+seasalt+dust+sulfate`, `gocart+seasalt+dust+sulfate+nitrate`)
  `GOCART2G_GridComp.yaml` to reflect the real status (code ported and
  runnable). `CA` itself is still left as `[]` (disabled) in both, now
  purely to keep those two cases' scope and published baselines unchanged;
  CA coverage lives in `carbon-bc-sa` and in the
  `gocart+seasalt+dust+sulfate+nitrate+carbon-bc` case below.

### 9. `gocart+seasalt+dust+sulfate+nitrate+carbon-bc` combined-case harness

Added (in a follow-up pass) a new combined-case directory,
`regression/gocart+seasalt+dust+sulfate+nitrate+carbon-bc/`, extending
`gocart+seasalt+dust+sulfate+nitrate` by activating all three CA instances
(`CA_oc`, `CA_bc`, `CA_br`) through the real GOCART2G parent, matching the
production `ACTIVE_INSTANCES: CA: [CA_oc, CA_bc, CA_br]` shape -- unlike
`carbon-bc-sa`, which only exercises one instance standalone via `FakeGOCART`.

- `CA2G_instance_CA.{oc,bc,br}.yaml` (top-level): copies of the production
  `CA2G_GridComp/CA2G_instance_CA.{oc,bc,br}.yaml` with the same
  regression-vs-production trim DU/SS/SU/NI already use in this directory
  (drop the "TODO: remove this section..." comment, drop
  `activate_all_imports: true`, drop `restart: import: true`; keep
  `restart: internal: true`). No `setServices:` block is added (or needed),
  since these three instances are added via the parent's `add_children__`
  procedure call, the same as DU/SS/SU/NI here -- not via `dso`/DSO loading
  as in the standalone `carbon-bc-sa` harness.
- `data-for-gocart2g.yaml`: extended with the 34 CA-specific
  emissions/aircraft/VOC-SOA fields from `carbon-bc-sa`'s `data-for-ca.yaml`
  (`OC_*`/`BC_*`/`BRC_*`, `pSOA_ANTHRO_VOC`, `pSOA_BIOB_VOC`); all met/
  surface fields CA needs were already present (shared with DU/SS/SU/NI).
- `GOCART2G_GridComp.yaml`: `CA: [CA_oc, CA_bc, CA_br]` enabled (vs `[]` in
  the parent `+nitrate` directory), with a comment noting that the
  instance names use `_` rather than the retired dotted spelling.
- `run{1,2}/root.yaml`: the existing 92-field connection list (from
  `+nitrate`) extended with the 34 new CA fields, same
  `DataForGOCART2G -> GOCART2G` pattern.
- `run{1,2}/CA_instance_CA.{oc,bc,br}.yaml`: symlinks to
  `../CA2G_instance_CA.{oc,bc,br}.yaml`, following the naming convention
  `add_children__` expects: `config_file = species%name // "_instance_" //
  child_name // ".yaml"`, i.e. `species%name="CA"`, `child_name="CA_oc"` ->
  `CA_instance_CA_oc.yaml` (confirmed by reading
  `GOCART2G_GridCompMod.F90`'s `add_children__`).
- `run{1,2}/cmpchk.py`: `GRIDCOMPS` dict extended with `"CA_oc"`, `"CA_bc"`,
  `"CA_br"` entries (alongside `DU`/`SS`/`SU`/`NI`/`GOCART2G`).
- Registered `gocart+seasalt+dust+sulfate+nitrate+carbon-bc/run1`/`run2` in
  `regression/cases.txt`.
- **Confirmed via manual `GEOS.x mapl.yaml`** (same recipe as `carbon-bc-sa`,
  build/run environment: `module load ifx-stack`) that all three CA
  instances are added, SetServices'd and Initialize'd by the *real*
  `GOCART2G_GridCompMod.F90` `add_children__` call site, not just by the
  standalone `carbon-bc-sa` scaffolding:
  ```
  CA_oc: SetServices:: ...completed
  CA_bc: SetServices:: ...completed
  CA_br: SetServices:: ...completed
  CA_oc: Initialize:: ...completed
  CA_bc: Initialize:: ...completed
  CA_br: Initialize:: ...completed
  ```
  The case still cannot be driven to completion by `ctest -L REGRESSION`
  because no baseline/restart data has been published for it yet (it
  reports "Regression data not found ... skipping", same as every case
  before its baselines exist); the run above was done against a scratch
  expdir with the restart reads disabled purely to exercise
  SetServices/Initialize.

## Run-testing update (manual `GEOS.x`, not yet `ctest -L REGRESSION`)

Build/run environment: `module load ifx-stack` (ifx/2025.3,
intelmpi/2021.17, Baselibs/8.33.0 -- the stack the build tree was
configured against).

Built `CA2G_GridComp`/`GOCART2G_GridComp` and relinked `GEOS.x`: **all
build and link cleanly, zero errors/warnings**. Inspected the generated
`GOCART2G_DeclarePointer___.h`/`GetPointer___.h`: the ALIAS-column pointers
(`bcexttau`, `ocexttau`, `brexttau`, etc.) are generated correctly and map
to the ESMF field names as intended.

**`carbon-bc-sa/run1` and `carbon-bc-sa/run2` both run to completion**
(`mpiexec -n 6 GEOS.x mapl.yaml`, scratch expdir, `ExtData` symlinked from
the shared regression-data tree, restarts from the partial baseline tree):
exit status 0, zero `FAIL` lines, and the expected checkpoints are written:

```
checkpoints/last/CA_bc_export.nc
checkpoints/last/CA_bc_import.nc
checkpoints/last/CA_bc_internal.nc
```

### Instance naming: `CA_oc`/`CA_bc`/`CA_br` (formerly `CA.oc`/`CA.bc`/`CA.br`)

The carbon instances were originally named with a period separator
(`CA.oc`, `CA.bc`, `CA.br`), inherited from the MAPL2-era GOCART2G design.
An earlier pass of this port recorded a hard blocker: MAPL3's
`MAPL_GridCompAddChild` (in
`src/Shared/@MAPL/superstructure/generic/MAPL_Generic.F90`) asserts
`is_valid_name(child_name)` on every child, and `is_valid_name`
(`src/Shared/@MAPL/utils/Validation.F90`) rejected `.`, so the dotted
names failed with *"Child name <CA.bc> does not conform to GEOS
standards."* before `Initialize` was ever reached. That was worked around
by relaxing `Validation.F90` to admit `.`:

```fortran
is_valid = (verify(name(2:), LOWER // UPPER // DIGITS // '_'//'.') == 0)
```

**The instances have since been renamed to use `_`**, for three reasons:

- MAPL3 History uses `.` as the `<gridcomp>.<variable>` dive-down
  delimiter (`gridcomps/history/HistoryUtilities.F90`,
  `replace_delimiter`, which substitutes only the *first* `.`).
  `HistoryGridComp_private.F90` states the contract explicitly: *"Collection
  names are permitted to include period ('.') characters, but gridcomps are
  not."* A dotted child name contradicts that, and composed names such as
  `CA_oc_data` are unparseable under a first-dot-wins grammar.
- The `.` carried no information: nothing in GOCART ever split on it. The
  names were flat, opaque tokens everywhere they appeared.
- `CA_bcEXTTAU` is a legal Fortran identifier, so the dot no longer forces
  the `GOCART2G_StateSpecs.rc` ALIAS column (kept for now only so the
  existing pointer names in `GOCART2G_GridCompMod.F90` are undisturbed).

The rename touched `GOCART2G_GridCompMod.F90`, `CA2G_GridCompMod.F90`'s
`GCsuffix` comparisons, `Chem_AeroGeneric.F90`,
`GOCART2G_AerosolActivationMod.F90`, `g2g_rst_generator.py`, every instance
yaml/`.rc` filename (and the symlinks in `AMIP/`, `AMIP.20C/` and the
regression `run{1,2}/` dirs), the `GOCART2G_StateSpecs.rc` CA rows, and the
`carbon-bc-sa` / `...+carbon-bc` harnesses.

Note `Chem_AeroGeneric.F90`'s `varname(6:)` offset is unaffected: `CA_bc`
is the same 5 characters wide as `CA.bc` was.

Remaining consequences:

- `Validation.F90`'s `PUNCTUATION = '_'//'.'` relaxation is no longer
  needed by GOCART and could be reverted to `'_'` (its doc comment still
  claims only letters/digits/underscores are allowed).
- **Regression baselines must be regenerated or renamed.** Checkpoint and
  restart filenames derive from the child name, so
  `CA.bc_{import,internal,export}.nc` under `${REGRESSION_DATA_DIR}`
  becomes `CA_bc_*.nc`. Likewise, existing GEOS internal restarts carry
  variables named `CA.bcphobic`/`CA.bcphilic`, which become
  `CA_bcphobic`/`CA_bcphilic`.
- The commented-out data-driven instance variants (`CA_oc_data` etc., and
  `DU_data`/`SS_data`/`SU_data`/`NI_data`) were renamed to underscores
  alongside the rest; they remain commented out and unwired.
- The aerosol-optics ExtData previously feared missing
  (`opticsBands_{BC,OC,BR}.v2.1.0.RRTMGP.nc4` /
  `optics_{BC,OC,BR}.v2.1.0.nc4`) resolves fine from the regression ExtData
  tree; all three instances log their optics files and complete
  `Initialize`.

## Baseline comparison (`cmpchk.py`) and the bug it exposed

With `carbon-bc-sa/run{1,2}` running, `cmpchk.py` (plain `baseline - current`,
no NaN masking) was used against the MAPL2-era baselines. Three distinct
classes of difference showed up; only one was a real bug.

### 1. REAL BUG, now fixed: `diurnal_bb` was never initialized

`run1` (Run1 phase = emissions) showed non-zero diffs in `CA_bcphobic` /
`CA_bcphilic`:

```
CA_bcphilic   diff (min, max, norm): -1.13e-11, 2.17e-11, 6.35e-11
CA_bcphobic   diff (min, max, norm): -4.53e-11, 8.68e-11, 2.54e-10
```

These are small in absolute terms only because the whole field is ~1e-10
kg/kg; locally the discrepancy reached **3200%**. The fingerprint pinned
down the cause precisely:

- only 2051 of 156384 points differed, in just **129 of 864 columns**
  (sparse -> a biomass-burning source, not a global one);
- diffs were confined to levels 175-180 (the near-surface injection
  layers) and were a *constant* increment per level;
- `diff(philic) / diff(phobic)` was exactly **0.25** = `(1 -
  fHydrophobic) / fHydrophobic` with `fHydrophobic = 0.8`, i.e. the
  missing mass went through the emission hydrophobic/hydrophilic split;
- imports were bitwise identical, so it was not an input/ExtData issue.

Root cause: `Run1` guards the biomass-burning diurnal cycle with
`if (self%diurnal_bb)`, but **`self%diurnal_bb` was never assigned
anywhere in the port** -- the MAPL2 original reads it in `Initialize`

```fortran
call ESMF_ConfigGetAttribute(universal_cfg, diurnal_bb, &
     label='DIURNAL_BIOMASS_BURNING:', default='YES', __RC__)
```

and that block was dropped during porting. The field was therefore an
uninitialized logical, which in practice evaluated `.false.`, silently
**disabling diurnal biomass burning** (MAPL2 default is *on*).
`SU2G_GridCompMod.F90` had been ported correctly and still had the
equivalent block, which is what the fix was modelled on:

```fortran
! Check whether to de-activate diurnal biomass burning (default is *on*)
call MAPL_GridCompGetResource(gc, "DIURNAL_BIOMASS_BURNING", diurnal_bb, default='YES', _RC)
diurnal_bb = ESMF_UtilStringUpperCase(diurnal_bb, _RC)
self%diurnal_bb = (trim(diurnal_bb) == 'YES')
```

**After the fix, `carbon-bc-sa/run1` and `run2` reproduce the baseline
bitwise** (`diff (min, max, norm): 0.0, 0.0, 0.0`) for `CA_bcphobic`,
`CA_bcphilic`, and for `CA_bcDP`/`CA_bcSD`/`CA_bcWT` in `run2`.

This class of bug (a resource read silently dropped, leaving an
uninitialized `logical` that defaults to the *non-default* behaviour) is
worth checking for in the other GOCART2G species ports.

### 2. Expected: phase-1 checkpoint has NaN for Run2-only exports

In `run1` the harness runs `run_phases: [Run1]` only, while `DP`
(dry deposition), `SD` (sedimentation) and `WT` (wet deposition) are
computed exclusively in `Run2`. MAPL3 initializes export fields to NaN
("undefined"), so they are correctly still NaN in the `run1` checkpoint.
The MAPL2 baseline shows values there only because MAPL2 zero-filled /
carried over export memory. Not a defect -- MAPL3 is simply making
"never written in this phase" visible. `cmpchk.py` has no NaN masking, so
`NaN - value = NaN` propagates to all three printed statistics.

### 3. Pre-existing: `SV` is declared but never computed

`CA_bcSV` ("Aerosol Convective Scavenging") is NaN in *both* runs. `SV` is
declared in `CA2G_StateSpecs.rc`, has a pointer generated into
`CA2G_DeclarePointer___.h`, is fetched in `CA2G_GetPointer___.h` and is
appended to the `_AERO_DP` bundle -- but **it is never assigned a value**,
in the MAPL3 port *or* in the MAPL2 original (verified against
`v11.11.0`). So this is an inherited dead export, not a porting
regression. It shows as NaN now only because MAPL3 no longer zero-fills
unwritten exports. Either CA2G should fill it or the spec should drop it.

Note that `ctest`'s own comparison uses `NANS_ARE_EQUAL`, which tolerates
NaN on both sides; `cmpchk.py` does not, which is why these show up here
and not in the `ctest` runs.

## Open items (follow-ups, not blocking the code port)

1. **Publish regression baselines.** `carbon-bc-sa/run{1,2}` and
   `gocart+seasalt+dust+sulfate+nitrate+carbon-bc/run{1,2}` are registered in
   `regression/cases.txt` and verified by hand, but no baseline/restart
   data exists for them under
   `/ford1/share/gmao_SIteam/regression-data/.../GOCART2G_GridComp/`, so
   `ctest -L REGRESSION` reports "skipping" for all four. Publishing the
   baselines is what turns the manual verification above into real CI
   coverage. (The `carbon-bc-sa` runs above used a partial baseline tree that
   is present locally but not published.)
2. **CA left disabled in the two older combined cases.**
   `gocart+seasalt+dust+sulfate` and `gocart+seasalt+dust+sulfate+nitrate`
   still carry `CA: []`. This is now a deliberate scope/baseline-stability
   choice, not a technical blocker: enabling CA there would invalidate
   their published baselines. CA coverage lives in `carbon-bc-sa` and
   `gocart+seasalt+dust+sulfate+nitrate+carbon-bc`.
3. **Dead `SV` export**: `CA_bcSV` etc. are declared and bundled but never
   assigned, in MAPL3 *and* MAPL2. Either compute it or remove it from
   `CA2G_StateSpecs.rc`. See "Baseline comparison" above.
4. **Audit other species for dropped resource reads**: the `diurnal_bb`
   bug (see above) was a silently-omitted `Initialize` resource read
   leaving an uninitialized `logical`. DU/SS/NI should be checked for the
   same class of omission.
5. **GMI/point-emissions paths untested**: `self%doing_point_emissions`
   defaults to `false` (`/dev/null`) in all three regression/production
   instance configs; the point-emissions code path itself is unexercised,
   same caveat class as SU2G's/NI2G's untested branches.
6. **Data-driven (`_data`) instances**: still guarded by
   `_FAIL("data driven section has not been activated yet")`, unchanged
   from the MAPL2 source; `CA_oc_data`/`CA_bc_data` are not wired into
   `ACTIVE_INSTANCES`/`PASSIVE_INSTANCES` anywhere, consistent with
   DU/SS/SU's equivalents.
