# CA2G_GridComp MAPL3 porting plan

Status: **code-complete, built, run-tested, blocked at Initialize** by a
MAPL3 framework constraint (see "Run-testing update" below). Ported on top
of `release/MAPL-v3`, following the SU2G port
(`../SU2G_GridComp/mapl3-porting-notes.md`) as the primary template for the
single-state-bundle `_AERO`/`_AERO_DP` export pattern, with NI2G
(`../NI2G_GridComp/mapl3-porting-notes.md`) as the most recent precedent for
multi-instance parent wiring and the "GOCART2G already anticipated this
species" starting point.

Scope of this pass, agreed up front:
- Full port: code + GOCART2G parent wiring + a new standalone `carbon-sa`
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
  (`ACTIVE_INSTANCES: CA: [CA.oc, CA.bc, CA.br]` already present in
  `GOCART2G_GridComp.yaml`, `get_mixRatioSum`/`serialize_bundle`/
  `run_aerosol_optics`/`getAerosolSum` already branching on `"CA.oc"`/
  `"CA.bc"`/`"CA.br"`), but `create_instances_` was unconditionally calling a
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
GridComp instances** (`CA.oc`, `CA.bc`, `CA.br` -- organic carbon, black
carbon, brown carbon), each with its own instance YAML/`.rc` (own
`nbins=2`, `fscav`/`fwet`/particle properties, aerosol-optics LUT file,
etc.), unlike DU/SS/SU/NI which each have exactly one active computational
instance. This drives most of the differences from the SU2G template:

- **`internal_variable_name`**: set to the scalar `trim(comp_name)` (e.g.
  `"CA.bc"`), used only as an identifying token for
  `get_mixR`/`get_mixRatioSum`, not an actual field name -- same idiom NI2G
  introduced for its `"NI"` token.
- **`GCsuffix`**: a local `character(len=2)` derived by comparing
  `trim(comp_name)` against the literal strings `'CA.oc'`/`'CA.bc'`/
  `'CA.br'` (`if/else if` chain, three places: `Run1`, the AERO-optics
  routine, and `Run_data`), mapping to `'OC'`/`'BC'`/`'BR'` respectively.
  This drives which emission import fields (`OC_*`/`BC_*`/`BRC_*`), which
  process-library behavior (`GOCART2G_Process`'s settling/wet-removal calls
  take a `GCsuffix` argument), and which climatology import names
  (`climOCphobic`, `climBCphobic`, ... in the still-inactive data-driven
  path) apply to a given instance.
- **`determine_data_driven(comp_name, data_driven, _RC)`**: the same shared
  helper NI2G uses (`Chem_AeroGeneric.F90`), called identically in
  `SetServices`/`Initialize`/`Run` to decide the `.data`-suffix branch.
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

- Three instance YAMLs: `CA2G_instance_CA.oc.yaml`, `CA2G_instance_CA.bc.yaml`,
  `CA2G_instance_CA.br.yaml` (plus matching `.rc` files kept alongside, as
  DU/SS/SU/NI do), each with its own `aerosol_radBands_optics_file`/
  `aerosol_monochromatic_optics_file` pointing at
  `ExtData/chemistry/AerosolOptics/v2.x.x/x/opticsBands_{BC,OC,BR}...`/
  `optics_{BC,OC,BR}...`, own `fscav`/`fwet`/particle-property lists sized
  to `nbins: 2`, and `point_emissions_srcfilen: /dev/null` (disabling point
  emissions by default, matching SU2G's regression default).
- Each instance yaml's `mapl:` block carries `misc: activate_all_imports:
  true, activate_all_exports: true, restart: {import: true, internal:
  true}, checkpoint: {export/import/internal: true}` -- used as-is for the
  `carbon-sa` regression copy except with `activate_all_imports` dropped
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
    `self%CA%instances`, using `child%name` (e.g. `"CA.bc"`) as the literal
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
  `CA.bcEXTTAU`/`CA.ocEXTTAU`/`CA.brEXTTAU`, etc.). The ALIAS column exists
  **specifically because of CA**: the ACG (`MAPL_GridCompSpecs_ACG.py`)
  normally derives a field's Fortran pointer variable name from its `NAME`
  column, but `CA.bcEXTTAU` etc. contain a literal `.`, which is not a
  legal Fortran identifier character. Adding `ALIAS` as the last column
  (e.g. `bcexttau` for `CA.bcEXTTAU`) makes the generated
  `GOCART2G_DeclarePointer___.h`/`GetPointer___.h` use the alias verbatim
  as the Fortran variable name while keeping `CA.bcEXTTAU` as the actual
  ESMF field/state name used for connections. Verified the existing
  DU/SS/SU/NI rows are unaffected (the ACG's row parser truncates/ignores
  trailing columns that are absent, so their un-aliased rows parse
  unchanged). A single IMPORT header for the whole file is required for
  this to work: `MAPL_GridCompSpecs_ACG.py` reads one column header per
  `category: IMPORT` block and multiple separate `category: IMPORT` blocks
  in the same file would not safely append/merge an extra column across
  blocks.
- `../CMakeLists.txt`: `CA2G_GridComp` enabled in `alldirs`.
- `GOCART2G_GridComp.yaml`: `ACTIVE_INSTANCES: CA: [CA.oc, CA.bc, CA.br]`
  was already present (pre-dating this port, left disabled only by the
  parent's `deactivate_instances__` stub); no yaml change needed in the
  production config.

### 8. `carbon-sa` regression harness

Added a new standalone case, `regression/carbon-sa/`, mirroring
`sulfate-sa`'s single-instance structure (one `DataForCA` +
`FakeGOCART(CA.bc)` pair) rather than NI2G's combined-case structure,
because -- unlike NI2G, which unconditionally reads `DU`/`SS`/`SO4` import
fields and their custom `radius`/`fnum` `Info` attributes at `Initialize`
-- a `grep` of `CA2G_GridCompMod.F90` confirmed **no cross-species
`ESMF_StateGet` calls on `import` at all**: CA2G only needs generic
met/emissions fields, the same standalone-friendly shape as SU2G. Picked
the `CA.bc` (black carbon) instance arbitrarily as the one to exercise,
since all three instances share one code path and only differ by
`GCsuffix`/instance-config values.

- `run{1,2}/data-for-ca.yaml`: full 58-field import set from
  `CA2G_StateSpecs.rc`'s IMPORT section (met/surface fields, emissions for
  OC/BC/BRC, aircraft/VOC-SOA fields), as a `configurable_gridcomp`
  `internal` state.
- `CA2G_instance_CA.bc.yaml` (top-level, regression copy of
  `CA2G_GridComp/CA2G_instance_CA.bc.yaml`): same resource values, `mapl:`
  block adapted to the SU2G regression pattern (`setServices: sharedObj:
  libCA2G_GridComp, userRoutine: ca2g_setservices_`;
  `activate_all_imports` dropped since imports are wired explicitly by
  `root.yaml`; `restart: internal: true` only).
- `run{1,2}/fake-gocart.yaml`: one child, `CA.bc:` (dso
  `libCA2G_GridComp.so`, `setservices: ca2g_setservices_`, `config_file:
  CA_instance_CA.bc.yaml`, `run_phases: [Run1]`/`[Run2]`).
- `run{1,2}/root.yaml`: `DataForCA` + `FakeGOCART` children, `im_world: 12`/
  `num_levels: 181` geometry (matching SU2G/NI2G), and an explicit
  `connections:` list for all 58 fields (`src_comp: DataForCA, src_intent:
  internal` -> `dst_comp: FakeGOCART`, same name both sides).
- `run{1,2}/cap_gridcomp.yaml`, `cap_driver.yaml`, `cap_restart.yaml`,
  `mapl.yaml`, `cmpchk.py`: near-verbatim copies of `sulfate-sa`'s,
  substituting `CA.bc` for `SU`.
- `run{1,2}/CA_instance_CA.bc.yaml`: symlinks to `../CA2G_instance_CA.bc.yaml`,
  matching the `SU_instance_SU.yaml -> ../SU2G_instance_SU.yaml` convention.
- Registered `carbon-sa/run1`/`carbon-sa/run2` in `regression/cases.txt`.
- Updated the stale `CA: []  # [CA.oc, CA.bc, CA.br] not ported to MAPL3
  yet` comment in both existing combined cases'
  (`gocart+seasalt+dust+sulfate`, `gocart+seasalt+dust+sulfate+nitrate`)
  `GOCART2G_GridComp.yaml` to reflect the real, more specific status (code
  ported, blocked by the `is_valid_name` issue below) -- `CA` itself is
  still left as `[]` (disabled) in both, since enabling it would hit the
  same blocker described below.

### 9. `gocart+seasalt+dust+sulfate+nitrate+carbon` combined-case harness

Added (in a follow-up pass, after the blocker in "Run-testing update" below
was already known and documented) a new combined-case directory,
`regression/gocart+seasalt+dust+sulfate+nitrate+carbon/`, extending
`gocart+seasalt+dust+sulfate+nitrate` by activating all three CA instances
(`CA.oc`, `CA.bc`, `CA.br`) through the real GOCART2G parent, matching the
production `ACTIVE_INSTANCES: CA: [CA.oc, CA.bc, CA.br]` shape -- unlike
`carbon-sa`, which only exercises one instance standalone via `FakeGOCART`.

- `CA2G_instance_CA.{oc,bc,br}.yaml` (top-level): copies of the production
  `CA2G_GridComp/CA2G_instance_CA.{oc,bc,br}.yaml` with the same
  regression-vs-production trim DU/SS/SU/NI already use in this directory
  (drop the "TODO: remove this section..." comment, drop
  `activate_all_imports: true`, drop `restart: import: true`; keep
  `restart: internal: true`). No `setServices:` block is added (or needed),
  since these three instances are added via the parent's `add_children__`
  procedure call, the same as DU/SS/SU/NI here -- not via `dso`/DSO loading
  as in the standalone `carbon-sa` harness.
- `data-for-gocart2g.yaml`: extended with the 34 CA-specific
  emissions/aircraft/VOC-SOA fields from `carbon-sa`'s `data-for-ca.yaml`
  (`OC_*`/`BC_*`/`BRC_*`, `pSOA_ANTHRO_VOC`, `pSOA_BIOB_VOC`); all met/
  surface fields CA needs were already present (shared with DU/SS/SU/NI).
- `GOCART2G_GridComp.yaml`: `CA: [CA.oc, CA.bc, CA.br]` enabled (vs `[]` in
  the parent `+nitrate` directory), with a comment documenting the expected
  `is_valid_name` failure inline (mirroring `nitrate-sa`'s
  "KNOWN LIMITATION" comment style).
- `run{1,2}/root.yaml`: the existing 92-field connection list (from
  `+nitrate`) extended with the 34 new CA fields, same
  `DataForGOCART2G -> GOCART2G` pattern.
- `run{1,2}/CA_instance_CA.{oc,bc,br}.yaml`: symlinks to
  `../CA2G_instance_CA.{oc,bc,br}.yaml`, following the naming convention
  `add_children__` expects: `config_file = species%name // "_instance_" //
  child_name // ".yaml"`, i.e. `species%name="CA"`, `child_name="CA.oc"` ->
  `CA_instance_CA.oc.yaml` (confirmed by reading
  `GOCART2G_GridCompMod.F90`'s `add_children__`).
- `run{1,2}/cmpchk.py`: `GRIDCOMPS` dict extended with `"CA.oc"`, `"CA.bc"`,
  `"CA.br"` entries (alongside `DU`/`SS`/`SU`/`NI`/`GOCART2G`).
- Registered `gocart+seasalt+dust+sulfate+nitrate+carbon/run1`/`run2` in
  `regression/cases.txt`.
- **Confirmed via manual `GEOS.x mapl.yaml`** (same recipe as `carbon-sa`)
  that this case hits the identical, already-documented `is_valid_name`
  blocker, but this time inside the *real* `GOCART2G_GridCompMod.F90`
  itself rather than only the standalone test scaffolding:
  ```
  FAIL at line=00442    MAPL_Generic.F90            <Child name <CA.oc> does not conform to GEOS standards.>
  FAIL at line=01269    GOCART2G_GridCompMod.F90    <status=1>
  FAIL at line=01250    GOCART2G_GridCompMod.F90    <status=1>
  ```
  This is the first time the production parent's own `add_children__` call
  site (not just `carbon-sa`'s copy of the same MAPL3 API) has been
  exercised end to end with CA active, and it fails exactly as predicted.
  `ctest -L REGRESSION` still passes (29/29) with this case included, since
  it reports "skipping" (no baseline data) rather than actually invoking
  `GEOS.x`.

## Run-testing update (manual `GEOS.x`, not yet `ctest -L REGRESSION`)

Built `CA2G_GridComp` and `GOCART2G_GridComp` (`cmake --build . --target
GOCART2G_GridComp -j4`): **both build cleanly, zero errors/warnings**.
Inspected the generated `GOCART2G_DeclarePointer___.h`/`GetPointer___.h`:
the ALIAS-column pointers (`bcexttau`, `ocexttau`, `brexttau`, etc.) are
generated correctly and map to the dotted ESMF field names as intended.

Reconfigured (`cmake .`, without `GEOSpyD` loaded -- loading it breaks MPI
detection in this environment) to pick up the new `carbon-sa/run{1,2}`
cases in `cases.txt`; confirmed via `ctest -N` and `ctest -R carbon-sa -V`
that the harness is wired correctly and reports "Regression data not found
... skipping" exactly like `nitrate-sa` and every other case did before its
baseline data existed -- i.e. no YAML/registration errors.

**Manually ran `GEOS.x mapl.yaml` against `carbon-sa/run1` and `run2`**
(`mpiexec -n 6`, scratch expdir, `ExtData` symlinked from the shared
regression-data tree, no restart data) to check CA2G's `Initialize`
specifically, since no `ctest`-level baseline exists yet. **Found a real,
previously-undiscovered blocker, not a missing-data issue:**

- Both runs fail identically, at `SetServices` (before `Initialize` is even
  reached):
  ```
  FAIL at line=00035 ...uterMetaComponent_SMOD_add_child_by_spec
  <Child name <CA.bc> does not conform to GEOS standards.>
  ```
- Root cause: MAPL3's `MAPL_GridCompAddChild` (all overloads, in
  `src/Shared/@MAPL/superstructure/generic/MAPL_Generic.F90`) asserts
  `is_valid_name(child_name)` before adding any child. `is_valid_name`
  (`src/Shared/@MAPL/utils/Validation.F90`) requires names to start with a
  letter and contain only letters/digits/underscores -- **no `.`
  allowed**. This is a hard, repo-wide MAPL3 constraint (added in a
  "Rename mapl3g_ module namespace to mapl_" / "enforce gridcomps use MAPL
  umbrella" pass per `git log -- utils/Validation.F90`), not something
  specific to the `carbon-sa` harness.
- **This is not limited to `carbon-sa`.** `GOCART2G_GridCompMod.F90`'s own
  `add_children__` calls the exact same `MAPL_GridCompAddChild(gc,
  child_name, setservices, hconfig, _RC)` with `child_name =
  "CA.oc"`/`"CA.bc"`/`"CA.br"` (from `ACTIVE_INSTANCES`). If the production
  `GOCART2G_GridComp.yaml`'s `CA: [CA.oc, CA.bc, CA.br]` were enabled and
  `GEOS.x` actually run, it would hit this identical assertion. This was
  never caught before this session because all prior verification of the
  CA2G parent wiring was compile-only (`cmake --build`), never an actual
  `GEOS.x` execution with CA wired in -- the dotted-instance-name
  convention is inherited unmodified from the MAPL2-era GOCART2G design,
  which did not have this validation.
- Per explicit direction, **not fixing this in this pass** (renaming the
  instance-naming convention from dotted (`CA.oc`/`CA.bc`/`CA.br`) to
  underscore-based (`CA_oc`/`CA_bc`/`CA_br`) would touch
  `GOCART2G_GridCompMod.F90`, `CA2G_GridCompMod.F90`'s `GCsuffix`
  comparisons, every instance yaml/`.rc` filename, the
  `GOCART2G_StateSpecs.rc` ALIAS rows (which would likely become
  unnecessary, since e.g. `CA_bcEXTTAU` is already a legal Fortran
  identifier without an alias), and the `carbon-sa` harness -- a
  nontrivial, cross-cutting rename, left for a follow-up). **Documenting as
  a known blocker instead**, same spirit as NI2G's missing aerosol-optics
  ExtData finding: the port is code-complete and builds, but cannot
  actually run end to end yet.
- `ctest -L REGRESSION` passes cleanly today with `carbon-sa/run{1,2}`
  included (both report "skipping", same as before baseline data exists);
  this is orthogonal to the `is_valid_name` blocker above, which only
  surfaces when `GEOS.x` is actually invoked (manually, since no baseline
  exists yet to drive a real `ctest` run).

## Open items (follow-ups, not blocking the code port)

1. **Blocking for actual execution: dotted CA instance names
   (`CA.oc`/`CA.bc`/`CA.br`) violate MAPL3's `is_valid_name` child-name
   check.** See "Run-testing update" above. Until this is resolved (most
   likely by renaming the convention to use underscores, and simplifying
   away the `GOCART2G_StateSpecs.rc` ALIAS-column workaround at the same
   time), neither `carbon-sa` nor the production GOCART2G parent can
   actually initialize with CA active.
2. **Missing ExtData, not yet confirmed**: even once item 1 is fixed, the
   `opticsBands_{BC,OC,BR}.v2.1.0.RRTMGP.nc4`/`optics_{BC,OC,BR}.v2.1.0.nc4`
   files referenced by each instance yaml are not present anywhere on this
   filesystem (checked via `find / -xdev`), matching the NI2G precedent for
   newly-ported species' aerosol-optics data. Expect `Initialize` to get
   further than `SetServices` once item 1 is fixed, then fail on these
   missing files, the same way NI2G's combined-case test did for
   `optics_NI`.
3. **No `carbon-sa` baseline data published** under
   `/ford1/share/gmao_SIteam/regression-data/.../GOCART2G_GridComp/`; until
   it is (and item 1 is fixed), `ctest -L REGRESSION` will keep reporting
   "skipping" for `carbon-sa/run{1,2}`.
4. **GMI/point-emissions paths untested**: `self%doing_point_emissions`
   defaults to `false` (`/dev/null`) in all three regression/production
   instance configs; the point-emissions code path itself is unexercised,
   same caveat class as SU2G's/NI2G's untested branches.
5. **Data-driven (`.data`) instances**: still guarded by
   `_FAIL("data driven section has not been activated yet")`, unchanged
   from the MAPL2 source; `CA.oc.data`/`CA.bc.data` are not wired into
   `ACTIVE_INSTANCES`/`PASSIVE_INSTANCES` anywhere, consistent with
   DU/SS/SU's equivalents.
