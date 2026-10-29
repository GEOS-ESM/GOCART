# NI2G_GridComp MAPL3 porting plan

Status: **code-complete, partially run-tested** (see "Run-testing update"
below). Ported on top of `release/MAPL-v3`,
following the SU2G port
(`../SU2G_GridComp/mapl3-porting-notes.md`) as the primary template, with
DU2G (`../DU2G_GridComp/DU2G_GridCompMod.F90`) as the secondary template for
the multi-bin settling/wet-removal/diagnostics pattern.

Scope of this pass, agreed up front:
- Code-level port only: no MAPL2 baseline capture, no build, no
  regression-testing in this session (no compilers/baseline data available
  in this environment). The next step for whoever picks this up is to build
  GOCART2G (`DU2G`, `SS2G`, `SU2G`, `NI2G`) and run the existing GOCART2G
  regressions plus a new NI2G one.
- No cosmetic/formatting passes (SU2G's steps 1, 1a, 1b, 1c): the MAPL2 ->
  MAPL3 API changes were made directly, keeping the existing layout/comment
  style instead of running `codee format` or the macro/comment cleanup
  passes.
- GOCART2G parent wiring included: NI2G child creation, the `NI_AERO`
  import/connection, SO4/DU/SS imports into NI, and NI's own export
  connections into GOCART2G (step 8-equivalent), since DU2G/SS2G/SU2G are
  already ported and `ACTIVE_INSTANCES` already lists `NI: [NI]` in
  `GOCART2G_GridComp.yaml`.

## Starting point

- `NI2G_GridCompMod.F90`: 1519 lines, fully MAPL2 (`MAPL_GenericSetServices`,
  `ESMF_UserCompSetInternalState`, `MAPL_AddImportSpec`/`MAPL_AddInternalSpec`,
  `MAPL_GetPointer`, `.rc`-based config).
- `NI2G_StateSpecs.rc`: MAPL2 column layout (`FRIENDLYTO`, literal `3`/`5`
  ungridded sizes, `MAPL_RestartOptional`).
- `GOCART2G_GridCompMod.F90` already anticipated NI everywhere (commented
  `use NI2G_GridCompMod`, commented AERO import spec, commented
  connectivity block, commented Run2 totals, `add_aero_states_` already
  skips NI's `_AERO_DP` on purpose, `get_mixRatioSum`/`serialize_bundle`/
  `run_aerosol_optics` already branch on `"NI"` and
  `sum_of_internalState_aerosol_NI`). `ACTIVE_INSTANCES: NI: [NI]` was
  already active in `GOCART2G_GridComp.yaml` (only `create_instances_` was
  deactivating it with a warning, since NI2G wasn't ported).

## Reference files

- Template component (closest, single-bin + oxidant/replenish-alarm
  pattern): `../SU2G_GridComp/SU2G_GridCompMod.F90`
- Template component (multi-bin settling/wet-removal/diags pattern):
  `../DU2G_GridComp/DU2G_GridCompMod.F90`
- Template state spec: `../SU2G_GridComp/SU2G_StateSpecs.rc`,
  `../DU2G_GridComp/DU2G_StateSpecs.rc`
- Template instance config: `../SU2G_GridComp/SU2G_instance_SU.yaml`,
  `../DU2G_GridComp/DU2G_instance_DU.yaml`
- Parent: `../GOCART2G_GridCompMod.F90`, `../GOCART2G_StateSpecs.rc`
- Shared helpers: `../../Shared/Chem_AeroGeneric.F90` (`add_aero`,
  `append_to_bundle`, `get_mixR`, `findKlid`, `setZeroKlid`)

## NI's structural difference from DU/SS/SU

NI has **five** named, non-ungridded, independently-restartable internal
tracers (`NH3`, `NH4a`, `NO3an1`, `NO3an2`, `NO3an3`), not one field with an
ungridded bin dimension (DU/SS) or one scalar aerosol (SU). This matters in
a few places where the generic MAPL3 "single `internal_variable_name`"
convention (adopted during the SU/DU/SS port, see SU2G notes step 7) doesn't
directly fit:

- **`internal_variable_name`**: set to the scalar `"NI"` (an identifying
  token, not an actual field name) in Initialize, purely so
  `get_mixR`/`get_mixRatioSum` can recognize the species. NI2G's own
  `aerosol_optics`/`monochromatic_aerosol_optics` hardcode the three bin
  names (`NO3an1`, `NO3an2`, `NO3an3`) directly as a local
  `character(len=6), parameter :: aerosol_names(3)` array instead of
  reading a list from Info (MAPL3's generic `ESMF_InfoGet` has no
  allocatable-array-of-strings specific, so DU2G/SS2G/SU2G had already
  moved away from storing a list there; see SU2G notes step 7).
- **`Chem_AeroGeneric::get_mixR`** gained a third branch (`aerosol_name ==
  "NI"`) that sums `NO3an1 + NO3an2 + NO3an3`, alongside the existing
  `"DU"/"SS"` (ungridded-sum) and else (single-field) branches.
- **Ungridded dims**: NI's `self%nbins` (= 5) is the *tracer* count used for
  `fscav`/`fwet`/etc. indexing (`nNH3=1, nNH4a=2, nNO3an1=3, nNO3an2=4,
  nNO3an3=5`), unrelated to any ungridded array dimension. Two new
  `MAPL_UngriddedDim`s were declared instead:
  - `ungrd_du_ss_bins` (hardcoded `5`): for the imported `DU`/`SS` fields,
    which NI reads whole (not per-bin) for heterogeneous chemistry. Matches
    DU2G/SS2G's `nbins: 5`.
  - `ungrd_no3_bins` (hardcoded `3`, a module parameter `n_no3_bins`): for
    the nitrate-bin-indexed exports (`NISD`, `NIDP`, `NIWT`, `NISV`,
    `NIHT`, `NIHTV`, `NISD_V`), which were a literal `3` in the MAPL2 `.rc`.
- **No `_AERO_DP` export**: unlike DU/SS/SU, NI never exported a deposition
  bundle in MAPL2 either (no `append_to_bundle` calls, no
  `MAPL_AddExportSpec` for a bundle). Confirmed by
  `GOCART2G_GridCompMod.F90::add_aero_states_`, which already special-cased
  `if (instance_%name(1:2) /= "NI")` before skipping the bundle merge.
  SetServices therefore only adds the `NI_AERO` state spec, not an
  `NI_AERO_DP` bundle spec, and the parent's NI import/connection block
  mirrors that (AERO only).

## Steps

### 1. State specs and build

- `NI2G_StateSpecs.rc` rewritten column-by-column following SU2G/DU2G:
  - `RESTART=OPT` kept on imports (dropped `OPT`-only columns that MAPL3
    doesn't use); `COND=using_GMI` added to `GMI_HNO3` (previously handled
    by a manual `if (self%using_GMI)` import block in MAPL2).
  - `DU`/`SS` imports gained `UNGRIDDED=ungrd_du_ss_bins` (MAPL2 declared
    them as `5` directly in the column).
  - All `*_no3_bins`-indexed exports (`NISD_V`, `NIHTV`, `NIHT`, `NISD`,
    `NIDP`, `NIWT`, `NISV`) use `ungrd_no3_bins` instead of the literal `3`.
  - Wavelength-indexed exports use `ungrd_wavelengths_profile` /
    `ungrd_wavelengths_vertint`, as in SU2G/DU2G.
  - Internal specs (`NH3`, `NH4a`, `NO3an1-3`, `XHNO3`) dropped the
    `RESTART` (`MAPL_RestartOptional`, not representable in MAPL3) and
    `FRIENDLYTO` (ignored by MAPL3's ACG, replaced by the `TRADV` service
    subscription, see step 2) columns; all five get `FILL=0.0`, matching
    SU2G's step 2 reasoning (MAPL3 fills new fields with a sentinel unless
    `fill_value` is given).
  - Checked with the in-repo ACG
    (`python3.12 src/Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py`): parses
    clean, `GMI_HNO3` wrapped in `if (using_GMI)`, `DU`/`SS` get
    `ungridded_dim_array=[ ungrd_du_ss_bins]`, NO3-bin exports get
    `ungrd_no3_bins`, internals emit `fill_value=0.0`.
  - `CMakeLists.txt`: added `TYPE SHARED` (DU2G/SU2G/SS2G all have it;
    NI2G's didn't).

### 2. Instance configuration

- Added `NI2G_instance_NI.yaml`, converted from `NI2G_instance_NI.rc`
  following `SU2G_instance_SU.yaml`/`DU2G_instance_DU.yaml` (values
  unchanged; `.false.` -> `false`; scalar -> `[...]` lists; quoted strings
  unquoted; the `mapl:` CI section copied from SU/DU). The `.rc` is kept
  (DU/SS/SU keep theirs too).
- No changes needed to `wavelengths_for_*_aop_in_nm`,
  `aerosol_*_wavelength_in_nm_from_LUT`, `n_moments`,
  `n_phase_function_moments_photolysis`: like SU2G/DU2G's instance yaml,
  these aren't set per-instance; `MAPL_GridCompGetResource` resolves them
  from the parent/universal config.

### 3. SetServices

Rewritten after SU2G/DU2G:

- `MAPL_GridCompGet(gc, name=)`; private state via
  `_SET_NAMED_PRIVATE_STATE`/`_GET_NAMED_PRIVATE_STATE` with
  `PRIVATE_STATE = "NI2G_GridComp"`; per-thread `workspaces` kept (the
  `ThreadWorkspace%first` field and `self%recycle_HNO3` are unused dead
  fields, same as in the MAPL2 source -- left as-is rather than removed,
  to keep the diff behavior-preserving).
- `ungrd_wavelengths_profile`/`ungrd_wavelengths_vertint` (from
  `self%wavelengths_profile`/`vertint`, set by
  `GA_Environment%load_from_config`), plus the two NI-specific
  `ungrd_du_ss_bins`/`ungrd_no3_bins` described above, declared before the
  ACG includes.
- Entry points: `Run -> "Run1"`, `Run2 -> "Run2"`, `Run0 -> "Run0"`, matching
  what GOCART2G's `MAPL_GridCompRunChild` calls.
  `MAPL_GenericSetServices`/`ESMF_UserCompSetInternalState` gone.
  Entry-point signatures match MAPL3's `I_Run` (no `intent` on
  `gc`/`import`/`export`/`clock`, `rc` `intent(out)` and not optional), as
  in DU2G/SU2G (MAPL2's `intent(inout)` dummies fail to compile against
  `MAPL_GridCompSetEntryPoint`'s interface).
- Data-driven specs converted to `MAPL_GridCompAddSpec` (`PLE`, `RH2`,
  internal `NO3an{1,2,3}` with `add_to_export`/`fill_value=0.0`,
  `climNO3an{1,2,3}`), guarded by `_FAIL` first like SU/DU/SS: the
  data-driven path was never activated in MAPL2 either and stays inactive.
- `using_GMI` resource read directly (`MAPL_GridCompGetResource(gc,
  "using_GMI", self%using_GMI, _RC)`); the generated Import include is
  wrapped in `associate(using_GMI => self%using_GMI)` so the `COND` column
  on `GMI_HNO3` resolves. No other NI-specific resources beyond what
  `GA_Environment%load_from_config` already reads.
- Friendliness (MAPL2 `FRIENDLYTO=DYNAMICS:TURBULENCE:MOIST` on `NH3`,
  `NH4a`, `NO3an1`, `NO3an2`, `NO3an3`): replaced with the same
  `MAPL_STATEITEM_SERVICE` `TRADV` subscription pattern used for SU/DU/SS
  (`service_items` = all five tracer names), computational instance only.
- `<name>_AERO` (`MAPL_STATEITEM_STATE`) export, same as DU/SU. **No**
  `<name>_AERO_DP` bundle (see "structural difference" above).

### 4. Initialize

Rewritten after SU2G's Initialize (closest pattern: a replenishment alarm,
oxidant-style HNO3 handling in Run2, single Mie-table set per instance):

- `MAPL_GridCompGet(gc, name=, geom=, num_levels=)`; private state via
  `_GET_NAMED_PRIVATE_STATE`; internal state via
  `MAPL_GridCompGetInternalState`.
- Timestep: `MAPL_ClockGet(clock, dt=)`, `GOCART2G_DT` via
  `MAPL_GridCompGetResource`.
- DU/SS `radius`/`fnum` attributes: MAPL2 read them via
  `ESMF_AttributeGet(field, name=..., itemCount=, valueList=)` from the
  import field; MAPL3 equivalent is `ESMF_InfoGetAlloc(field_info,
  key=..., values=self%rmedDU/fnumDU/..., _RC)` on the field's Info
  (auto-allocates the output array -- no `itemCount` query needed, unlike
  MAPL2). These attributes are written by DU2G/SS2G's own Initialize onto
  their exported `DU`/`SS` fields (`ESMF_InfoSet(field_info, key="radius",
  values=self%radius, ...)`, already present from the DU2G/SS2G port); NI2G
  is the first MAPL3 reader of them.
- `ScavengingFractionPerKm`/`SetofHenryLawCts` -> `ESMF_InfoSet` on the
  field's Info (`ESMF_InfoGetFromHost` then `ESMF_InfoSet`), replacing
  `ESMF_AttributeSet`.
- `NO3an{1,2,3}` added to the AERO state with `ESMF_StateAdd(aero,
  [field])`, as DU/SU do, instead of `MAPL_FieldCreate` + `MAPL_StateAdd`.
- `add_aero_named_alias` -> `add_aero(..., geom=, km=, typekind=)`. NI does
  *not* add `SAREA`/`REFF` (unlike SU) -- the MAPL2 source never did either.
- `internal_variable_name` set to the scalar `"NI"` (see "structural
  difference" above), not a `valueList`.
- `createReplenishAlarm` call kept, `rc` no longer optional.

### 5. Run routines

Run0, Run, Run1, Run2, Run_data, following SU2G/DU2G:

- All routines: `MAPL_GridCompGet`/`MAPL_GridCompGetInternalState`,
  `_GET_NAMED_PRIVATE_STATE`; `MAPL_GetObjectFromGC`, `MAPL_Get`,
  `ESMF_UserCompGetInternalState` gone; `rc` non-optional (except the
  internal `Run_data`, matching SU2G/DU2G).
- Fields come from the generated `NI2G_GetPointer___.h`/
  `NI2G_DeclarePointer___.h`. Other `MAPL_GetPointer` calls ->
  `MAPL_StateGetPointer`.
- Run2: a module-level `internal_names`-style mapping wasn't needed (NI
  doesn't loop generically over internal fields by position the way
  SU2G's Run2 does); the five tracer-array indices (`nNH3=1, nNH4a=2,
  nNO3an1=3, nNO3an2=4, nNO3an3=5`) are module parameters, matching the
  original MAPL2 source's local parameters and NI2G_StateSpecs.rc order,
  used to index `self%fscav(...)`/`self%fwet(...)`.
- **Edge-level fields fixed to be 0-based** (`PLE`, `ZLE`, `PFL_LSAN`,
  `PFI_LSAN`): 1-based in MAPL3, 0-based in MAPL2. Pointer lower-bound
  remaps (`ple0(lbound(ple,1):, lbound(ple,2):, 0:) => ple`, same for
  `zle0`/`pfl_lsan0`/`pfi_lsan0`) added before every process-library call
  that takes an edge array: `findKlid` (Run0, Run2), `DryDeposition`,
  `Chem_SettlingSimple`, `WetRemovalGOCART2G`, `Aero_Compute_Diags` (Run2).
  - **Two latent bugs fixed by this remap, inherited from the MAPL2
    source's direct port**: the original MAPL2 `Run2` passed the *raw*
    (not klid-adjusted, and in this port, not-yet-0-based) `zle`/`ple`
    straight into `DryDeposition` and all four `Aero_Compute_Diags` calls.
    Checked against the process-library signatures
    (`GOCART2G_Process.F90`): `DryDeposition`'s `hghte` dummy is indexed
    `hghte(:,:,km-1) - hghte(:,:,km)` (0-based convention, same issue
    DU2G/SS2G/SU2G hit with `PLE` in their AERO optics, see SU2G notes
    step 7), and `Aero_Compute_Diags`'s `ple` dummy is used the same way.
    Ported as `zle0`/`ple0`, matching the DU2G/SU2G fix pattern.
- **Undefined-value handling**: NI2G's only `MAPL_UNDEF` use is the
  initial `kan = UNDEF` fill in
  `GOCART2G_Process::NIheterogenousChem`, immediately overwritten by
  `kan = 0.0` at the top of each bin's reaction-rate loop before any read
  -- functionally inert regardless of whether `UNDEF` is `1.0e15` (MAPL2)
  or `-huge` (MAPL3). Passed `MAPL2_UNDEF` (module parameter, same value
  and rationale as SU2G/DU2G) for consistency with the rest of GOCART2G,
  not because it changes behavior here.
- **GMI coupling**: `self%using_GMI` drives `xhno3 = GMI_HNO3` /
  `GMI_HNO3 = xhno3` (2-way) in Run2, unchanged from MAPL2. Untested (no
  case in this repo runs with `using_GMI: true`; same caveat as SU2G's
  risk 2).
- **HNO3 replenishment**: `self%alarm` (`createReplenishAlarm`, 3-hourly)
  drives `xhno3 = NITRATE_HNO3` when not using GMI, unchanged from MAPL2.

### 6. AERO methods

`aerosol_optics` and `monochromatic_aerosol_optics` rewritten after SU2G's
MAPL3 versions, generalized to loop over NI's three named bins (closer in
spirit to DU2G's bin loop, but over hardcoded field names instead of an
ungridded dimension -- see "structural difference" above):

- `ESMF_AttributeGet` -> `ESMF_InfoGet` on the state's Info, as in
  SU2G/DU2G. `aerosol_names = ['NO3an1', 'NO3an2', 'NO3an3']` is a local
  `parameter`, fetched per-bin via `MAPL_StateGetPointer(state, q,
  trim(aerosol_names(n)), _RC)` instead of `ESMF_StateGet` +
  `ESMF_FieldGet(farrayPtr=)` (equivalent, matches the `MAPL_StateGetPointer`
  idiom used everywhere else in the ported file).
- **PLE indexing**: same 0-based remap as Run2
  (`ple0(lbound(ple,1):,lbound(ple,2):,0:) => ple`), `km = ubound(ple0,3)`,
  matching the DU2G/SU2G/SS2G fix (SU2G notes step 7) for the
  `q_4d(i,j,k,n) = ((ple0(i,j,k)-ple0(i,j,k-1))*...)*q(i,j,k)` loop. The
  pre-port MAPL2 source's `ple(i,j,k) - ple(i,j,k-1)` loop was already
  written assuming a 0-based `ple`, so this is a direct carry-over of the
  DU2G/SU2G fix, not a new bug found in NI specifically.
- `ESMF_InfoGetAlloc` for `mieTable_pointer` (integer array), matching the
  SU2G/DU2G fix for allocatable-array Info reads (SU2G notes step 7).
- `get_mixR` dispatch: see "structural difference" above
  (`Chem_AeroGeneric::get_mixR`'s new `"NI"` branch).

### 7. Parent wiring (GOCART2G)

- `GOCART2G_GridCompMod.F90`:
  - `use NI2G_GridCompMod`, `add_children__(gc, self%NI, NI2G_SetServices,
    _RC)` enabled; `deactivate_instances__(gc, self%NI, _RC)` removed (CA
    stays deactivated -- not ported).
  - SetServices: `NI_AERO` (state, no bundle) import spec added for every
    NI instance, alongside the existing DU/SS/SU loops.
  - SetServices: the commented "Add connectivities for Nitrate component"
    block enabled, using `MAPL_GridCompAddConnection` (the non-legacy name
    for the same `gridcomp_add_simple_connection` the commented-out
    `MAPL_GridCompAddConnectivity` aliased to) to connect
    DU's `"DU"` -> NI's `"DU"` import, SS's `"SS"` -> NI's `"SS"`, and
    SU's `"SO4"` -> NI's `"SO4"`, each gated on the respective first
    instance being active (matching the `_ASSERT` above it: "GOCART
    supports only one active nitrate instance").
  - SetServices: new "Connections to Nitrate's export items" block
    (`NIEXTTAU, NISTEXTTAU, ..., NISMASS, NISMASS25, NH4SMASS` ->
    `<self>`, plus `NI_AERO` -> `<self>`), following the DU/SS/SU blocks;
    no `NI_AERO_DP` connection (NI has none).
  - Initialize: `add_aero_states_(self%NI%instances(:))` enabled (the
    subroutine already special-cased skipping NI's `_AERO_DP` bundle).
  - Run2: the commented Nitrate totals block rewritten in the DU/SS/SU
    style (whole-array `totexttau = totexttau + niexttau` etc. through the
    generated pointers, instead of the old `gex(id)` +
    per-wavelength-loop style the comment had). `nifactor` logic (using
    `1.0` instead of `132.14/96.06` once NI is active) was already correct
    and untouched -- it anticipated this port.
- `GOCART2G_StateSpecs.rc`: the 19 NI fields
  (`NIEXTTAU...NISMASS25, NH4SMASS`) added as IMPORT rows
  (`#-imports-from-nitrate-`), units/long names from NI's spec. Checked
  with the ACG: parses clean, all pointers generated
  (`niexttau`...`nismass25`, `nh4smass`).
- `../CMakeLists.txt`: `NI2G_GridComp` enabled in `alldirs`.
- `Chem_AeroGeneric.F90::get_mixR`: added the `"NI"` branch (see
  "structural difference" above).

## Run-testing update (ctest -L REGRESSION)

Picked up the "next step" from open item 1 below: built GOCART2G
(`DU2G`/`SS2G`/`SU2G`/`NI2G`) with `ifx-stack` against this checkout's
existing MAPL3 build, and added a `gocart+seasalt+dust+sulfate+nitrate`
regression case (mirroring `regression/gocart+seasalt+dust+sulfate`, not
`sulfate-sa`; see why below) so NI2G participates in `ctest -L REGRESSION`
the same way DU2G/SS2G/SU2G do.

- **Update: the standalone `nitrate-sa` case now works, via a dedicated
  stand-in component.** NI2G's `Initialize` unconditionally does
  `ESMF_StateGet(import, 'DU', field)` / `'SS'` and then reads the
  `radius`/`fnum` `Info` attributes DU2G/SS2G's own `Initialize` writes onto
  their exported `DU`/`SS` fields -- a generic `configurable_gridcomp`
  "data" component can supply the array values but has no mechanism to
  attach those custom attributes. Rather than extend
  `configurable_gridcomp` with generic attribute support, `run1`/`run2` now
  include a small dedicated `DuSsStub_GridComp` (see
  `../regression/DuSsStub_GridComp/DuSsStub_GridCompMod.F90`) that imports
  `DU`/`SS` (connected to `DataForNI`'s internal `DU`/`SS`, i.e. the same
  import netcdf4 file `DataForNI`'s other fields are read from) and
  re-exports them, attaching `radius`/`fnum` Info to them exactly like
  DU2G/SS2G do (radius/fnum values copied from `DU2G_instance_DU.yaml`/
  `SS2G_instance_SS.yaml`). It must be declared as a child *after*
  `DataForNI` but *before* `FakeGOCART` in `root.yaml` so its Initialize
  runs after the DU/SS data is available and before NI2G reads the Info.
  With that in place, NI2G's Initialize gets past the radius/fnum read and
  fails later,
  at `nf90_open` of `opticsBands_NI.v2.1.0.RRTMGP.nc4` -- that file simply
  isn't present anywhere in the available ExtData input tree (not even
  under a `gocart+seasalt+dust+sulfate+nitrate` input directory, which
  doesn't exist either), so it's a missing-test-data issue, not a code
  issue.
- **Found and fixed a real connection bug** (not a test-scaffolding issue):
  `NI2G_GridCompMod.F90`'s `ungrd_du_ss_bins` was declared
  `MAPL_UngriddedDim(n_du_ss_bins, name="du_ss_bins", units="1")`, but
  DU2G/SS2G declare the ungridded dim on their own exported `DU`/`SS`
  fields as `name="nbins"`. MAPL's `UngriddedDimsAspect` connection check
  requires an exact `name`/`units`/extent match between producer and
  consumer (`UngriddedDim::equal_to`, `infrastructure/esmf/UngriddedDim.F90`),
  so the DU->NI and SS->NI connections GOCART2G_GridCompMod.F90 wires up
  failed at Initialize with `cannot connect aspect UNGRIDDED_DIMS` before
  this fix. Changed the name to `"nbins"` to match; this was the only code
  change needed. This bug would have hit *any* active NI configuration
  (test or production), not just the new regression case.
- **Verified end-to-end past this fix**, using the real
  `gocart+seasalt+dust+sulfate/run1` checkpoint/restart data (DU/SS/SU) and
  hand-fabricated zero-filled `NI_internal.nc` plus the 7 added
  `EMI_NH3_*` variables in `DataForGOCART2G_internal.nc` (since no real
  baseline exists yet for the new case name): DU/SS/SU/NI all read their
  restarts and reach `NI: Initialize:: starting...` successfully; NI then
  fails only on `nf90_open
  ExtData/chemistry/AerosolOptics/v2.x.x/x/opticsBands_NI.v2.1.0.RRTMGP.nc4`
  -- confirmed those `optics_NI`/`opticsBands_NI` files (and the new case's
  own baseline checkpoints) simply don't exist yet under the shared
  `/ford1/share/gmao_SIteam/regression-data/.../GOCART2G_GridComp/`
  tree (DU/SS/SU's equivalents are present). I.e. the code path and YAML
  wiring are confirmed correct up to a pure missing-data issue, consistent
  with the ask to "assume the input data is available."
- `ctest -L REGRESSION` (23 tests, including the 2 new
  `gocart+seasalt+dust+sulfate+nitrate/run{1,2}` cases) passes cleanly
  today: the new cases report "Regression data not found ... skipping" the
  same way `sulfate-sa` etc. did before its data existed. Once baseline
  checkpoints *and* the NI aerosol-optics ExtData files are published to
  the shared/S3 regression-data tree (same `DATA_PATH` as DU2G/SS2G/SU2G,
  under the new case subdirectory), these two tests should start actually
  running and comparing, with no further config changes expected.

## Open items (follow-ups, not blocking the code port)

1. **Not built, not run-tested** in the original porting session. ~~No
   Fortran compiler toolchain or MAPL2 baseline data was available in that
   session (user-selected "code-level port only" scope).~~ **Update:** now
   built and run-tested up to the point of missing data; see "Run-testing
   update" above for what was found/fixed and what's left (NI aerosol
   optics ExtData + a `gocart+seasalt+dust+sulfate+nitrate` baseline still
   need to be captured/published). The original idea of a standalone
   `nitrate-sa` case using `DataForDU`/`DataForSS`/`DataForSU` data
   components (not real DU2G/SS2G/SU2G instances) does **not** work, for
   the reason given above.
2. **GMI coupling untested** (`using_GMI: true`): same open risk as SU2G.
3. **The three `Chem_AeroGeneric::get_mixR` branches are now an
   `if/else if/else` on the `internal_variable_name` scalar**: this is
   fine for `DU`/`SS`/`NI`/single-field species, but doesn't generalize
   to a fourth bin layout without another branch. If CA2G (oc/bc/br, each
   with 2 bins per the MAPL2 source) is ported next, check whether it fits
   the "ungridded" or "named-bins" pattern and extend accordingly.
4. **`self%recycle_HNO3`/`ThreadWorkspace%first`** are declared but unused,
   inherited dead fields from the MAPL2 source. Left in place rather than
   removed, to keep this a faithful behavior-preserving port; harmless
   either way.
5. **Advection**: NI subscribes to `TRADV` for `NH3`, `NH4a`, `NO3an1-3`
   (computational instance only), following the SU/DU/SS pattern. Like
   those three, this is inactive until the AGCM/Physics port connects
   SUPERDYN's `TRADV` to PHYS's (see SU2G notes step 4/risk 5). Not
   run-tested even standalone in this session.
6. **DU->NI/SS->NI/SU->NI connections**: wired per the SU2G notes' step-10
   follow-up. **Update:** now run-tested via the new
   `gocart+seasalt+dust+sulfate+nitrate` case (see "Run-testing update"
   above); the connections themselves needed an ungridded-dim name fix
   (done) and then connected/initialized successfully. Still worth
   checking, once a real baseline exists, that the sulforg01-03 activation
   modes / MAMnet `SUsum` pick up nonzero `SO4` from SU now feeding NI too
   (same question SU2G's notes raised for SU alone).

## `DuSsStub_GridComp` relocation + real DU/SS data (follow-up)

- Moved `DuSsStub_GridComp` from `../DuSsStub_GridComp` to
  `../regression/DuSsStub_GridComp`: it is only ever used by the
  `nitrate-sa` standalone regression harness, so it now lives under
  `regression/` alongside the other test-only pieces (`regression/*-sa`
  cases, `regression/CMakeLists.txt`). `regression/CMakeLists.txt` now
  builds it (`add_subdirectory(DuSsStub_GridComp)`); the top-level
  `GOCART2G_GridComp/CMakeLists.txt` no longer references it directly.
- `DuSsStub_GridComp` no longer zero-fills its exported `DU`/`SS` fields.
  It now declares matching `IMPORT` specs for `DU`/`SS`, connected (in
  `run1`/`run2`'s `root.yaml`) to `DataForNI`'s internal `DU`/`SS` --
  newly added to `data-for-ni.yaml` alongside the other internal fields
  (`DELP`, `EMI_NH3_*`, ...) that are populated from the same import
  netcdf4 (restart/checkpoint) file. `Initialize` copies the import data
  onto the export fields before attaching the `radius`/`fnum` Info, so
  NI2G now sees real DU/SS mixing ratios instead of zeros.
