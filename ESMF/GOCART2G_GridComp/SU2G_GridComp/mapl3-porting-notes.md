# SU2G_GridComp MAPL3 porting plan

Status: in progress on `feature/pchakrab/port-su2g-to-mapl3`. Steps 1, 1a, 1b, 1c, 2, 3, 4, 5, 6, 7 and 8 done.

Derived from the Radiation porting notes
(`@GEOSradiation_GridComp/GEOSsolar_GridComp/mapl3-porting-notes.md` and
`GEOSirrad_GridComp/mapl3-porting-notes.md`) and from the in-repo DU2G/SS2G
ports, which are the closest templates. DU2G was ported routine-by-routine
(`07a4c2b` ... `51ba8bf`), and SU2G has the same routine set: SetServices,
Initialize, Run0, Run, Run1, Run2, Run_data, aerosol_optics,
monochromatic_aerosol_optics.

## Starting point
- `SU2G_GridCompMod.F90`: 1791 lines, about 236 calls still on the MAPL2 API.
- `SU2G_StateSpecs.rc`: 17 rows whose sizes use the MAPL2 `self%...` form.
- `GOCART2G_GridComp.yaml` already lists `SU: [SU]` as active. The parent
  doesn't create it: `add_children__` and `use SU2G_GridCompMod` are
  commented out.

## Reference files
- Template component: `../DU2G_GridComp/DU2G_GridCompMod.F90`
- Template state spec: `../DU2G_GridComp/DU2G_StateSpecs.rc`
- Template build file: `../DU2G_GridComp/CMakeLists.txt`
- Template instance config: `../DU2G_GridComp/DU2G_instance_DU.yaml`
- Parent: `../GOCART2G_GridCompMod.F90`

## Steps

0. **MAPL2 baseline (do this first).**
   - Capture an SU-only MAPL2 run at C12L181, the way the Solar and GWD
     baselines were captured: `COMPONENT_TO_RECORD: SU`, with before and
     after Run-phase files for import, internal and export.
   - The baseline must match the standalone's configuration, or it will
     never reproduce bit-for-bit (Solar "Comparison status"; GWD grid name,
     timestep and efficiency settings).
   - Record every SU and GOCART2G setting the baseline used: `using_GMI`,
     `disable_emissions`, the volcano and point-emission files,
     `DIURNAL_BIOMASS_BURNING`, `GOCART2G_DT`, `MAXCLEAN`, `CCNTUNING`,
     `USE_MAMNET`.

1. **DONE - codee format pass** (before any other change).
   `codee format --verbose ~/config/codee-format SU2G_GridCompMod.F90`,
   following `.vscode/instructions/fortran-formatting.instructions.md`.
   - Five `allocate(x, mold=y, __STAT__)` lines (L1105, L1247-1249, L1298
     before formatting) were syntax errors to codee; masked with
     `!CODEE_MASK!`, formatted, unmasked. codee re-indented the masked
     comment lines, so unmasking kept codee's indentation; each statement
     text is present unchanged.
   - Verified layout-only: code token stream (lowercased outside strings,
     whitespace/`&`/`;` removed, codee's `'a'//'b'` literal splits joined),
     comments and cpp lines all identical; subroutine/end-subroutine,
     `contains`, `#if`/`#endif` counts unchanged. 1791 -> 1800 lines.

1a. **DONE - Style-only pass (no logic changes)**, as in Solar step 1.
   codee (step 1) had already lowercased keywords/`intent` and indented
   every continuation line +5, so this step only changed macro lines
   (1800 -> 1811 lines, from the `__Iam__` expansions); every changed line
   is a macro rename.
   - Error-handling macros: `__RC__`/`__STAT__`/`VERIFY_`/`RETURN_` ->
     `_RC`/`_STAT`/`_VERIFY`/`_RETURN`; each `__Iam__('x')` expanded to
     its two declarations (`integer :: status`,
     `character(len=255) :: Iam = 'x'`); `#include "MAPL_Generic.h"` ->
     `#include "MAPL.h"`. All expand identically (`MAPL_ErrLog.h`,
     `MAPL_Exceptions.h`; `MAPL_Generic.h` only includes `MAPL.h`).
   - Lowercase keywords and dot-operators in code; `intent( IN )` ->
     `intent(in)`; continuation lines indented 5 spaces from the
     statement.
   - Verified by normalizing both versions (cpp lines and comments as-is,
     code lowercased with whitespace removed outside strings): identical.
     Macro spellings inside comments were left unchanged.

1b. **DONE - Comment cleanup** (1811 -> 1691 lines).
   - Removed 116 decorative comment lines: `!` followed only by `-`/`=`/
     `*` runs, bare `!`, and ProTeX `!BOP`/`!EOP`.
   - `!   text` -> `! text` and `!   ! text` -> `! text` (full-line and
     trailing comments); trailing whitespace stripped; blank-line runs
     collapsed. ProTeX headings `! !ROUTINE:` -> `! ROUTINE:` (44 lines).
     `!$` OpenMP sentinels kept.
   - Verified: code and cpp lines identical; the 286 remaining comments
     have the same text as before modulo whitespace.

1c. **DONE - codee format re-run** after 1a/1b (same masking of the five
   `allocate(..., mold=..., _STAT)` lines). codee fixed 161 issues, mostly
   identifier case (e.g. `GC`/`RC` -> `gc`/`rc`); 232 changed lines, line
   count unchanged (1674). Verified layout/case-only with the same
   token-stream check as step 1.

2. **DONE - State specs and build.**
   - `SU2G_StateSpecs.rc`:
     - `self%nbins` -> `ungrd_nbins`.
     - `size(self%wavelengths_*)` -> `ungrd_wavelengths_profile` /
       `ungrd_wavelengths_vertint`, matching DU and SS. SetServices must
       declare and set these before the ACG includes (step 4).
     - All five internals get `FILL=0.0`. MAPL2 zero-filled internals;
       MAPL3 fills new fields with a sentinel (`FieldFill` in
       `field_empty_complete`) unless `fill_value` is given.
     - Kept `ADD2EXPORT=T`; removed the `FRIENDLYTO` column (ignored by
       MAPL3's ACG).
     - `RESTART` column removed from the internals (was
       `MAPL_RestartOptional`). MAPL3 has only `RESTART_REQUIRED` (the
       default) and `RESTART_SKIP` (`infrastructure/field/RestartModes.F90`);
       the ACG maps `OPTIONAL` to `MAPL_RESTART_OPTIONAL`, which doesn't
       exist and fails to compile. Leaving it blank gives the default, the
       same as DU2G (whose `MAPL_RestartOptional` is silently dropped).
   - Checked with the in-repo ACG
     (`python3.12 src/Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py`): 17
     `ungridded_dim_array=[ ungrd_*]`, and internals emit
     `fill_value=0.0`.
   - `CMakeLists.txt`: added `TYPE SHARED` and `NetCDF::NetCDF_Fortran`, as
     DU2G has.

3. **DONE - Instance configuration.**
   - Added `SU2G_instance_SU.yaml`, converted from `SU2G_instance_SU.rc`
     following `DU2G_instance_DU.yaml` (values unchanged; `.false.` ->
     `false`; lists in `[...]`; the `mapl:` CI section copied from DU).
     The `.rc` is kept, as DU and SS keep theirs.
     - `fwet_ice`/`fwet_snow`/`fwet_rain` were a single `1.0` in the `.rc`;
       MAPL2 read them as `nbins`-arrays with `default=1.0`, so every bin
       was 1.0. MAPL3's `load_from_config` reads an array, so the YAML
       lists four `1.0` values.
     - Integer-looking reals (`-1`, `1700`) written as `-1.0`, `1700.0`.
   - No `.data` variant: the parent loads
     `<species>_instance_<child>.yaml`, so `SU2G_instance_SU.data.yaml` is
     only needed if `SU.data` is made active. DU and SS don't have one
     either.
   - SetServices: the `.rc` loading (with its `SU2G_instance_SU.rc`
     fallback) is gone; `GA_Environment%load_from_config(gc)` and
     `MAPL_GridCompGetResource(gc, "label", ...)` read the child's YAML.
     `volcano_srcfilen_*` and `point_emissions_srcfilen` became
     `character(len=:), allocatable`; `aviation_layers` became an
     allocatable array with an `_ASSERT` that it has 4 values.
   - Initialize still loads the `.rc` (optics files, `n_moments`, the
     universal config); that moves in step 5.

4. **DONE (except friendliness) - SetServices.** Rewritten after DU/SS:
   - `MAPL_GridCompGet(gc, name=, logger=)`; private state via
     `_SET_NAMED_PRIVATE_STATE`/`_GET_NAMED_PRIVATE_STATE` with
     `PRIVATE_STATE = "SU2G_GridComp"`; per-thread `workspaces` kept.
     The module still has `type wrap_`, used by the other routines until
     they are ported (steps 5-7); they won't find the state until then.
   - `ungrd_nbins` / `ungrd_wavelengths_profile` /
     `ungrd_wavelengths_vertint` (`MAPL_UngriddedDim`) declared and set
     before the ACG includes.
   - Entry points: `Run` -> `"Run1"`, `Run2` -> `"Run2"`, `Run0` ->
     `"Run0"` (the names GOCART2G's `MAPL_GridCompRunChild` calls use).
     `MAPL_GenericSetServices` and `ESMF_UserCompSetInternalState` gone.
   - Entry-point signatures: every routine passed to
     `MAPL_GridCompSetEntryPoint` (Initialize, Run, Run2, Run0) must match
     MAPL3's `I_Run` (`infrastructure/esmf/ESMF_Interfaces.F90`): `gc`,
     `import`, `export`, `clock` declared with no `intent`, `rc`
     `intent(out)` and not optional, as in DU2G. MAPL2's
     `intent(inout)` dummies fail with ifx errors #7061-7063. (Missed in
     the original plan and fixed at the first build; see step 9.)
   - Data-driven specs converted to `MAPL_GridCompAddSpec` (`PLE`, `RH2`,
     internal `SO4` with `add_to_export`/`fill_value=0.0`, `climSO4`,
     per-bin `climSU{DP,WT,SD,SV}nnn`), but guarded by `_FAIL` first, as
     SS does: the data-driven path is not active. The per-bin fields were
     `xy` + `VLocationCenter` in MAPL2; now `xy` + `NONE`.
   - Manual `SU_H2O2/OH/NO3` imports dropped: they duplicated the spec
     file's rows (MAPL2 tolerated the duplicates). The manual
     `if (using_GMI)` `GMI_*` imports moved into the spec file instead:
     new `COND` column (`using_GMI`), units `1` -> `mol/mol` (what Run2
     assumes, and what the manual specs said). The Import include is
     wrapped in `associate(using_GMI => self%using_GMI)`. ACG checked: 40
     import specs as before, the three `GMI_*` inside `if (using_GMI)`;
     export and internal output unchanged.
   - `<name>_AERO` (`MAPL_STATEITEM_STATE`) and `<name>_AERO_DP`
     (`MAPL_STATEITEM_FIELDBUNDLE`) as in DU.
   - `MAPL_AM_I_ROOT` prints -> `logger%info`; `Iam` and `DEFVAL` dropped.
   - Import `RESTART` cells say `OPT`, which the ACG doesn't map
     (only `OPTIONAL`, `SKIP`, ...), so they're silently dropped; `SKIP`
     works. Left as is; pre-existing. (MAPL2 `FRIENDLYTO=DYNAMICS:TURBULENCE:MOIST` on DMS,
     SO2, SO4, MSA). The MAPL3 ACG ignores `FRIENDLYTO` (DU2G's spec still
     carries it as a no-op), so nothing currently puts SU's tracers into
     the advection/diffusion/convection bundles. The MAPL3 replacement is
     `MAPL_STATEITEM_SERVICE`:
     - Subscribe to advection: an import `TRADV` with
       `itemtype=MAPL_STATEITEM_SERVICE` and
       `service_items=['DMS','SO2','SO4','MSA']` (`MAPL_GridCompAddSpec`
       argument; `items:` in YAML). DynCore/AdvCore already provide
       `TRADV`, and Superdyn re-exports it; no component subscribes yet.
     - Turbulence and moist have no MAPL3 service providers yet; subscribe
       once they exist.
     - Check whether the subscription belongs in SU or in a GOCART2G/PHYS
       ancestor, and do the same for DU and SS.
     - DONE (SU only): SetServices subscribes to `TRADV` with
       `service_items` DMS/SO2/SO4/MSA (a `StringVector`), for the
       computational instance only. Placement settled: `items` must be the
       subscriber's own internals, so it has to be in SU, not an ancestor.
       Builds; not run-tested (first Fortran use of `service_items`; the
       FV3 `adv-dyn`/`dyn-sa` regressions subscribe via YAML, with sibling
       connections only). Open: whether a nested subscriber's `TRADV`
       propagates up to a DYN/ADV connection in the parent, and duplicate
       names if two SU instances subscribe. DU and SS still to do.

5. **DONE - Initialize.** Rewritten after DU2G's MAPL3 Initialize:
   - `MAPL_GridCompGet(gc, name=, geom=, num_levels=)`; private state via
     `_GET_NAMED_PRIVATE_STATE`; internal state via
     `MAPL_GridCompGetInternalState`. `MAPL_GetObjectFromGC`, `MAPL_Get`,
     `MAPL_GridGet` and `MAPL_GenericInitialize` gone. The unused
     `lons`/`lats` fetch was dropped.
   - Timestep: `MAPL_ClockGet(clock, dt=)`; `GOCART2G_DT`,
     `DIURNAL_BIOMASS_BURNING` (default `YES`) and all optics settings
     (`aerosol_*_optics_file`, `n_moments`,
     `n_phase_function_moments_photolysis`, the LUT wavelengths) via
     `MAPL_GridCompGetResource`. The `.rc` loading and the universal
     `ESMF_Config` are gone.
   - Attributes -> `ESMF_InfoSet`: `ScavengingFractionPerKm` on
     DMS/SO2/SO4/MSA, `SetofHenryLawCts` (`values=`), and on the AERO
     state `use_photolysis_table`, `mie_table_instance`,
     `band_for_aerosol_optics`, `wavelength_for_aerosol_optics`,
     `mieTable_pointer`, `internal_variable_name`. The MSA check now uses
     `ESMF_StateGet(..., itemType)` instead of testing the `rc` of a
     failed get.
   - SO4 goes into the AERO state with `ESMF_StateAdd(aero, [field])`, as
     DU does, instead of `MAPL_FieldCreate` + `MAPL_StateAdd`.
   - `add_aero_named_alias` -> `add_aero(..., geom=, km=, typekind=)`,
     including `SAREA` and `REFF` (R4, 3D).
   - `createReplenishAlarm` call kept. `rc` is no longer optional, as in
     DU.

6. **DONE - Run routines** (Run0, Run, Run1, Run2, Run_data, in one
   change rather than one per commit).
   - All routines: `MAPL_GridCompGet`, `MAPL_GridCompGetInternalState`,
     `_GET_NAMED_PRIVATE_STATE`; `MAPL_GetObjectFromGC`, `MAPL_Get`,
     `ESMF_UserCompGetInternalState`, `Iam` gone; `rc` non-optional (except
     the internal `Run_data`). The now-unused `type wrap_` was removed.
   - Fields come from the generated `SU2G_GetPointer___.h`, wrapped in
     `associate(using_GMI => self%using_GMI)` because the GMI imports are
     `COND`-ed (step 4). Other `MAPL_GetPointer` -> `MAPL_StateGetPointer`.
     `dummyMSA` now just points at the generated `MSA` pointer.
   - Run2: `MAPL_VarSpecGet(InternalSpec(n))` (no MAPL3 equivalent) ->
     a local `internal_names = ['DMS','SO2','SO4','MSA']` in the spec-file
     order (`nDMS..nMSA`).
   - Run2: `SAREA`/`REFF` written into the AERO state by looking up the
     field name from the AERO Info key (`surface_area_density`,
     `effective_radius_in_microns`), as `add_aero` stores them.
   - `lons`/`lats` via `MAPL_GridGetCoordinates(grid, ...)` (radians, as
     MAPL2's `MAPL_Get(lons=, lats=)`); Run1 fetches them only when
     `diurnal_bb` is on.
   - `MAPL_GetHorzIJIndex(n, i, j, grid=, lon=, lat=)` ->
     `MAPL_GeomGetHorzIJIndex(geom, lon, lat, i, j, local_indices=.true.)`.
     MAPL2 returned local indices (-1 off-PE); without `local_indices`
     MAPL3 returns global ones. DU2G's Run1 omits `local_indices`, which
     looks like a bug there (only hit with point emissions on).
   - `StrTemplate` -> `MAPL_StrTemplate`; root-only prints ->
     `logger%info` (the "not found" one in an `omp critical`, as DU).
     No timers in SU's Run routines.
   - **Undefined-value masks:** about 10 lines like
     `where(1.01*x > MAPL_UNDEF) x = 0.` in Run1 (around L930-L964).
     MAPL2's undefined value was `+huge`; MAPL3's is `-huge`
     (`InternalConstants.F90` L18-19). Ported as-is the condition is
     always true and every emission is zeroed.
     - DONE: module parameter `MAPL2_UNDEF = huge(1.)`, used in the ten
       Run1 masks and passed as `undefval` to `SulfateUpdateOxidants`
       (which applies the same `1.01*x > undefval` test to the oxidant
       climatologies, Process_Library L7158-7164). This keeps the MAPL2
       behaviour exactly. Still to check: what MAPL3 ExtData writes for
       missing data (if `-huge`, these masks no longer catch it).
     - DU2G (`where (1.01*du_src > MAPL_UNDEFINED_REAL)`, around L648 on
       `release/MAPL-v3`) has the same pattern. Check whether DU's
       emissions are being zeroed too.
   - **Edge-level fields** (`PLE`, `ZLE`, `PFL_LSAN`, `PFI_LSAN`): 1-based
     in MAPL3, 0-based in MAPL2. DONE: pointer lower-bound remaps
     (`ple0(lbound(ple,1):, lbound(ple,2):, 0:) => ple`), passed to every
     process-library call: `findKlid` (Run0, Run2), `SUvolcanicEmissions`,
     `SulfateDistributeEmissions`, `updatePointwiseEmissions` (Run1),
     `Chem_SettlingSimple`, `SulfateChemDriver`, `SU_Wet_Removal`,
     `WetRemovalUFS`, `SU_Compute_Diags` (Run2). Remaps rather than DU's
     allocated copies because most of these dummies are
     `real, pointer` (some without `intent(in)`), which an allocatable
     actual can't match. `findKlid` takes assumed shape, so it was
     unaffected either way.
   - **Oxidants:** `SulfateUpdateOxidants` writes into its `oh_clim`/
     `no3_clim`/`h2o2_clim` arguments, i.e. the `SU_OH`/`SU_NO3`/
     `SU_H2O2` imports: it zeroes negative and undefined values in
     place (L7158-7164). It does not write chemistry results back; the
     `recycle_H2O2` path only copies `h2o2_clim` into `xh2o2`. Unchanged
     from MAPL2; check at run time that MAPL3 allows the in-place
     cleanup (it modifies the provider's field).

7. **DONE - AERO methods.** `aerosol_optics` and
   `monochromatic_aerosol_optics` rewritten after DU's MAPL3 versions:
   - `ESMF_AttributeGet` -> `ESMF_InfoGet` on the state's Info; every
     field (PLE, RH, SO4, EXT/SSA/ASY/MOM, monochromatic_EXT) is found
     through the name stored under its Info key (`add_aero`), and outputs
     are skipped when that name is empty. `MAPL_GetPointer` ->
     `MAPL_StateGetPointer`.
   - The aerosol field name comes from the scalar
     `internal_variable_name` (`SO4`, set in Initialize); `nbins = 1`,
     as before (the MAPL2 list had one entry).
   - **PLE indexing:** `add_aero` creates PLE with `num_levels=km+1`, and
     MAPL3 field bounds start at 1 (`make_bounds_from_args`, `LU_Bound(1,
     num_levels)`). The layer loops use `ple(k) - ple(k-1)` for
     `k = 1..ubound(ple,3)`, so they now go through a 0-based remap
     `ple0`, as in step 6 (`km = ubound(ple0,3)` = number of layers).
     DU2G's MAPL3 optics read `ple(:,:,k-1)` at `k = 1` and take
     `km = ubound(ple,3)` (= km+1): out of bounds and one level too many.
     Same fix needed there (see step 10).
   - `get_mixR` lives in `ESMF/Shared/Chem_AeroGeneric.F90` (shared, not
     touched). It reads `internal_variable_name` with `values=` (a list)
     while SU/DU/SS store a scalar, still calls `ESMF_AttributeGet`, and
     never sets `rc`. Check when it's first exercised (step 9).

8. **DONE - Wrapper and parent wiring** (not built yet).
   - External `SetServices` wrapper after `end module SU2G_GridCompMod`, as
     in DU.
   - The standalone wrappers now have unique names: `DU2G_SetServices`,
     `SS2G_SetServices`, `SU2G_SetServices`, `GOCART2G_SetServices` (were
     all `SetServices`). Regression YAMLs that load them by symbol now say
     `setservices: <comp>_setservices_` (children `dso:` entries) or
     `userRoutine: <comp>_setservices_` (`mapl: setServices: sharedObj:`).
     Fortran callers (GOCART2G, Aerosol_GridComp, GEOS_ChemGridComp) use
     the module procedures and are unaffected.
   - `../GOCART2G_GridCompMod.F90`:
     - `use SU2G_GridCompMod`, `add_children__(gc, self%SU, ...)`,
       `add_aero_states_(self%SU%instances(:))` enabled.
     - SetServices: `<name>_AERO` (state) and `<name>_AERO_DP` (bundle)
       import specs for every SU instance, and connections from each
       active computational SU instance for `_AERO`, `_AERO_DP` and
       `SUEXTTAU, SUSTEXTTAU, SUSCATAU, SUSTSCATAU, SUEXTCOEF,
       SUEXTCOEFRH20/80, SUSCACOEF, SUSCACOEFRH20/80, SUBCKCOEF, SUANGSTR,
       SO4SMASS, PSO4`, following the SS/DU blocks.
     - `PSO4`: `MAPL_AddExportSpec(..., CHILD_ID=)` ->
       `MAPL_GridCompReexport(gc, src_comp=<first SU>, src_name="PSO4")`
       when the first SU instance is active and computational. `PSO4` is
       also connected into GOCART2G's import for `pso4tot`; check at run
       time that re-exporting and connecting the same child export works.
     - Run2: SU totals restored, after the DU/SS blocks (whole-array
       sums through the generated pointers instead of `gex(id)` and
       per-wavelength loops). As in MAPL2, `totextt25`/`totexttfm` use
       `suexttau` and `totscat25`/`totscatfm` use `suscatau` (SU has no
       PM2.5/PM1 optics). `nifactor` = 132.14/96.06 unless an active NI
       instance exists (none on MAPL3 yet). The commented SU pointer
       declarations were dropped.
   - `../GOCART2G_StateSpecs.rc`: the 14 SU fields above added as IMPORT
     rows (`#-imports-from-sulfate-`), units/long names from SU's spec.
     ACG checked: 54 import specs, pointers generated.
   - `../CMakeLists.txt`: `SU2G_GridComp` enabled in `alldirs`.
   - Instance config name: `add_children__` loads
     `<species>_instance_<child>.yaml` with species `SU`, i.e.
     `SU_instance_SU.yaml`, not `SU2G_instance_SU.yaml` (step 3). DU/SS
     regressions ship `DU_instance_DU.yaml`/`SS_instance_SS.yaml` in their
     run dirs; the `su-sa` case needs `SU_instance_SU.yaml` too.
   - Noticed, not changed: the DU `_AERO_DP` import-spec loop in
     GOCART2G's SetServices iterates over `size(self%SS%instances)`.

9. **Run-testing.**
   - First build (ifx Debug) of SU2G + GOCART2G succeeds after three
     fixes: entry-point dummies lose their `intent(inout)` (step 4 should
     have covered this); no
     `MAPL_RESTART_OPTIONAL` (step 2, and the data-driven `SO4`/`climSO4`
     specs); `MAPL_MaxMin` is now a function `MAPL_MaxMin(array, comm)`
     returning `[max, min]`, so the GMI diagnostics in Run2 log it with
     the VM's MPI communicator.
   - Add a `regression/su-sa` case: `root.yaml` (current vertical-grid
     schema), `cap_driver.yaml`, `cap_restart.yaml`, and SU's ExtData
     needs (anthropogenic/ship/biomass/aircraft/DMS emissions, the
     `SU_OH/H2O2/NO3` oxidants, volcano tables).
   - Hand runs: `export OMP_NUM_THREADS=1`, set `LD_LIBRARY_PATH` to the
     install `lib`, reset `cap_restart.yaml` before each run. Compare with
     a `cmpchk.py`-style script that counts both `+huge` and `-huge` as
     undefined.
   - Target: bit-for-bit against the step-0 baseline in all three states.
     Publish data under `regression-data/.../GOCART2G_GridComp/su-sa` and
     add `su-sa` to `regression/cases.txt`.

10. **Follow-ups unlocked by SU.**
    - Check that the `sulforg01-03` activation modes and MAMnet's `SUsum`
      input (`ESMF/Shared/GOCART2G_AerosolActivationMod.F90`) now pick up
      non-zero `SO4`.
    - NI2G's `SO4` import from SU can be connected when NI is ported; wire
      DU->NI and SS->NI at the same time.
    - DU2G fixes found while porting SU (separate PR, not on this branch):
      - Undefined-value mask: `where (1.01 * du_src > MAPL_UNDEFINED_REAL)`
        in `run1` is always true under MAPL3 (`-huge`), so `du_src` is
        zeroed everywhere and Ginoux/K14 emissions are lost. Use the same
        `MAPL2_UNDEF = huge(1.)` approach as SU (step 6).
      - Point emissions: `run1` calls `MAPL_GeomGetHorzIJIndex` without
        `local_indices=.true.`, so it gets global indices where
        `updatePointwiseEmissions` expects local ones (-1 off-PE). Only hit
        when `point_emissions_srcfilen` is set.
      - AERO optics: `aerosol_optics` and `monochromatic_aerosol_optics`
        index the 1-based PLE as if 0-based (`ple(:,:,k-1)` at `k = 1`,
        `km = ubound(ple,3)`). Use a 0-based remap as SU does (step 7).
        SS2G likely has the same pattern.

## Risks to settle early
Nothing below has been checked at run time yet (no step-0 baseline, no
`su-sa` case).

1. The undefined-value comparisons (step 6). This also affects DU on the
   current branch.
   - SU: fixed in code (`MAPL2_UNDEF = +huge` in the Run1 masks and
     `SulfateUpdateOxidants`).
   - Open: what MAPL3 ExtData writes for missing data; if `-huge`, the
     masks don't catch it.
   - DU: not fixed (step 10, separate PR).
2. GMI coupling: `using_GMI` changes which fields SU imports. Test with
   GMI off first, and confirm the baseline used the same setting.
   - Code done: the `GMI_*` imports are `COND=using_GMI` (step 4).
   - Open: no GMI-off run yet; baseline setting unknown until step 0.
3. Whether SU writes back into the `SU_H2O2` import (step 6).
   - Settled by reading the code: no chemistry write-back.
     `SulfateUpdateOxidants` only zeroes negative/undefined values in the
     `SU_OH`/`SU_NO3`/`SU_H2O2` imports in place.
   - Open: confirm at run time that MAPL3 allows modifying the import.
4. Volcano and point-emission files are read from disk, so the regression
   data must include them.
   - Not addressed: the `su-sa` case and its data don't exist yet (step 9).
5. Tracer transport (step 4, friendliness): until SU subscribes to the
   `TRADV` service (and turbulence/moist equivalents exist), DMS, SO2,
   SO4 and MSA are not advected, diffused or convected in a coupled run.
   The SU-only standalone regression is unaffected. DU and SS have the
   same gap.
   - Not addressed: no subscription, and placement (SU vs. a parent
     component) is undecided.
