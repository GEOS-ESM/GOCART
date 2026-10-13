# SU2G_GridComp MAPL3 porting plan

Status: in progress on `feature/pchakrab/port-su2g-to-mapl3`. Steps 1, 1a, 1b, 1c, 2, 3 and 4 done.

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
     - `RESTART`: `MAPL_RestartOptional` -> `OPTIONAL`. The MAPL2 spelling
       isn't in the ACG's mapping and was silently dropped (DU2G's spec
       still has this problem).
   - Checked with the in-repo ACG
     (`python3.12 src/Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py`): 17
     `ungridded_dim_array=[ ungrd_*]`, and internals emit
     `restart_mode=MAPL_RESTART_OPTIONAL` and `fill_value=0.0`.
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
     - NOT DONE: no MAPL3 component subscribes to `TRADV` yet, and the
       placement question above is open.

5. **Initialize.**
   - Get grid, levels and coordinates from MAPL (Solar step 18: never ask
     ESMF for the grid): `MAPL_GridCompGet(gc, geom=, grid=, num_levels=,
     logger=)`, `MAPL_GridGetCoordinates` for `LONS`/`LATS`,
     `MAPL_GridCompGetInternalState`.
   - Timestep from the clock; `GOCART2G_DT` and `DIURNAL_BIOMASS_BURNING`
     via `MAPL_GridCompGetResource`.
   - Field and state attributes become `ESMF_InfoSet`:
     - `ScavengingFractionPerKm` on DMS/SO2/SO4/MSA
     - the `SetofHenryLawCts` vector, as `values=`
     - on the AERO state: `use_photolysis_table`, `mie_table_instance`,
       `band_for_aerosol_optics`, `wavelength_for_aerosol_optics`,
       `mieTable_pointer`, `internal_variable_name`
   - Replace `add_aero_named_alias` (develop's old-API helper) with the
     MAPL3 `add_aero(..., geom=, km=, ...)`, including the SU-only
     `SAREA` and `REFF` fields used by chemistry.
   - Mie tables: radiation bands, photolysis and diagnostics, as DU does.
   - `createReplenishAlarm` is already on MAPL3; keep the call.

6. **Run routines (one routine per commit, as DU was done).**
   - Fetch fields through the generated `#include "SU2G_GetPointer___.h"`;
     replace any other `MAPL_GetPointer` with `MAPL_StateGetPointer`.
   - **Undefined-value masks (must fix):** about 10 lines like
     `where(1.01*x > MAPL_UNDEF) x = 0.` in Run1 (around L930-L964).
     MAPL2's undefined value was `+huge`; MAPL3's is `-huge`
     (`InternalConstants.F90` L18-19). Ported as-is the condition is
     always true and every emission is zeroed. Rewrite the test (e.g.
     `abs(x) >= abs(MAPL_UNDEF)/1.01`) and run-test it.
     - DU2G (`where (1.01*du_src > MAPL_UNDEFINED_REAL)`, around L648 on
       `release/MAPL-v3`) has the same pattern. Check whether DU's
       emissions are being zeroed too.
   - **Edge-level fields** (`PLE`, `ZLE`): their vertical index starts at
     1 in MAPL3 (0 in MAPL2). Remap code that assumes `0:LM`
     (`p1d => PLE; PLE(0:LM) => p1d`), or confirm it only uses `lbound`.
   - **Oxidants:** check whether `SulfateUpdateOxidants` /
     `SulfateChemDriver` write back into the `SU_H2O2` import (the
     `recycle_H2O2` path). MAPL3 may not let a component write into an
     import.
   - Timers -> `MAPL_GridCompTimerStart/Stop`. Threading keeps
     `MAPL_get_current_thread`.
   - File-name templating -> `MAPL_StrTemplate`. Date/time packing is
     already ported.

7. **AERO methods.** `aerosol_optics`, `monochromatic_aerosol_optics` and
   `get_mixR` follow DU's MAPL3 versions (`ESMF_InfoGet`, field lookup
   through Info keys). The current copies use develop's long-name lookup,
   which needs reverting to the Info-key style.

8. **Wrapper and parent wiring.**
   - Add an external `SU2G_SetServices` wrapper after `end module`, as in
     DU and Solar step 13.
   - In `../GOCART2G_GridCompMod.F90`:
     - uncomment `use SU2G_GridCompMod`, `add_children__(gc, self%SU, ...)`
       and `add_aero_states_(self%SU%instances(:))`
     - add connections for SU's export diagnostics and
       `<name>_AERO`/`_AERO_DP`, following the SS/DU blocks
     - re-export `PSO4` (replacing the commented
       `MAPL_AddExportSpec(CHILD_ID=...)` around L230-233)
     - restore the SU parts of the Run2 totals (the commented block around
       L795-850, including `pso4tot`)
   - Uncomment `SU2G_GridComp` in `../CMakeLists.txt`.

9. **Run-testing.**
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

## Risks to settle early
1. The undefined-value comparisons (step 6). This also affects DU on the
   current branch.
2. GMI coupling: `using_GMI` changes which fields SU imports. Test with
   GMI off first, and confirm the baseline used the same setting.
3. Whether SU writes back into the `SU_H2O2` import (step 6).
4. Volcano and point-emission files are read from disk, so the regression
   data must include them.
5. Tracer transport (step 4, friendliness): until SU subscribes to the
   `TRADV` service (and turbulence/moist equivalents exist), DMS, SO2,
   SO4 and MSA are not advected, diffused or convected in a coupled run.
   The SU-only standalone regression is unaffected. DU and SS have the
   same gap.
