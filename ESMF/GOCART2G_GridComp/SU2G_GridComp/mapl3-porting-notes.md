# SU2G_GridComp MAPL3 porting plan

Status: in progress on `feature/pchakrab/port-su2g-to-mapl3`. Steps 1, 1a, 1b, 1c and 2 done.

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

3. **Instance configuration.**
   - Convert `SU2G_instance_SU.rc` to `SU2G_instance_SU.yaml` (plus a
     `.data` variant), following `DU2G_instance_DU.yaml`.
   - Generic settings already go through
     `GA_Environment%load_from_config(gc)`, which is on MAPL3.
   - SU-specific settings (`using_GMI`, `disable_emissions`, the volcano
     files, `aircraft_fuel_emission_factor`, `so4_anthropogenic_fraction`,
     `aviation_vertical_layers(4)`, `point_emissions_srcfilen`) move to
     `MAPL_GridCompGetResource(gc, "label", ...)`, with the trailing `:`
     dropped from labels.

4. **SetServices.**
   - Store SU's state with `_SET_NAMED_PRIVATE_STATE` /
     `_GET_NAMED_PRIVATE_STATE` in place of the `wrap_` +
     `ESMF_UserComp*InternalState` pattern. Keep the per-thread
     `workspaces` (`MAPL_get_num_threads`).
   - Register entry points with `phase_name` (`"Run1"`, `"Run2"`,
     `"Run0"`), as DU does. Drop `MAPL_GenericSetServices`.
   - Convert the manually declared fields to `MAPL_GridCompAddSpec`,
     following DU2G's manual `MAPL_GridCompAddSpec` calls in SetServices:
     - data-driven: `PLE`, `RH2`, `SO4` (internal with add-to-export),
       `climSO4`, and the per-bin `climSU{DP,WT,SD,SV}nnn`
     - `SU_H2O2/OH/NO3` (computational instance only), and
       `GMI_OH/H2O2/NO3` (only `if (using_GMI)`)
     - the `<name>_AERO` state and `<name>_AERO_DP` bundle, as in DU
   - Friendliness (MAPL2 `FRIENDLYTO=DYNAMICS:TURBULENCE:MOIST` on DMS,
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
   - Replace the `MAPL_AM_I_ROOT` prints with the component logger.

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
