# GOCART: previous architectural decisions

**Audience:** coding agents and developers working on the GOCART repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GOCART` @ `96b3c9a` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`; repository `GOCART`, mounted at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GOCART`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Each entry records a pattern seen in the code, the likely reason, and where useful a question. The patterns are verified; the reasons are the author's reading, not confirmed by maintainers (author decision, for review). The entries below were written for the MAPL2 line; each carries a "MAPL3 line" note where the code on `release/MAPL-v3` differs. The last entries are specific to the MAPL3 line. Overview: [ARCHITECTURE.md](ARCHITECTURE.md). Component detail: [GOCART2G](../ESMF/GOCART2G_GridComp/ARCHITECTURE.md).

### One MAPL child per aerosol instance
- Pattern: the parent creates one child component per instance name listed in `ACTIVE_INSTANCES_*` and `PASSIVE_INSTANCES_*`, so several carbonaceous instances (`CA.oc`, `CA.bc`, `CA.br`) are separate children with their own rc files (`ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `getInstances_`, `createInstances_`; `ESMF/GOCART2G_GridComp/GOCART2G_GridComp.rc`). [Verified]
- MAPL3 line: the pattern holds; instance lists are read in `setup_constituents_` and children are made in `create_instances_`; each child is created with `MAPL_GridCompAddChild(gc, child_name, setservices, hconfig)`, where `hconfig` is read from `<species>_instance_<name>.yaml`; instance lists come from `ACTIVE_INSTANCES` and `PASSIVE_INSTANCES` in `GOCART2G_GridComp.yaml` (`GOCART2G_GridCompMod.F90`, `setup_constituents_`, `add_children__`). [Verified]
- Likely reason: a species can be configured, run, restarted and made data-driven independently, by name. [Inferred]
- Question for maintainers: is the single-NI-instance limit (nitrate reads the first DU, SS and SU instance only) intended to be permanent?

### Data-driven instances selected by name
- Pattern: an instance whose name contains `data` registers only a `Run` phase that copies climatologies from ExtData; a non-data instance registers three (`Run`, `Run2`, `Run0`) (`ESMF/Shared/Chem_AeroGeneric.F90`, `determine_data_driven`; each child's `SetServices`; parent `Run2`, `index(child_name, 'data')`). [Verified]
- MAPL3 line: the name test (`determine_data_driven`) remains, but the data-driven branch of each ported child's `SetServices` calls `_FAIL("... section has not been activated yet")` (DU2G says "data driver"); a non-data child registers phases `Run1`, `Run2`, `Run0` by name, and the parent `run2` skips names containing `data` (`DU2G_GridCompMod.F90`, `SS2G_GridCompMod.F90`, `SU2G_GridCompMod.F90`). [Verified]
- Likely reason: the same AERO interface can be supplied by a prescribed-aerosol configuration without changing the consumers. [Inferred]
- Question for maintainers: is selection by substring of the instance name the intended long-term mechanism?

### `AERO` is a state of methods, not only data
- Pattern: consumers obtain optics, mixing ratio sums and activation properties by calling ESMF methods attached to the `AERO` state (`ESMF_MethodAdd` in `GOCART2G_GridCompMod.F90`, `Initialize`; `ESMF/Shared/GOCART2G_AerosolActivationMod.F90`, `setup_aerosol_activation`). [Verified]
- Likely reason: Radiation, Moist, GAAS and GMI can use aerosol properties without linking to or knowing the GOCART2G species list. [Inferred]
- MAPL3 line: unchanged in kind (`ESMF_MethodAdd` in `Initialize`; `setup_aerosol_activation`). [Verified]
- Question for maintainers: is the method interface the contract other aerosol providers (MAM, the `none` option) must meet?

### Chemistry chooses the aerosol provider by name
- Pattern (MAPL2 code in `../GEOS_ChemGridComp.F90`, not compiled on the MAPL3 line): Chemistry promotes `AERO` and `AERO_DP` from the child named by `AERO_PROVIDER`; for `none` it creates empty states and marks the optics and activation methods as not implemented (`../GEOS_ChemGridComp.F90`, `GetProvider_`, `SetServices`). [Verified]
- Likely reason: aerosol providers can be swapped without changing Physics connections. [Inferred]
- Question for maintainers: see the `AERO_PROVIDER` default observation in [ARCHITECTURE.md](ARCHITECTURE.md), section 9.

### Shared environment type for all species
- Pattern: every child's derived type extends `GA_Environment`, which holds Mie tables, per-bin parameters, scheme names and wavelength lists, and reads them with one `load_from_config` (`ESMF/GOCART2G_GridComp/GA_Environment/GA_EnvironmentMod.F90`). [Verified]
- Likely reason: common configuration and optics set-up are written once for five species. [Inferred]

### Process code in a separate library
- Pattern: emission, settling, removal, chemistry and diagnostics routines are in `Process_Library` and take plain arrays and parameters; the child modules gather pointers and call them (`Process_Library/GOCART2G_Process.F90`; `Process_Library/CMakeLists.txt`; `DU2G_GridCompMod.F90`, `run1`). [Verified]
- Likely reason: the science code can be reused outside the MAPL component (`@GMI` calls `ReadPointEmissions`; the UFS cap reuses the component; `@GMI` and `ESMF/UFS` are not built on this line) and tested apart from the framework. [Inferred]
- Question for maintainers: `GOCART2G_Process.F90`'s header says "only uses intrinsic Fortran functions" while the library depends on NetCDF Fortran; is the statement current?

### State specifications generated by ACG
- Pattern: each component's imports, exports and internal fields are listed in a `*_StateSpecs.rc` that `mapl_acg` turns into generated headers; hand-written `MAPL_GridCompAddSpec` calls (MAPL3; `MAPL_AddExportSpec` on the MAPL2 line) add what the rc cannot, such as `AERO`, `AERO_DP`, `<name>_AERO` and `TRADV` (`ESMF/GOCART2G_GridComp/*/CMakeLists.txt`; `GOCART2G_StateSpecs.rc`; `GOCART2G_GridCompMod.F90`, `SetServices`). [Verified]
- Likely reason: uniform field declarations across many species and bins with less repeated code. [Inferred]

### Parent runs children in two Chemistry phases plus a pre-alarm phase
- Pattern: Parent: two run phases (`Run1`, `Run2`); each non-data child: three. On the MAPL3 line the phases are registered by name (`phase_name="Run1"`, `"Run2"`, `"Run0"`) and the parent runs them with `MAPL_GridCompRunChild`: `run1` runs every child's `Run1` (emissions); `run2` runs each non-data child's `Run0` (klid zeroing) and then its `Run2` (removal, chemistry, diagnostics) (`GOCART2G_GridCompMod.F90`, `SetServices`, `run1`, `run2`; `DU2G_GridCompMod.F90`, `SetServices`). The MAPL2 code also had a phase-number protocol driven by the uncompiled `../GEOS_ChemGridComp.F90`, `Run1`, `Run2`. [Verified]
- MAPL3 line: the parent has no alarm check in `run1` or `run2`; `run2` runs `Run0` then `Run2` of each non-data child every call (`GOCART2G_GridCompMod.F90`). [Verified]
- Likely reason: emissions are needed before turbulence and the rest of the step is needed before radiation; the klid zeroing keeps fields above the pressure lid at zero whether or not the alarm rings. [Inferred]

### Replenish alarm for prescribed oxidants
- Pattern (MAPL2 line): SU2G and NI2G create an alarm with `createReplenishAlarm(gc, clock, 30000, ...)` whose consistency with the run-alarm convention is asserted against `RUN_AT_INTERVAL_START:` (`ESMF/Shared/ReplenishAlarm.F90`; `SU2G_GridCompMod.F90`, `NI2G_GridCompMod.F90`, `Initialize`). [Verified]
- MAPL3 line: `ESMF/Shared/ReplenishAlarm.F90` reads `RUN_AT_INTERVAL_START` with `MAPL_GridCompGetResource` only to offset `ringTime` (`if (.not.run_at_interval_start) ringTime = ringTime-TSTEP`); there is no assert. Only SU2G creates the alarm on this line (`SU2G_GridCompMod.F90`, `Initialize`); NI2G is not built. [Verified]
- Likely reason: `H2O2_INIT` and `XHNO3` are periodically reset to their prescribed values. [Inferred]
- Question for maintainers: is the three-hour period a tuned value or a convention?

### Two component families kept side by side
- Pattern: legacy `GOCART_GridComp` (`Chem_Registry`-driven trace gases and aerosol) sits beside `GOCART2G_GridComp`; on the MAPL3 line it is commented out of the build as "not yet ported to MAPL3" (`ESMF/CMakeLists.txt`), and the Chemistry file that adds it is not compiled (`../CMakeLists.txt`). On the MAPL2 line it was built and added by Chemistry. [Verified]
- Likely reason: trace-gas children (CO, CO2, CH4, O3, CFC, Rn) have not been moved to the newer pattern. [Inferred]
- Question for maintainers: which parts of `GOCART_GridComp` are still used by supported configurations?

### Framework variants held as placeholders
- Pattern: `CCPP/GOCART_CCPP.F90` and `ESMF/NUOPC/NUOPC_Cap.F90` are tracked but empty; the UFS NUOPC cap is in `ESMF/UFS` and built only with `UFS_GOCART` on the MAPL2 line; on the MAPL3 line `UFS` is commented out of `ESMF/CMakeLists.txt` (`README.md`). [Verified]
- Likely reason: directories reserved for further framework interfaces. [Inferred]

### Unported species are ignored, not failed (MAPL3 line)
- Pattern: the parent's `add_children__` calls for CA and NI are commented out; `deactivate_instances__` marks every CA and NI instance inactive, logs a warning when any was listed active, and the source comment says anything that checks `is_active` "must not see them as active" (`ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `create_instances_`, `deactivate_instances__`). [Verified]
- Likely reason: the shipped `GOCART2G_GridComp.yaml` lists CA and NI instances as active while those children are not ported, so the parent must still run. [Inferred]
- Question for maintainers: when the CA and NI children are ported, are `deactivate_instances__` and the commented calls to be removed together?

### Advection by service subscription (MAPL3 line)
- Pattern: each ported child lists its advected tracers in an import of item type `MAPL_STATEITEM_SERVICE` named `TRADV`; a source comment says it "Replaces MAPL2's FRIENDLYTO=DYNAMICS" (`DU2G_GridCompMod.F90`, `SS2G_GridCompMod.F90`, `SU2G_GridCompMod.F90`, `SetServices`). [Verified]
- Likely reason: in MAPL3 a tracer is declared to its consumer by a service subscription instead of a friendly-to list. [Inferred]
- Question for maintainers: `FRIENDLYTO` values for `TURBULENCE` and `MOIST` remain in the `*_StateSpecs.rc` files; do they have an effect on this line?
