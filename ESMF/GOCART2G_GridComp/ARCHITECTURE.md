# GOCART2G_GridComp: architecture

**Audience:** coding agents and developers working on the GOCART2G aerosol gridded component and its children, `Process_Library` and `GA_Environment`.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GOCART` @ `0d19ada` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`; repository `GOCART`, mounted at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GOCART`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths are relative to the GOCART repository root. Up: [repository architecture](../../docs/ARCHITECTURE.md), [agent rules](../../AGENTS.md), [decisions](../../docs/decisions.md), parent repository `GEOSchem_GridComp` (`../../../docs/ARCHITECTURE.md`, written separately), Physics (`../../../../ARCHITECTURE.md` in the fixture).
"Mie table" means a precomputed table of aerosol optical properties as a function of wavelength, relative humidity and size bin. "Klid" is the model level at the pressure lid above which aerosol is zeroed.

## 1. Purpose
GOCART2G simulates dust, sea salt, sulfate, carbonaceous and nitrate aerosol: emission, settling, dry and wet removal, chemistry, and diagnostics. It also builds the aerosol optical, microphysical and mixing-ratio services that Radiation, Moist, Surface, GAAS and GMI consume through the exported `AERO` state (`ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `Initialize`). Each species is a separate child component; the parent only creates them, sums their exports and assembles the combined `AERO` state.

## 2. Place in the hierarchy
- Parent: Chemistry (`GEOS_ChemGridComp.F90` in `GEOSchem_GridComp`) adds `GOCART2G` via `MAPL_AddChild` (see [repository architecture](../../docs/ARCHITECTURE.md), section 2).
- Children: `createInstances_` adds one MAPL child per instance name, in the order DU, SS, CA, SU, NI, active instances before passive ones. The source comment says active-first ordering is needed for filling the `AERO` states (`GOCART2G_GridCompMod.F90`, `createInstances_`).
- Instance names come from `ACTIVE_INSTANCES_{DU,SS,SU,CA,NI}` and `PASSIVE_INSTANCES_*` in `ESMF/GOCART2G_GridComp/GOCART2G_GridComp.rc` (`getInstances_`). Shipped values: one instance each for DU, SS, SU, NI; CA has `CA.oc CA.bc CA.br`; passive lists are empty.
- A name containing `data` makes a data-driven instance (`ESMF/Shared/Chem_AeroGeneric.F90`, via `determine_data_driven` in each child's `SetServices`). Data instances read climatologies instead of computing. The parent comment says "Default is the data component", but the shipped rc lists computational names (`GOCART2G_GridCompMod.F90`, `getInstances_`; `GOCART2G_GridComp.rc`). The fixture's `gcm_setup` rewrites the names to `*.data` when `DATA_DRIVEN` is true (`Applications/@GEOSgcm_App/gcm_setup`, outside this repository).

| Child | Directory | Species | Bins (`nbins` in rc) | Internal fields (state spec) |
|---|---|---|---|---|
| DU2G | `ESMF/GOCART2G_GridComp/DU2G_GridComp` | Dust | 5 | `DU` (`DU2G_StateSpecs.rc`) |
| SS2G | `.../SS2G_GridComp` | Sea salt | 5 | `SS`, `DEEP_LAKES_MASK` |
| SU2G | `.../SU2G_GridComp` | Sulfate and precursors | 4 | `DMS`, `SO2`, `SO4`, `MSA`, `H2O2_INIT` |
| CA2G | `.../CA2G_GridComp` | Carbonaceous (instances `CA.oc`, `CA.bc`, `CA.br`) | 2 (hydrophobic, hydrophilic) | `<name>phobic`, `<name>philic` |
| NI2G | `.../NI2G_GridComp` | Nitrate and ammonium | 5 | `NH3`, `NH4a`, `NO3an1`, `NO3an2`, `NO3an3`, `XHNO3` |

Bin counts: each `*_instance_*.rc`, key `nbins`. Internal specs: the `*_StateSpecs.rc` in each child directory (processed by `mapl_acg`).

## 3. Entry points
### Parent (`GOCART2G_GridCompMod.F90`)
- `SetServices`: registers Initialize and two run phases, `Run1` and `Run2` (`MAPL_GridCompSetEntryPoint`); Parent: two run phases (`Run1`, `Run2`); each non-data child: three (`Run`, `Run2`, `Run0`); each data child: one (`Run`); no Finalize is registered. Reads `GOCART2G_GridComp.rc` into the component configuration, copies the wavelength lists into the universal configuration with `MAPL_ConfigSetAttribute`, reads `use_threads:` (default false, calls `MAPL%set_use_threads`), builds child instances, forwards `RUN_AT_INTERVAL_START:` to children (`alarmResourcesToChildren`), declares exports `AERO` (state) and `AERO_DP` (bundle), and adds connectivities among children (section 4).
- `Initialize`: `MAPL_GenericInitialize`, then loads `AGCM.rc`; `add_aero_states_` adds each active child's `<name>_AERO` state into the parent `AERO` state and each child's `<name>_AERO_DP` fields into `AERO_DP` (except NI); adds named aliases (`PLE`, `RH`, `EXT`, `SSA`, `ASY`, optionally `MOM`, `monochromatic_EXT`, and the `sum_of_internalState_aerosol_*` fields); sets state attributes (`band_for_aerosol_optics`, `wavelength_for_aerosol_optics`, `use_photolysis_table`, `n_phase_function_moments`, `implements_aerosol_optics_method`, among others); attaches ESMF methods (section 5); calls `setup_aerosol_activation`.
- `Run1`: returns if the run alarm is not ringing; otherwise runs phase 1 of every child (`ESMF_GridCompRun`, `phase=1`).
- `Run2`: first runs phase 3 (klid zeroing) of every non-data child, before the alarm check; then returns if the alarm is not ringing; zeroes the total exports; runs phase 2 of non-data children; then sums per-instance exports into `TOT*`, `PM*`, `PSO4TOT`, `TOTANGSTR`, `TOTABCKTOA`, `TOTABCKSFC` (`Run2`). `TOTANGSTR` requires 550 nm in the vertically integrated wavelength list, else it returns error code 100 (`VERIFY_(100)`). The backscatter exports require 532 nm and the requested exports `TOTEXTCOEF` and `TOTBCKCOEF`.
- Chemistry calls phase 1 for components with more than one run phase and phase 2 afterwards, so GOCART2G `Run1` and `Run2` map onto Physics' two Chemistry calls (`../GEOS_ChemGridComp.F90`, `Run1`, `Run2`).

### Child (same shape in all five)
Each child registers, in this order: Initialize, `Run` (phase 1), and, only if the instance is not data-driven, `Run2` (phase 2) and `Run0` (phase 3) (for example `DU2G_GridComp/DU2G_GridCompMod.F90`, `SetServices`).
- `Run` is a dispatcher: it calls `Run1` for a computational instance and `Run_data` for a data instance.
- `Run1`: emissions (and for NI, ammonia emissions only).
- `Run2`: settling, dry deposition, wet removal, chemistry (SU, NI, CA) and diagnostics.
- `Run0`: zeroes internal fields above the pressure lid (`Chem_AeroGeneric.F90`, `setZeroKlid`, `setZeroKlid4d`, `findKlid`).
- `Initialize`: loads `<X>2G_instance_<name>.rc`, reads `RUN_DT:` and `GOCART2G_DT:` (default `RUN_DT`) to get the chemistry time step, loads the Mie tables, builds the `<name>_AERO` state and `<name>_AERO_DP` bundle, attaches `aerosol_optics`, `monochromatic_aerosol_optics` and `get_mixR` methods and the `mieTable_pointer` attribute.
- Instance rc fallback when the file for a name is missing: DU, SS, NI, SU fall back to `<X>2G_instance_<X>.rc`; CA uses `CA2G_instance_<first five characters of the name>.rc` (for example `CA.oc`) (each child, `SetServices`).

## 4. State
- Parent imports `DELP`, `AIRDENS`, `T`, `PLE`; exports `PSO4TOT`, `TOT*` and `PM*` fields; no internal state (`ESMF/GOCART2G_GridComp/GOCART2G_StateSpecs.rc`). Generated headers `GOCART2G_Export___.h`, `GOCART2G_Import___.h` are included in the module.
- Additional parent exports: `AERO`, `AERO_DP` (added in `SetServices`); `PSO4` promoted from the first SU instance unless its name contains `data` (`MAPL_AddExportSpec` with child id).
- Child state: ACG `*_StateSpecs.rc` for imports, exports and internal fields; hand-written additions in each `SetServices` (data-instance branches, `MAPL_AddExportSpec`). Internal fields are `MAPL_RestartOptional` and, for most species, `friendlyto` `DYNAMICS:TURBULENCE:MOIST` so Dynamics advects them and Turbulence and Moist treat them as tracers. In the DU data-driven branch the `friendlyto` setting is commented out (`DU2G_GridCompMod.F90`, `SetServices`). Fields `DEEP_LAKES_MASK` (SS), `H2O2_INIT` (SU) and `XHNO3` (NI) are not exported.
- Restart file names are set in the fixture's `AGCM.rc.tmpl` (`@MODELATM{DU,SS,CA.oc,CA.bc,CA.br,NI,SU}_INTERNAL_RESTART_FILE`).
- Couplings among siblings (parent `SetServices`, `MAPL_AddConnectivity`): when any NI instance exists, NI instance 1 receives `DU` from the first DU instance, `SS` from SS instance 1, and `SO4` from SU instance 1, when those are active. These are same-spec sibling connections, so the field is shared between components (no copy). The source comment says NI "only uses first active dust and sea salt instance".
- Couplings with other Chemistry children (imports from `CHEMENV`, `ACHEM`, `HEMCO`, `GMICHEM`; `AERO` to `GAAS`, `GMICHEM`, `CARMA`) are made in `GEOS_ChemGridComp.F90`; table in [repository architecture](../../docs/ARCHITECTURE.md), section 5.3.
- Per-thread workspaces: DU, SU, CA, NI hold `workspaces(0:num_threads-1)` indexed by `MAPL_get_current_thread`; SS has none (each child module, `Initialize`/`Run`). Active when `use_threads` is true.
- Data instances take climatology imports `clim*` through ExtData (the `*_GridComp_ExtData.yaml` in each child) and copy them in `Run_data`.

## 5. The `AERO` state and its methods
Consumers call methods stored on `AERO`; they are not ordinary data. All are module procedures in `GOCART2G_GridCompMod.F90`, attached in `Initialize` unless noted.
| Method | Does | Consumer (outside this repository) |
|---|---|---|
| `run_aerosol_optics` | Sums over the child states the extinction, extinction-weighted single-scattering albedo and asymmetry (and optionally phase function moments) that each child's `aerosol_optics` produced | Radiation (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`); GMI photolysis (`@GMI/GmiPhotolysis_GridCompClassMod.F90`) |
| `get_mixRatioSum` | Sums aerosol mixing ratios over the child states (header comment: sum of aerosol mixing ratio for GAAS) | `@GAAS/GAAS_GridCompMod.F90` |
| `serialize_bundle` | Returns a serialized field bundle of the aerosol fields of the child states (header comment) | `@GAAS/GAAS_GridCompMod.F90` |
| `get_monochromatic_aop` | Monochromatic aerosol optical properties summed over children (wavelengths from `aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:`) | `@GAAS/GAAS_GridCompMod.F90` |
| `aerosol_activation_properties` | Mode-wise number, size, hygroscopicity for cloud-droplet activation; attached by `setup_aerosol_activation` in `ESMF/Shared/GOCART2G_AerosolActivationMod.F90`, not in the parent | Moist (`GEOSmoist_GridComp/aer_actv_single_moment.F90`; see the Moist document (`../../../../GEOSmoist_GridComp/ARCHITECTURE.md` in the fixture)) |
- `run_aerosol_optics` returns sums; the Solar consumer then divides to get single-scattering albedo and asymmetry (`SSAA=SSAA/TAUA`, `ASYA=ASYA/SSAA`) (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, in the aerosol optics block). *Unverified:* that the infrared consumer does the same.
- Activation modes: 14 GOCART modes (`du001`-`du005`, `ss001`-`ss003`, `sulforg01`-`sulforg03`, `bcphilic`, `ocphilic`, `brcphilic`), or 7 MAMnet modes (`NUM_A_*`) when `USE_MAMNET:` is true in `AGCM.rc`; `MAXCLEAN:` (default 1e-9) and `CCNTUNING:` (default 1.8) are read from the same configuration. The routine excludes NI children (`GOCART2G_AerosolActivationMod.F90`, `setup_aerosol_activation`).
- Child-level Mie lookups: `Process_Library/GOCART2G_MieMod.F90` (type `GOCART2G_Mie`, `QueryByWavelength_*`, `QueryByChannel_*`). Table files come from `aerosol_radBands_optics_file:` and `aerosol_monochromatic_optics_file:` in the instance rc; wavelength lists come from the parent rc (`wavelengths_for_profile_aop_in_nm:`, `wavelengths_for_vertically_integrated_aop_in_nm:`, `aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:`, `aerosol_photolysis_wavelength_in_nm_from_LUT:`, `n_phase_function_moments_photolysis:`).

## 6. Configuration
Only keys that change structure or scheme. File in parentheses.
- Instances: `ACTIVE_INSTANCES_*`, `PASSIVE_INSTANCES_*` (`GOCART2G_GridComp.rc`); `use_threads:` (same file); `RUN_AT_INTERVAL_START:` (copied to children; `ReplenishAlarm` checks consistency).
- Time step: `RUN_DT:`, `GOCART2G_DT:` (`AGCM.rc`).
- Common to every child, read by `GA_Environment::load_from_config` (`ESMF/GOCART2G_GridComp/GA_Environment/GA_EnvironmentMod.F90`) from `<X>2G_instance_<name>.rc`: `nbins`, `particle_radius_microns`, `particle_density`, `fscav`, `fwet`, `molecular_weight`, `fnum`, `pressure_lid_in_hPa`, `fwet_ice`, `fwet_snow`, `fwet_rain`, `wet_radius_thr`, `washout_tuning`, `wet_removal_scheme` (`gocart` or `ufs`; selects `WetRemovalGOCART2G` or `WetRemovalUFS`), `settling_scheme` (`gocart` or `ufs`; selects `SettlingSolver` or `SettlingSolverUFS`). `validate_config` asserts valid scheme names and rainout efficiency.
- DU2G: `emission_scheme:` `ginoux` (default), `k14`, `fengsha`, plus scheme-specific keys; `point_emissions_srcfilen:` (point emissions read once per day) (`DU2G_GridCompMod.F90`, `SetServices`).
- SS2G: `emission_scheme` (an integer in the shipped rc), `sstEmisFlag`, `weibullFlag`, `hoppelFlag`, `emission_scale`, `radius_lower`, `radius_upper`.
- SU2G: `using_GMI`, `disable_emissions`, `volcano_srcfilen_degassing`, `volcano_srcfilen_explosive`, `aircraft_fuel_emission_factor`, `so4_anthropogenic_fraction`.
- CA2G: `hydrophobic_fraction`, `time_days_hydrophobic_to_hydrophilic` (default 2.5; negative means no transfer), `pom_ca_ratio`, `monoterpenes_emission_fraction`, `isoprene_emission_fraction`, `my_day_of_week`, point emissions. `Chem_BiomassDiurnal` reads `DIURNAL_BIOMASS_BURNING:` from the universal configuration.
- NI2G: `using_GMI` selects `GMI_HNO3` instead of the ExtData `NITRATE_HNO3`.
- `AMIP/` and `AMIP.20C/` subdirectories of each child hold variants of the rc and yaml files (for example volcano source files in SU2G). The fixture copies `AMIP/` files when `AMIP_EMISSIONS` is set (`Applications/@GEOSgcm_App/gcm_run.j`, outside this repository).
- Resource lookup: `MAPL_GetResource` tries `<component>_<label>`, `<type>_<label>` and `<label>` (`MAPL` `base/MAPL_Resource.F90`, `get_labels_with_prefix`); a child's `p:` prefixed key (such as `p:RUN_AT_INTERVAL_START:`) is a per-child override.

## 7. Dependencies
- MAPL (state, resource, ACG, threading), ESMF, NetCDF Fortran, OpenMP (not with NAG) (`ESMF/GOCART2G_GridComp/CMakeLists.txt`).
- `Process_Library` (`Process_Library/GOCART2G_Process.F90`, about 90 procedures, and `GOCART2G_MieMod.F90`). Process routines called per child: DU `DustEmissionGOCART2G`, `DustEmissionK14`, `DustEmissionFENGSHA`, `Chem_SettlingSimple`, `DryDeposition`, `WetRemovalGOCART2G`, `WetRemovalUFS`, `Aero_Compute_Diags`; SS `jeagleSSTcorrection`, `weibullDistribution`, `SeasaltEmission`; SU `SulfateDistributeEmissions`, `DMSemission`, `SUvolcanicEmissions`, `SulfateUpdateOxidants`, `SulfateChemDriver`, `SU_Wet_Removal`, `SU_Compute_Diags`; CA `CAEmission`, `Chem_BiomassDiurnal`, `phobicTophilic`, `carbonChemLoss`; NI `NIthermo` (uses `RPMARES`), `NIheterogenousChem`.
- `GA_Environment` (`ESMF/GOCART2G_GridComp/GA_Environment/GA_EnvironmentMod.F90`): type `GA_Environment` holds the Mie tables (`rad_Mie`, `diag_Mie`, `phot_Mie`), per-bin parameters, `nbins`, `km`, `CDT`, `klid`, scheme names, wavelength lists. Each child's derived type `extends(GA_Environment)`. `load_from_config` and `validate_config` are the only type-bound procedures.
- `Chem_Shared2G` (`ESMF/Shared`): `Chem_AeroGeneric.F90` (`add_aero`, `add_aero_named_alias`, `append_to_bundle`, `setZeroKlid`, `findKlid`, `get_mixR`, `determine_data_driven`), `ReplenishAlarm.F90` (`createReplenishAlarm`), `GOCART2G_AerosolActivationMod.F90`.
- `createReplenishAlarm(gc, clock, 30000, ...)` sets an alarm that recycles prescribed oxidant fields: SU2G for `H2O2`, NI2G for `HNO3`, every three hours (the argument 30000 is unpacked as hhmmss by `MAPL_UnpackTime`) (`ESMF/Shared/ReplenishAlarm.F90`, `createReplenishAlarm`; `SU2G_GridCompMod.F90`, `NI2G_GridCompMod.F90`, `Initialize`). *Unverified:* that the alarm's only consumer is the oxidant recycling described in the source comments.
- Outside the repository: `GEOS_ChemGridComp.F90` (parent); `AGCM.rc` (read in parent `Initialize` and for `USE_MAMNET`); ExtData yaml files (listed in `Applications/@GEOSgcm_App/construct_extdata_yaml_list.py`).
- `Process_Library` is also linked by `@GMI` (`ReadPointEmissions`).

## 8. Known traps
- **Run phase numbers.** Parent: two run phases (`Run1`, `Run2`); each non-data child: three (`Run`, `Run2`, `Run0`); each data child: one (`Run`). The parent runs child phase 3 at the start of its own `Run2`, before the alarm check, so klid zeroing happens at every Chemistry phase 2 call, not only when the alarm rings.
- **Data instances have one phase.** A data instance registers only `Run`. Code that loops over children with phase 2 or 3 must skip them (the parent tests `index(child_name, 'data')`).
- **Single-instance assumptions.** The parent dereferences `instances(1)` for DU, SS and SU when connecting to NI and warns when NI has more than one active instance (`getInstances_`). The dust, sea-salt and sulfate fields reaching NI are those of the first instance only.
- **Order of instances.** Active before passive, and DU, SS, CA, SU, NI; `AERO` and `AERO_DP` contents follow this order (`createInstances_`, `add_aero_states_`).
- **`AERO_DP` bin counts** must match the Surface component's expectation or Surface stops (see [repository architecture](../../docs/ARCHITECTURE.md), section 5.6). The SU2G data instance appends only bin `003` (sulfate) to `AERO_DP` (`SU2G_GridCompMod.F90`, `append_to_bundle`).
- **Fields in `AERO` are allocated by the child.** `add_aero` allocates the `AERO` fields for a child with `MAPL_FieldAllocCommit` and adds them to the child's `<name>_AERO` state (`ESMF/Shared/Chem_AeroGeneric.F90`, `add_aero`); the parent then adds each child state into its own `AERO` state (`GOCART2G_GridCompMod.F90`, `add_aero_states_`). Consumers reach species data by iterating the nested child states, not by a flat field list.
- **Wavelength lists.** `TOTANGSTR` needs 550 nm, backscatter exports need 532 nm in the lists in `GOCART2G_GridComp.rc`; removing them fails `Run2`.
- **`nifactor` in the parent.** Sulfate particulate-matter sums in parent `Run2` use `nifactor` (132.14/96.06) unless the first NI instance is active and not data-driven, then 1.0 (`GOCART2G_GridCompMod.F90`, `Run2`). Changing the NI instances changes the `PM*` totals.
- **Resource variants.** Files under `AMIP/`, `AMIP.20C/` and the base directory define the same keys with different values; the run directory receives one set (fixture scripts).
- **Threads.** Workspaces are per thread (`workspaces(0:num_threads-1)`); *Inferred:* shared module-level state added to a child would be unsafe when `use_threads` is true.

## 9. Decisions
Recorded in [decisions.md](../../docs/decisions.md): one child per instance; `AERO` as a state of methods; shared environment type `GA_Environment`; `Process_Library` separated from MAPL state; hand-ordered coupling of NI to sibling instances; data instances by name; dual run-phase protocol; ACG for state specs. All are author decisions for review.

## 10. Observations (noted, not fixed) and open questions
- `pi = 3.141529265` in `GOCART2G_GridCompMod.F90` (`Run2`, backscatter) and in `Process_Library/GOCART2G_Process.F90` (`DustAerosolDistributionMeng2022`, `Chem_CalcVsettle`, `Chem_CalcVsettle2Gorig`). The standard value is 3.14159265. Question for maintainers: is the value intentional?
- `SS2G_GridCompMod.F90` (`Run1`, emission section): the call to `hoppelCorrection` is commented out and `fhoppel` stays 1.0, although `hoppelFlag` is read in `load_from_config`. Question for maintainers: is `hoppelFlag` intended to have no effect?
- `GOCART2G_GridCompMod.F90` (`Run2`): the sulfate contribution to `TOTEXTT25` and `TOTEXTTFM` uses `suexttau`, the same field as for `TOTEXTTAU`; no size-specific sulfate field is used. Question for maintainers: is this intended?
- `DU2G_GridCompMod.F90` (`SetServices`, `Initialize`): a message says "Loading DU2G_GridComp_DU.rc instead", while the file loaded is `DU2G_instance_DU.rc`.
- `GOCART2G_GridCompMod.F90` (`SetServices`): comment "Default is the data component" disagrees with the shipped `GOCART2G_GridComp.rc`.
- `Process_Library/GOCART2G_Process.F90` header says the library uses only intrinsic Fortran functions; `GOCART2G_MieMod.F90` uses `netcdf`. Question for maintainers: is the header current?
- `GOCART2G_GridCompMod.F90` (`SetServices`): `instances(1)` of DU, SS, SU is accessed for the NI connectivity; no guard against an empty list was seen. Question for maintainers: is NI with an empty DU, SS or SU list supported?
- `GEOSsurface_GridComp/Shared/StieglitzSnow.F90` (GEOSgcm_GridComp repository) comments refer to the legacy path `GOCART_GridComp/GOCART_GridCompMod.F90` for the bin counts it must match (outside this repository). Question for maintainers: what is the intended reference for the bin counts that Surface checks?
- Question for maintainers: is the single-NI-instance limit permanent?
