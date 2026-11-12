# GOCART architecture

A map of the repository: what lives where, how the pieces fit, and what must not break. Paths are relative to the repository root. It names modules and routines rather than linking to lines; use them to search the code. It describes only things that change rarely; what is ported and built right now, rules for changes, known traps and to-dos are in [`AGENTS.md`](../AGENTS.md).

## Bird's-eye view

GOCART is the aerosol model of GEOS: emission, settling, dry and wet removal, chemistry and diagnostics for dust, sea salt, sulfate, nitrate and carbonaceous aerosol, plus the aerosol optical and microphysical properties that Radiation and Moist need. In GEOS it is a child of Chemistry (`GEOS_ChemGridComp.F90`, in `GEOSchem_GridComp`); it also builds standalone (`components.yaml`, `CMakeLists.txt`).

**On the MAPL3 line, GOCART2G (`ESMF/GOCART2G_GridComp`) is the only GOCART available.**

## Codemap

- `ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`: the parent component. Creates one child per aerosol instance, connects the children, runs their phases, sums their exports into totals (`TOT*`, `PM*`), and assembles the exported `AERO` state and `AERO_DP` bundle. Its specs: `GOCART2G_StateSpecs.rc`; its configuration: `GOCART2G_GridComp.yaml` (instance lists, wavelength lists).
- `ESMF/GOCART2G_GridComp/<X>2G_GridComp/`: one directory per species: `DU2G` (dust), `SS2G` (sea salt), `SU2G` (sulfate and precursors), `NI2G` (nitrate and ammonium), `CA2G` (carbonaceous). Each holds the component module, `<X>2G_StateSpecs.rc` and per-instance resource files `<X>2G_instance_<name>.yaml` (with `.rc` counterparts).
- `ESMF/GOCART2G_GridComp/GA_Environment/`: `GA_Environment`, the base type of every child's private state: Mie tables, bin parameters, scheme choices, read from the instance configuration (`load_from_config`).
- `ESMF/Shared/` (`Chem_Shared2G`): helpers shared by the children: `Chem_AeroGeneric.F90` (the per-child `AERO` states and bundles, zeroing above the pressure lid, data-instance detection), `GOCART2G_AerosolActivationMod.F90` (the aerosol-activation callback for Moist), `ReplenishAlarm.F90`.
- `Process_Library/`: the aerosol physics (`GOCART2G_Process.F90`: emission, settling, deposition, wet removal, chemistry, diagnostics) and Mie table lookups (`GOCART2G_MieMod.F90`).
- `ESMF/GOCART2G_GridComp/regression/`: standalone regression cases, and `DuSsStub_GridComp`, a stand-in used by `nitrate-sa`.
- `ESMF/GOCART_GridComp` (legacy GOCART), `ESMF/Aerosol_GridComp`, `ESMF/Apps` (offline aerosol optics), `ESMF/UFS` (UFS cap), `CCPP`, `ESMF/NUOPC`: other GOCART code, not part of GOCART2G.

## How GOCART2G is put together

- **One child per instance.** `ACTIVE_INSTANCES` / `PASSIVE_INSTANCES` in `GOCART2G_GridComp.yaml` name the instances (for example `DU`, `SU`); the parent adds each as a child with its instance YAML as `hconfig` (`create_instances_`, `add_children__`), in the order DU, SS, SU, NI.
- **Run phases by name.** The parent registers `Run1` (emissions) and `Run2`; its `run2` runs every child's `Run0` (zero above the pressure lid) and then its `Run2` (settling, deposition, removal, chemistry, diagnostics), then sums the exports.
- **Children talk through connections made by the parent:** each child's diagnostics and `<name>_AERO` state, and the `<name>_AERO_DP` bundle of DU, SS and SU (NI has none), go to the parent.

## Contracts with the rest of GEOS

### Imports

Meteorology and emissions, from ExtData and sibling components.

### Exports

- **`AERO` state:** nests one state per instance (`<name>_AERO`, created by `add_aero` in `ESMF/Shared/Chem_AeroGeneric.F90` and added by `add_aero_states_`); consumers reach species data by iterating those states. It also carries callbacks that Radiation, GAAS and Moist run by label to get aerosol optics, mixing ratios and activation properties (Appendix: `AERO` callbacks).
- **`AERO_DP` bundle:** an ordinary field bundle, without callbacks: the deposition, settling and wet-removal fields per species and bin of DU, SS and SU (NI contributes none), copied in by `add_aero_states_` and read by Surface, which expects fixed counts (DU 5, SS 5, SU 1, BC 2, OC 2; `NUM_*` in `StieglitzSnow.F90`). SU has 4 bins but contributes only sulfate (bin `003`).
- **Totals and `PSO4`:** `TOT*`, `PM*` exports; `PSO4` re-exported from the first SU instance for other Chemistry children.

### Services

GOCART2G provides no service; its children subscribe. Each computational instance is meant to subscribe its transported tracers (the `service_items`: DU `DU`; SS `SS`; SU `DMS`, `SO2`, `SO4`, `MSA`; NI `NH3`, `NH4a`, `NO3an1`-`NO3an3`) to three services, which replace MAPL2's `FRIENDLYTO=DYNAMICS:TURBULENCE:MOIST`: advection (`TRADV`, provided by FV3), convective transport and scavenging (`MTR`, provided by Moist's GF), and turbulent diffusion (`TRDIFF`, provided by Turbulence). The providers update the subscribed fields in place.

Scavenging parameters travel with the fields: the children set `ScavengingFractionPerKm` (and, for the gases of SU and NI, the Henry's-law constants `SetofHenryLawCts`) in the ESMF Info of their tracer fields; Moist reads them when it convects the tracers (`CNV_Tracers_Init` in GEOSmoist_GridComp's `Process_Library.F90`).

## Invariants

- `Process_Library` depends only on NetCDF and its own Mie module: no MAPL or ESMF. Physics goes there; state handling stays in the components.
- Every child's private state extends `GA_Environment`, and the configuration common to all species is read in one place (`load_from_config`).
- The order of instances is the order of the child states in `AERO` and of the fields in `AERO_DP` (which has DU, SS and SU only).
- At most one NI instance is active, and it imports `DU`, `SS` and `SO4` from the first DU, SS and SU instances.

## Cross-cutting concerns

- **Configuration:** structure comes from `GOCART2G_GridComp.yaml`, per-instance settings from `<X>2G_instance_<name>.yaml`. Each resource `.yaml` has a `.rc` counterpart with the same settings, and `AMIP/` and `AMIP.20C/` hold variants chosen outside this repository.
- **Generated specs:** imports, exports and internals are declared in the `*_StateSpecs.rc` files and turned into `*_Import___.h`, `*_Export___.h`, `*_Internal___.h` (and the pointer helpers `*_GetPointer___.h`, `*_DeclarePointer___.h`) by MAPL's ACG at build time; only state that the table cannot express is added in `SetServices`.
- **Computational and data-driven instances:** an instance whose name contains `data` is data-driven. The distinction runs through every child's `SetServices` and `run_data` (`Run_data` in SU2G and NI2G), the parent's loops (which skip data instances) and `AERO_DP`.
- **Threading:** each child keeps per-thread workspaces (`workspaces(0:num_threads-1)`); module-level scratch data would be shared between threads.
- **Undefined values:** DU2G, SU2G and NI2G compare inputs against `MAPL2_UNDEF = 1.0e15`, MAPL2's value, not MAPL3's `MAPL_UNDEF`.
- **Testing:** each regression case (`ESMF/GOCART2G_GridComp/regression/`) runs one run phase (`Run1` or `Run2`) of GOCART2G or of a single child for one step, with every import supplied by a data component restarted from a MAPL2 capture, and compares the checkpoints exactly against a stored MAPL3 baseline that was validated against the MAPL2 capture.

## Where to go next

- [`AGENTS.md`](../AGENTS.md): rules, status on this line, to-dos, known traps, build and test.
- `ESMF/GOCART2G_GridComp/SU2G_GridComp/mapl3-porting-notes.md`, `ESMF/GOCART2G_GridComp/NI2G_GridComp/mapl3-porting-notes.md`: how the port was done and what is left.

## Appendix: `AERO` callbacks

GOCART2G exposes some computations as callbacks on `AERO` and on its nested `<name>_AERO` states rather than as fields: a callback is registered on a state with `ESMF_MethodAdd` and run by label with `ESMF_MethodExecute`. Inputs and results go through the state (its fields and its ESMF Info keys), not through the call, and the callback's own errors come back in `userRC`.

The parent attaches these callbacks to `AERO` in its `Initialize` unless noted. A consumer outside this repository depends on each one:

| Callback (label) | Does | Consumer |
|---|---|---|
| `run_aerosol_optics` | Sums extinction, extinction-weighted single-scattering albedo and asymmetry (optionally phase-function moments) over the child states | Radiation: `GEOS_SolarGridComp.F90`, `GEOS_IrradGridComp.F90` (Solar divides the sums to get albedo and asymmetry) |
| `get_mixRatioSum`, `serialize_bundle`, `get_monochromatic_aop` | Mixing-ratio sum, serialized aerosol bundle, monochromatic optical properties | GAAS (`@GAAS/GAAS_GridCompMod.F90`) |
| `aerosol_activation_properties` | Mode-wise number, size and hygroscopicity for cloud-droplet activation; attached by `setup_aerosol_activation` (`ESMF/Shared/GOCART2G_AerosolActivationMod.F90`), which skips NI | Moist (`aer_actv_single_moment.F90`) |

Each child also attaches `aerosol_optics`, `monochromatic_aerosol_optics` and `get_mixR` to its own `<name>_AERO` state in its `Initialize`; only the parent calls them, which is why its callbacks loop over the child states.
