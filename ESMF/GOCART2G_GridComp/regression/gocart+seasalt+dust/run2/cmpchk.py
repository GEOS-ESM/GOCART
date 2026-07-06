#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

def alias(name):
    # # SS exports
    # if name == "SS_AERO_DP_SSDP001": return "SS_AERO_DP_1"
    # if name == "SS_AERO_DP_SSDP002": return "SS_AERO_DP_2"
    # if name == "SS_AERO_DP_SSDP003": return "SS_AERO_DP_3"
    # if name == "SS_AERO_DP_SSDP004": return "SS_AERO_DP_4"
    # if name == "SS_AERO_DP_SSDP005": return "SS_AERO_DP_5"
    # if name == "SS_AERO_DP_SSSV001": return "SS_AERO_DP_6"
    # if name == "SS_AERO_DP_SSSV002": return "SS_AERO_DP_7"
    # if name == "SS_AERO_DP_SSSV003": return "SS_AERO_DP_8"
    # if name == "SS_AERO_DP_SSSV004": return "SS_AERO_DP_9"
    # if name == "SS_AERO_DP_SSSV005": return "SS_AERO_DP_10"
    # if name == "SS_AERO_DP_SSWT001": return "SS_AERO_DP_11"
    # if name == "SS_AERO_DP_SSWT002": return "SS_AERO_DP_12"
    # if name == "SS_AERO_DP_SSWT003": return "SS_AERO_DP_13"
    # if name == "SS_AERO_DP_SSWT004": return "SS_AERO_DP_14"
    # if name == "SS_AERO_DP_SSWT005": return "SS_AERO_DP_15"
    # if name == "SS_AERO_DP_SSSD001": return "SS_AERO_DP_16"
    # if name == "SS_AERO_DP_SSSD002": return "SS_AERO_DP_17"
    # if name == "SS_AERO_DP_SSSD003": return "SS_AERO_DP_18"
    # if name == "SS_AERO_DP_SSSD004": return "SS_AERO_DP_19"
    # if name == "SS_AERO_DP_SSSD005": return "SS_AERO_DP_20"
    # # DU exports
    # if name == "DU_AERO_DP_DUDP001": return "DU_AERO_DP_1"
    # if name == "DU_AERO_DP_DUDP002": return "DU_AERO_DP_2"
    # if name == "DU_AERO_DP_DUDP003": return "DU_AERO_DP_3"
    # if name == "DU_AERO_DP_DUDP004": return "DU_AERO_DP_4"
    # if name == "DU_AERO_DP_DUDP005": return "DU_AERO_DP_5"
    # if name == "DU_AERO_DP_DUSV001": return "DU_AERO_DP_6"
    # if name == "DU_AERO_DP_DUSV002": return "DU_AERO_DP_7"
    # if name == "DU_AERO_DP_DUSV003": return "DU_AERO_DP_8"
    # if name == "DU_AERO_DP_DUSV004": return "DU_AERO_DP_9"
    # if name == "DU_AERO_DP_DUSV005": return "DU_AERO_DP_10"
    # if name == "DU_AERO_DP_DUWT001": return "DU_AERO_DP_11"
    # if name == "DU_AERO_DP_DUWT002": return "DU_AERO_DP_12"
    # if name == "DU_AERO_DP_DUWT003": return "DU_AERO_DP_13"
    # if name == "DU_AERO_DP_DUWT004": return "DU_AERO_DP_14"
    # if name == "DU_AERO_DP_DUWT005": return "DU_AERO_DP_15"
    # if name == "DU_AERO_DP_DUSD001": return "DU_AERO_DP_16"
    # if name == "DU_AERO_DP_DUSD002": return "DU_AERO_DP_17"
    # if name == "DU_AERO_DP_DUSD003": return "DU_AERO_DP_18"
    # if name == "DU_AERO_DP_DUSD004": return "DU_AERO_DP_19"
    # if name == "DU_AERO_DP_DUSD005": return "DU_AERO_DP_20"
    # # GOCART2G exports
    # if name == "AERO_DP_DUDP001": return "AERO_DP_1"
    # if name == "AERO_DP_DUDP002": return "AERO_DP_2"
    # if name == "AERO_DP_DUDP003": return "AERO_DP_3"
    # if name == "AERO_DP_DUDP004": return "AERO_DP_4"
    # if name == "AERO_DP_DUDP005": return "AERO_DP_5"
    # if name == "AERO_DP_DUSD001": return "AERO_DP_6"
    # if name == "AERO_DP_DUSD002": return "AERO_DP_7"
    # if name == "AERO_DP_DUSD003": return "AERO_DP_8"
    # if name == "AERO_DP_DUSD004": return "AERO_DP_9"
    # if name == "AERO_DP_DUSD005": return "AERO_DP_10"

    # if name == "AERO_DP_DUSV001": return "AERO_DP_11"
    # if name == "AERO_DP_DUSV002": return "AERO_DP_12"
    # if name == "AERO_DP_DUSV003": return "AERO_DP_13"
    # if name == "AERO_DP_DUSV004": return "AERO_DP_14"
    # if name == "AERO_DP_DUSV005": return "AERO_DP_15"
    # if name == "AERO_DP_DUWT001": return "AERO_DP_16"
    # if name == "AERO_DP_DUWT002": return "AERO_DP_17"
    # if name == "AERO_DP_DUWT003": return "AERO_DP_18"
    # if name == "AERO_DP_DUWT004": return "AERO_DP_19"
    # if name == "AERO_DP_DUWT005": return "AERO_DP_20"

    # if name == "AERO_DP_SSDP001": return "AERO_DP_21"
    # if name == "AERO_DP_SSDP002": return "AERO_DP_22"
    # if name == "AERO_DP_SSDP003": return "AERO_DP_23"
    # if name == "AERO_DP_SSDP004": return "AERO_DP_24"
    # if name == "AERO_DP_SSDP005": return "AERO_DP_25"
    # if name == "AERO_DP_SSSD001": return "AERO_DP_26"
    # if name == "AERO_DP_SSSD002": return "AERO_DP_27"
    # if name == "AERO_DP_SSSD003": return "AERO_DP_28"
    # if name == "AERO_DP_SSSD004": return "AERO_DP_29"
    # if name == "AERO_DP_SSSD005": return "AERO_DP_30"

    # if name == "AERO_DP_SSSV001": return "AERO_DP_31"
    # if name == "AERO_DP_SSSV002": return "AERO_DP_32"
    # if name == "AERO_DP_SSSV003": return "AERO_DP_33"
    # if name == "AERO_DP_SSSV004": return "AERO_DP_34"
    # if name == "AERO_DP_SSSV005": return "AERO_DP_35"
    # if name == "AERO_DP_SSWT001": return "AERO_DP_36"
    # if name == "AERO_DP_SSWT002": return "AERO_DP_37"
    # if name == "AERO_DP_SSWT003": return "AERO_DP_38"
    # if name == "AERO_DP_SSWT004": return "AERO_DP_39"
    # if name == "AERO_DP_SSWT005": return "AERO_DP_40"

    return name

PHASE = "2"
GRIDCOMPS = {
    "DU": ["import", "internal", "export"],
    "SS": ["import", "internal", "export"],
    "GOCART2G": ["import", "export"]
}

BASELINE_DIR = "/home/pchakrab/input/gocart+seasalt+dust/C12-L181/after/new-style/"
for gridcomp in GRIDCOMPS:
    print("\nGRIDCOMP:", gridcomp)
    for state in GRIDCOMPS[gridcomp]:
        print("STATE:", state)
        baseline = f"{BASELINE_DIR}/{gridcomp}_{state}_after_runPhase{PHASE}.nc"
        current = f"checkpoints/last/{gridcomp}_{state}.nc"
        with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
            for var in bas.variables:
                if alias(var) in cur.variables:
                    bas_var = bas[var][:]
                    cur_var = cur[alias(var)][0] # remove the time=1 dimension
                    diff = bas_var - cur_var
                    diffnorm = np.linalg.norm(diff)
                    # print()
                    # print(bas_var.shape, cur_var.shape)
                    # print(f"BAS_{var:<20} (min, max): {np.min(bas_var)}, {np.max(bas_var)}")
                    # print(f"CUR_{var:<20} (min, max): {np.min(cur_var)}, {np.max(cur_var)}")
                    print(f"   {var:<20} diff (min, max, norm): {np.min(diff)}, {np.max(diff)}, {diffnorm}")

