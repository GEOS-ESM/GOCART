#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

PHASE = "1"
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
                if var in cur.variables:
                    bas_var = bas[var][:]
                    cur_var = cur[var][0] # remove the time=1 dimension
                    diff = bas_var - cur_var
                    diffnorm = np.linalg.norm(diff)
                    # print()
                    # print(bas_var.shape, cur_var.shape)
                    # print(f"BAS_{var:<20} (min, max): {np.min(bas_var)}, {np.max(bas_var)}")
                    # print(f"CUR_{var:<20} (min, max): {np.min(cur_var)}, {np.max(cur_var)}")
                    print(f"   {var:<20} diff (min, max, norm): {np.min(diff)}, {np.max(diff)}, {diffnorm}")

