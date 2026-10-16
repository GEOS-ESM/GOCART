#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

PHASE = "1"
GRIDCOMPS = {
    "DU": ["import", "internal", "export"],
    "SS": ["import", "internal", "export"],
    "SU": ["import", "internal", "export"],
    "GOCART2G": ["import", "export"],
}

for gridcomp, states in GRIDCOMPS.items():
    print("\nGRIDCOMP:", gridcomp)
    for state in states:
        print("STATE:", state)
        baseline = f"baseline/{gridcomp}_{state}_after_runPhase{PHASE}.nc"
        current = f"checkpoints/last/{gridcomp}_{state}.nc"
        with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
            for var in bas.variables:
                if var in cur.variables:
                    bas_var = bas[var][:]
                    cur_var = cur[var][0]  # remove the time=1 dimension
                    diff = bas_var - cur_var
                    diffnorm = np.linalg.norm(diff)
                    print(f"   {var:<20} diff (min, max, norm): {np.min(diff)}, {np.max(diff)}, {diffnorm}")
