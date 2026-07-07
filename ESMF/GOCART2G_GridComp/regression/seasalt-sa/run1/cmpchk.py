#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

PHASE = "1"
for state in ["import", "internal", "export"]:
    print("\nSTATE:", state)
    baseline = f"/home/pchakrab/input/gocart+seasalt+dust/C12-L181/after/new-style/SS_{state}_after_runPhase{PHASE}.nc"
    current = f"checkpoints/last/SS_{state}.nc"
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        for var in bas.variables:
            if var in cur.variables:
                bas_var = bas[var][:]
                cur_var = cur[var][0]
                diff = bas_var - cur_var
                diffnorm = np.linalg.norm(diff)
                print(f" {var:<20} diff (min, max, norm): {np.min(diff)}, {np.max(diff)}, {diffnorm}")
