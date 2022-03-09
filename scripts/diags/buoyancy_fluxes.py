import os
from pathlib import Path

import xarray as xr
import xgcm
import numpy as np
from xbasin import eos, surface_flux

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)

if __name__ == "__main__":
    ds, grid, old_vars = io.open_all(snakemake.input)

    ds = ds.merge(surface_flux.nemo_wrap(ds), compat="override")

    ds.drop_vars(old_vars + ["gdept_1d", "gdept_0"], errors="ignore").to_netcdf(
        snakemake.output[0]
    )
