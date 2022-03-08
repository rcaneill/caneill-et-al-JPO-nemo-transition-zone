import xarray as xr
import xgcm
import numpy as np

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)

if __name__ == "__main__":
    ds, grid, old_vars = io.open_all(snakemake.input)

    ds["psi_baro"] = grid.cumint((ds.vo * ds.e3v).sum(dim="z_c"), "X") * 1e-6

    ds.drop_vars(old_vars, errors="ignore").to_netcdf(snakemake.output[0])
