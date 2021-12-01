import xarray as xr
import xgcm
import numpy as np

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)

def tmean(da):
    return da.mean('t')

if __name__ == '__main__':
    ds, grid, old_vars = io.open_all(snakemake.input)

    # Mean T-S

    ds["mean_T"] = grid.average(ds.thetao.where(ds.tmask), ["X", "Y", "Z"])
    ds["mean_S"] = grid.average(ds.so.where(ds.tmask), ["X", "Y", "Z"])

    ds["abyss_T"] = grid.average(ds.thetao.where(ds.tmask).isel(z_c=-2), ["X", "Y"])
    ds["abyss_S"] = grid.average(ds.so.where(ds.tmask).isel(z_c=-2), ["X", "Y"])

    # the next variables make only sense for annual averages
    mld = ds.mldr10_1 #.groupby('t.year').map(tmean)
    tos = ds.tos #.groupby('t.year').map(tmean)
    alpha = ds.alpha #.groupby('t.year').map(tmean)
    ds["T_deep_ML"] = grid.average(tos.where(mld > 800), ["X", "Y"])
    ds["alpha_deep_ML"] = grid.average(
        alpha.isel(z_c=0).where(mld > 800), ["X", "Y"]
    )
    
    ds["alpha_north"] = ds.alpha.isel(z_c=0, y_c=-2).mean("x_c")

    ds.drop_vars(old_vars, errors='ignore').to_netcdf(snakemake.output[0])
