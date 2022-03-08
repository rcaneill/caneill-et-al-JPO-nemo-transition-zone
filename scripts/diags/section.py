import xarray as xr
import xgcm
import numpy as np

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)


if __name__ == "__main__":
    ds, grid, old_vars = io.open_all(snakemake.input["data"])

    # We want to interpolate the sci and the mld so we only take them
    ds = ds[["sci", "mldr10_1", "delta_t"]]

    # we change the x_c coordinate to glamt
    ds = ds.swap_dims({"x_c": "glamt"})

    if ds.delta_t.values[0] == "1m":
        out = ds.interp(glamt=snakemake.params["lon_sections"]).sel(
            month=snakemake.params["month_sections"]
        )
    elif ds.delta_t.values[0] == "1y":
        out = ds.interp(glamt=snakemake.params["lon_sections"])

    out.to_netcdf(snakemake.output[0])
