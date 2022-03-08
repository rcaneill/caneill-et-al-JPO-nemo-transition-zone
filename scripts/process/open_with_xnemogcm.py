import xnemogcm
import f90nml
from pathlib import Path
import xarray as xr
import numpy as np
import os

xr.set_options(keep_attrs=True)

if __name__ == "__main__":

    remove_outputs = [
        "depthu_bounds",
        "t_bounds",
        "deptht_bounds",
        "depthv_bounds",
        "depthw_bounds",
    ]

    coordinates = [
        "glamt",
        "glamu",
        "glamv",
        "glamf",
        "gphit",
        "gphiu",
        "gphiv",
        "gphif",
        "gdept_0",
        "gdepw_0",
    ]

    filepath = Path(snakemake.input[0])

    t_file = xr.open_dataset(snakemake.input["delta_t_file"])
    delta_t = t_file.delta_t.values[0]
    t_file.close()

    nemo_files = list((filepath / "VARS").glob(f"*{delta_t}*grid*"))
    domcfg_filepath = filepath

    ds = xnemogcm.open_nemo_and_domain_cfg(
        nemo_files=nemo_files, domcfg_files=domcfg_filepath
    ).chunk({"t": 10})

    # Take the time average
    if delta_t == "1y":
        ds = ds.mean("t")
    elif delta_t == "1m":
        time_vars = [i for i in ds if "t" in ds[i].dims]
        ds = ds.drop_vars(time_vars).merge(ds[time_vars].groupby("t.month").mean("t"))
    else:
        raise (NotImplementedError)

    ds = ds.drop_vars(["t"], errors="ignore")

    nam = xnemogcm.open_namelist(domcfg_filepath, files=["namelist_cfg"], ref=False)
    ds = ds.merge(nam)

    # We remove unnecessary outputs
    ds = ds.drop_vars(remove_outputs, errors="ignore")

    # We add a dimension for the experience
    for i in ds:
        ds[i] = ds[i].expand_dims({"exp": [snakemake.params["expname"]]})

    # Here all the phi, and lambda can be 1d
    for i in ds.variables:
        if "gphi" in i:
            for dim in ds[i].dims:
                if "x" in dim:
                    ds[i] = ds[i].isel({dim: 0}).drop_vars(dim)
        if "glam" in i:
            for dim in ds[i].dims:
                if "y" in dim:
                    ds[i] = ds[i].isel({dim: 0}).drop_vars(dim)

    # We add coordinates for the plots
    for coord in coordinates:
        ds.coords[coord] = ds[coord]

    # Add cf coordinate
    for i in ds.variables:
        if "glam" in i:
            ds[i].attrs.update({"units": "degrees_east"})
        if "gphi" in i:
            ds[i].attrs.update({"units": "degrees_north"})
    # Here glam{T,U,V,F} are 1D, we can drop half of them
    # same for gphi
    ds = ds.drop_vars(["glamv", "glamf", "gphiu", "gphif"])

    for eiv in [i for i in ["uoce_eiv", "voce_eiv"] if i in ds]:
        ds[eiv] = xr.where(np.isnan(ds[eiv]), 0, ds[eiv])

    ds["tmask_surf"] = ds.tmask.isel(z_c=0)
    ds["Cb"] = ds.rn_a0 * ds.rn_lambda1

    # equator metrics
    for i in "tuvf":
        ds[f"e2{i}"] = ds[f"e2{i}"].copy().load()
        # needed to load a copy to get a numpy array and not dask
        # + get an error if not using a copy:
        # ValueError: assignment destination is read-only
    ds.e2t[{"y_c": 0}] = 0
    ds.e2u[{"y_c": 0}] = 0
    ds.e2t[{"y_c": 1}] = ds.e2t[{"y_c": 1}] / 2
    ds.e2u[{"y_c": 1}] = ds.e2u[{"y_c": 1}] / 2
    ds.e2v[{"y_f": 0}] = 0
    ds.e2f[{"y_f": 0}] = 0

    ds.drop_vars(
        [
            "time_instant_bounds",
            "time_instant",
            "time_centered_bounds",
            "time_centered",
        ],
        errors="ignore",
    ).to_netcdf(snakemake.output[0])
