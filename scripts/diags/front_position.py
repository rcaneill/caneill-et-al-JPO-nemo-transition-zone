import xarray as xr
import xgcm
import numpy as np

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)


def tmean(da):
    return da.mean("t")


if __name__ == "__main__":
    ds, grid, old_vars = io.open_all(snakemake.input)

    ds["beta_volume"] = grid.integrate(xr.where(ds.sci <= -1, 1, 0), ["X", "Y", "Z"])

    # what makes sense is to use the year average of the buoyancy flux
    try:
        F_b = ds.F_b.groupby("t.year").map(tmean)
    except KeyError:
        # this means that the data are already time averaged
        F_b = ds.F_b
    # same for mld
    try:
        mld = ds.mldr10_1.groupby("t.year").map(tmean)
    except KeyError:
        # this means that the data are already time averaged
        mld = ds.mldr10_1

    ds["b_gain_surf_north"] = grid.integrate(
        xr.where((F_b.where(ds.tmask_surf) >= 0) & (ds.gphit >= 48), 1, 0), ["X", "Y"]
    )
    ds["fresh_pool_surf_using_mld"] = grid.integrate(
        xr.where((mld.where(ds.tmask_surf) <= 100) & (ds.gphit >= 48), 1, 0),
        ["X", "Y"],
    )
    sci_under_ml = ds.sci_under_ml
    ds["fresh_pool_surf"] = grid.integrate(
        xr.where((sci_under_ml <= -1), 1, 0), ["X", "Y"]
    )

    phi2 = ds.gphiv.isel(y_f=-2)
    lam1 = ds.glamu.isel(x_f=0)
    lam2 = ds.glamu.isel(x_f=-2)

    ds["lat_buoyancy_front"] = np.rad2deg(
        np.arcsin(
            np.sin(np.deg2rad(phi2))
            - ds["b_gain_surf_north"] / (40 * 111198.92344855**2 * 180 / np.pi)
        )
    )
    ds["lat_fresh_pool"] = np.rad2deg(
        np.arcsin(
            np.sin(np.deg2rad(phi2))
            - ds["fresh_pool_surf"] / (40 * 111198.92344855**2 * 180 / np.pi)
        )
    )
    ds["lat_fresh_pool_using_mld"] = np.rad2deg(
        np.arcsin(
            np.sin(np.deg2rad(phi2))
            - ds["fresh_pool_surf_using_mld"]
            / (40 * 111198.92344855**2 * 180 / np.pi)
        )
    )

    ds.drop_vars(old_vars, errors="ignore").to_netcdf(snakemake.output[0])
