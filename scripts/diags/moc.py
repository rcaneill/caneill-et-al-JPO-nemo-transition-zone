import xarray as xr
import xgcm
import numpy as np

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)

if __name__ == "__main__":
    ds, grid, old_vars = io.open_all(snakemake.input)

    e3 = grid.interp(ds.e3w_0, "Y")
    ds.coords["gdepv_0"] = (
        grid.cumsum(e3, axis="Z", boundary="fill", fill_value=0)
        - e3.isel({"z_f": 0}).drop_vars("z_f") / 2
    )

    e3 = grid.interp(ds.e3w_0, "X")
    ds.coords["gdepu_0"] = (
        grid.cumsum(e3, axis="Z", boundary="fill", fill_value=0)
        - e3.isel({"z_f": 0}).drop_vars("z_f") / 2
    )

    ds = ds.isel(z_c=slice(None, -1))

    grid = xgcm.Grid(
        ds,
        metrics={
            ("X",): ["e1t", "e1u", "e1v", "e1f"],
            ("Y",): ["e2t", "e2u", "e2v", "e2f"],
            ("Z",): ["e3t", "e3u", "e3v", "e3w"],
        },
        periodic=False,
    )

    ds["Fv"] = (ds.vo + ds.voce_eiv) * ds.e3v
    ds.Fv.attrs = {
        "standard_name": "velocity flux along j-axis",
        "long_name": "velocity_flux_along_j-axis",
        "units": "m^2/s",
    }
    ds["depth_W_uniform"] = xr.DataArray(
        np.append(
            np.linspace(0, 1000, 101),
            np.linspace(1100, 4000, int((4000 - 1100) / 100 + 1)),
        ),
        dims=["depth_W_uniform"],
    )

    # target: new horizontally uniform depths of W points
    # target_data: original depth of W points
    ds["Fv_transformed"] = grid.transform(
        ds.Fv,
        "Z",
        target=xr.DataArray(ds["depth_W_uniform"].values, dims=["depth_T_uniform"]),
        method="conservative",
        target_data=grid.interp(ds.gdepw_0, "Y", boundary="extend"),
    )

    ds["depth_W_uniform"].attrs = {"axis": "Zu", "c_grid_axis_shift": +0.5}
    ds["depth_T_uniform"].attrs = {"axis": "Zu"}

    ds["thickness"] = xgcm.Grid(ds, periodic=False, boundary="extend").diff(
        ds.depth_W_uniform, "Zu"
    )
    ds["vo_transformed"] = ds.Fv_transformed / ds["thickness"]

    metrics = {
        ("X",): ["e1t", "e1u", "e1v", "e1f"],
        ("Y",): ["e2t", "e2u", "e2v", "e2f"],
        ("Z",): ["e3t", "e3u", "e3v", "e3w"],
        ("Zu",): ["thickness"],
    }
    gridu = xgcm.Grid(ds, periodic=False, metrics=metrics, boundary="extend")

    ds["psi_moc"] = (
        gridu.cumsum(ds.thickness * gridu.integrate(ds.vo_transformed, ["X"]), "Zu")
        * 1e-6
    )
    ds["moc_depth_strength"] = ds.psi_moc.where(ds.gphiv > 20).max(
        ["depth_W_uniform", "y_f"]
    )

    # MOC on density levels
    ds["v_transport_sigma"] = grid.transform(
        ds.Fv,
        "Z",
        xr.DataArray(np.arange(22, 28, 0.1), dims=["sigma_uniform"]),
        target_data=grid.interp(ds.sigma, "Y", boundary="extend"),
        method="conservative",
    )
    ds["sigma_thickness"] = np.diff(ds["sigma_uniform"]).mean()
    ds["psi_moc_sigma"] = (grid.integrate(ds.v_transport_sigma, "X")).cumsum(
        dim="sigma_uniform"
    ) * 1e-6
    ds["psi_moc_sigma"].attrs = {"long_name": "meridional overturning streamfunction"}
    ds["moc_density_strength"] = ds.psi_moc_sigma.where(ds.sigma_uniform > 25.5).max(
        ["sigma_uniform", "y_f"]
    )

    ds.drop_vars(old_vars, errors="ignore").to_netcdf(snakemake.output[0])
