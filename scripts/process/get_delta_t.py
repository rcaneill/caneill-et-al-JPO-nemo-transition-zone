from pathlib import Path
import xarray as xr

if __name__ == "__main__":
    filepath = Path(snakemake.input[0])
    # 2 cases: 1) we have yearly outputs + monthly surface or 2) monthly outputs
    # In case 1), we filter out the surface fields
    nemo_files = (filepath / "VARS").glob("*1y*grid*")

    try:
        next(nemo_files)
        delta_t = "1y"
    except StopIteration:
        delta_t = "1m"

    ds = xr.Dataset()
    exp = [snakemake.params["expname"]]
    ds["delta_t"] = xr.DataArray(delta_t, coords=[exp], dims=["exp"])
    ds.to_netcdf(snakemake.output[0])
