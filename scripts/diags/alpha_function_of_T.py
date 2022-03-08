import xarray as xr
import numpy as np
from xbasin.eos import compute_alpha


xr.set_options(keep_attrs=True)

if __name__ == "__main__":
    ds0 = xr.open_dataset(snakemake.input[0])

    ds = xr.Dataset()

    T = xr.DataArray(np.linspace(0, 10, 6), name="Theta")
    alpha = compute_alpha(T, 35, 0, rn_a0=ds0.rn_a0, rn_lambda1=ds0.rn_lambda1)
    alpha.coords["T"] = T

    ds["alpha"] = alpha

    ds.to_netcdf(snakemake.output[0])
