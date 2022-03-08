import xarray as xr
import csv

if __name__ == "__main__":
    # inefficient, we open all the file when only exp would be needed
    ds0 = xr.open_dataset(snakemake.input["xnemo"], decode_times=False)

    with open(snakemake.input["expnames"], "r") as file:
        expnames = {i[0]: i[1] for i in csv.reader(file)}
        exp_short_name = expnames[ds0.exp.values[0]]

        ds = xr.Dataset()
        ds["short_name"] = xr.DataArray(exp_short_name, coords=[ds0.exp])
        ds.to_netcdf(snakemake.output[0])
