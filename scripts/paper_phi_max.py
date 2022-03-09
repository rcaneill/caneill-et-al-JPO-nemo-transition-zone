import xarray as xr

ds = xr.open_mfdataset(snakemake.input)

with open(snakemake.output[0], "w") as f:
    f.write(ds.phi_max.to_pandas().to_csv())
    f.close()
