import xarray as xr
import xgcm
import numpy as np

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)

if __name__ == '__main__':
    ds, grid, old_vars = io.open_all(snakemake.input)

    zeta = 1/(ds.e1f*ds.e2f) * (grid.diff(ds.vo*ds.e2v, 'X') - grid.diff(ds.uo*ds.e1u, 'Y'))
    g = 9.81
    rho0 = 1026
    N2 = -g/rho0 * (-grid.derivative(ds.sigma, 'Z'))
    # We have an extra minus sign in the derivative as the Z axis is oriented downward in NEMO
    ds['qv'] = (zeta + ds.ff_f) * grid.interp(N2, ['X', 'Y', 'Z'])
    ds.coords['gdepf_0'] = grid.interp(ds.gdept_0, ['X', 'Y'])

    ds.drop_vars(old_vars, errors='ignore').to_netcdf(snakemake.output[0])
