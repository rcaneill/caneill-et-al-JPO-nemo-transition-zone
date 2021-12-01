import xarray as xr
import xgcm
import numpy as np

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)

"""
# Ekman velocity in spherical coordinates

$w_{ek}= 
- \frac{1}{\rho_0} \pdv{}{y}\qty(\frac{\tau_\lambda}{f}) + \frac{\tau_\lambda}{\rho_0 \beta R^2}$
"""


if __name__ == '__main__':
    ds, grid, old_vars = io.open_all(snakemake.input)
    
    R = 6.3781e6
    beta = grid.interp(grid.derivative(ds.ff_t, 'Y'), 'X') # F point
    
    w1 = - 1/1026 * grid.derivative(ds.utau / grid.interp(ds.ff_t, 'X'), 'Y').where(ds.fmask.isel(z_c=0)).mean(['exp', 'x_f']) # F point
    w2 = (grid.interp(ds.utau, 'Y') / (1026 * beta * R**2)).where(ds.fmask.isel(z_c=0)).mean(['exp', 'x_f']) # F point
    ds['w_ek'] = w1 + w2

    ds.drop_vars(old_vars, errors='ignore').to_netcdf(snakemake.output[0])
