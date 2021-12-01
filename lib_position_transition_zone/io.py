"""Functions for the input / output"""


import xarray as xr
import xgcm

def open_all(files):
    """
    Open all provided nc files, and create the grid.
    
    Parameters
    ----------
        files: list or generator of file paths
    
    Returns
    -------
        (ds, grid, old_vars)
    """
    ds = xr.open_mfdataset(files)
    
    old_vars = [i for i in ds.variables if i not in ds.dims]
    
    metrics = {
        ('X',): ['e1t', 'e1u', 'e1v', 'e1f'],
        ('Y',): ['e2t', 'e2u', 'e2v', 'e2f'],
        ('Z',): ['e3t', 'e3u', 'e3v', 'e3w']
    }
    grid = xgcm.Grid(ds, periodic=False, metrics=metrics, boundary='extend')

    return (ds, grid, old_vars)
