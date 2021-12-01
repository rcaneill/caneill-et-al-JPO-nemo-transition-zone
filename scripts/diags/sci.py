import xarray as xr
import xgcm
import numpy as np
from xbasin import stratification, eos

from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)


def sci_under_ml(ds):
    return ds.sci_under_ml_at_every_time.where(
        (ds.mldr10_1 == ds.mldr10_1.max('month'))
        &
        (ds.mldr10_1.max('month')<ds.gdepw_0.isel(z_f=-2))
    ).mean('month')
    
    


if __name__ == '__main__':
    ds, grid, old_vars = io.open_all(snakemake.input)

    ds['sci'] = stratification.compute_strati_ratio(ds.thetao, ds.so, grid.interp(ds.alpha, 'Z'), grid.interp(ds.beta, 'Z'), grid)

    # SCI under the ML
    """
    We follow the method imagined by Fabien Roquet:
    we compute a new variable phi=depth - MLD and we remap the SCI
    (or equivalent, T and S profiles) onto this new variable phi.
    We then get the value at e.g. 0 (bottom of the ML), 20 (20 meters under the ML), etc.
    """

    rn_a0 = 1.655e-1
    rn_b0 = 7.6554e-1
    rn_lambda1 = 5.9520e-2
    rn_lambda2 = 5.4914e-4
    rn_nu = 2.4341e-3
    rn_mu1 = 1.4970e-4
    rn_mu2 = 1.1090e-5
    nameos = dict(
        rn_lambda1=ds.get("rn_lambda1", default=rn_lambda1),
        rn_lambda2=ds.get("rn_lambda2", default=rn_lambda2),
        rn_a0=ds.get("rn_a0", default=rn_a0),
        rn_b0=ds.get("rn_b0", default=rn_b0),
        rn_mu1=ds.get("rn_mu1", default=rn_mu1),
        rn_mu2=ds.get("rn_mu2", default=rn_mu2),
        rn_nu=ds.get("rn_nu", default=rn_nu),
    )
    ds['phi_dep_mld_t'] = ds.gdept_0 - ds.mldr10_1
    target = np.array([10, 30])

    T = grid.transform(ds.thetao, 'Z', target=target, target_data=ds.phi_dep_mld_t, method='linear')
    S = grid.transform(ds.so, 'Z', target=target, target_data=ds.phi_dep_mld_t, method='linear')
    alpha = eos.compute_alpha(T.mean('phi_dep_mld_t'), S.mean('phi_dep_mld_t'), 0, **nameos)
    beta = eos.compute_beta(T.mean('phi_dep_mld_t'), S.mean('phi_dep_mld_t'), 0, **nameos)
    
    g = 9.81
    N2_T = g * alpha * T.diff('phi_dep_mld_t').squeeze('phi_dep_mld_t') / (target[1] - target[0])
    N2_S = - g * beta * S.diff('phi_dep_mld_t').squeeze('phi_dep_mld_t') / (target[1] - target[0])
    
    N2 = N2_T + N2_S
    SCI = (N2_T - N2_S) / N2

    ds['sci_under_ml_at_every_time'] = SCI.drop_vars('phi_dep_mld_t')

    # Now if we have monthly data we can compute the sci_under_ml at
    # the time of deepest ML over the year

    if ds.delta_t.values[0] == '1m':
        ds['sci_under_ml'] = sci_under_ml(ds)
    elif ds.delta_t.values[0] == '1y':
        ds['sci_under_ml'] = ds['sci_under_ml_at_every_time']
    
    ds = ds.drop_vars(['phi_dep_mld_w', 'phi_dep_mld_t'], errors='ignore')

    ds.drop_vars(old_vars, errors='ignore').to_netcdf(snakemake.output[0])
