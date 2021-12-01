import xarray as xr
import xgcm
import numpy as np

from lib_position_transition_zone import io

from xbasin import eos

xr.set_options(keep_attrs=True)

# variables for the forcings (from usrdef_sbc.F90)

zL = 61       # [degrees] Approximative meridional extend of the basin       
# Wind stress 
ztau0  =   0.1    # [Pa] 
zatau  =   0.8    # [no unit]   
zdeltatau     =   5.77   # [deg North] 
# T star and qns     
#ztrp   = -40.     # [W/m2/K] retroaction term on heat fluxes      
zts_eq =  25.     # [deg C] Temperature at the equator      
zts_n  =   0.     # [deg C] Temperature in the north 
zdeltasst     =  16.22   # [deg North] 
z1_2L  =   1. / (2. * zL)    
zdeltaT = 2   # [deg C] half difference of temperature during winter and summer in the north (magnitude of the cos) ##rc TODO set in namelist      
# EMP  
zconv  =   1. / ( 86400.)   # convertion factor: 1 mm/day => 1/(3600*24) mm/s 
##rc TODO put a1, a2 and a3 in namelist   
za1 = -3.24       # [mm/day] Set the amplitude of EMP at the equator      
za2 = 4.15 # [mm/day] Set the amplitude of EMP at mid latitude     
za3 = -1.59       # [mm/day] Set the amplitude of EMP at the northern part       
zphi01 = 0.       # [deg North] 
zphi02 = 20.      # [deg North] 
zphi03 = 50.      # [deg North] 
z1_d1 = 1. / 8.       # [1 / deg North]    
zd2 = 30.  # [deg North] 
z1_d2 = 1. / zd2  # [1 / deg North]    
zd3 = 40.  # [deg North] 
z1_d3 = 1. / zd3  # [1 / deg North]    
z1_s = 1. / 10.       # streching of the tanh function (i.e. smoothness of the filter)      
za1 = za1 * zconv    # [mm/s] after  conversion  
za2 = za2 * zconv    # [mm/s] after  conversion  
za3 = za3 * zconv    # [mm/s] after  conversion

if __name__ == '__main__':
    ds, grid, old_vars = io.open_all(snakemake.input)

    ds['zsstar'] = 37.12 * np.exp( - ds.gphit**2 / 260.**2 ) - 1.1 * np.exp( - ds.gphit**2 / 7.5**2)

    ds['day'] = xr.DataArray(np.arange(360), coords=[("day", np.arange(360))])
    if ds.ln_ann_cyc:
        zcos_sais2 = np.cos((ds.day - 21 + 180)/180 * np.pi)
        zcos_sais1 = np.cos((ds.day + 19 + 180)/180 * np.pi)
    else:
        zcos_sais2 = 0 * ds.day
        zcos_sais1 = 0 * ds.day
    ds['ztstar'] = (zts_eq - (zts_n + zcos_sais2 * zdeltaT) ) * np.cos((np.pi * ds.gphit) * z1_2L)**2 + (zts_n + zcos_sais2 * zdeltaT)
    ds['qsolar_analytic'] = 230. * np.cos(np.pi * (ds.gphit - 23.5 * zcos_sais1 ) / ( 0.9 * 180 ) )

    ds.ztstar.attrs['long_name'] = 'T*'
    ds.ztstar.attrs['units'] = 'degC'
    ds.zsstar.attrs['long_name'] = 'S*'
    ds.zsstar.attrs['units'] = '1e-3'
    ds.qsr.attrs['long_name'] = 'Q solar'

    ds['sigma0_star'] = eos.compute_sigma0(
        ds.ztstar, ds.zsstar,
        rn_lambda1=ds.rn_lambda1,
        rn_lambda2=ds.rn_lambda2,
        rn_a0=ds.rn_a0,
        rn_b0=ds.rn_b0,
        rn_mu1=ds.rn_mu1,
        rn_mu2=ds.rn_mu2,
        rn_nu=ds.rn_nu
    )
    ds['sigma0_star'].attrs['long_name'] = 'sigma0*'
    ds['sigma0_star'].attrs['units'] = 'kg/m3'


    # We use a higher resolution for gphit to compute phi_max
    gphit = xr.DataArray(np.linspace(0,60,1000), dims=['y_high_res'])
    tstar = (zts_eq - (zts_n + zcos_sais2 * zdeltaT) ) * np.cos((np.pi * gphit) * z1_2L)**2 + (zts_n + zcos_sais2 * zdeltaT)
    sstar = 37.12 * np.exp( - gphit**2 / 260.**2 ) - 1.1 * np.exp( - gphit**2 / 7.5**2)
    rho = eos.compute_sigma0(
        tstar, sstar,
        rn_lambda1=ds.rn_lambda1,
        rn_lambda2=ds.rn_lambda2,
        rn_a0=ds.rn_a0,
        rn_b0=ds.rn_b0,
        rn_mu1=ds.rn_mu1,
        rn_mu2=ds.rn_mu2,
        rn_nu=ds.rn_nu
    )
    ds['phi_max_day'] = gphit.where(rho == rho.max('y_high_res')).mean('y_high_res')
    ds['phi_max_day'].attrs['long_name'] = 'latitude_of_maximum_sigma0*'
    ds['phi_max_day'].attrs['units'] = 'degrees_north'
    ds['phi_max'] = ds['phi_max_day'].max('day')
    
    ds.drop_vars(old_vars, errors='ignore').to_netcdf(snakemake.output[0])
