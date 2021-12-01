import xarray as xr
import cf_xarray
import matplotlib.pyplot as plt
import cmocean.cm as cmo

if __name__ == '__main__':
    ds = xr.open_mfdataset(snakemake.input)

    fig, ax = plt.subplots(1, 1, figsize=(4,4))
    ax.grid(True)
    
    ax.scatter(ds.sigma0_star.cf.max('Y').max('day'), ds.sigma.isel(z_c=0).cf.max(['X','Y']).max('month'), label='time max')
    ax.scatter(ds.sigma0_star.cf.max('Y').mean('day'), ds.sigma.isel(z_c=0).cf.max(['X','Y']).mean('month'), label='time mean')

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.plot(*([[xlim[0]-1, xlim[1]+1]]*2), label='x=y')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.legend()
    
    ax.set_xlabel('$\mathrm{max}(\sigma_0^*)$ [kg$\,$m$^{-3}$]')
    ax.set_ylabel('$\mathrm{max}(\sigma_0^{surface})$ [kg$\,$m$^{-3}$]')
    
    fig.tight_layout()
    fig.savefig(snakemake.output[0])
