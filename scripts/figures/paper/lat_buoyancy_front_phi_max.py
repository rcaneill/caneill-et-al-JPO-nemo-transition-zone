import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    ds = time_mean(xr.open_mfdataset(snakemake.input))
    
    fig, ax = plt.subplots(2, 1, figsize=(pc19, 4))
    ds.plot.scatter("phi_max", "lat_buoyancy_front", ax=ax[0], label='Latitude of buoyancy flux inversion')
    ds.plot.scatter("phi_max", "lat_fresh_pool", ax=ax[0], label='Southern latitude of beta ocean')
    ax[0].legend()
    ax[0].set_xlabel(r"$\varphi_M$ [$^\circ$N]")
    ax[0].set_ylabel("[$^\circ$N]")
    ax[0].set_title("a)")

    ds.plot.scatter("lat_buoyancy_front", "lat_fresh_pool", ax=ax[1])
    ax[1].set_xlabel("Latitude of the buoyancy flux inversion [$^\circ$N]")
    ax[1].set_ylabel("Southern latitude of the beta ocean [$^\circ$N]")
    ax[1].set_title("b)")


    fig.tight_layout()
    fig.savefig(snakemake.output[0])
