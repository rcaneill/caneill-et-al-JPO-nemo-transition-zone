import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    ds = time_mean(xr.open_mfdataset(snakemake.input))
    
    fig, axe = plt.subplots(1, 1, figsize=(pc19, 3))
    ds.plot.scatter("lat_buoyancy_front", "lat_fresh_pool", ax=axe)
    axe.set_xlabel("Latitude of the buoyancy flux inversion [$^\circ$N]")
    axe.set_ylabel("Southern latitude of the beta ocean [$^\circ$N]")
    axe.set_title("")
    #axe.set_xticks([52, 54, 56, 58])
    #axe.set_yticks([50, 52, 54, 56, 58])
    #axe.set_xlim(50,58)
    #axe.set_ylim(50,58)

    fig.tight_layout()
    fig.savefig(snakemake.output[0])
