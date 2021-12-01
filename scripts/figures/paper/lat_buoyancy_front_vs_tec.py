import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    ds = time_mean(xr.open_mfdataset(snakemake.input))

    fig, axe = plt.subplots(1, 1, figsize=(pc19, 3))
    ds.plot.scatter("alpha_deep_ML", "lat_buoyancy_front", ax=axe)
    axe.set_ylabel("Latitude of the buoyancy flux inversion [$^\circ$N]")
    axe.set_xlabel(r"TEC in the convective zone [$^\circ$C$^{-1}$]")
    axe.set_title("")

    axe.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
    
    fig.tight_layout()
    fig.savefig(snakemake.output[0])
