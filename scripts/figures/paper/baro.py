import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == "__main__":
    with sns.axes_style("ticks"):
        ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0))

        fig, axe = plt.subplots(1, 1, figsize=(pc19, 3.5))

        cf = ds.psi_baro.cf.plot.contourf(
            x="longitude", y="latitude", levels=21, cmap=cmo.curl, ax=axe
        )
        cf.colorbar.set_label(r"$\Psi_{baro}$ [Sv]")
        axe.set_title("Barotropic streamfunction")
        axe.set_xlabel(r"$\lambda$ [$^\circ$E]")
        axe.set_ylabel(r"$\varphi$ [$^\circ$N]")
        axe.set_ylim(0, 60)
        axe.set_xlim(0, 40)
        fig.tight_layout()
        fig.savefig(snakemake.output[0])
