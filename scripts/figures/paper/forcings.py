import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    with sns.color_palette("dark"):
        ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0), day=False)

        fig, ax = plt.subplots(4, 1, figsize=(pc19, 4.5), sharex=True)
        
        ax[0].fill_between(
            ds.gphit, ds.ztstar.min("day"), ds.ztstar.max("day"), color="silver"
        )
        ds.ztstar.isel({"day": 21 + 90}).cf.plot.line(x="latitude", ax=ax[0])
        
        ds.zsstar.cf.plot(x="latitude", ax=ax[1])
        
        ds.utau.isel({"x_f": 10}).cf.plot(x="latitude", ax=ax[2])
        
        ax[3].fill_between(
            ds.gphit,
            ds.qsolar_analytic.min("day"),
            ds.qsolar_analytic.max("day"),
            color="silver",
        )
        ds.qsr.isel({"x_c": 10}).cf.plot.line(x="latitude", ax=ax[3])
        
        for axe, abcd in zip(ax.flatten(), "abcd"):
            axe.grid(True)
            axe.set_title(f"{abcd})")
            axe.set_xlim(0, 60)
            axe.set_xlabel("")
            
        ax[0].set_ylabel(r"$T^*$ [$^\circ$C]")
        ax[1].set_ylabel(r"$S^*$ [g$\,$kg$^{-1}$]")
        ax[2].set_ylabel(r"$\tau_i$ [Pa]")
        ax[3].set_ylabel(r"$Q_{solar}$ [W$\,$m$^{-2}$]")
        
        ax[-1].set_xlabel(r"$\varphi$ [$^\circ$N]")
        
        fig.tight_layout()
        fig.savefig(snakemake.output[0])
