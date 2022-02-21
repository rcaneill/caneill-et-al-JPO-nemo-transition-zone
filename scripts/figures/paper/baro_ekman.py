import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    with sns.axes_style("ticks"):
        ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0))
        
        fig, ax = plt.subplots(
            1, 2, figsize=(pc33, 3.5), sharey=True, gridspec_kw={"width_ratios": [0.5, 1]}
        )
        
        axe = ax[1]
        cf = ds.psi_baro.cf.plot.contourf(
            x="longitude", y="latitude", levels=21, cmap=cmo.curl, ax=axe
        )
        cf.colorbar.set_label(r"$\Psi_{baro}$ [Sv]")
        axe.set_title("b) Barotropic streamfunction")
        axe.set_xlabel(r"$\lambda$ [$^\circ$E]")
        axe.set_ylabel("")
        axe.set_ylim(0, 60)
        axe.set_xlim(0, 40)
        
        axe = ax[0]
        phi = ds.gphiv
        w_ek = ds.w_ek.where(ds.gphiv > 2)
        axe.fill_betweenx(
            phi, w_ek, where=(w_ek > 0), color="C1", interpolate=True, alpha=0.5
        )
        axe.fill_betweenx(
            phi, w_ek, where=(w_ek < 0), color="C0", interpolate=True, alpha=0.5
        )
        axe.plot(w_ek, phi, color="k", linewidth=1)
        xlim = axe.get_xlim()
        Y = 34.8
        axe.plot(xlim, [Y, Y], "k--")
        axe.set_xlim(xlim)
        axe.text(-4.2e-6, 50, "Upwelling", size="small")
        axe.text(-4.2e-6, 20, "Downwelling", size="small")
        axe.set_title("a) Ekman vertical velocity")
        axe.set_ylabel(r"$\varphi$ [$^\circ$N]")
        axe.set_xlabel(r"Velocity [m$\,$s$^{-1}$]")
        
        fig.tight_layout()
        fig.savefig(snakemake.output[0])
