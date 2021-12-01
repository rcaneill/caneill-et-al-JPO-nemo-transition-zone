import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    with sns.axes_style("ticks"):
        ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0))
        fig, axe = plt.subplots(1, 1, figsize=(pc19, 2))

        im = ds.psi_moc_sigma.cf.plot.contourf(
            ax=axe,
            x="latitude",
            y="sigma_uniform",
            yincrease=False,
            levels=19,
            cmap=cmo.amp,
            vmin=0,
            vmax=9,
        )
        ds.sigma.where(ds.tmask_surf).isel(z_c=0).max("x_c").cf.plot(
            x="latitude", ax=axe
        )
        ds.sigma.where(ds.tmask_surf).isel(z_c=0).min("x_c").cf.plot(
            x="latitude", ax=axe, color="C2"
        )

        axe.set_title(r"MOC in density coordinates")
        axe.set_xlabel(r"$\varphi$ [$^\circ$N]")
        axe.set_ylabel(r"$\sigma_0$ [kg m$^{-3}$]")
        axe.set_ylim(28, 23.5)
        axe.set_xlim(0, 60)
        im.colorbar.set_label("[Sv]")

        fig.tight_layout()
        fig.savefig(snakemake.output[0])
        
