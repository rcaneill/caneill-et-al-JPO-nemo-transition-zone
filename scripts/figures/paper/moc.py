import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == "__main__":
    with sns.axes_style("ticks"):
        ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0))
        fig, ax = plt.subplots(2, 1, figsize=(pc19, 4), sharex=True)

        im = ds.psi_moc.cf.plot.contourf(
            ax=ax[0],
            x="latitude",
            y="depth_W_uniform",
            yincrease=False,
            levels=13,
            vmin=-6,
            vmax=6,
            cmap=cmo.balance,
        )
        ds.mldr10_1.max("x_c").cf.plot(ax=ax[0], x="latitude", color="C0")

        ax[0].set_xlim(0, 60)
        ax[0].set_title("a) Depth coordinates")
        ax[0].set_ylabel("Depth [m]")
        ax[0].set_xlabel("")
        im.colorbar.set_label("[Sv]")

        im = ds.psi_moc_sigma.cf.plot.contourf(
            ax=ax[1],
            x="latitude",
            y="sigma_uniform",
            yincrease=False,
            levels=19,
            cmap=cmo.amp,
            vmin=0,
            vmax=9,
        )
        ds.sigma.where(ds.tmask_surf).isel(z_c=0).max("x_c").cf.plot(
            x="latitude", ax=ax[1]
        )
        ds.sigma.where(ds.tmask_surf).isel(z_c=0).min("x_c").cf.plot(
            x="latitude", ax=ax[1], color="C2"
        )

        ax[1].set_title(r"b) Density coordinates")
        ax[1].set_xlabel(r"Latitude [$^\circ$N]")
        ax[1].set_ylabel(r"$\sigma_0$ [kg m$^{-3}$]")
        ax[1].set_ylim(28, 23.5)
        im.colorbar.set_label("[Sv]")

        fig.tight_layout()
        fig.savefig(snakemake.output[0])
