import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    with sns.axes_style("ticks"):
        ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0))
        
        fig, ax = plt.subplots(1, 3, figsize=(pc39, 5), sharey=True, sharex=True)

        T = ds.thetao.where(ds.tmask).isel({"z_c": 0})
        S = ds.so.where(ds.tmask).isel({"z_c": 0})
        sigma0 = ds.sigma.where(ds.tmask).isel({"z_c": 0})

        # TEMPERATURE
        plot1 = T.cf.plot.contourf(
            ax=ax[0],
            cmap=cmo.thermal,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            levels=29,
        )
        plt.clabel(
            T.cf.plot.contour(
                ax=ax[0],
                x="longitude",
                y="latitude",
                colors="w",
                levels=np.arange(2, 25, 4),
            ),
            inline=1,
            fontsize=10,
            fmt="%i",
        )

        # SALINITY
        plot2 = S.cf.plot.contourf(
            ax=ax[1],
            cmap=cmo.haline,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            levels=18,
            vmin=35.3,
            vmax=37,
        )
        plt.clabel(
            S.cf.plot.contour(
                ax=ax[1],
                x="longitude",
                y="latitude",
                colors="w",
                levels=np.arange(35, 37, 0.4),
            ),
            inline=1,
            fontsize=10,
            fmt="%.1f",
        )

        # DENSITY
        plot3 = sigma0.cf.plot.contourf(
            ax=ax[2],
            cmap=cmo.dense,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            levels=22,
            vmin=23.6,
            vmax=27.8,
        )
        plt.clabel(
            sigma0.cf.plot.contour(
                ax=ax[2],
                x="longitude",
                y="latitude",
                colors="w",
                levels=np.arange(23.4, 27.8, 0.6),
            ),
            inline=1,
            fontsize=10,
            fmt="%.1f",
        )
        ax[2].plot([0,40], [ds.phi_max, ds.phi_max], color='C3')

        ax[0].set_ylim(0, 60)
        ax[0].set_xlim(0, 40)
        ax[0].set_title("a) SST")
        ax[1].set_title("b) SSS")
        ax[2].set_title("c) $\sigma_0$")

        ax[0].set_ylabel(r"$\varphi$ [$^\circ$N]")
        for axe in ax[1:]:
            axe.set_ylabel("")
        for axe in ax:
            axe.set_xlabel(r"$\lambda$ [$^\circ$E]")

        plot1.colorbar.set_label(r"$\Theta$ [$^\circ$C]")
        plot2.colorbar.set_label(r"$S_A$ [g$\,$kg$^{-1}$]")
        plot3.colorbar.set_label(r"$\sigma_0$ [kg$\,$m$^{-3}$]")

        fig.tight_layout()
        fig.savefig(snakemake.output[0])
        
