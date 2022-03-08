import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean


if __name__ == "__main__":
    with sns.axes_style("ticks"):
        ds = time_mean(xr.open_mfdataset(snakemake.input["data"]).isel(exp=0))

        fig, ax = plt.subplots(2, 3, figsize=(pc39, 6.1), sharey=True, sharex=False)

        T = ds.thetao.where(ds.tmask).isel({"z_c": 0})
        S = ds.so.where(ds.tmask).isel({"z_c": 0})
        sigma0 = ds.sigma.where(ds.tmask).isel({"z_c": 0})

        # TEMPERATURE
        plot1 = T.cf.plot.contourf(
            ax=ax[0, 0],
            cmap=cmo.thermal,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            levels=29,
        )
        plt.clabel(
            T.cf.plot.contour(
                ax=ax[0, 0],
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
            ax=ax[0, 1],
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
                ax=ax[0, 1],
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
            ax=ax[0, 2],
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
                ax=ax[0, 2],
                x="longitude",
                y="latitude",
                colors="w",
                levels=np.arange(23.4, 27.8, 0.6),
            ),
            inline=1,
            fontsize=10,
            fmt="%.1f",
        )
        ax[0, 2].plot([0, 40], [ds.phi_max, ds.phi_max], color="C3")

        plot1.colorbar.set_label(r"$\Theta$ [$^\circ$C]")
        plot2.colorbar.set_label(r"$S_A$ [g$\,$kg$^{-1}$]")
        plot3.colorbar.set_label(r"$\sigma_0$ [kg$\,$m$^{-3}$]")

        ############ Lower plot #################
        ds = time_mean(
            xr.open_mfdataset(snakemake.input["data"]).isel(exp=0), month=False
        )
        try:
            mld = ds.mldr10_1.max("month")
        except ValueError:
            mld = ds.mldr10_1
        ds = time_mean(ds)

        plot1 = mld.where(ds.tmask.isel(z_c=0)).cf.plot.contourf(
            ax=ax[1, 1],
            cmap=cmo.deep,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            vmin=0,
            vmax=3600,
            levels=19,
        )

        strat = ds.sci_under_ml.cf.plot.contourf(
            ax=ax[1, 2],
            cmap=cm_sci,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            vmin=-2,
            vmax=2,
            levels=41,
        )
        mld_levels = [1800, 3600]
        for axe, c in zip(ax[1, 1:], ["w", "k"]):
            mld.where(ds.tmask.isel(z_c=0)).cf.plot.contour(
                ax=axe, colors=[c], x="longitude", y="latitude", levels=mld_levels
            )
        cb = plot1.colorbar
        ylim = cb.ax.get_ylim()
        for i in mld_levels:
            cb.ax.plot([i, i], ylim, color="w")

        plot2 = strat
        plot3 = np.abs(ds.sci_under_ml).cf.plot.contourf(
            ax=ax[1, 2],
            x="longitude",
            y="latitude",
            colors="none",
            levels=[-1, 1],
            linewidths=2,
            hatches=["/", None],
            add_colorbar=False,
        )
        for collection in plot3.collections:
            collection.set_edgecolor("w")
            collection.set_linewidth(2)
        strat.colorbar.ax.vlines(
            [-1, 1], -10, 10, colors="w", linestyle="-", linewidth=2
        )
        plot4 = strat.colorbar.ax.contourf(
            np.linspace(-1, 1, 10),
            np.linspace(-10, 10, 10),
            np.meshgrid(np.linspace(-1, 1, 10), np.linspace(-10, 10, 10))[0],
            colors="none",
            levels=[-1, 1],
            hatches=["/", None],
        )
        for collection in plot4.collections:
            collection.set_edgecolor("w")
            collection.set_linewidth(2)

        plot_ssh = ds.zos.cf.plot.contourf(
            ax=ax[1, 0],
            cmap=cmo.balance,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            levels=21,
        )

        for axe in ax[1]:
            for lon in snakemake.params["lon_sections"]:
                axe.plot([lon, lon], [0, 60], "k-.")

        for axe in ax.flatten():
            axe.set_ylabel("")
            axe.set_xlabel("")
            axe.set_ylim(0, 60)
            axe.set_xlim(0, 40)
        ax[0, 0].set_title("a) SST")
        ax[0, 1].set_title("b) SSS")
        ax[0, 2].set_title("c) $\sigma_0$")
        ax[1, 1].set_title("e) Mixed Layer Depth")
        ax[1, 2].set_title(f"f) SCI under the ML")
        ax[1, 0].set_title(f"d) SSH")

        for axe in ax[:, 0]:
            axe.set_ylabel(r"$\varphi$ [$^\circ$N]")

        for axe in ax.flatten():
            axe.set_xlabel(r"$\lambda$ [$^\circ$E]")

        plot1.colorbar.set_label(r"Depth [m]")
        plot1.colorbar.set_ticks(np.arange(0, 3600 + 1, 600))
        plot2.colorbar.set_label("SCI")
        plot_ssh.colorbar.set_label("SSH [m]")

        fig.tight_layout()
        fig.savefig(snakemake.output[0])
