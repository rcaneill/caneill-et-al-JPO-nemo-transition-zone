import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean


if __name__ == "__main__":
    with sns.axes_style("ticks"):
        ds = time_mean(
            xr.open_mfdataset(snakemake.input["data"]).isel(exp=0), month=False
        )

        fig, ax = plt.subplots(1, 3, figsize=(pc39, 5), sharey=True)

        try:
            mld = ds.mldr10_1.max("month")
        except ValueError:
            mld = ds.mldr10_1
        ds = time_mean(ds)

        # To render nice plots
        mld.load()
        mld[{"x_c": 0}] = mld[{"x_c": 1}]
        mld[{"x_c": -1}] = mld[{"x_c": -2}]
        plot1 = mld.where(ds.tmask.isel(z_c=0)).cf.plot.contourf(
            ax=ax[1],
            cmap=cmo.deep,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            vmin=0,
            vmax=3600,
            levels=19,
        )

        # To render nice plots
        ds.sci_under_ml.load()
        ds.sci_under_ml[{"x_c": 0}] = ds.sci_under_ml[{"x_c": 1}]
        ds.sci_under_ml[{"x_c": -1}] = ds.sci_under_ml[{"x_c": -2}]
        strat = ds.sci_under_ml.cf.plot.contourf(
            ax=ax[2],
            cmap=cm_sci,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            vmin=-2,
            vmax=2,
            levels=41,
        )
        mld_levels = [1800, 3600]
        for axe, c in zip(ax[1:], ["w", "k"]):
            mld.where(ds.tmask.isel(z_c=0)).cf.plot.contour(
                ax=axe, colors=[c], x="longitude", y="latitude", levels=mld_levels
            )
        cb = plot1.colorbar
        ylim = cb.ax.get_ylim()
        for i in mld_levels:
            cb.ax.plot([i, i], ylim, color="w")

        plot2 = strat
        plot3 = np.abs(ds.sci_under_ml).cf.plot.contourf(
            ax=ax[2],
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
            ax=ax[0],
            cmap=cmo.balance,
            x="longitude",
            y="latitude",
            cbar_kwargs=dict(orientation="horizontal"),
            levels=21,
        )

        for axe in ax:
            for lon in snakemake.params["lon_sections"]:
                axe.plot([lon, lon], [0, 60], "k-.")
            axe.set_ylim(0, 60)
            axe.set_xlim(0, 40)
        ax[1].set_title("b) Mixed Layer Depth")
        ax[2].set_title(f"c) SCI under the ML")
        ax[0].set_title(f"a) SSH")

        ax[0].set_ylabel(r"$\varphi$ [$^\circ$N]")
        for axe in ax[1:]:
            axe.set_ylabel("")
        for axe in ax:
            axe.set_xlabel(r"$\lambda$ [$^\circ$E]")

        plot1.colorbar.set_label(r"Depth [m]")
        plot1.colorbar.set_ticks(np.arange(0, 3600 + 1, 600))
        plot2.colorbar.set_label("SCI")
        plot_ssh.colorbar.set_label("SSH [m]")

        fig.tight_layout()
        fig.savefig(snakemake.output[0])
