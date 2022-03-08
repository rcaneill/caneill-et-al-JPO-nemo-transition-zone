import numpy as np
import xarray as xr
import calendar
from lib_position_transition_zone.figures.paper import *


def plot_sec(ds, cm_sci, axe):
    sci_plot = ds.sci.where(ds.mldr10_1 < ds.gdepw_0).plot.contourf(
        ax=axe,
        cmap=cm_sci,
        x="gphit",
        y="gdepw_0",
        vmin=-2,
        vmax=2,
        levels=41,
        yincrease=False,
        add_colorbar=False,
        add_labels=False,
    )
    ds.sci.where(ds.mldr10_1 < ds.gdepw_0).cf.plot.contour(
        ax=axe,
        colors="w",
        x="gphit",
        y="gdepw_0",
        levels=[-1, 1],
        yincrease=False,
        add_colorbar=False,
        add_labels=False,
    )

    plot3 = np.abs(ds.sci).plot.contourf(
        ax=axe,
        x="gphit",
        y="gdepw_0",
        colors="none",
        levels=[-1, 1],
        linewidths=2,
        hatches=["/", None],
        add_colorbar=False,
        add_labels=False,
        yincrease=False,
    )
    for collection in plot3.collections:
        collection.set_edgecolor("w")
        collection.set_linewidth(2)

    ds.mldr10_1.cf.plot(ax=axe, x="gphit", color="k", _labels=False)
    return sci_plot


if __name__ == "__main__":
    with sns.axes_style("ticks"):
        ds = xr.open_mfdataset(snakemake.input).isel(exp=0)
        N = 2 if ("month" in ds) else 1
        fig, axes = plt.subplots(
            N, 2, figsize=(pc39, 5 / 2 * N), sharex=True, sharey=True
        )
        if N == 1:
            axes = np.array([axes])

        if "month" in ds.dims:
            for i, t in enumerate(ds.month):
                for j, x in enumerate(ds.glamt):
                    sci_plot = plot_sec(
                        ds.sel(month=t, glamt=x), cm_sci, axe=axes[i, j]
                    )
        else:
            # yearly averaged data
            for j, x in enumerate(ds.glamt):
                sci_plot = plot_sec(ds.sel(glamt=x), cm_sci, axe=axes[0, j])

        cb = fig.colorbar(
            sci_plot, ax=list(axes.flatten()), location="bottom", fraction=0.05
        )

        plot4 = cb.ax.contourf(
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

        for axe in axes.flatten():
            axe.set_title("")
            axe.set_ylim(4000, 0)

        for i in range(2):
            l1 = ds.glamt.isel(glamt=i).values
            axes[0, i].set_title(f"Section {i+1}, $\lambda = {l1} \,^\circ$E")

        if N == 2:
            axes[0, 1].yaxis.set_label_position("right")
            axes[1, 1].yaxis.set_label_position("right")
            axes[0, 1].set_ylabel(calendar.month_name[ds.month.values[0]])
            axes[1, 1].set_ylabel(calendar.month_name[ds.month.values[1]])

        cb.set_label("SCI")

        if N == 2:
            axes[1, 0].set_xlabel(r"$\varphi$ [$^\circ$N]")
            axes[1, 1].set_xlabel(r"$\varphi$ [$^\circ$N]")
            axes[0, 0].set_ylabel("Depth [m]")
            axes[1, 0].set_ylabel("Depth [m]")
        else:
            axes[0, 0].set_xlabel(r"$\varphi$ [$^\circ$N]")
            axes[0, 1].set_xlabel(r"$\varphi$ [$^\circ$N]")
            axes[0, 0].set_ylabel("Depth [m]")

        fig.savefig(snakemake.output[0])
