import xarray as xr
from lib_position_transition_zone.figures.paper import *

sns.axes_style("ticks")

if __name__ == "__main__":
    ds = xr.open_mfdataset(snakemake.input)

    ds.coords["lon_mercatoru"] = np.cos(np.deg2rad(ds.gphit)) * (ds.glamu - 20) + 20
    ds.coords["lon_mercatort"] = np.cos(np.deg2rad(ds.gphit)) * (ds.glamt - 20) + 20

    fig, axe = plt.subplots(1, 1, figsize=(pc19, 3))

    p = ds.gdepw_0.isel({"z_f": -1}).cf.plot.contourf(
        x="lon_mercatort",
        y="latitude",
        levels=11,
        cmap=cmo.deep_r,
        ax=axe,
        vmin=2000,
        vmax=4000,
        extend="neither",
    )
    p.colorbar.set_label("Depth [m]")

    # We only add the lines for the 1 degree case
    if ds.rn_e1_deg.values == 1.0:
        # We add the lines
        for i in axe.get_yticks():
            axe.plot([0, 40], [i, i], "silver")
        # We hide the horizontal lines out of the plot
        axe.fill_betweenx(ds.gphit, ds.lon_mercatoru.isel(x_f=0), color="w", zorder=2.1)
        axe.fill_betweenx(
            ds.gphit, ds.lon_mercatoru.isel(x_f=40), 50, color="w", zorder=2.1
        )
        for i in range(0, 41, 10):
            ds.lon_mercatoru.isel(x_f=i).cf.plot.line(
                y="latitude", color="silver", zorder=2.2
            )

    axe.set_ylim(0, 60)
    axe.set_xlim(0, 40)

    axe.set_title("")
    axe.set_xlabel(r"$\lambda$ [$^\circ$E]")
    axe.set_ylabel(r"$\varphi$ [$^\circ$N]")

    fig.tight_layout()
    fig.savefig(snakemake.output[0])
