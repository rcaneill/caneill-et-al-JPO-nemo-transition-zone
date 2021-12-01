import xarray as xr
from lib_position_transition_zone.figures.paper import *

sns.axes_style("ticks")

if __name__ == '__main__':
    ds = xr.open_mfdataset(snakemake.input)

    ds.coords["lon_mercatoru"] = (
        np.cos(np.deg2rad(ds.gphit)) * (ds.glamu - 20) + 20
    )
    ds.coords["lon_mercatort"] = (
        np.cos(np.deg2rad(ds.gphit)) * (ds.glamt - 20) + 20
    )

    fig, ax = plt.subplots(
        2, 1, figsize=(pc19, 4), sharex=True,
        gridspec_kw={"height_ratios":[4, 1]}
    )

    bathy = ds.gdepw_0.isel({"z_f": -1}) / 1000

    hmax = ds.rn_h.values[0] / 1000
    hmin = ds.rn_hborder.values[0] / 1000

    p = bathy.cf.plot.contourf(
        x="lon_mercatort", y="latitude", levels=11, cmap=cmo.deep_r, ax=ax[0],
        vmin=hmin, vmax=hmax, extend='neither', add_colorbar=False,
    )
    #p.colorbar.set_label("Depth [m]")

    # We only add the lines for the 1 degree case
    if ds.rn_e1_deg.values == 1.0:
        # We add the lines
        for i in ax[0].get_yticks():
            ax[0].plot([0, 40], [i, i], "silver")
        # We hide the horizontal lines out of the plot
        ax[0].fill_betweenx(
            ds.gphit, ds.lon_mercatoru.isel(x_f=0), color="w", zorder=2.1
        )
        ax[0].fill_betweenx(
            ds.gphit, ds.lon_mercatoru.isel(x_f=40), 50, color="w", zorder=2.1
        )
        for i in range(0, 41, 10):
            ds.lon_mercatoru.isel(x_f=i).cf.plot.line(
                y="latitude", color="silver", zorder=2.2, ax=ax[0]
            )

    # We plot the section
    # plot of the length scale of the exponential
    #ax[1].plot([0,ds.rn_distlam.values[0]], [hmin, hmax], color='gray')
    #ax[1].plot([40,40-ds.rn_distlam.values[0]], [hmin, hmax], color='gray')
    # we change the y_c coordinate to gphit
    bathy.swap_dims({'y_c':'gphit'}).interp(gphit=30).cf.plot.line(
        x="longitude", ax=ax[1], yincrease=False
    )

    ax[0].set_ylim(0, 60)
    ax[0].set_xlim(0, 40)
    ax[1].set_ylim(hmax+0.02, hmin)
    
    for axe in ax:
        axe.set_xlabel("")
        axe.set_ylabel("")
    ax[1].set_xlabel(r"$\lambda$ [$^\circ$E]")
    ax[0].set_ylabel(r"$\varphi$ [$^\circ$N]")
    ax[1].set_ylabel('Depth [km]')

    ax[0].set_title('a) Bathymetry map')
    ax[1].set_title('b) Section at $30^\circ$N')

    cb = fig.colorbar(p, ax=ax)
    cb.set_label("Depth [km]")

    #fig.tight_layout()
    fig.savefig(snakemake.output[0])
