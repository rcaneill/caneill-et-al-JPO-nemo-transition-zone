import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

sns.axes_style("ticks")

if __name__ == '__main__':
    ds = xr.open_mfdataset(snakemake.input)

    ds.coords["lon_mercatoru"] = (
        np.cos(np.deg2rad(ds.gphit)) * (ds.glamu - 20) + 20
    )
    ds.coords["lon_mercatort"] = (
        np.cos(np.deg2rad(ds.gphit)) * (ds.glamt - 20) + 20
    )

    fig, axes = plt.subplots(
        5, 2, figsize=(pc39, 5.5), sharex='col',
    )

    ##############  Bathymetry ###################
    gs = axes[0,0].get_gridspec()
    # remove the underlying axes
    for ax in axes[:-1, 0]:
        ax.remove()
    axbig = fig.add_subplot(gs[:-1, 0], sharex=axes[-1,0])
    ax = [axbig, axes[-1,0]]
    bathy = ds.gdepw_0.isel({"z_f": -1}) / 1000

    hmax = ds.rn_h.values[0] / 1000
    hmin = ds.rn_hborder.values[0] / 1000

    p = bathy.cf.plot.contourf(
        x="lon_mercatort", y="latitude", levels=11, cmap=cmo.deep_r, ax=ax[0],
        vmin=hmin, vmax=hmax, extend='neither', cbar_kwargs={'location':'bottom'}#add_colorbar=False,
    )
    p.colorbar.set_label("Depth [m]")

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
        axe.set_xticks(np.arange(0,41,10))
    ax[1].set_xlabel(r"$\lambda$ [$^\circ$E]")
    ax[0].set_ylabel(r"$\varphi$ [$^\circ$N]")
    ax[1].set_ylabel('Depth [km]')

    ax[0].set_title('a) Bathymetry map')
    ax[1].set_title('b) Section at $30^\circ$N')

    #cb = fig.colorbar(p, ax=ax)
    #cb.set_label("Depth [km]")


    ###################### Forcings ######################
    ax = axes[:,1].flatten()
    ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0), day=False)
    ax[0].fill_between(
        ds.gphit, ds.ztstar.min("day"), ds.ztstar.max("day"), color="silver"
    )
    ds.ztstar.isel({"day": 21 + 90}).cf.plot.line(x="latitude", ax=ax[0], color='k')
    
    ds.zsstar.cf.plot(x="latitude", ax=ax[1], color='k')
    
    ds.utau.isel({"x_f": 10}).cf.plot(x="latitude", ax=ax[2], color='k')
    
    ax[3].fill_between(
        ds.gphit,
        ds.qsolar_analytic.min("day"),
        ds.qsolar_analytic.max("day"),
        color="silver",
    )
    ds.qsr.isel({"x_c": 10}).cf.plot.line(x="latitude", ax=ax[3], color='k')
    
    ax[4].fill_between(
        ds.gphit, ds.sigma0_star.min("day"), ds.sigma0_star.max("day"), color="silver"
    )
    ds.sigma0_star.isel({"day": 21 + 90}).cf.plot(x="latitude", ax=ax[4], color='k')
    ax[4].plot([ds.phi_max, ds.phi_max], [20,30], color='gray')
    
    for axe, a in zip(ax.flatten(),
        [r'c) $T^*$',
         r'd) $S^*$',
         r'e) $\tau_i$',
         r'f) $Q_{solar}$',
         r'g) $\sigma_0^*$']
    ):
        axe.grid(True)
        axe.set_title(f"{a}")
        axe.set_xlim(0, 60)
        axe.set_xlabel("")
        axe.yaxis.set_label_position("right")
        #axe.yaxis.tick_right()
        #axe.tick_params(right=False)
    ax[1].set_ylim(35,37)
    # ax[0].set_ylabel(r"$T^*$ [$^\circ$C]")
    # ax[1].set_ylabel(r"$S^*$ [g$\,$kg$^{-1}$]")
    # ax[2].set_ylabel(r"$\tau_i$ [Pa]")
    # ax[3].set_ylabel(r"$Q_{solar}$ [W$\,$m$^{-2}$]")
    # ax[4].set_ylabel(r"$\sigma_0^*$ [kg$\,$m$^{-3}$]")
    ax[0].set_ylabel(r"[$^\circ$C]")
    ax[1].set_ylabel(r"[$\,$kg$^{-1}$]")
    ax[2].set_ylabel(r"[Pa]")
    ax[3].set_ylabel(r"[W$\,$m$^{-2}$]")
    ax[4].set_ylabel(r"[kg$\,$m$^{-3}$]")
    
    ax[4].set_ylim(23,28)
    ax[4].set_yticks(np.arange(23,28+1))
    ax[4].set_yticklabels(['23','','25','','27',''])
    
    ax[-1].set_xlabel(r"$\varphi$ [$^\circ$N]")
    

    fig.tight_layout()
    fig.savefig(snakemake.output[0])
