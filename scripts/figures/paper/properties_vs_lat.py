import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    ds = time_mean(xr.open_mfdataset(snakemake.input))

    fig = plt.figure(figsize=(pc39, 6))
    gs = fig.add_gridspec(4, 2)

    axes = [None,]*5
    axes[0] = fig.add_subplot(gs[0:2, 0])
    axes[1] = fig.add_subplot(gs[0:2, 1], sharey=axes[0])
    axes[2] = fig.add_subplot(gs[2:4, 0])
    axes[3] = fig.add_subplot(gs[2, 1], sharex=axes[1])
    axes[4] = fig.add_subplot(gs[3, 1], sharex=axes[3])

    axe = axes[0]
    ds = ds.sortby('phi_max')
    # ds.plot.scatter("phi_max", "lat_buoyancy_front", ax=axe, label=r'Latitude of buoyancy flux inversion $\varphi_b$')
    # ds.plot.scatter("phi_max", "lat_fresh_pool", ax=axe, label=r'Southern latitude of beta ocean $\varphi_\beta$')
    axe.plot(ds["phi_max"], ds["lat_buoyancy_front"], 'o-', label=r'Latitude of buoyancy flux inversion $\varphi_b$')
    axe.plot(ds["phi_max"], ds["lat_fresh_pool"], 'o-', label=r'Southern latitude of beta ocean $\varphi_\beta$')
    axe.legend()
    axe.set_xlabel(r"$\varphi_M$ [$^\circ$N]")
    axe.set_ylabel("[$^\circ$N]")
    axe.set_title(r"a) $\varphi_\beta$ and $\varphi_b$ versus $\varphi_M$")

    axe = axes[1]
    ds = ds.sortby("lat_buoyancy_front")
    axe.plot(ds["lat_buoyancy_front"], ds["lat_fresh_pool"], 'o-', color='C1')
    #ds.plot.scatter("lat_buoyancy_front", "lat_fresh_pool", ax=axe, color='C1')
    axe.set_xlabel(r"$\varphi_b$ [$^\circ$N]")
    axe.set_ylabel(r"$\varphi_\beta$ [$^\circ$N]")
    axe.set_title(r"b) $\varphi_\beta$ versus $\varphi_b$")

    axe = axes[2]
    ds = ds.sortby("alpha_deep_ML")
    axe.plot(ds["alpha_deep_ML"], ds["lat_buoyancy_front"], 'o-')
    #ds.plot.scatter("alpha_deep_ML", "lat_buoyancy_front", ax=axe)
    axe.set_ylabel(r"$\varphi_b$ [$^\circ$N]")
    axe.set_xlabel(r"TEC in the convective zone [$^\circ$C$^{-1}$]")
    axe.set_title(r"c) $\varphi_b$ versus TEC")
    axe.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
    axe.set_xlim(0.8e-4, 1.1e-4)

    axe = axes[3]
    ds = ds.sortby("lat_buoyancy_front")
    axe.plot(
        ds["lat_buoyancy_front"], ds["mean_T"], 'v-', label=r"$\overline{\Theta}$"
    )
    axe.plot(
        ds["lat_buoyancy_front"], ds["abyss_T"], 'o-', label=r"$\Theta$ abyssal"
    )
    # ds.plot.scatter(
    #     "lat_buoyancy_front", "mean_T", ax=axe, label=r"$\overline{\Theta}$", marker="v"
    # )
    # ds.plot.scatter(
    #     "lat_buoyancy_front", "abyss_T", ax=axe, label=r"$\Theta$ abyssal", marker="o"
    # )
    axe.legend()
    axe.set_ylim(1, 5)
    axe.set_title("d) Thermal properties")
    axe.set_ylabel(r"Temperature [$^\circ$C]")
    #axe.set_xlabel("")
    axe.set_xlabel(r"$\varphi_b$ [$^\circ$N]")

    axe = axes[4]
    ds = ds.sortby("lat_buoyancy_front")
    axe.plot(
        ds["lat_buoyancy_front"], ds["mean_S"], 'v-', label=r"$\overline{S}$"
    )
    axe.plot(
        ds["lat_buoyancy_front"], ds["abyss_S"], 'o-', label=r"$S$ abyssal"
    )
    # ds.plot.scatter(
    #     "lat_buoyancy_front", "mean_S", ax=axe, label=r"$\overline{S}$", marker="v"
    # )
    # ds.plot.scatter(
    #     "lat_buoyancy_front", "abyss_S", ax=axe, label=r"$S$ abyssal", marker="o"
    # )
    axe.legend()
    axe.set_title("e) Haline properties")
    axe.set_ylim(35.7, 36.1)
    axe.set_ylabel(r"Salinity [g$\,$kg$^{-1}$]")
    axe.set_xlabel(r"$\varphi_b$ [$^\circ$N]")

    for axe in axes[1:2] + axes[3:]:
        axe.yaxis.set_label_position("right")
        axe.yaxis.tick_right()
        axe.tick_params(right=False)
    
    fig.tight_layout()
    fig.savefig(snakemake.output[0])
