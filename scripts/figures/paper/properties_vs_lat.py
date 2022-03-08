import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    ds = time_mean(xr.open_mfdataset(snakemake.input))

    
    fig, ax = plt.subplots(3,2,figsize=(pc39, 6))
    
    axes = list(ax.flatten())

    axe = axes[0]
    ds = ds.sortby('phi_max')
    axe.plot(ds["phi_max"], ds["lat_buoyancy_front"], 'o-', label=r'Latitude of buoyancy flux inversion $\varphi_b$')
    axe.plot(ds["phi_max"], ds["lat_fresh_pool"], 'o-', label=r'Southern latitude of beta ocean $\varphi_\beta$')
    axe.legend()
    axe.set_xlabel(r"$\varphi_M$ [$^\circ$N]")
    axe.set_ylabel("[$^\circ$N]")
    axe.set_title(r"a) $\varphi_\beta$ and $\varphi_b$ versus $\varphi_M$")

    axe = axes[1]
    ds = ds.sortby("lat_buoyancy_front")
    axe.plot(ds["lat_buoyancy_front"], ds["lat_fresh_pool"], 'o-', color='C2')
    axe.set_xlabel(r"$\varphi_b$ [$^\circ$N]")
    axe.set_ylabel(r"$\varphi_\beta$ [$^\circ$N]")
    axe.set_title(r"b) $\varphi_\beta$ versus $\varphi_b$")

    axe = axes[2]
    ds = ds.sortby("alpha_deep_ML")
    axe.plot(ds["alpha_deep_ML"], ds["lat_buoyancy_front"], 'o-', label=r'$\varphi_b$')
    axe.plot(ds["alpha_deep_ML"], ds["lat_fresh_pool"], 'o-', label=r'$\varphi_\beta$')
    axe.legend()
    axe.set_ylabel(r"[$^\circ$N]")
    axe.set_xlabel(r"TEC in the convective zone [$^\circ$C$^{-1}$]")
    axe.set_title(r"c) $\varphi_\beta$ and $\varphi_b$ versus TEC")
    axe.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
    # if EXP_main
    if len(ds.exp) == 5:
        axe.set_xlim(0.8e-4, 1.1e-4)

    
    axe = axes[3]
    ds = ds.sortby('Cb')
    axe.plot(ds["Cb"], ds["lat_buoyancy_front"], 'o-', label=r'$\varphi_b$')
    axe.plot(ds["Cb"], ds["lat_fresh_pool"], 'o-', label=r'$\varphi_\beta$')
    axe.legend()
    axe.set_ylabel(r"[$^\circ$N]")
    axe.set_xlabel("$C_b$ [$^\circ$C$^{-2}\,$kg$\,$m$^{-3}$]")
    axe.set_title(r"d) $\varphi_\beta$ and $\varphi_b$ versus $C_b$")

    axe = axes[4]
    ds = ds.sortby("lat_buoyancy_front")
    axe.plot(
        ds["lat_buoyancy_front"], ds["mean_T"], 'v-', label=r"$\overline{\Theta}$", color='C3'
    )
    axe.plot(
        ds["lat_buoyancy_front"], ds["abyss_T"], 'o-', label=r"$\Theta$ abyssal", color='C7'
    )
    axe.legend()
    # if EXP_main
    if len(ds.exp) == 5:
        axe.set_ylim(1, 5)
    axe.set_title("e) Thermal properties")
    axe.set_ylabel(r"Temperature [$^\circ$C]")
    axe.set_xlabel(r"$\varphi_b$ [$^\circ$N]")

    axe = axes[5]
    ds = ds.sortby("lat_buoyancy_front")
    axe.plot(
        ds["lat_buoyancy_front"], ds["mean_S"], 'v-', label=r"$\overline{S}$", color='C3'
    )
    axe.plot(
        ds["lat_buoyancy_front"], ds["abyss_S"], 'o-', label=r"$S$ abyssal", color='C7'
    )
    axe.legend()
    axe.set_title("f) Haline properties")
    # if EXP_main
    if len(ds.exp) == 5:
        axe.set_ylim(35.7, 36.1)
    axe.set_ylabel(r"Salinity [g$\,$kg$^{-1}$]")
    axe.set_xlabel(r"$\varphi_b$ [$^\circ$N]")

    for axe in axes[1::2]:
        axe.yaxis.set_label_position("right")
        axe.yaxis.tick_right()
        axe.tick_params(right=False)
    
    fig.tight_layout()
    fig.savefig(snakemake.output[0])
