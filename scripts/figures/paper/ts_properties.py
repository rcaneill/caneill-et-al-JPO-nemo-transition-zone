import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == "__main__":
    ds = time_mean(xr.open_mfdataset(snakemake.input))

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(pc19, 3))
    ax[0].grid(True)
    ax[1].grid(True)

    ds.plot.scatter(
        "lat_buoyancy_front",
        "mean_T",
        ax=ax[0],
        label=r"$\overline{\Theta}$",
        marker="v",
    )
    ds.plot.scatter(
        "lat_buoyancy_front", "abyss_T", ax=ax[0], label=r"$\Theta$ abyssal", marker="o"
    )
    ax[0].legend()

    ds.plot.scatter(
        "lat_buoyancy_front", "mean_S", ax=ax[1], label=r"$\overline{S}$", marker="v"
    )
    ds.plot.scatter(
        "lat_buoyancy_front", "abyss_S", ax=ax[1], label=r"$S$ abyssal", marker="o"
    )
    ax[1].legend()
    ax[0].set_ylim(1, 5)

    ax[0].set_title("a) Thermal properties")
    ax[1].set_title("b) Haline properties")
    ax[1].set_ylim(35.7, 36.1)

    ax[0].set_ylabel(r"Temperature [$^\circ$C]")
    ax[1].set_ylabel(r"Salinity [g$\,$kg$^{-1}$]")
    ax[0].set_xlabel("")
    ax[1].set_xlabel(r"Latitude of the buoyancy inversion [$^\circ$N]")

    fig.tight_layout()
    fig.savefig(snakemake.output[0])
