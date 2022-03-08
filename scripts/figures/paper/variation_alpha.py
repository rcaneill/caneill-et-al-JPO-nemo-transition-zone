import xarray as xr
from lib_position_transition_zone.figures.paper import *

sns.axes_style("ticks")

if __name__ == "__main__":
    ds = xr.open_mfdataset(snakemake.input)

    fig, axe = plt.subplots(1, 1, figsize=(pc19, 3))

    linestyles = ["--", "-.", "solid", "dotted", (0, (5, 0.5))]
    for exp, linestyle in zip(ds.exp, linestyles * int(len(ds.exp) / 5 + 5)):
        expname = ds.short_name.sel(exp=exp).values
        if "ref" in str(expname):
            marker = "o"
        else:
            marker = ""
        ds.alpha.sel(exp=exp).plot(
            ax=axe, label=expname, marker=marker, x="T", linestyle=linestyle
        )
    axe.legend()
    axe.set_title("")
    axe.set_xlabel(r"$\Theta$ [$^\circ$C]")
    axe.set_ylabel(r"$\alpha$ [$^\circ$C$^{-1}$]")
    axe.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    fig.tight_layout()
    fig.savefig(snakemake.output[0])
