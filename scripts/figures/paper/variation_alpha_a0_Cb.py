import xarray as xr
from lib_position_transition_zone.figures.paper import *

sns.axes_style("ticks")

if __name__ == "__main__":
    ds = xr.open_mfdataset(snakemake.input)

    fig, ax = plt.subplots(1, 2, figsize=(pc33, 2.5))

    # Plot variation alpha
    axe = ax[1]
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
    axe.set_title(r"b) $\alpha_{surface} = f(\Theta)$")

    # scatter plot a0 vs Cb
    axe = ax[0]
    if len(ds.exp) == 5:
        markers = ["x", "v", "o", "*", "^"]
    else:
        markers = [
            "o",
        ] * len(ds.exp)
    for exp, m in zip(ds.exp, markers):
        expname = ds.short_name.sel(exp=exp).values
        axe.scatter(
            ds.sel(exp=exp)["rn_a0"], ds.sel(exp=exp)["Cb"], label=expname, marker=m
        )
    axe.legend()
    axe.ticklabel_format(axis="x", style="sci", scilimits=(-1, 1))
    axe.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3))
    axe.set_xlabel(r"$a_0$ [$^\circ$C$^{-1}\,$kg$\,$m$^{-3}$]")
    axe.set_ylabel(r"$C_b$ [$^\circ$C$^{-2}\,$kg$\,$m$^{-3}$]")
    axe.set_title("a) $C_b$ vs $a_0$")

    fig.tight_layout()
    fig.savefig(snakemake.output[0])
