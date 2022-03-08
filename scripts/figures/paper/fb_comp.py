import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

"""
The figure is made to look good for 5 experiences
"""


if __name__ == '__main__':
    ds = time_mean(xr.open_mfdataset(snakemake.input))
    ds['F_b'] *= 1e8
    N = len(ds.exp)
    with sns.axes_style("ticks"):
        fig, ax = plt.subplots(
            1,
            N+1,
            figsize=(pc39*N/5, 3.5),
            sharex=False,
            sharey=False,
            gridspec_kw={"width_ratios": [1,]*N + [0.1,]},
        )
        fig.subplots_adjust(wspace=0.2)

        ds = ds.sortby(ds.alpha_deep_ML)
        ims = []

        for i in range(N):
            sub_ds = ds.isel(exp=i)
            #print(f'{i} / {N} : {sub_ds.exp.values}')
            im = (
                sub_ds.where(sub_ds.tmask_surf)
                .F_b.cf.plot.contourf(
                    levels=11,
                    vmax=1.5,
                    vmin=-1.5,
                    cmap=cmo.balance,
                    x="longitude",
                    y="latitude",
                    ax=ax[i],
                    add_colorbar=False,
                )
            )
            sub_ds.where(sub_ds.tmask_surf).F_b.cf.plot.contour(
                levels=np.linspace(-9, 9, 13),
                colors="gray",
                x="longitude",
                y="latitude",
                ax=ax[i],
                add_colorbar=False,
            )
            sub_ds.where(sub_ds.tmask_surf).F_b.cf.plot.contour(
                levels=[0],
                colors="w",
                x="longitude",
                y="latitude",
                ax=ax[i],
                add_colorbar=False,
                linewidths=2,
            )
            ims += [im]
            sci_ml = sub_ds.where(sub_ds.tmask_surf).sci_under_ml
            sci_ml = xr.where(sub_ds.sci_under_ml <= -1, -2, 0).where(sub_ds.tmask_surf)
            if sci_ml.min() <= -1:
                sci_ml.cf.plot.contour(
                    levels=[-1],
                    colors="k",
                    x="longitude",
                    y="latitude",
                    ax=ax[i],
                    linewidths=3,
                    linestyles="solid",
                )
                sci_ml.cf.plot.contour(
                    levels=[-1],
                    colors="w",
                    x="longitude",
                    y="latitude",
                    ax=ax[i],
                    linewidths=1,
                    linestyles="solid",
                )
            name = sub_ds.short_name.values
            ax[i].plot([0,40], [sub_ds.phi_max, sub_ds.phi_max], color='C3')
            ax[i].set_title(name)
            ax[i].set_xlabel("")
            ax[i].set_ylabel("")
            ax[i].set_yticks(np.arange(40, 61, 5))
            ax[i].set_ylim(40, 60)
            ax[i].set_xlim(0, 40)
            if i != 0:
                ax[i].yaxis.set_ticklabels("")

        ax[0].set_ylabel(r"$\varphi$ [$^\circ$N]")
        fig.supxlabel(r"$\lambda$ [$^\circ$E]")

        fig.suptitle(
            r"$\longleftarrow$ Smaller polar TEC"
            + "~" * 20
            + "Larger polar TEC $\longrightarrow$"
        )

        cb = fig.colorbar(ims[0], cax=ax[-1])
        cb.set_label(r"Buoyancy flux [$\times 10^{-8}\,$m$^2\,$s$^{-3}$]")
        fig.savefig(snakemake.output[0])
