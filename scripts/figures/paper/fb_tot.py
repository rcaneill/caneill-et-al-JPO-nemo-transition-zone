import xarray as xr
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean

if __name__ == '__main__':
    with sns.axes_style("ticks"):
        ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0))
        fig, ax = plt.subplots(
            1,
            4,
            figsize=(pc39, 4),
            sharex=False,
            sharey=False,
            gridspec_kw={"width_ratios": [1, 1, 1, 0.05]},
        )

        fig.subplots_adjust(wspace=0.3)

        ims = []
        f = [ds.F_b_heat, ds.F_b_salt, ds.F_b]

        for i in range(3):
            im = f[i].plot.contourf(
                levels=11,
                vmax=1.5e-8,
                vmin=-1.5e-8,
                cmap=cmo.balance,
                x="glamt",
                y="gphit",
                ax=ax[i],
                add_colorbar=False,
            )
            f[i].plot.contour(
                levels=np.linspace(-9e-8, 9e-8, 13),
                colors="gray",
                x="glamt",
                y="gphit",
                ax=ax[i],
                add_colorbar=False,
            )
            f[i].plot.contour(
                levels=[-1, 0],
                colors="w",
                x="glamt",
                y="gphit",
                ax=ax[i],
                add_colorbar=False,
                linewidths=2,
            )
            ims += [im]
            ax[i].set_title("")
            ax[i].set_xlabel("")
            ax[i].set_ylabel("")
            ax[i].set_ylim(0, 60)
            ax[i].set_xlim(0, 40)
            if i != 0:
                ax[i].yaxis.set_ticklabels("")

        ax[0].set_ylabel(r"$\varphi$ [$^\circ$N]")
        fig.supxlabel(r"$\lambda$ [$^\circ$E]")

        ax[0].set_title(r"a) Heat contribution")
        ax[1].set_title(r"b) Salt contribution")
        ax[2].set_title(r"c) Total")
        # fig.suptitle('Buoyancy flux, and its heat and salinity components')

        # fig.colorbar(im, ax=ax.ravel().tolist(), orientation='vertical')
        cb = fig.colorbar(ims[0], cax=ax[-1])
        cb.set_label(r"Buoyancy flux [m$^2\,$s$^{-3}$]")
        fig.tight_layout()
        fig.savefig(snakemake.output[0])
        
