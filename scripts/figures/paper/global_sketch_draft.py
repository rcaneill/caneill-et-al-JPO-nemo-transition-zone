import xarray as xr
import xgcm
from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean
from lib_position_transition_zone import io

xr.set_options(keep_attrs=True)

if __name__ == "__main__":
    with sns.axes_style("ticks"):
        ds = time_mean(xr.open_mfdataset(snakemake.input).isel(exp=0))
        metrics = {
            ("X",): ["e1t", "e1u", "e1v", "e1f"],
            ("Y",): ["e2t", "e2u", "e2v", "e2f"],
            ("Z",): ["e3t", "e3u", "e3v", "e3w"],
        }
        grid = xgcm.Grid(ds, periodic=False, metrics=metrics, boundary="extend")
        fig, ax = plt.subplots(
            4,
            1,
            figsize=(5, 5),
            gridspec_kw={"height_ratios": [0.3, 1, 0.3, 0.3]},
            sharex=True,
        )

        X = np.arange(5) + 33
        # X = np.arange(20) + 10

        for f in ["F_b", "F_b_heat", "F_b_salt"]:
            ds[f].where(ds.tmask_surf).cf.isel(X=X).cf.mean("X").cf.plot.line(
                x="latitude", ax=ax[0]
            )
        ax[0].set_xlim(20, 60)
        ax[0].set_ylim([-2.5e-8, 1e-8])

        ds.sci.where(ds.tmask_surf).cf.isel(X=X).cf.mean("X").cf.plot.contourf(
            x="latitude",
            y="gdepw_1d",
            yincrease=False,
            cmap=cm_sci,
            vmin=-1,
            vmax=1,
            levels=2,
            add_colorbar=False,
            ax=ax[1],
        )
        ds.mldr10_1.cf.isel(X=X).cf.mean("X").cf.plot.line(x="latitude", ax=ax[1])
        # ds.sigma.where(ds.tmask).cf.isel(X=X).cf.plot.contour(
        #    x="latitude", y='gdept_0', yincrease=False, cmap=cmo.dense, levels=1001
        # )
        ds.thetao.where(ds.tmask).cf.isel(X=X).cf.mean("X").cf.plot.contour(
            x="latitude",
            y="gdept_1d",
            yincrease=False,
            cmap=cmo.thermal,
            levels=np.linspace(0, 20, 21),
            ax=ax[1],
        )
        ax[1].set_ylim(2000, 0)

        ds.thetao.where(ds.tmask).isel(z_c=0).cf.isel(X=X).cf.mean("X").cf.plot.line(
            x="latitude", ax=ax[2]
        )
        ds.ztstar.cf.plot.line(x="latitude", ax=ax[2], color="gray")

        for i in ["qt", "qsr", "qns"]:
            ds[i].cf.isel(X=X).cf.mean("X").cf.plot.line(
                x="latitude", ax=ax[3], label=i
            )
        ax[3].legend()

        for axe in list(ax):
            axe.set_title("")
            axe.set_xlabel("")
            axe.set_ylabel("")
            axe.grid()

        fig.tight_layout()
        fig.savefig(snakemake.output[0])
