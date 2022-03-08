import calendar
import xarray as xr
import cf_xarray
import matplotlib.pyplot as plt
import cmocean.cm as cmo

if __name__ == "__main__":
    ds = (
        xr.open_mfdataset(snakemake.input)
        .groupby("t.month")
        .map(lambda x: x.mean("t"))
        .isel(exp=0)
    )
    N = len(ds.month)
    fig, ax = plt.subplots(
        1,
        N * 2 + 1,
        figsize=(N * 3, 4),
        gridspec_kw={
            "width_ratios": [
                1,
            ]
            * N
            * 2
            + [
                0.05,
            ]
        },
    )
    fig.subplots_adjust(wspace=0.3)
    for i in range(N):
        sub_ds = ds.isel(month=i)
        im = sub_ds.F_b.cf.plot.contourf(
            levels=11,
            vmax=1.5e-8,
            vmin=-1.5e-8,
            cmap=cmo.balance,
            x="longitude",
            y="latitude",
            ax=ax[i * 2],
            add_colorbar=False,
        )
        sub_ds.F_b.cf.plot.contour(
            levels=[0],
            colors="w",
            x="longitude",
            y="latitude",
            ax=ax[i * 2],
        )
        sub_ds.mldr10_1.cf.plot.contourf(
            levels=21,
            vmax=3800,
            vmin=0,
            cmap=cmo.deep,
            x="longitude",
            y="latitude",
            ax=ax[i * 2 + 1],
            add_colorbar=False,
        )
        ax[i * 2].set_title(f"{calendar.month_name[sub_ds.month.values]}")
        ax[i * 2].set_xlabel("")
        ax[i * 2].set_ylabel("")
        ax[i * 2].set_ylim(0, 60)
        ax[i * 2].set_xlim(0, 40)
        ax[i * 2 + 1].set_title(f"{calendar.month_name[sub_ds.month.values]}")
        ax[i * 2 + 1].set_xlabel("")
        ax[i * 2 + 1].set_ylabel("")
        ax[i * 2 + 1].set_ylim(0, 60)
        ax[i * 2 + 1].set_xlim(0, 40)
        if i != 0:
            ax[i * 2].yaxis.set_ticklabels("")

    cb = fig.colorbar(im, cax=ax[-1])
    cb.set_label(r"Buoyancy flux [m$^2\,$s$^{-2}$]")

    fig.tight_layout()
    fig.savefig(snakemake.output[0])
