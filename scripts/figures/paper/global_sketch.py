from lib_position_transition_zone.figures.paper import *
from lib_position_transition_zone.tools import time_mean
import numpy as np

red_alpha = "#8c3e52"


def isotherms_func(mld, x):
    T = np.array([18, 14.5, 10, 5]).reshape(4, 1)
    mld = mld.reshape(1, mld.shape[0])

    v0 = 200 - 2 * T
    v1 = 0.1
    v2 = 45 - 1.5 * T
    v3 = -1250 + 50 * T

    isotherm = v0 * (1 - np.exp(-v1 * (-x.reshape(1, x.shape[0]) + v2))) ** 4 + v3
    isotherm0 = 200 * (1 - np.exp(0.3 * (-x + 54))) ** 4 - 300

    isotherm = -np.append(isotherm, isotherm0.reshape(1, x.shape[0]), axis=0)
    isotherm[isotherm < mld.reshape(1, x.shape[0])] = 0
    return isotherm


def plot_sci_bezier(axe):
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches

    # alpha ocean
    Path = mpath.Path
    verts = np.array([(20, 0), (60, 0), (60, 2000), (20, 2000), (0, 0)])
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor=red_alpha, alpha=1, zorder=0)
    axe.add_patch(patch)
    # The vertices points come from inkscape export to latex
    # The inkscape size of the figure was set to 400 x 200
    # hence the scaling needed to convert back to lat / depth coordinates
    # transition zone
    Path = mpath.Path
    path_data = np.array(
        [
            (400.39006999, 199.12410169),
            (400.28646999, 166.47676917),
            (400.50326999, 154.58939614),
            (400.50326999, 122.06260362),
            (399.84028334, 123.98463823),
            (384.26537706, 164.30952923),
            (360.32899099, 140.57862315),
            (343.8170314, 124.20841289),
            (308.62384562, 43.74657091),
            (296.71208592, 56.73808791),
            (287.25889949, 67.04816366),
            (304.25603239, 168.13546246),
            (306.95547233, 183.8603554),
            (309.60603226, 199.30050168),
            (310.12445891, 196.64411508),
            (310.11695225, 204.34724822),
            (313.1602855, 210.21987474),
            (396.16572343, 199.19924835),
            (400.39006999, 199.12410169),
            (400.39006999, 199.12410169),
        ]
    )
    path_data[:, 0] /= 10
    path_data[:, 0] += 20
    path_data[:, 1] *= -10
    path_data[:, 1] += 2000
    verts = path_data
    codes = [Path.MOVETO] + [Path.CURVE4] * (path_data.shape[0] - 2) + [Path.CLOSEPOLY]
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor="#f3efee", alpha=1, zorder=0)
    axe.add_patch(patch)

    # beta ocean
    Path = mpath.Path
    path_data = np.array(
        [
            (400.11428333, 198.56003504),
            (400.11428333, 122.43323294),
            (394.6503768, 135.75424994),
            (391.72213687, 150.35928957),
            (374.49719064, 150.98631623),
            (349.16429794, 151.90848954),
            (330.64708507, 110.68723323),
            (318.47992537, 119.17657169),
            (313.62296549, 122.56539694),
            (312.88011218, 154.24626281),
            (313.76212549, 174.53567564),
            (314.15096548, 183.48024875),
            (315.26457879, 196.22870176),
            (315.29349878, 199.10400836),
            (315.35389878, 205.11382154),
            (332.09608503, 198.34524837),
            (354.34612447, 198.98387503),
            (369.49365743, 199.41866168),
            (387.14348365, 198.5472217),
            (400.11428333, 198.56000837),
            (400.11428333, 198.56003504),
        ]
    )
    path_data[:, 0] /= 10
    path_data[:, 0] += 20
    path_data[:, 1] *= -10
    path_data[:, 1] += 2000
    verts = path_data
    codes = (
        [Path.MOVETO, Path.LINETO]
        + [Path.CURVE4] * (path_data.shape[0] - 3)
        + [Path.CLOSEPOLY]
    )
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor="#505c99", alpha=1, zorder=0)
    axe.add_patch(patch)


terms = [
    -4.0000001032503008e0,
    8.7828572927567916e-1,
    -6.9645238928795206e-2,
    2.4061904984305315e-3,
    -3.7047619333308569e-5,
    2.0952381094281399e-7,
]


def f_b_salt_func(x):
    result = 0 * x
    for n, t in enumerate(terms):
        result += t * x**n
    return result


def f_b_theta_func(x):
    b = 45
    d = 50.6
    c = -2.4
    a1 = 0.017
    a2 = 0.05
    f = 2 * a2 * (d - b)

    r1 = (x < b) * (a1 * (x - b) ** 2 + c)
    r2 = ((b <= x) & (x < d)) * (a2 * (x - b) ** 2 + c)
    r3 = (d <= x) * (-f * (np.exp(-(x - d)) - 1) + (a2 * (d - b) ** 2 + c))
    return r1 + r2 + r3


def mld_func(x):
    h1 = 51
    h2 = 1200
    h3 = 60

    r1 = (x < h1) * (h3 + h2 * np.exp(-(((x - h1) / 9) ** 2)))
    r2 = (x > h1) * (h3 + h2 * np.exp(-(((x - h1) / 0.8) ** 2)))
    return r1 + r2


def remove_ticks(axe):
    axe.set_yticks([])
    axe.set_xticks([20, x1, x2, 60])
    axe.set_xticklabels("")


def curly_arrow(x, axe, up=True, n=5, col="gray", linew=0.8, width=0.1):
    # adapted from https://stackoverflow.com/questions/45365158/matplotlib-wavy-arrow
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches
    import matplotlib as mpl
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch

    if up:
        ymin = 0
        ymax = 0.8
    else:
        ymin = 1
        ymax = 0.2
    dist = np.sqrt((ymin - ymax) ** 2)
    n0 = dist / (2 * np.pi)

    y_plot = np.linspace(ymin, ymax, 50)
    x_plot = width * np.sin(n * y_plot / n0) + x
    axe.plot(x_plot, y_plot, color=col, linewidth=linew)

    y_tail = 0.2
    y_head = 0
    if up:
        y_tail = 1 - y_tail
        y_head = 1
    codes = [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]
    l = 0.3
    vertices = [(x - l, y_tail), (x + l, y_tail), (x, y_head), (0, 0)]
    path = Path(vertices, codes)
    pathpatch = PathPatch(path, facecolor=col, edgecolor=col)
    axe.add_patch(pathpatch)


if __name__ == "__main__":
    x = np.linspace(20, 60, 1000)
    f_b_theta = f_b_theta_func(x)
    f_b_salt = f_b_salt_func(x)
    f_b = f_b_theta + f_b_salt

    (x1, x2) = (min(x[f_b < 0]), max(x[f_b < 0]))

    fig, ax = plt.subplots(
        3,
        1,
        figsize=(pc39, 4),
        gridspec_kw={"height_ratios": [0.3, 0.08, 1]},
        sharex=True,
    )
    ax = [ax[0], ax[2], ax[1]]

    # Upper plot with buoyancy fluxes
    ax[0].plot(x, f_b_salt)
    ax[0].plot(x, f_b_theta)
    ax[0].plot(x, f_b)
    ax[0].set_ylim(-2.5, 1)
    ax[0].set_xlim(20, 60)
    ax[0].set_yticks([])

    # Lower plot with alpha -- beta -- transition
    ax[1].set_xlabel("")
    ax[1].set_ylabel("")
    ax[1].set_title("")
    # sci
    plot_sci_bezier(ax[1])
    # mld
    mld = mld_func(x)
    ax[1].plot(x, mld, color="k")
    ax[1].text(28, 140, "MLD", size="x-large", color="k")
    # isotherms
    isotherm = isotherms_func(mld, x)
    ax[1].plot(x, isotherm.T, linestyle="solid", color="silver", linewidth=1)
    for i, y in enumerate(isotherms_func(np.array([0, 0]), np.array([22, 58]))):
        _x = 21 if (i != 4) else 58
        _y = y[0] if (i != 4) else y[1] - 10
        ax[1].text(_x, _y - 40, f"$\Theta_{i}$", size="x-large", color="silver")

    # yticks
    ax[1].set_ylim(2000, 20)
    ax[1].set_yticks([50, 1700], minor=True)
    ax[1].set_yticklabels(["surface", "mid-depth"], minor=True)
    ax[1].tick_params(axis="y", rotation=90, which="both")
    # xticks
    ax[1].set_xticks([23, 56], minor=True)
    # ax[1].set_xticklabels([u'20$^\circ$N',u'60$^\circ$N'], minor=True)
    ax[1].set_xticklabels(["subtropics", "(sub)polar region"], minor=True)
    # remove extra ticks
    remove_ticks(ax[1])

    # text alpha beta transition
    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    ax[1].text(37, 800, "Alpha", size="xx-large", bbox=bbox_props)
    ax[1].text(53.5, 200, "Beta", size="xx-large", bbox=bbox_props)
    ax[1].text(52, 1000, "Transition", size="xx-large", bbox=bbox_props)

    # text buoyancy fluxes
    ax[0].text(60, f_b_salt[-1], "$\mathcal{B}^S$", size="xx-large", color="C0")
    ax[0].text(60, f_b_theta[-1] - 0.2, "$\mathcal{B}^\Theta$", size="xx-large", color="C1")
    ax[0].text(40, -1.3, "$\mathcal{B}$", size="xx-large", color="C2")

    # Arrow buoyancy fluxes
    ax[2].axis("off")
    ax[2].set_yticks([])
    ax[2].set_ylim(-0.15, 1.15)
    l = 4
    for (x, up) in zip(
        [
            20 + (x1 - 20) / l,
            20 + (x1 - 20) / l * (l - 1),
            x1 + (x2 - x1) / l,
            x1 + (x2 - x1) / l * (l - 1),
            x2 + (60 - x2) / l,
            x2 + (60 - x2) / l * (l - 1),
        ],
        [False, False, True, True, False, False],
    ):
        curly_arrow(x, up=up, n=3, axe=ax[2], col="C2")
    for (x, i) in zip([(20 + x1) / 2, (x1 + x2) / 2, (x2 + 60) / 2], "><>"):
        ax[2].text(x, 0.2, f"$\mathcal{{B}}{i}0$", ha="center", size='xx-large')

    # Upper plot background
    xlim = ax[0].get_xlim()
    ylim = ax[0].get_ylim()
    ax[0].fill_between(
        xlim, [ylim[1]] * 2, 0, color=ax[0].get_xgridlines()[0].get_color(), zorder=0
    )
    ax[0].set_yticks([-1, 0.8], minor=True)
    ax[0].set_yticklabels(["loss", "gain"], minor=True)
    ax[0].tick_params(axis="y", rotation=90, which="both")
    # ax[0].text(21, -1, 'buoyancy loss', size="large", color='k')
    # ax[0].text(21, .5, 'buoyancy gain', size="large", color='k')

    # Upper plot phi_b
    #ax[2].text(x2 - 0.5, 0.8, r"$\varphi_b$", size="large")

    # Intermediate water creation
    ax[1].annotate(
        "",
        xy=(40, 1500),
        xycoords="data",
        xytext=(50, 1000),
        textcoords="data",
        arrowprops=dict(
            facecolor="w",
            edgecolor="k",
            shrinkB=5,
            connectionstyle="arc3,rad=-0.2",
            width=4,
        ),
    )
    ax[1].text(38, 1200, "Intermediate water", color="w", size="xx-large")

    # Surface adv
    ax[1].annotate(
        "",
        xy=(50, 100),
        xycoords="data",
        xytext=(40, 100),
        textcoords="data",
        arrowprops=dict(facecolor="w", edgecolor="k", shrinkB=5, width=2),
    )
    ax[1].text(40.8, 210, "Surface advection", color="w", size="xx-large")

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.0)
    fig.savefig(snakemake.output[0])
