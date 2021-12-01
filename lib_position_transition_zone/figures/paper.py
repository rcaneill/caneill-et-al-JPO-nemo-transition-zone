import os
from pathlib import Path

import cf_xarray as cfxr
import cmocean
import cmocean.cm as cmo
import numpy as np
import seaborn as sns
import xarray as xr

sns.set_style("whitegrid")
sns.set_context("paper")

import matplotlib.pyplot as plt

SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 10

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
}

plt.rcParams.update(tex_fonts)

# Allowed sizes for the figures:
# 19pc, 27pc, 33pc, 39pc
# 1pc = 0.16666666666667 in
conversion = 0.16666666666667
pc19 = 19 * conversion
pc27 = 27 * conversion
pc33 = 33 * conversion
pc39 = 39 * conversion

# Not sure the 2 next lines are needed
plt.rcParams["figure.facecolor"] = "w"
plt.rcParams["figure.dpi"] = 300


palette = "colorblind"
sns.set_palette(palette)


cm_sci = cmocean.tools.crop_by_percent(cmo.balance, 15, which="both", N=None)
