from os import listdir
import os.path
from pathlib import Path


DATA = Path("data")
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
ALL_EXP = [i for i in listdir(RAW) if (RAW / i).is_dir()]
FIGURES = Path("figures")
FIGURES_PAPER = FIGURES / "paper"
CONFIGURATION = Path("configuration")

# Organization of the folders
""" For each experiment in data/raw
- data/raw
  |- EXP_blabla
    |- Experiments
      |- EXP_bla1
        |- VARS
          |- BASIN_blabla_grid_X.nc
        |- domcfg_blabla
        |- namelist_cfg
      |- EXP_bla2
        |- ...

- data/processed
  |- EXP_blabla
    |- EXP_bla1
      |- delta_t.txt
      |- xnemogcm.nemo.nc
      |- etc
    |- EXP_bla2
      |- ...
"""

# Naming conventions
"""
An experiment (exp) contains one of multiple configurations (config) that each contain a run.

Experiments are the main folder containing there configurations runs.
The name of an experiment should be explicit of what is done
(e.g. EXP_1_degree_seasonal_variation to study the effect of seasonality)

The configurations are called from the name of the namelist parameters that change
between all the configurations of a same experiment
(e.g. EXP_ln_ann_cyc_.false. and EXP_ln_ann_cyc_.true.)
"""

exp_configs = []
exp_configs_dict = {}

for exp in ALL_EXP:
    configs = [i for i in listdir(RAW / exp / "Experiments") if "EXP_" in i]
    exp_configs_dict[exp] = configs
    for config in configs:
        exp_configs.append(f"{exp}/{config}")


SECTION_PARAM_LON = CONFIGURATION / "section_longitudes.txt"
with open(SECTION_PARAM_LON) as f:
    lon_sections = [float(i) for i in f.read().splitlines()]
SECTION_PARAM_TIME = CONFIGURATION / "section_times.txt"
with open(SECTION_PARAM_TIME) as f:
    month_sections = [int(i) for i in f.read().splitlines()]

##############################
# RULES                      #
##############################


include: "rules/process.skm"
include: "rules/diags.skm"
include: "rules/plots_config.skm"
include: "rules/plots_exp.skm"
include: "rules/paper.skm"
include: "rules/plots_all.skm"
