# The polar transition from alpha to beta regions set by a surface buoyancy flux inversion.
*Authors: Romain Caneill, Fabien Roquet, Gurvan Madec, Jonas Nycander

Sources for the paper `Caneill, R., Roquet, F., Madec, G., and Nycander, J. (2022). The polar transition from alpha to beta regions
set by a surface buoyancy flux inversion.`, published by the Journal of Physical Oceanography.

DOI: https://doi.org/10.1175/JPO-D-21-0295.1  
URL: https://journals.ametsoc.org/view/journals/phoc/52/8/JPO-D-21-0295.1.xml

This repository along with the zenodo data zip can be used to reproduce all the analyzes of the study.
All the study, from the model runs to their output analyze, is entirely reproducible from a few commands.

Released on Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5747663.svg)](https://doi.org/10.5281/zenodo.5747663)


## Content

* BASIN/  
  This folder contains the sources for the NEMO configuration used in this study
* configuration/  
  Some configuration files for the analyzes
* data/  
  Location for the raw and processed data
* figures/  
  Location for the figures
* lib_position_transition_zone/  
  Python library containing some functions for the data processing
* rules/  
  Snakemake rules
* scripts/  
  Python scripts for the data processing and analyzes
* Pipfile*  
  Environment files for python's pipenv
* setupy.py  
  Installation script for lib_position_transition_zone
* Snakefile  
  Main file used by Snakemake

## Analyze organization and tools

The analyze is based on python scripts. The main used libraries are xarray, xgcm, cf-xarray, xnemogcm, dask, and matplotlib.

Snakemake is used as workflow management tool.
The mentality is that snakemake will process every intermediate step needed to achieve the rule it is asked to compute,
[more to read on the official documentation](https://snakemake.readthedocs.io).

The svg images present in `figures/paper_DAG` and `figures/paper_RULEGRAPH` represent the
directed acyclic graph or rule graph for each figure of the paper. These file are
generated by snakemake (which means that you can regenerate them as any other figure),
but we include them here. They can be used to see all the steps needed to compute each
figure of the paper, in case you would like to look at the code behind.

## Reproduce the analyze

### Clone this directory

```
git clone https://github.com/rcaneill/caneill-et-al-JPO-nemo-transition-zone
cd caneill-et-al-JPO-nemo-transition-zone
```

Or use the release from zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5747663.svg)](https://doi.org/10.5281/zenodo.5747663)

### Install the virtual environment

If you haven't already installed pipenv:
```
pip3 install --user pipenv
```

And then (in the base directory):
```
pipenv install
```

### Download the raw data from zenodo

You first need to download the raw data from zenodo, and unzip the content into the folder `data/raw`.
These data are the output files of all the NEMO runs used in this study.

You need to have zip and unzip installed (e.g. on a debian like linux
`sudo apt install zip unzip`).

You can use the `download_and_extract_data_from_zenodo.sh` script:
```
pipenv install
./download_and_extract_data_from_zenodo.sh
```

### Run the analyze

As mentionned, the analyzes are organized by the use of snakemake.
One can think of it as a python friendly GNU make equivalent.

A first good step is to reproduce the figures used in the paper:
```
pipenv run snakemake --cores 8 paper_all_fig
```
You can adapt the number of cores depending on you computer.
With 8 cores in a modern laptop this takes few minutes to run.

It is then possible to generate other figures, like e.g. the MOC in another configuration:
```
pipenv run snakemake --cores 8 figures/EXP_main/EXP_rn_lambda1_0.06_rn_a0_0.145/moc.pdf
```

To generate all the figures at once (about 5 minutes to run for one experiment):
```
pipenv run snakemake --cores 8 figures/plot_all.done
```

You can also generate the figures giving the figure number
from the paper, e.g.:
```
pipenv run snakemake --cores 8 figures/paper_by_number/figure_1.pdf
```

The rules are described in the files located in rules/, and each rule uses a python script located in scripts/.
The names are self-describing.

## Reproduce the runs

### Compile and run

Follow the [BASIN/README.md](BASIN/README.md) instructions.

### Download data

e.g. `rsync -avz user@supercomputer:/path/to/nemo/NEMOGCM/tests/BASIN_compiled_v0.1.0/EXP_* data/raw/`

### Re-run the analyze

Use snakemake to run the analyze by asking it to produce the desired figures.
