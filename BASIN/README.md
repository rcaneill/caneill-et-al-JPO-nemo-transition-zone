# BASIN configuration

## Reproducible research

This folder (along with the nemo sources) make possible to reproduce all the runs used in this study.

As this is not in the dev / research period but the production period,
all the runs should be able to be re-ran from scratch with a single command / script,
but only for the number of runs needed.
It is however possible to add new experiments in the `EXPERIMENTS_PARAMS` folder,
in a new csv file, following the examples.

### Workflow management
We use [snakemake](https://snakemake.readthedocs.io).

### Specifications
* 1 experiment contains multiple configurations
* In an experiment, between the configs we need to change 1 or 2 parameters of the namelist
* Run the model for 2000 years, with no output. Then run 50 more years with monthly outputs (done in the run_nemo_restart.sh file)
* The timestep must be variable between 1 degree and 2 degrees runs. Same in the run_nemo_restart script (number time step per year)
* Depending on 1 or 2 deg resolution, choose the proper amount of processors and running time

### Folder organisation

Before the run, after compilation and copying the proper folders (see next section),
the folder are organized as follow:
```
BASIN_compiled
|-- MY_SRC    ‾| <= These folders are created when compiling nemo
|-- EXP00      |
|-- WORK       |
|-- BLD       _|
|
|-- EXPERIMENTS_PARAMS       ‾| <= Input files
|   |-- EXP_bathymetry.csv    |
|   |-- EXP_flat_bottom.csv  _|
|
|-- CREATE_RUNS              ‾| <= Scripts
|   |-- Snakefile             |
|   |-- change_nam.sh         |
|   |-- create_exp.sh         |
|   |-- create_exp_config.sh  |
|   |-- EXP_REF              _|
```

Creating the configurations will create these folder (this only is an example,
it of course depends on the experiments parameters located in `EXPERIMENTS_PARAMS`):
```
|-- EXP_bathymetry                          ‾| <= Created automatically
|   |-- experiments.csv                      |    using the input
|   |-- Experiments                          |    csv files in
|       |-- EXP00                            |    EXPERIMENTS_PARAMS
|       |-- EXP_nn_botcase_0                 |
|       |-- EXP_nn_botcase_1                 |
|-- EXP_flat_bottom                          |
|   |-- experiments.csv                      |
|   |-- Experiments                          |
|       |-- EXP00                            |
|       |-- EXP_rn_lambda1_0.06_rn_a0_0.165  |
|       |-- EXP_rn_lambda1_0.04_rn_a0_0.145 _|
```

Here is an example of a csv parameter file:
```
$ cat EXPERIMENTS_PARAMS/EXP_flat_bottom.csv
short_name, nn_botcase
ref, 0
A, 1
```
All the fields from the namelist that we want to change are added as header,
their values is then added for each experiment.
The column `short_name` corresponds to the name of the configuration (used in
the analyze in certain figures).

## How to use

1. This configuration runs with NEMO release 4.0  
   You can download this version using snv  
   `svn co https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 NEMOGCM`
2. Copy the folder BASIN in NEMOGCM/tests/BASIN
2. run the following command in the BASIN folder:  
    `echo "BASIN  OCE" >> ../work_cfgs.txt`
3. That's done, you can now compile, in the NEMO root directory, e.g.:  
    `./makenemo -n 'BASIN_compiled_v0.1.0' -a 'BASIN' -m 'X64_TETRALITH' -j 32`  
    replacing `X64_TETRALITH` by your arch file
4. Copy the folders from BASIN:  
    ```
    cd tests/BASIN_compiled_v0.1.0
    cp -r ../BASIN/CREATE_RUNS .
    cp -r ../BASIN/EXPERIMENTS_PARAMS .
    ```
5. Copy the xios_server.exe into the `CREATE_RUNS/EXP_REF` folder
6. Install the pipenv virtual environment to get snakemake and f90nml:  
    ```
    cd CREATE_RUNS
    pip3 install --user pipenv
    pipenv install
    ```
7. The project is make to be run in an environment using sbatch (common on supercomputers).
  Change your sbatch project name in the config.yaml file (or give `--config sbatch_proj=yourProjName` as parameter
  when executing snakemake).
8. You can relaunch the experiments using snakemake:  
    ```
    pipenv run snakemake --cores 1 launch_all
    ```
    Remark: the number of cores given to snakemake only corresponds to the number of cores to create all the batches.
    The number of cores used for the nemo runs are set in the `create_exp_config.py` file.
    Remark 2: For one configuration at 1 degree resolution, the 2000 years integration roughly corresponds to
    using 1500 hours (1 day and 22 hours of computation with 32 cores). For one configuration at 2 degrees,
    he integration roughly corresponds to 150 hours (16 hours of computation with 9 cores).
