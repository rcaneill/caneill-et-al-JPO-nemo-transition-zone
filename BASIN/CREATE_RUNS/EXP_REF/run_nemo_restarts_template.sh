#!/bin/bash
#
#SBATCH -J BASIN
#SBATCH -t {TIME}
#SBATCH --tasks={N_PROCS_TOT}
#SBATCH -A {SBATCH_PROJ}
#

# Loop on the experiment years
yearInit=0
nYear=2000
# for 2 degrees 5760
# for 1 degree 11520
nTimestepPerYear={N_TIME_STEP_PER_YEAR}
nYearPerLoop=50

echo '' > years.log

nrest=''

for((year = $yearInit ; year <= $nYear ; year=$(( $year+$nYearPerLoop )) ))
do
    # set variables
    tinit=$(( 1 + $year * $nTimestepPerYear ))
    tend=$(( $tinit + $nTimestepPerYear * $nYearPerLoop - 1 ))
    echo tinit $tinit
    echo tend $tend
    # Change namelist
    # restart file or not
    if [ $year -eq 0 ]
    then
	sed -e "s/ln_rstart\s*=.*\!/ln_rstart = .false. \!/" namelist_cfg > namelist_cfg.tmp
	mv namelist_cfg.tmp namelist_cfg
    else
	nrest=$(printf "%08d" $(( $tinit - 1 )) )
	sed -e "s/ln_rstart\s*=.*\!/ln_rstart = .true. \!/" namelist_cfg > namelist_cfg.tmp
	mv namelist_cfg.tmp namelist_cfg
	sed -e "s/cn_ocerst_in\s*=.*\!/cn_ocerst_in = \"BASIN_${{nrest}}_restart\" \!/" namelist_cfg > namelist_cfg.tmp
	mv namelist_cfg.tmp namelist_cfg
	sed -e "s/nn_rstctl\s*=.*\!/nn_rstctl   =    2 \!/" namelist_cfg > namelist_cfg.tmp
	mv namelist_cfg.tmp namelist_cfg
    fi
    # time step
    sed -e "s/nn_it000\s*=.*\!/nn_it000 = $tinit \!/" -e "s/nn_itend\s*=.*\!/nn_itend = $tend \!/" namelist_cfg > namelist_cfg.tmp
    mv namelist_cfg.tmp namelist_cfg
    # If last year, we need to change the xios file-def to monthly outpus
    if [ $year -eq $nYear ]
    then
	mv file_def_nemo-oce.monthly.xml file_def_nemo-oce.xml
    fi
    # log the year
    echo 'starting year ' $year >> years.log
    # run NEMO
    mpprun -np {N_PROCS_NEMO} ./nemo : -n {N_PROCS_XIOS} ./xios_server.exe
    if grep 'E R R O R' ocean.output
    then
	echo "errors in ocean.output => stopping jobs"
	break
    else
	echo "no error in ocean.output => removing/saving previous restart files and continue"
    fi
done

# clean the run from useless files, removes nemo, xios_server, etc
rm output.init*
rm slope_ratio*
rm xios_server.exe
rm nemo
rm time.step
rm slurm*
rm communication_report.txt
# clean all restart files except the last one
mkdir RESTARTS
mv BASIN_${{nrest}}_restart* RESTARTS
rm BASIN*restart*.nc
# move output to VARS
mkdir VARS
mv BASIN_*grid*.nc VARS
