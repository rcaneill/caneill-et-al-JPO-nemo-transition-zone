#!/bin/bash

cd $1

for i in $(ls .)
do
    cd $i
    # make sure script is executable
    chmod +x run_nemo_restart.sh
    # submit to sbatch
    sbatch run_nemo_restart.sh
    cd ..
done
