#!/bin/bash

# Create an experiment
# $1 EXP_expname/experiments.csv
# $2 EXP_REF
# $3 EXP_expname/
# $4 nemo
# $5 EXP_expname/EXP00

mkdir -p $3 #/Experiments
cp -r $2 $5
cp $4 $5
echo "full_exp_name,short_name" > $1
