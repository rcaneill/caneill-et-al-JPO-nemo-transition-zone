#!/bin/bash

DATA=data
RAW=${DATA}/raw
TMP=${DATA}/tmp

mkdir -p $RAW $TMP
cd $TMP
pipenv run zenodo_get 5607673
zip -s 0 EXP_main.zip --out out.zip
unzip out.zip
mv EXP_main ../../${RAW}/EXP_main
