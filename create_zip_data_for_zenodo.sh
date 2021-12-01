#!/bin/bash
65;6003;1c
# This script creates zip files ready to be uploaded to zenodo.
# The file will contain the data/raw/EXP_experiment folder

cd data/raw
for i in $(ls)
do
    zip -r -1 ../${i}.zip $i
done
