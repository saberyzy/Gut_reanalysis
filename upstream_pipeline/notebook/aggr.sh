#!/bin/bash
export PATH=/domino/edv/id-td-virology/Aligner/cellranger/cellranger-7.0.0:$PATH
cd /domino/edv/id-td-virology/Public_dataset/2022_Gut/aggr/
cellranger aggr \
--id=Gut_merge \
--csv=/domino/edv/id-td-virology/Public_dataset/2022_Gut/aggr/aggr_csv.csv \
--normalize=none \
--nosecondary
