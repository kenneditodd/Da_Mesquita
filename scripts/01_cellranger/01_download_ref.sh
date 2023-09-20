#!/bin/bash

# go to refs folder
cd ../../refs

# From cellranger
# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

# Extract files
tar -xf refdata-gex-mm10-2020-A.tar.gz
