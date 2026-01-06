#!/bin/bash

# Edited by Riley Porter 12/31/25

# Load required modules
module load shared rc-base
module load GCCcore/13.2.0
module load R
module load Anaconda3

# Create software directory
mkdir -p /home/$(whoami)/software/homer

# Download and install HOMER
wget http://homer.ucsd.edu/homer/configureHomer.pl -P /home/$(whoami)/software/homer
cd /home/$(whoami)/software/homer
perl configureHomer.pl -install
perl configureHomer.pl -install mm10
perl configureHomer.pl -install rn6
perl configureHomer.pl -install hg38
perl configureHomer.pl -update
cd -

# Download and install A.C.Rsuite
cd /home/$(whoami)/software/
git clone https://github.com/rikutakei/A.C.Rsuite.git
cd -

# Download and install IDR (A.C.Rsuite dependency)
wget https://github.com/nboley/idr/archive/2.0.2.zip -P /home/$(whoami)/software/
cd /home/$(whoami)/software
unzip 2.0.2.zip && rm 2.0.2.zip
cd /home/$(whoami)/software/idr-2.0.2/
cp /home/$(whoami)/software/A.C.Rsuite/idr.py /home/$(whoami)/software/idr-2.0.2/idr/
pip install --user blosc
pip install --user Cython
python3 setup.py install --user
pip3 install --user matplotlib
cd -

# Update PATH
echo 'export PATH=/home/$(whoami)/software/homer/bin/:/home/$(whoami)/software/A.C.Rsuite/:$PATH' >> ~/.bashrc

# Install R packages
Rscript --vanilla -e '
  pkgs <- c("optparse", "ggplot2", "data.table", "dplyr", "BiocManager")
  user_lib <- Sys.getenv("R_LIBS_USER")
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  install.packages(pkgs, lib = user_lib, repos = "https://cloud.r-project.org")

  if (!requireNamespace("BiocManager", quietly=TRUE, lib.loc=user_lib)) {
    install.packages("BiocManager", lib=user_lib, repos="https://cloud.r-project.org")
  }
  BiocManager::install("DESeq2", ask=FALSE, update=FALSE)
'

# Define R_LIBS_USER
echo 'export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.4' >> ~/.bash_profile

# Source updated bash config
source ~/.bashrc

