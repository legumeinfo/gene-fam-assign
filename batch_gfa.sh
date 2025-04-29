#!/bin/bash
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="hmmsrch"
## #SBATCH --mail-user=YOU@DOMAIN  # email address
## #SBATCH --mail-type=BEGIN
## #SBATCH --mail-type=END
## #SBATCH --mail-type=FAIL

set -o errexit
set -o nounset
set -o xtrace

date   # print timestamp

# If using conda environment for dependencies:
ml miniconda
source activate hmmer

PATH=$PWD/bin:$PATH

gfa.sh -l data/lis.protein_files -c config/gfa.conf

date   # print timestamp

