#!/bin/bash
#SBATCH --account=XXX
#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=96
#SBATCH --partition=ceres    # standard node(s)
#SBATCH --job-name="sample"
#SBATCH --mail-user=XXX   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -o errexit
set -o nounset
set -o xtrace

date   # print timestamp

# If using conda environment for dependencies:
ml miniconda
source activate pandagma 

PATH=$PWD/bin:$PATH

prot_fam_path="/project/legume_project/common_data/genefamilies/legume.fam3.VLMQ.sup1B_proteomes"

make_reduced_fam_sample.sh -f "$prot_fam_path" -w work_dir -m hybrid -n 100

date   # print timestamp

