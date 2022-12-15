#!/bin/bash
#
#SBATCH --job-name=all_combine
#SBATCH --time=72:00:00
#SBATCH --mem=90000
#SBATCH -n 1
#SBATCH --chdir=/proj/hs_micro_div_072022
#SBATCH --output=./Reports/output_%j.out
#SBATCH --error=./Reports/error_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=swamini.khurana@natgeo.su.se
#
# Run a single task in the foreground.
module load buildtool-easybuild/4.5.3-nsce8837e7
module load foss/2020b
module load Anaconda/2021.05-nsc1
conda activate ds-envsci-env
python "/home/x_swakh/tools/HoliSoils/Scripts/Linux/transient/Simulation_processing/Combine_all_datasets.py"
# Scripts ends here