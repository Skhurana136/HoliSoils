#!/bin/bash
#
#SBATCH -J test_carbon_null
#SBATCH -t 00:60:00
#SBATCH --mem=2000
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=swamini.khurana@natgeo.su.se
#
# Run a single task in the foreground.
module load buildtool-easybuild/4.5.3-nsce8837e7
module load foss/2020b
module load Anaconda/2021.05-nsc1
conda activate ds-envsci-env
python "/home/x_swakh/tools/HoliSoils/Scripts/Linux/transient/Run_scripts/activity_loss_-02/carbon_switch_off_null.py" "null"
#
# Scripts ends here