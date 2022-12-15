#!/bin/bash
#
#SBATCH --job-name=cue_01
#SBATCH --time=48:00:00
#SBATCH --mem=2000
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
#python "/home/x_swakh/tools/HoliSoils/Scripts/Linux/transient/Simulation_processing/calculate_cue_s.py" --sim_label "competition_adaptation" --scenario "gen_spec_lognorm_0_1"
python "/home/x_swakh/tools/HoliSoils/Scripts/Linux/transient/Simulation_processing/calc_cue_os_fd_wremainingc.py" "gen_spec_lognorm_0_1"
# Scripts ends here