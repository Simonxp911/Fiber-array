#!/bin/bash
#SBATCH -J imperf_atom_cav
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=1GB
#SBATCH --array=1-6
####SBATCH --account=p70623


echo "========= started at `date` =========="


# Load Julia
module purge
module load julia/1.8.2-gcc-12.2.0-zn5hofz
# module load julia/1.10.4-gcc-12.2.0-b3e3xer   

# Input filename
INPUT=input_ff_0.9,0.95,6_pur_0.1_w0r_5.0_L_6.2_rad_6.txt

# Print some info
echo "  "
echo "  INPUT = $INPUT"
echo "  mem   = $SLURM_MEM_PER_NODE MB"
echo "  array = $SLURM_ARRAY_TASK_MIN - $SLURM_ARRAY_TASK_MAX"


# Run code 
julia main.jl $INPUT $DATA $SLURM_ARRAY_TASK_ID


echo "========= job finished at `date` =========="

