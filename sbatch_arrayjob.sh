#!/bin/bash
#SBATCH -J imperf_atom_cav
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:10:00
#SBATCH --mem=1MB


echo "========= started at `date` =========="


# Load Julia
module purge
# module load julia/1.8.2-gcc-12.2.0-zn5hofz
module load julia/1.10.4-gcc-12.2.0-b3e3xer   


# Print some info
echo "  "
echo "  mem   = $SLURM_MEM_PER_NODE MB"
echo "  #proc = $SLURM_CPUS_PER_TASK"


# Run code 
# julia main.jl $INPUT $DATA $SLURM_ARRAY_TASK_ID
mpiexec -n $SLURM_CPUS_PER_TASK julia --project test.jl

echo "========= job finished at `date` =========="

