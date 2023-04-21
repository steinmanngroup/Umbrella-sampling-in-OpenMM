#!/usr/bin/env bash
#SBATCH --job-name F0_benzene
#SBATCH --partition teach
#SBATCH --time 02:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0

export PATH=/usr/local/cuda-11.1/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.1/lib64:$LD_LIBRARY_PATH
export OPENMM_CPU_THREADS=$SLURM_CPUS_PER_TASK

source activate pyth39

ID=$SLURM_ARRAY_TASK_ID
L=20
MAIN_PATH=$1
RUN_ON_GPU=$2

# HIGH_INTERVAL: $3: if yes the interval is [0:26], else interval is [abs(-26):0]

# low interval
python F0_3.py $MAIN_PATH $RUN_ON_GPU no ${L[$ID]}