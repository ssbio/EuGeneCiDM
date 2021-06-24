#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16000
#SBATCH --time=168:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
module load gams
gams buffer_Cu/EuGeneCiDM_buffer_Cu.gms