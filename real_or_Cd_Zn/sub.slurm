#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16000
#SBATCH --time=168:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
module load gams
gams real_or_Cd_Zn/EuGeneCiDM_or_Cd_Zn.gms