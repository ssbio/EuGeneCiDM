#usr/bin/perl

#Written by: Wheaton Schroeder
#Latest Version: 06/24/2021

#Written to more quickly submit all "real" circuit solvers for EuGeneCiD/M quickly

#real buffers
printf "buffer_Cd: ";
system("sbatch buffer_Cd/sub.slurm"); #1
printf "buffer_Cu: ";
system("sbatch buffer_Cu/sub.slurm"); #2
printf "buffer_Zn: ";
system("sbatch buffer_Zn/sub.slurm"); #3

#test circuits
printf "test_adder: ";
system("sbatch test_adder/sub.slurm"); #4
printf "test_and: ";
system("sbatch test_and/sub.slurm"); #5
printf "test_buffer: ";
system("sbatch test_buffer/sub.slurm"); #6
printf "test_half_adder: ";
system("sbatch test_half_adder/sub.slurm"); #7
printf "test_nand: ";
system("sbatch test_nand/sub.slurm"); #8
printf "test_nor: ";
system("sbatch test_nor/sub.slurm"); #9
printf "test_not: ";
system("sbatch test_not/sub.slurm"); #10
printf "test_or: ";
system("sbatch test_or/sub.slurm"); #11
printf "test_xnor: ";
system("sbatch test_xnor/sub.slurm"); #12
printf "test_xor: ";
system("sbatch test_xor/sub.slurm"); #13