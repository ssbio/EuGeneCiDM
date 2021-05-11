#usr/bin/perl

#Written by: Wheaton Schroeder
#Latest Version: 04/26/2021

#Written to more quickly submit all "real" circuit solvers for EuGeneCiD/S quickly

#AND
printf "real_and_Cd_Cu: ";
system("sbatch real_and_Cd_Cu/sub.slurm"); #1
printf "real_and_Cd_Zn: ";
system("sbatch real_and_Cd_Zn/sub.slurm"); #2
printf "real_and_Cu_Zn: ";
system("sbatch real_and_Cu_Zn/sub.slurm"); #3


#NIMPLY
printf "real_Cd_nimply_Cu: ";
system("sbatch real_Cd_nimply_Cu/sub.slurm"); #4
printf "real_Cd_nimply_Zn: ";
system("sbatch real_Cd_nimply_Zn/sub.slurm"); #5
printf "real_Cu_nimply_Zn: ";
system("sbatch real_Cu_nimply_Zn/sub.slurm"); #6


#CDI
printf "real_Cu_nimply_Cd: ";
system("sbatch real_Cu_nimply_Cd/sub.slurm"); #7
printf "real_Zn_nimply_Cd: ";
system("sbatch real_Zn_nimply_Cd/sub.slurm"); #8
printf "real_Zn_nimply_Cu: ";
system("sbatch real_Zn_nimply_Cu/sub.slurm"); #9


#HALF ADDER
printf "real_half_adder_Cd_Cu: ";
system("sbatch real_half_adder_Cd_Cu/sub.slurm"); #10
printf "real_half_adder_Cd_Zn: ";
system("sbatch real_half_adder_Cd_Zn/sub.slurm"); #11
printf "real_half_adder_Cu_Zn: ";
system("sbatch real_half_adder_Cu_Zn/sub.slurm"); #12


#NAND
printf "real_nand_Cd_Cu: ";
system("sbatch real_nand_Cd_Cu/sub.slurm"); #13
printf "real_nand_Cd_Zn: ";
system("sbatch real_nand_Cd_Zn/sub.slurm"); #14
printf "real_nand_Cu_Zn: ";
system("sbatch real_nand_Cu_Zn/sub.slurm"); #15


#NOR
printf "real_nor_Cd_Cu: ";
system("sbatch real_nor_Cd_Cu/sub.slurm"); #16
printf "real_nor_Cd_Zn: ";
system("sbatch real_nor_Cd_Zn/sub.slurm"); #17
printf "real_nor_Cu_Zn: ";
system("sbatch real_nor_Cu_Zn/sub.slurm"); #18


#OR
printf "real_or_Cd_Cu: ";
system("sbatch real_or_Cd_Cu/sub.slurm"); #19
printf "real_or_Cd_Zn: ";
system("sbatch real_or_Cd_Zn/sub.slurm"); #20
printf "real_or_Cu_Zn: ";
system("sbatch real_or_Cu_Zn/sub.slurm"); #21


#XNOR
printf "real_xnor_Cd_Cu: ";
system("sbatch real_xnor_Cd_Cu/sub.slurm"); #22
printf "real_xnor_Cd_Zn: ";
system("sbatch real_xnor_Cd_Zn/sub.slurm"); #23
printf "real_xnor_Cu_Zn: ";
system("sbatch real_xnor_Cu_Zn/sub.slurm"); #24


#XOR
printf "real_xor_Cd_Cu: ";
system("sbatch real_xor_Cd_Cu/sub.slurm"); #25
printf "real_xor_Cd_Zn: ";
system("sbatch real_xor_Cd_Zn/sub.slurm"); #26
printf "real_xor_Cu_Zn: ";
system("sbatch real_xor_Cu_Zn/sub.slurm"); #27


close LOG;