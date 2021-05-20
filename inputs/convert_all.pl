#usr/bin/perl -w

#Written by: Wheaton Schroeder
#Latest Version: 5/20/21

#written to parallelize conversion of the unified inputs

use strict;

#and conceptualizations
system("perl convert_input_and_Cd_Cu.pl");
system("perl convert_input_and_Cd_Zn.pl");
system("perl convert_input_and_Cu_Zn.pl");

#half-adder conceptualizations
system("perl convert_input_half_adder_Cd_Cu.pl");
system("perl convert_input_half_adder_Cd_Zn.pl");
system("perl convert_input_half_adder_Cu_Zn.pl");

#nand conceptualizations
system("perl convert_input_nand_Cd_Cu.pl");
system("perl convert_input_nand_Cd_Zn.pl");
system("perl convert_input_nand_Cu_Zn.pl");

#nor conceptualizations
system("perl convert_input_nor_Cd_Cu.pl");
system("perl convert_input_nor_Cd_Zn.pl");
system("perl convert_input_nor_Cu_Zn.pl");

#or conceptualizations
system("perl convert_input_or_Cd_Cu.pl");
system("perl convert_input_or_Cd_Zn.pl");
system("perl convert_input_or_Cu_Zn.pl");

#xnor conceptualizations
system("perl convert_input_xnor_Cd_Cu.pl");
system("perl convert_input_xnor_Cd_Zn.pl");
system("perl convert_input_xnor_Cu_Zn.pl");

#xor conceptualizations
system("perl convert_input_xor_Cd_Cu.pl");
system("perl convert_input_xor_Cd_Zn.pl");
system("perl convert_input_xor_Cu_Zn.pl");

#nimply conceptualizations
system("perl convert_input_Cd_nimply_Cu.pl");
system("perl convert_input_Cd_nimply_Zn.pl");
system("perl convert_input_Cu_nimply_Zn.pl");
system("perl convert_input_Cu_nimply_Cd.pl");
system("perl convert_input_Zn_nimply_Cd.pl");
system("perl convert_input_Zn_nimply_Cu.pl");



