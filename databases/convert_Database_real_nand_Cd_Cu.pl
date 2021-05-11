#usr/bin/perl -w

#Written by: Wheaton Schroeder
#Latest version: 01/14/2020

#Written to convert the Microsoft Excel version
#of the database for EuGeneCiD/S to the requisite
#input files

use strict;
use Spreadsheet::Read qw(ReadData);
use Algorithm::Combinatorics qw(variations);

my $gate = "real_nand_Cd_Cu";

#create a file to write the output to
#create a directory for all of the output files, if it doesn't already exist

my $directory = $gate."_inputs";

#checks if directory exists, if not, tries to make it
unless(-e $directory or mkdir $directory) {
	
	#if the directory does not exist and could not make the directory
	die "unable to create directory $directory\n";
	
}

open(LOG, ">".$directory."/log.txt") or die "could not write to or create log file, reason: $!\n";
select LOG;

#create arrays for all the set files
my @pro = ( );				#will be the identifier of the promotor 
my @tran = ( );				#will be the transcript identifier
my @term = ( );				#will be the terminator identifer
my @enz = ( );				#will be the enzyme indentifier
my @lig = ( );				#list of ligans
my @des_lig = ( );			#desired ligands, those ligands that are used in the logic table
my @des_enz = ( ); 			#desired enzymes, those enzymes that are used in the logic table
my @all = ( );				#basically set of all molecules in the system used to identify self-self iteractions

#create hashes for all the requisite output files
#promotor-related parameters
my %pro_name = ( );			#promotor name
my %pro_state = ( );		#promotor normal state
my %pro_str = ( );			#promotor strength
my %pro_leak = ( );			#promotor leakiness
my %pro_ind = ( );			#promotor inducers
my %pro_ind_str = ( );		#promotor inducer strength 2-d hash
my %pro_rep = ( );			#promotor repressors
my %pro_rep_str = ( );		#promotor repressors strength 2-d hash

#transcript-related parameters
my %tran_name = ( );		#transcript names
my %tran_eff = ( ); 		#transcript translational efficiency
my %tran_enz = ( );			#enzyme encoded by transcript

#terminator-related parameters
my %term_name = ( ); 		#terminator names
my %term_half = ( );		#terminator half-life

#enzyme-related parameters
my %enz_name = ( );			#enzyme name
my %enz_state = ( );		#enzyme normal state
my %enz_thresh = ( );		#enzyme concentration threshold for activity
my %enz_half = ( );			#enzyme half-life
my %enz_ind = ( );			#promotor inducers
my %enz_ind_str = ( );		#promotor inducer strength 2-d hash
my %enz_rep = ( );			#promotor repressors
my %enz_rep_str = ( );		#promotor repressors strength 2-d hash

#ligand-associated names
my %lig_names = ( );		#ligand names

#logic table
my @logic_table = ( ); 		#will be strings specifying the logic table

#stores the number of ligands responding to in the logic table
my $num_resp = 0;

#read the workbook
my $data_book = Spreadsheet::Read->new('database_'.$gate.'.xlsx');

#read the worksheets
#note: sheet names are case-sensitive
my $pro_sheet = $data_book->sheet("Promotors");				#"Promotors"
my $trans_sheet = $data_book->sheet("Transcripts");			#"Transcripts"
my $term_sheet = $data_book->sheet("Terminators");			#"Terminators"
my $enz_sheet = $data_book->sheet("Enzymes&Proteins");		#"Enzymes&Proteins"
my $lig_sheet = $data_book->sheet("Ligands");				#"Ligands"
my $log_sheet = $data_book->sheet("LogicTable");			#"LogicTable"
my $other_sheet = $data_book->sheet("Other");				#"other"

#read the data from the promotors sheet

#read by row, get the largest number of rows
my $pro_row_num = $pro_sheet->maxrow;

printf "\nNumber of promotors: %s\n",($pro_row_num-1);

#note, will be skipping the first row since that is just labels
#start with 2 bacause the indexing is 1-based 
for(my $a = 2; $a <= $pro_row_num; $a++) {
	
	#should have a consistent number of columns
	#remove superfluous white space
	chomp(my $name = $pro_sheet->cell(1,$a));
	chomp(my $ident = $pro_sheet->cell(2,$a));
	chomp(my $state = $pro_sheet->cell(3,$a));
	chomp(my $str = $pro_sheet->cell(4,$a));
	chomp(my $leak = $pro_sheet->cell(5,$a));
	chomp(my $ind_pair = $pro_sheet->cell(6,$a));
	chomp(my $rep_pair = $pro_sheet->cell(7,$a));
	
	#add to the array and hashes
	
	#parse the inducer and repressor info
	#note that each inducer and repressor pair is seprated by ";"
	my @ind_list = split /;/, $ind_pair;
	
	my @ind_array = ( );
	
	#get the data for the inducers 
	for(my $b = 0; $b <= $#ind_list; $b++) {
		
		#spit, get the strength from the parentheses 
		(my $ind, my $ind_str) = split /\(/, $ind_list[$b];
		
		$ind =~ s/\s+$//;
		$ind =~ s/^\s+//;
		
		$ind_str =~ s/\)//g;
		
		push @ind_array, $ind;
		
		$pro_ind_str{$ident}{$ind} = $ind_str;
		
	}
	
	#the inducers are an array hash is array handle
	$pro_ind{$ident} = \@ind_array;
	
	#parse the inducer and repressor info
	#note that each inducer and repressor pair is seprated by ";"
	my @rep_list = split /;/, $rep_pair;
	
	my @rep_array = ( );
	
	#get the data for the inducers 
	for(my $c = 0; $c <= $#rep_list; $c++) {
		
		#spit, get the strength from the parentheses 
		(my $rep, my $rep_str) = split /\(/, $rep_list[$c];
		
		$rep =~ s/\s+$//;
		$rep =~ s/^\s+//;
		
		$rep_str =~ s/\)//g;
		
		push @rep_array, $rep;
		
		$pro_rep_str{$ident}{$rep} = $rep_str;
		
	}
	
	#the repressors are an array hash is array handle
	$pro_rep{$ident} = \@rep_array;
	
	push @pro, $ident;
	push @all, $ident;
	
	$pro_name{$ident} = $name;
	$pro_state{$ident} = $state;
	$pro_str{$ident} = $str;
	$pro_leak{$ident} = $leak;
	
}

#read the data from the transcripts sheet

#read by row, get the largest number of rows
my $tran_row_num = $trans_sheet->maxrow;

printf "Number of transcripts: %s\n",($tran_row_num-1);

#note, will be skipping the first row since that is just labels
#start with 2 bacause the indexing is 1-based 
for(my $d = 2; $d <= $tran_row_num; $d++) {
	
	#should have a consistent number of columns
	#remove superfluous white space
	chomp(my $name = $trans_sheet->cell(1,$d));
	chomp(my $ident = $trans_sheet->cell(2,$d));
	chomp(my $eff = $trans_sheet->cell(3,$d));
	chomp(my $encoded = $trans_sheet->cell(4,$d));
	
	push @tran, $ident;
	push @all, $ident;
	
	$tran_name{$ident} = $name;
	$tran_eff{$ident} = $eff;
	$tran_enz{$ident}{$encoded} = 1;
	
}

#read the data from the terminators sheet

#read by row, get the largest number of rows
my $term_row_num = $term_sheet->maxrow;

printf "Number of terminators: %s\n",($term_row_num-1);

#note, will be skipping the first row since that is just labels
#start with 2 bacause the indexing is 1-based 
for(my $e = 2; $e <= $term_row_num; $e++) {
	
	#should have a consistent number of columns
	#remove superfluous white space
	chomp(my $name = $term_sheet->cell(1,$e));
	chomp(my $ident = $term_sheet->cell(2,$e));
	chomp(my $half = $term_sheet->cell(3,$e));
	
	$term_name{$ident} = $name;
	$term_half{$ident} = $half;
	
	push @term, $ident;
	push @all, $ident;
	
}

#read the data from the enzyme sheet

#read by row, get the largest number of rows
my $enz_row_num = $enz_sheet->maxrow;

printf "Number of enzymes: %s\n",($enz_row_num-1);

#note, will be skipping the first row since that is just labels
#start with 2 bacause the indexing is 1-based 
for(my $f = 2; $f <= $enz_row_num; $f++) {
	
	#should have a consistent number of columns
	#remove superfluous white space
	chomp(my $name = $enz_sheet->cell(1,$f));
	chomp(my $ident = $enz_sheet->cell(2,$f));
	chomp(my $state = $enz_sheet->cell(3,$f));
	chomp(my $thresh = $enz_sheet->cell(4,$f));
	chomp(my $half = $enz_sheet->cell(5,$f));
	chomp(my $ind_pair = $enz_sheet->cell(6,$f));
	chomp(my $rep_pair = $enz_sheet->cell(7,$f));
	
	#add to the array and hashes
	
	#parse the inducer and repressor info
	#note that each inducer and repressor pair is seprated by ";"
	my @ind_list = split /;/, $ind_pair;
	
	my @ind_array = ( );
	
	#get the data for the inducers 
	for(my $g = 0; $g <= $#ind_list; $g++) {
		
		#spit, get the strength from the parentheses 
		(my $ind, my $ind_str) = split /\(/, $ind_list[$g];
		
		$ind =~ s/^\s+//;
		$ind =~ s/\s+$//;
		
		$ind_str =~ s/\)//g;
		
		push @ind_array, $ind;
		
		$enz_ind_str{$ident}{$ind} = $ind_str;
		
	}
	
	#the inducers are an array hash is array handle
	$enz_ind{$ident} = \@ind_array;
	
	#parse the inducer and repressor info
	#note that each inducer and repressor pair is seprated by ";"
	my @rep_list = split /;/, $rep_pair;
	
	my @rep_array = ( );
	
	#get the data for the inducers 
	for(my $h = 0; $h <= $#rep_list; $h++) {
		
		#spit, get the strength from the parentheses 
		(my $rep, my $rep_str) = split /\(/, $rep_list[$h];
		
		$rep =~ s/^\s+//;
		$rep =~ s/\s+$//;
		
		$rep_str =~ s/\)//g;
		
		push @rep_array, $rep;
		
		$enz_rep_str{$ident}{$rep} = $rep_str;
		
	}
	
	#the inducers are an array hash is array handle
	$enz_rep{$ident} = \@rep_array;
	
	push @enz, $ident;
	push @all, $ident;
	
	$enz_name{$ident} = $name;
	$enz_state{$ident} = $state;
	$enz_thresh{$ident} = $thresh;
	$enz_half{$ident} = $half;	
	
}

#read the ligands sheet

#read by row, get the largest number of rows
my $lig_row_num = $lig_sheet->maxrow;

printf "Number of ligands: %s\n",($lig_row_num-1);

#note, will be skipping the first row since that is just labels
#start with 2 bacause the indexing is 1-based 
for(my $i = 2; $i <= $lig_row_num; $i++) {
	
	chomp(my $name = $lig_sheet->cell(1,$i));
	chomp(my $ident = $lig_sheet->cell(2,$i));
	
	push @lig, $ident;
	push @all, $ident;
	
	$lig_names{$ident} = $name;
	
}

#need to add "none" as a ligand
push @lig, "none";
push @all, "none";
$lig_names{"none"} = "none";

#read the logic table
my $log_row_num = $log_sheet->maxrow;
my $log_col_num = $log_sheet->maxcol;

#count the number of 
my $num_enz = 0;

#the first row has the desired ligand conditions and enzyme responders
#this is reading just the first row
for(my $j = 0; $j <= $log_col_num; $j++) {
	
	#get a header item
	chomp(my $header = $log_sheet->cell($j,1));
	
	#check if header item is a ligand
	if ($header ~~ @lig) {
		
		#if here then the header is a ligand
		push @des_lig, $header;
		$num_resp++;
		
	#or check if header item is an enzyme
	} elsif ($header ~~ @enz) {
		
		#if here then the header is an enzyme
		push @des_enz, $header;
		$num_enz++;
		
	#or header isn't a ligand or an enzyme which is a problem
	} else {
		
		#if here, then the header is neither this is a problem
		#we've got a problem if we are here
		
	}
	
}

#need to add the "none" condition to the desired ligands
#the "none" will not be in the count however
push @des_lig, "none";

#after that, read the lines for the conditions and responses
#each row is a specific condition and response
#note that condition should be order unspecific and also should get pair with "none" if only
#one condition is specified
for(my $k = 2; $k <= $log_row_num; $k++) {
	
	#array of ligands specified to exist in the current condition
	#this array will contain all repeats necessary to find all combinations
	my @temp_cond = ( );
	
	#this array contains the unique conditions
	my @temp_cond_uni = ( );
	
	#save the number of ligands in the current condition
	my $cond_lig_num = 0;
	
	#array of desired respondants
	my @temp_resp = ( );
	
	printf "\nrow #%s present ligands: ", ($k-1);
	
	for(my $cc = 1; $cc <= $log_col_num; $cc++) {
	
		#get a cell from the sheet
		chomp(my $binary = $log_sheet->cell($cc,$k));
	
		#check if the binary is specifying the ligand condition or enzyme response
		if($cc <= $#des_lig) {
			
			#then is a ligand specification can get from index and binary
			
			#if the binary is 1 then use the column to return the ligand
			if($binary =~ /^\s*1\s*$/) {
				
				#note that arrays are 0-based whereas spreadsheet reading is 10based
				push @temp_cond, $des_lig[($cc-1)];
				push @temp_cond_uni, $des_lig[($cc-1)];
				printf "%s ", $des_lig[($cc-1)];
				
				$cond_lig_num++;
				
			} #if a zero will ignore only looking for present makes life easier
			
		} else {
			
			#if got to here with no desired ligands, then is a "none" condition
			if(!@temp_cond) {
				
				printf "none";
				
			}
			
			#otherwise is a specified enzyme response
			#again will only look at those where a response is desired
			if($binary =~ /^\s*1\s*$/) {
				
				#again issue of 1-based vs 0-based conversion
				push @temp_resp, $des_enz[($cc-$#des_lig-1)];
				printf "\nenzyme response: %s ", $des_enz[($cc-$#des_lig-1)];				
				
			} #note will assume no response if incorrectly formatted
			
		}
	
	}
	
	if(!@temp_resp) {
		
		printf "\ndesired to have no enzyme response\n";
		
	} else {
		
		printf "\n";
	
	}
	
	#if no conditions, specify none as the desired unique
	if(!@temp_cond) {
		
		push @temp_cond_uni, "none";
		
	}
	
	#okay now that we know the ligands and the responding enzymes need to specify strings to write to logic
	#table file
	#note this relies on zero being default in GAMS
	
	#size of desired ligands minus one is the number of 
	#note that this combination of ligand should not re-occur
	#each of the ligands in the @temp_cond array should be included
	
	#the number of "none" that should be included then will be the number of ligands conditions by which
	#the @temp_cond array is short
	
	#next figure out the number of duplicate entries for each current ligand
	#for each ligand defining the condition
	
	#do the calculation to find out how many copies
	my $num_lig_copies = $num_resp - $cond_lig_num;
	
	for(my $hh = 1; $hh <= $cond_lig_num; $hh++) {
		
		#add desired number of copies
		for(my $ii = 1; $ii <= $num_lig_copies; $ii++) {
			
			push @temp_cond, $temp_cond[($hh-1)];
			
		}
		
	}
	
	#find the number of none's that need to be added
	#should be equal to the number of ligand copies
	my $num_none = $num_lig_copies;
	
	#add this number of "nones" to the @temp_cond array
	for(my $dd = 1; $dd <= $num_none; $dd++) {
		
		push @temp_cond, "none";
		
	}
	
	#at this point the number of conditions that we have should be equal to the number of desired ligands
	#need to produce each possible shuffled version of the array @temp_cond, will deliberately exclude 
	#none-none cases.
	
	#let us add comment lines to the logic table
	my $comment = join " ", @temp_cond_uni;
	$comment =~ s/^/\n\/*/;
	$comment =~ s/$/ ligands present condition*\//;
	
	push @logic_table, $comment;
	
	#need to write for each of the responses
	for(my $ff = 0; $ff <= $#temp_resp; $ff++) {
	
		#get each permutation of the conditions
		my $iter = variations(\@temp_cond,$num_resp);
		
		while(my $ee = $iter->next) {
			
			#note that $ee will be an array handle
			my @temp_arr = @$ee;
			
			#combine the permutations
			my $cond_str = join "\'.\'", @temp_arr;
			$cond_str =~ s/^/\'/;
			$cond_str =~ s/$/\'/;
			
			#next, add the current respondant enzyme
			$cond_str .= ".\'".$temp_resp[$ff]."\' 1";
			
			#state the condition to allow us to check things
			printf "condition string: %s\n", $cond_str;
			
			#variable to determine if has all the necessary bits
			#easiest to start by assuming it does, change to 0 if missing one bit
			my $has_all_bits = 1;
			
			#check to see if has each unique ligand for the condition at least once
			for(my $kk = 0; $kk <= $#temp_cond_uni; $kk++) {
				
				#check if the condition has this bit
				if($cond_str =~ /\'$temp_cond_uni[$kk]\'/) {
					
					#has the bit
					printf "found %s\n", $&;
					
				} else {
					
					printf "did not find \'%s\'\n", $temp_cond_uni[$kk];
					$has_all_bits = 0;
					
				}
				
			}
			
			#if the string is not already in the logic table, add it otherwise don'table
			if($cond_str ~~ @logic_table) {
				
				#don't add already there
				
			} elsif($has_all_bits eq 1) { #need to make sure has each unique ligand for the condition at least once
				
				#select a present ligand to fill this position
				push @logic_table, $cond_str;
				
			} else {
				
				#not already in the logic table and doesn't have all the bits do nothing
				
			}
			
		}
		
	}
	
}

#finally, read what the maximum time should be
chomp(my $max_time = $other_sheet->cell(2,1));

printf "\nmax time: %s\n", $max_time;

#write all of the necessary output files to the directory to keep things neat

#write the files 
open(SET_A, ">".$directory."/all_molecules_".$gate.".txt") or die "Could not write/create the file for set A, reason: $!\n";
open(PARAM_SIGMA, ">".$directory."/self_".$gate.".txt") or die "Could not write/create the file for parameter sigma, reason: $!\n";

printf SET_A "\/\n";
printf PARAM_SIGMA "\/\n";

for(my $n = 0; $n <= $#all; $n++) {
	
	printf SET_A "\'".$all[$n]."\'\n";
	printf PARAM_SIGMA "\'".$all[$n]."\'.\'".$all[$n]."\' 1\n";
	
}

printf SET_A "\/";
printf PARAM_SIGMA "\/";

close SET_A;
close PARAM_SIGMA;

#write transcript-based files
open(SET_J, ">".$directory."/transcripts_".$gate.".txt") or die "Could not write/create the file for set T, reason: $!\n";
open(PARAM_RHO, ">".$directory."/trans_to_enzyme_".$gate.".txt") or die "Could not write/create the file for parameter rho, reason: $!\n";
open(PARAM_ETA, ">".$directory."/translation_efficiency_".$gate.".txt") or die "Could not write/create the file for parameter eta, reason: $!\n";

printf SET_J "\/\n";
printf PARAM_RHO "\/\n";
printf PARAM_ETA "\/\n";

for(my $p = 0; $p <= $#tran; $p++) {
	
	printf SET_J "\'".$tran[$p]."\'\n";
	printf PARAM_ETA "\'".$tran[$p]."\' ".$tran_eff{$tran[$p]}."\n";
	
	for(my $q = 0; $q <= $#enz; $q++) {
		
		if($tran_enz{$tran[$p]}{$enz[$q]} eq 1) {
			
			printf PARAM_RHO "\'".$tran[$p]."\'.\'".$enz[$q]."\' ".$tran_enz{$tran[$p]}{$enz[$q]}."\n";
			
		}
		
	}
	
}

printf SET_J "\/";
printf PARAM_RHO "\/";
printf PARAM_ETA "\/";

#write promotor-based files
open(SET_P, ">".$directory."/promoters_".$gate.".txt") or die "Could not write/create the file for set P, reason: $!\n";
open(PARAM_S, ">".$directory."/promoter_strength_".$gate.".txt") or die "Could not write/create the file for parameter S, reason: $!\n";
open(PARAM_Z, ">".$directory."/promoter_normal_state_".$gate.".txt") or die "Could not write/create the file for parameter Z, reason: $!\n";
open(PARAM_I, ">".$directory."/promoter_ligand_interactions_".$gate.".txt") or die "Could not write/create the file for parameter I, reason: $!\n";
open(PARAM_H, ">".$directory."/promoter_ligand_strength_".$gate.".txt") or die "Could not write/create the file for parameter H, reason: $!\n";
open(PARAM_F, ">".$directory."/promoter_leakiness_".$gate.".txt") or die "Could not write/create the file for parameter F, reason: $!\n";

printf SET_P "\/\n";
printf PARAM_S "\/\n";
printf PARAM_Z "\/\n";
printf PARAM_I "\/\n";
printf PARAM_H "\/\n";
printf PARAM_F "\/\n";

for(my $r = 0; $r <= $#pro; $r++) {
	
	printf SET_P "\'".$pro[$r]."\'\n";
	printf PARAM_S "\'".$pro[$r]."\' ".$pro_str{$pro[$r]}."\n";
	printf PARAM_Z "\'".$pro[$r]."\' ".$pro_state{$pro[$r]}."\n";
	printf PARAM_F "\'".$pro[$r]."\' ".$pro_leak{$pro[$r]}."\n";
	
	#get the array of inducers for this promotor
	my $ind_handle = $pro_ind{$pro[$r]};
	
	#turn the array handle back into an array
	my @ind_list_temp = @$ind_handle;
	
	for(my $s = 0; $s <= $#ind_list_temp; $s++) {
		
		#state that it is an inducer
		printf PARAM_I "\'".$pro[$r]."\'.\'".$ind_list_temp[$s]."\' 1\n";
		
		#write the strength of the inducer interaction
		printf PARAM_H "\'".$pro[$r]."\'.\'".$ind_list_temp[$s]."\' ".$pro_ind_str{$pro[$r]}{$ind_list_temp[$s]}."\n";
		
	}
	
	#get the array of repressors for this promotor
	my $rep_handle = $pro_rep{$pro[$r]};
	
	#turn the array handle back into an array
	my @rep_list_temp = @$rep_handle;
	
	for(my $t = 0; $t <= $#rep_list_temp; $t++) {
		
		#state that it is a repressor
		printf PARAM_I "\'".$pro[$r]."\'.\'".$rep_list_temp[$t]."\' -1\n";
		
		#write the strength of the inducer interaction
		printf PARAM_H "\'".$pro[$r]."\'.\'".$rep_list_temp[$t]."\' ".$pro_rep_str{$pro[$r]}{$rep_list_temp[$t]}."\n";
		
	}
		
}

printf SET_P "\/";
printf PARAM_S "\/";
printf PARAM_Z "\/";
printf PARAM_I "\/";
printf PARAM_H "\/";
printf PARAM_F "\/";

#write ligand-based files
open(SET_L1, ">".$directory."/Ligands_".$gate.".txt") or die "Could not write/create the file for set L1, reason: $!\n";
open(SET_LD, ">".$directory."/desired_ligands_".$gate.".txt") or die "Could not write/create the file for set LD, reason: $!\n";
open(PARAM_LD_VAL, ">".$directory."/desired_ligands_".$gate."_vals.txt") or die "Could not write/create the file for parameter Ld_val, reason: $!\n";

printf SET_L1 "\/\n";
printf SET_LD "\/\n";
printf PARAM_LD_VAL "\/\n";

for(my $u = 0; $u <= $#lig; $u++) {
	
	printf SET_L1 "\'".$lig[$u]."\'\n";
	
}

for(my $v = 0; $v <= $#des_lig; $v++) {
	
	printf SET_LD "\'".$des_lig[$v]."\'\n";
	printf PARAM_LD_VAL "\'".$des_lig[$v]."\' 1\n";
	
}

printf SET_L1 "\/";
printf SET_LD "\/";
printf PARAM_LD_VAL "\/";

#write terminator-based files
open(SET_T, ">".$directory."/terminators_".$gate.".txt") or die "Could not write/create the file for set T, reason: $!\n";
open(PARAM_G, ">".$directory."/terminator_strength_".$gate.".txt") or die "Could not write/create the file for parameter G, reason: $!\n";

printf SET_T "\/\n";
printf PARAM_G "\/\n";

for(my $w = 0; $w <= $#term; $w++) {
	
	printf SET_T "\'".$term[$w]."\'\n";
	printf PARAM_G "\'".$term[$w]."\' ".$term_half{$term[$w]}."\n";
	
}

printf SET_T "\/";
printf PARAM_G "\/";

#write enzyme/protein-based files
open(SET_E, ">".$directory."/enzymes_".$gate.".txt") or die "Could not write/create the file for set E, reason: $!\n";
open(PARAM_Q, ">".$directory."/protein_ligand_strength_".$gate.".txt") or die "Could not write/create the file for parameter Q, reason: $!\n";
open(PARAM_ZETA, ">".$directory."/protein_normal_state_".$gate.".txt") or die "Could not write/create the file for parameter Zeta, reason: $!\n";
open(PARAM_TAU, ">".$directory."/protein_degradation_".$gate.".txt") or die "Could not write/create the file for parameter tau, reason: $!\n";
open(PARAM_THETA, ">".$directory."/expression_threshold_".$gate.".txt") or die "Could not write/create the file for parameter theta, reason: $!\n";
open(SET_ED, ">".$directory."/desired_enzymes_".$gate.".txt") or die "Could not write/create the file for set P, reason: $!\n";
open(PARAM_B, ">".$directory."/protein_ligand_interactions_".$gate.".txt") or die "Could not write/create the file for parameter B, reason: $!\n";
open(PARAM_ED_VAL, ">".$directory."/desired_enzymes_".$gate."_vals.txt") or die "Could not write/create the file for parameter Ed_val, reason: $!\n";

printf SET_E "\/\n";
printf PARAM_Q "\/\n";
printf PARAM_ZETA "\/\n";
printf PARAM_TAU "\/\n";
printf PARAM_THETA "\/\n";
printf SET_ED "\/\n";
printf PARAM_B "\/\n";
printf PARAM_ED_VAL "\/\n";

for(my $x = 0; $x <= $#enz; $x++) {
	
	printf SET_E "\'".$enz[$x]."\'\n";
	printf PARAM_ZETA "\'".$enz[$x]."\' ".$enz_state{$enz[$x]}."\n";
	printf PARAM_TAU "\'".$enz[$x]."\' ".$enz_half{$enz[$x]}."\n";
	printf PARAM_THETA "\'".$enz[$x]."\' ".$enz_thresh{$enz[$x]}."\n";
	
	#get the array of inducers for this promotor
	my $ind_handle = $enz_ind{$enz[$x]};
	
	#turn the array handle back into an array
	my @ind_list_temp = @$ind_handle;
	
	for(my $y = 0; $y <= $#ind_list_temp; $y++) {
		
		#state that it is an inducer
		printf PARAM_B "\'".$enz[$x]."\'.\'".$ind_list_temp[$y]."\' 1\n";
		
		#write the strength of the inducer interaction
		printf PARAM_Q "\'".$enz[$x]."\'.\'".$ind_list_temp[$y]."\' ".$enz_ind_str{$enz[$x]}{$ind_list_temp[$y]}."\n";
		
	}
	
	#get the array of repressors for this promotor
	my $rep_handle = $enz_rep{$enz[$x]};
	
	#turn the array handle back into an array
	my @rep_list_temp = @$rep_handle;
	
	for(my $bb = 0; $bb <= $#rep_list_temp; $bb++) {
		
		#state that it is a repressor
		printf PARAM_B "\'".$enz[$x]."\'.\'".$rep_list_temp[$bb]."\' -1\n";
		
		#write the strength of the inducer interaction
		printf PARAM_Q "\'".$enz[$x]."\'.\'".$rep_list_temp[$bb]."\' ".$enz_rep_str{$enz[$x]}{$rep_list_temp[$bb]}."\n";
		
	}
	
}

for(my $z = 0; $z <= $#des_enz; $z++) {
	
	printf SET_ED "\'".$des_enz[$z]."\'\n";
	printf PARAM_ED_VAL "\'".$des_enz[$z]."\' 1\n";
	
}

printf SET_E "\/";
printf PARAM_Q "\/";
printf PARAM_ZETA "\/";
printf PARAM_TAU "\/";
printf PARAM_THETA "\/";
printf SET_ED "\/";
printf PARAM_B "\/";
printf PARAM_ED_VAL "\/";

#write the time sets
open(SET_BET, ">".$directory."/time_set_".$gate.".txt") or die "Could not write/create the file for set bet, reason: $!\n";

printf SET_BET "\/\n";

#basically need to enumerate all time points as set elements
for(my $aa = 0; $aa <= $max_time; $aa++) {
	
	 printf SET_BET "\'".$aa."\'\n";
	
}

printf SET_BET "\/";

#write the logic table
open(PARAM_LAMDBA, ">".$directory."/logic_table_".$gate.".txt") or die "Could not write/create the file for parameter lambda, reason: $!\n";

#start the file
printf PARAM_LAMDBA "\/\n";

#write each of the lines 
for(my $gg = 0; $gg <= $#logic_table; $gg++) {
	
	printf PARAM_LAMDBA "%s\n", $logic_table[$gg];
	
}

#finish the file
printf PARAM_LAMDBA "\/";