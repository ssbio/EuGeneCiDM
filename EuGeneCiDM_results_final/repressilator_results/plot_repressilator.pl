#usr/bin/perl -w

#Written by: Wheaton Schroeder
#Latest version: 04/29/2021

#Since the circuit state file takes quite a bit of manual rearrangement to
#easily plot some results, this code will essential, given an input solution
#number, plot the relative transcript and enzyme levels under each
#unique condition

use strict;
use Excel::Writer::XLSX;

#inputs of the desired logic gate and solution number

my $gate = "rep_CiM_only";
my $soln_num = 0;

#set up to read the information
open(STATEDATA, 'circuit_state_real_'.$gate.'.csv') or die "Could not open input file, reason: $!\n";

#create a multi-dimensional hash to store treanscript and enzymes levels
#and for transcript and enzyme production
my %trans_level = ( );
my %trans_prod = ( );
my %enz_level = ( );
my %enz_prod = ( );

#store the maximum number for setting the y-axis
my $max_trans_level = 0;
my $max_trans_prod = 0;
my $max_enz_level = 0;
my $max_enz_prod = 0;

#set of transcripts
my @trans = ( );

#set of enzymes
my @enz = ( );

#set of indeces for the levels
my @trans_level_ind = ( );
my @enz_level_ind = ( );

#set of indeces for the production
my @trans_prod_ind = ( );
my @enz_prod_ind = ( );

#hashes to store if a transcript or enyzme has some production at some point in the circuit
#will only bother to plot if has some production
my %enz_is_prod = ();
my %trans_is_prod = ();

#set of timepoints
my @time_points = ( );

#set of ligands
my @ligands = ( );

#hash to convert column number to letter
my %num_to_alpha = (

	0 => 'A',
	1 => 'B',
	2 => 'C',
	3 => 'D',
	4 => 'E',
	5 => 'F',
	6 => 'G',
	7 => 'H',
	8 => 'I',
	9 => 'J',
	10 => 'K',
	11 => 'L',
	12 => 'M',
	13 => 'N',
	14 => 'O',
	15 => 'P',
	16 => 'Q',
	17 => 'R',
	18 => 'S',
	19 => 'T',
	20 => 'U',
	21 => 'V',
	22 => 'W',
	23 => 'X',
	24 => 'Y',
	25 => 'Z'

);

#have a count to skip the header
my $line_count = 0;

#read just the lines for the current solution number
foreach my $line (<STATEDATA>) {
	
	#ignore the first 3 lines the header
	if($line_count <= 2) {
		
		#skip
		
	} elsif ($line_count == 3) {
		
		#split the lines into items
		my @items = split /,/, $line;
		
		#this is the header line that actually has the variable names
		for(my $a = 0; $a <= $#items; $a++) {
			
			#if labeled "phi_" then is the production of transcripts
			#capture the transcript name and the indices production
			if($items[$a] =~ /^phi_(.+)$/i) {
				
				my $trans_temp = $1;
				
				#remove excess spaces
				$trans_temp =~ s/^\s+//g;
				$trans_temp =~ s/\s+$//g;
				
				printf "identified transcript: %s\n", $trans_temp;
				
				push @trans_prod_ind, $a;
				push @trans, $trans_temp;
				
				#by default assume the transcript is not produced
				$trans_is_prod{$trans_temp} = 0;
				
			}
			
			#if labeled "RNA_carry_" then is is the transcript level 
			#capture the level
			if($items[$a] =~ /^RNA_carry_(.+)$/i) {
				
				push @trans_level_ind, $a;
				
			}
			
			#if labeled "Enzyme_carry_" then is enzyme level
			#capture enzyme name and level index
			if($items[$a] =~ /^enzyme_carry_(.+)$/i) {
				
				my $enz_temp = $1;
				
				#remove excess spaces
				$enz_temp =~ s/^\s+//g;
				$enz_temp =~ s/\s+$//g;
				
				printf "identified enzyme: %s\n", $enz_temp;
				
				push @enz_level_ind, $a;
				push @enz, $enz_temp;
				
				#by default assume the enzyme is not produced
				$enz_is_prod{$enz_temp} = 0;
				
			}
			
			#if labeled "C_" then is enzyme production
			#capture enzyme production index
			if($items[$a] =~ /^C_(.+)$/i) {
				
				push @enz_prod_ind, $a;
				
			}
			
			#this will be the limit of what we capture
			
		}
		
	} else {
		
		#split the lines into items
		my @items = split /,/, $line;
		
		#first item is solution number
		chomp(my $number = $items[0]);
		
		#remove the quotation marks
		$number =~ s/\"//g;
		
		#pull the time point
		chomp(my $timepoint = $items[1]);
		
		#remove the quotation marks
		$timepoint =~ s/\"//g;
		
		#check if already have the time point captured
		if(!($timepoint ~~ @time_points)) {
			
			printf "identified time point %s\n", $timepoint;
			push @time_points, $timepoint;
			
		}
		
		#get the ligand
		chomp(my $lig_1 = $items[2]);
		
		#remove the quotation marks
		$lig_1 =~ s/\"//g;
		
		#check if we already have the ligand captured
		if(!($lig_1 ~~ @ligands)) {
			
			printf "identified ligand %s\n", $lig_1;
			push @ligands, $lig_1;
			
		}
		
		#if it is not the correct solution number, skip
		if($number < $soln_num) {
			
			#do nothing wrong solution number, still haven't seen the correct one
			
		#solution number are always added in order, therefore if the
		#current solution numer is greater than our target no need to
		#continue to read
		} elsif ($number > $soln_num) {
			
			last;
			
		} else {
		
			#at the correct solution number then get the requisite information
			#get the transcript levels
			for(my $b = 0; $b <= $#trans_level_ind; $b++) {
				
				#get the item which will be the transcript level
				chomp(my $trans_temp = $items[$trans_level_ind[$b]]);
				
				if($trans_temp != 0) {
					
					$trans_is_prod{$trans[$b]} = 1;
					
				}
				
				#add to appropriate hash
				$trans_level{$trans[$b]}{$timepoint}{$lig_1} = $trans_temp;
				
				#update the max
				if($trans_temp > $max_trans_level) {
					
					$max_trans_level = $trans_temp;
					
				}
				
			}
			
			#basically repeat the above for each thing we need to capture and plot
			#get the transcript production
			#at the correct solution number then get the requisite information
			#get the transcript levels
			for(my $c = 0; $c <= $#trans_prod_ind; $c++) {
				
				#get the item which will be the transcript production
				chomp(my $trans_temp = $items[$trans_prod_ind[$c]]);
				
				if($trans_temp != 0) {
					
					$trans_is_prod{$trans[$c]} = 1;
					
				}
				
				#add to appropriate hash
				$trans_prod{$trans[$c]}{$timepoint}{$lig_1} = $trans_temp;
				
				#update the max
				if($trans_temp > $max_trans_prod) {
					
					$max_trans_prod = $trans_temp;
					
				}
				
			}
			
			#get the enzyme levels
			for(my $d = 0; $d <= $#enz_level_ind; $d++) {
				
				#get the item which will be the transcript level
				chomp(my $enz_temp = $items[$enz_level_ind[$d]]);
				
				if($enz_temp != 0) {
					
					$enz_is_prod{$enz[$d]} = 1;
					
				}
				
				#add to appropriate hash
				$enz_level{$enz[$d]}{$timepoint}{$lig_1} = $enz_temp;
				
				#update the max
				if($enz_temp > $max_enz_level) {
					
					$max_enz_level = $enz_temp;
					
				}
				
			}
			
			#get the transcript levels
			for(my $e = 0; $e <= $#enz_prod_ind; $e++) {
				
				#get the item which will be the transcript level
				chomp(my $enz_temp = $items[$enz_prod_ind[$e]]);
				
				if($enz_temp != 0) {
					
					$enz_is_prod{$enz[$e]} = 1;
					
				}
				
				#add to appropriate hash
				$enz_prod{$enz[$e]}{$timepoint}{$lig_1} = $enz_temp;
				
				#update the max
				if($enz_temp > $max_enz_prod) {
					
					$max_enz_prod = $enz_temp;
					
				}
				
			}
		
		}
		
	}
	
	#increment the line count
	$line_count++;
	
}

#at this point should have all the data we want to plot
#let's create the output file so we can write tables and plot them

#create the workbook
my $workbook = Excel::Writer::XLSX->new('solution_'.$soln_num.'_plots_'.$gate.'.xlsx');

#create the worksheet for the transcript plots
my $trans_prod_plots = $workbook->add_worksheet('Trans_Prod');
my $trans_level_plots = $workbook->add_worksheet('Trans_Level');

#create the worksheet for the enzyme plots
my $enz_prod_plots = $workbook->add_worksheet('Enz_Prod');
my $enz_level_plots = $workbook->add_worksheet('Enz_Level');

#lets first create the tables we will use for the plots
#will be a total of four tables per worksheet

#table number counter
#used to get the column numbers for formatting where things are written 
my $table_num = 0;

printf "\nTranscript production data...\n";

#transcript production plot
for(my $f = 0; $f <= $#ligands; $f++) {
	
	printf "getting results for ligand conditions %s(%s)...\n", $ligands[$f], $f;
	
	my $title = $ligands[$f].' Production';

	#titles will always go on row 0
	$trans_prod_plots->write_string( 0, ($table_num*4), $title);
	
	#write the headers
	#headers will always go on row 1
	$trans_prod_plots->write_string( 1, ($table_num * 4 + 0), "Biopart");
	$trans_prod_plots->write_string( 1, ($table_num * 4 + 1), "Time");
	$trans_prod_plots->write_string( 1, ($table_num * 4 + 2), "Production");
	
	#counter to get the number of bioparts
	#used for creating lines for the plot
	my $part_count = 0;
	
	#create a chart to show this data
	my $trans_prod_chart = $workbook->add_chart( type => 'scatter', embedded => 1 );
	$trans_prod_chart->set_title( name 	=>	$title );
	$trans_prod_chart->set_x_axis( 
									name	=>	'Time Point (relative)', 
									min		=>	0,
									max		=>	$time_points[$#time_points]
								);
	$trans_prod_chart->set_y_axis(
									name	=>	'Transcript Production Rate (relative)',
									min		=>	0,
									max		=> 	(int($max_trans_prod) + 1)
								);
	
	for(my $h = 0; $h <= $#trans; $h++) {
		
		#only write to table if the transcript is produced
		if($trans_is_prod{$trans[$h]} eq 1) {
	
			for(my $i = 0; $i <= $#time_points; $i++) {
				
				$trans_prod_plots->write_string(($part_count * $#time_points + $i + 2 + $part_count), ($table_num * 4 + 0), $trans[$h]);
				$trans_prod_plots->write_number(($part_count * $#time_points + $i + 2 + $part_count), ($table_num * 4 + 1), $time_points[$i]);
				$trans_prod_plots->write_number(($part_count * $#time_points + $i + 2 + $part_count), ($table_num * 4 + 2), $trans_prod{$trans[$h]}{$time_points[$i]}{$ligands[$f]});
				
			}
			
			#time column
			my $time_col = $num_to_alpha{($table_num * 4 + 1)};
			
			#production column
			my $prod_col = $num_to_alpha{($table_num * 4 + 2)};
			
			#minimum row
			my $min_row = $part_count * $#time_points + 3 + $part_count;
			
			#maximum row
			my $max_row = $part_count * $#time_points + $#time_points + 3 + $part_count;
			
			#add this data as a line to the chart
			$trans_prod_chart->add_series(
				name 		=> $trans[$h],
				categories 	=> '=Trans_Prod!$'.$time_col.'$'.$min_row.':$'.$time_col.'$'.$max_row,
				values		=> '=Trans_Prod!$'.$prod_col.'$'.$min_row.':$'.$prod_col.'$'.$max_row,
			);
			
			#increment the number of parts in the table
			$part_count++;
			
		} else {
			
			#otherwise not produced so don't write
		
		}
		
	}
	
	#add the chart
	$trans_prod_plots->insert_chart( ($part_count * $#time_points + $#time_points + 4), ($table_num * 4), $trans_prod_chart);
			
	$table_num++;

}

printf "\nTranscript level data...\n";

#reset table_num
$table_num = 0;

#transcript level plot
for(my $j = 0; $j <= $#ligands; $j++) {
	
	printf "getting results for ligand conditions %s(%s)...\n", $ligands[$j], $j;
	
	my $title = $ligands[$j].' Level';
	
	#titles will always go on row 0
	$trans_level_plots->write_string( 0, ($table_num*4), $title);
	
	#write the headers
	#headers will always go on row 1
	$trans_level_plots->write_string( 1, ($table_num * 4 + 0), "Biopart");
	$trans_level_plots->write_string( 1, ($table_num * 4 + 1), "Time");
	$trans_level_plots->write_string( 1, ($table_num * 4 + 2), "Level");
	
	#counter to get the number of bioparts
	#used for creating lines for the plot
	my $part_count = 0;
	
	#create a chart to show this data
	my $trans_level_chart = $workbook->add_chart( type => 'scatter', embedded => 1 );
	$trans_level_chart->set_title( name => $title );
	$trans_level_chart->set_x_axis( 
									name	=>	'Time Point (relative)', 
									min		=>	0,
									max		=>	$time_points[$#time_points]
								);
	$trans_level_chart->set_y_axis(
									name	=>	'Transcript Level (relative)',
									min		=>	0,
									max		=> 	(int($max_trans_level) + 1)
								);
	
	for(my $m = 0; $m <= $#trans; $m++) {
		
		#only write to table if the transcript has a non-zero level
		if($trans_is_prod{$trans[$m]} eq 1) {
			
			for(my $n = 0; $n <= $#time_points; $n++) {
				
				$trans_level_plots->write_string(($part_count * $#time_points + $n + 2 + $part_count), ($table_num * 4 + 0), $trans[$m]);
				$trans_level_plots->write_number(($part_count * $#time_points + $n + 2 + $part_count), ($table_num * 4 + 1), $time_points[$n]);
				$trans_level_plots->write_number(($part_count * $#time_points + $n + 2 + $part_count), ($table_num * 4 + 2), $trans_level{$trans[$m]}{$time_points[$n]}{$ligands[$j]});
				
			}
			
			#time column
			my $time_col = $num_to_alpha{($table_num * 4 + 1)};
			
			#level column
			my $level_col = $num_to_alpha{($table_num * 4 + 2)};
			
			#minimum row
			my $min_row = $part_count * $#time_points + 3 + $part_count;
			
			#maximum row
			my $max_row = $part_count * $#time_points + $#time_points + 3 + $part_count;
					
			#add this data as a line to the chart
			$trans_level_chart->add_series(
				name 		=> $trans[$m],
				categories 	=> '=Trans_Level!$'.$time_col.'$'.$min_row.':$'.$time_col.'$'.$max_row,
				values		=> '=Trans_Level!$'.$level_col.'$'.$min_row.':$'.$level_col.'$'.$max_row,
			);
			
			#increment the number of parts in the table
			$part_count++;
			
		} else {
			
			#otherwise not produced so don't write
	
		}
		
	}
	
	#add the chart
	$trans_level_plots->insert_chart( ($part_count * $#time_points + $#time_points + 4), ($table_num * 4), $trans_level_chart);
	
	$table_num++;
	
}

printf "\nEnzyme production data...\n";
#reset table_num
$table_num = 0;

#enzyme production plot
for(my $p = 0; $p <= $#ligands; $p++) {
	
	printf "getting results for ligand conditions %s(%s)...\n", $ligands[$p], $p;
	
	my $title = $ligands[$p].' Production';
	
	#titles will always go on row 0
	$enz_prod_plots->write_string( 0, ($table_num*4), $title);
	
	#write the headers
	#headers will always go on row 1
	$enz_prod_plots->write_string( 1, ($table_num * 4 + 0), "Biopart");
	$enz_prod_plots->write_string( 1, ($table_num * 4 + 1), "Time");
	$enz_prod_plots->write_string( 1, ($table_num * 4 + 2), "Production");
	
	#counter to get the number of bioparts
	#used for creating lines for the plot
	my $part_count = 0;
	
	#create a chart to show this data
	my $enz_prod_chart = $workbook->add_chart( type => 'scatter', embedded => 1 );
	$enz_prod_chart->set_title( name => $title );
	$enz_prod_chart->set_x_axis( 
									name	=>	'Time Point (relative)', 
									min		=>	0,
									max		=>	$time_points[$#time_points]
								);
	$enz_prod_chart->set_y_axis(
									name	=>	'Enzyme Production Rate (relative)',
									min		=>	0,
									max		=> 	(int($max_enz_prod) + 1)
								);
	
	for(my $r = 0; $r <= $#enz; $r++) {
		
		#only write to table if the transcript is produced
		if($enz_is_prod{$enz[$r]} eq 1) {
			
			for(my $s = 0; $s <= $#time_points; $s++) {
				
				$enz_prod_plots->write_string(($part_count * $#time_points + $s + 2 + $part_count), ($table_num * 4 + 0), $enz[$r]);
				$enz_prod_plots->write_number(($part_count * $#time_points + $s + 2 + $part_count), ($table_num * 4 + 1), $time_points[$s]);
				$enz_prod_plots->write_number(($part_count * $#time_points + $s + 2 + $part_count), ($table_num * 4 + 2), $enz_prod{$enz[$r]}{$time_points[$s]}{$ligands[$p]});
				
			}
			
			#time column
			my $time_col = $num_to_alpha{($table_num * 4 + 1)};
			
			#production column
			my $prod_col = $num_to_alpha{($table_num * 4 + 2)};
			
			#minimum row
			my $min_row = $part_count * $#time_points + 3 + $part_count;
			
			#maximum row
			my $max_row = $part_count * $#time_points + $#time_points + 3 + $part_count;
			
			#add this data as a line to the chart
			$enz_prod_chart->add_series(
				name 		=> $enz[$r],
				categories 	=> '=Enz_Prod!$'.$time_col.'$'.$min_row.':$'.$time_col.'$'.$max_row,
				values		=> '=Enz_Prod!$'.$prod_col.'$'.$min_row.':$'.$prod_col.'$'.$max_row,
			);
			
			#increment the number of parts in the table
			$part_count++;
			
		} else {
			
			#otherwise not produced so don't write
			
		}
		
	}
	
	#add the chart
	$enz_prod_plots->insert_chart( ($part_count * $#time_points + $#time_points + 4), ($table_num * 4), $enz_prod_chart);
	
	$table_num++;
	
}

printf "\nEnzyme level data...\n";
#reset table_num
$table_num = 0;

#transcript level plot
for(my $t = 0; $t <= $#ligands; $t++) {
	
	printf "getting results for ligand conditions %s(%s)...\n", $ligands[$t], $t;
	
	my $title = $ligands[$t].' Level';
	
	#titles will always go on row 0
	$enz_level_plots->write_string( 0, ($table_num*4), $title);
	
	#write the headers
	#headers will always go on row 1
	$enz_level_plots->write_string( 1, ($table_num * 4 + 0), "Biopart");
	$enz_level_plots->write_string( 1, ($table_num * 4 + 1), "Time");
	$enz_level_plots->write_string( 1, ($table_num * 4 + 2), "Level");
	
	#counter to get the number of bioparts
	#used for creating lines for the plot
	my $part_count = 0;
	
	#create a chart to show this data
	my $enz_level_chart = $workbook->add_chart( type => 'scatter', embedded => 1 );
	$enz_level_chart->set_title( name => $title );
	$enz_level_chart->set_x_axis( 
									name	=>	'Time Point (relative)', 
									min		=>	0,
									max		=>	$time_points[$#time_points]
								);
	$enz_level_chart->set_y_axis(
									name	=>	'Enzyme Level (relative)',
									min		=>	0,
									max		=> 	(int($max_enz_level) + 1)
								);
	
	for(my $v = 0; $v <= $#enz; $v++) {
		
		#only write to table if the transcript has a non-zero level
		if($enz_is_prod{$enz[$v]} eq 1) {
			
			for(my $w = 0; $w <= $#time_points; $w++) {
				
				$enz_level_plots->write_string(($part_count * $#time_points + $w + 2 + $part_count), ($table_num * 4 + 0), $enz[$v]);
				$enz_level_plots->write_number(($part_count * $#time_points + $w + 2 + $part_count), ($table_num * 4 + 1), $time_points[$w]);
				$enz_level_plots->write_number(($part_count * $#time_points + $w + 2 + $part_count), ($table_num * 4 + 2), $enz_level{$enz[$v]}{$time_points[$w]}{$ligands[$t]});
				
			}
			
			#time column
			my $time_col = $num_to_alpha{($table_num * 4 + 1)};
			
			#level column
			my $level_col = $num_to_alpha{($table_num * 4 + 2)};
			
			#minimum row
			my $min_row = $part_count * $#time_points + 3 + $part_count;
			
			#maximum row
			my $max_row = $part_count * $#time_points + $#time_points + 3 + $part_count;
					
			#add this data as a line to the chart
			$enz_level_chart->add_series(
				name 		=> $enz[$v],
				categories 	=> '=Enz_Level!$'.$time_col.'$'.$min_row.':$'.$time_col.'$'.$max_row,
				values		=> '=Enz_Level!$'.$level_col.'$'.$min_row.':$'.$level_col.'$'.$max_row,
			);
			
			#increment the number of parts in the table
			$part_count++;
			
		} else {
			
			#otherwise not produced so don't write
			
		}
		
	}
	
	#add the chart
	$enz_level_plots->insert_chart( ($part_count * $#time_points + $#time_points + 4), ($table_num * 4), $enz_level_chart);
	
	$table_num++;
	
}

printf "done!\n";