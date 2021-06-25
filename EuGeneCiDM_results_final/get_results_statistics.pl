#!usr/bin/perl -w 

#Written by: Wheaton Schroeder
#Latest Version: 04/27/2021

use strict;

#Written to read the circuit design result files to get some automated summaries
#Since these files were originally built for human reading, but very difficult
#when we get into dozens of solutions

my $gate = "buffer_Zn";

printf "\nresult for gate: %s\n", $gate;

#we are going to read the file first
open(HUMANFILE, "<circuit_designs_file_".$gate.".txt") or die "Could not open the human-readable file, reason: $!\n";
chomp(my @text_lines = <HUMANFILE>);

my $text = join "\n", @text_lines;

#replace sections of 6 or more newlines with only three
$text =~ s/\n{6}/\n\n\n/g;

#split the file into solutions
my @solns = split /\n\n\n/, $text;

#remove the file header
shift @solns;

#for each solution:
#1. Identify if it is a solution or just a note that there is no more solutions of that size.
#
#2. IF a solution get the objective value, size, status, time, and parts
#
#Things to track:
#	1. number of solutions of a particular size (hash %size_tallies)
my %size_tallies = ( );

#	2. minimum size ($min_size)
my $min_size = 10000000;

#	3. maximum size ($max_size)
my $max_size = 0;

#	4. array of values for time so can do statistics (%soln_times)
my %soln_times = ( );

#	5. maximum time ($max_time)
my $max_time = 0;

#	6. minimum time ($min_time)
my $min_time = 10000000;

#	7. tallies of different solver statuses (hash %solver_status)
my %solver_status = ( );

#	8. tallies of different model statuses (hash %model_status)
my %model_status = ( );

#	9. list of objective values (hash %obj_vals)
my %obj_vals;

#  10. Maximum objective value ($max_obj)
my $max_obj = -10000000;

#  11. Solution number with maximum value ($max_obj_soln)
my $max_obj_soln = 0;

#  12. Minimum objective value ($min_obj)
my $min_obj = 10000000;

#  13. Solution number with minimum value ($min_obj_soln)
my $min_obj_soln = 0;

#  14. store part tallies where the key is the part name (hash %part_tallies)
my %part_tallies = ( );

#  16. Set of solution numbers
my @soln_num = ( ); 

#  17. Solution size
my %soln_size = ( );

#  18. Solution design
my %soln_design = ( );

#  21. if a fatal error had occured
my $fatal_error = "no";

#  22. tallies of model statuses
my %model_tallies = ( );

#  23. tallies of solver statuses
my %solver_tallies = ( );

#4. IF a note that no more solutions of that size, get the time to determine the statistics on this

#	1. store times for determining no solution (array @det_times)
my @det_times = ( );

no warnings 'uninitialized';

for(my $a = 0; $a <= $#solns; $a++) {
	
	#split the solution into lines
	my @soln_lines = split /\n/, $solns[$a];
	
	#remove empty elements if present
	my $done = 0;
	
	while($done == 0) {
		
		if(! defined $soln_lines[0]) {
			
			shift @soln_lines;
			
		}
		
		if($soln_lines[0] =~ /^$/) {
			
			shift @soln_lines;
			
		} else {
			
			$done = 1;
			
		}
		
	}
	
	#can determine if a true solution or not by the first line
	if($soln_lines[0] =~ /Solution Number:\s+(\d+)/i) {
		
		#add the solution number
		my $num = $1;
		push @soln_num, $num;
		
		#get the objective value
		if($soln_lines[1] =~ /Objective Value:\s+((-)?\d+\.\d+)/i) {
			
			my $obj = $1;
			$obj_vals{$num} = $obj;
			
			#update min objective value if nencessary
			if($obj < $min_obj) {
				
				$min_obj = $obj;
				$min_obj_soln = $num;
				
			}
			
			#update max objective value if necessary
			if($obj > $max_obj) {
				
				$max_obj = $obj;
				$max_obj_soln = $num;
				
			}
			
		}
		
		#get the circuit size
		if($soln_lines[2] =~ /Circuit Size:\s+(\d+)(\.\d+)?/i) {
			
			my $size = $1;
			$soln_size{$num} = $size;
			
			#update min size if nencessary
			if($size < $min_size) {
				
				$min_size = $size;
				
			}
			
			#update max size if necessary
			if($size > $max_size) {
				
				$max_size = $size;
				
			}
			
			#update the size tallies
			#if already tallied
				if(exists $size_tallies{$size}) {
					
					#update the tally
					$size_tallies{$size}++;
				
				} else {
					
					#otherwise call it as the first instance of that part
					$size_tallies{$size} = 1;
					
				}
			
		}
		
		#get the model status
		if($soln_lines[3] =~ /model status:\s+(\d+)(\.\d+)?/i) {
			
			my $status = $1;
			$model_status{$num} = $status;
			
			#update the tallies
			#if already tallied
			if(exists $model_tallies{$status}) {
				
				#update the tally
				$model_tallies{$status}++;
			
			} else {
				
				#otherwise call it as the first instance of that part
				$model_tallies{$status} = 1;
				
			}
			
		}
		
		#get the solver status
		if($soln_lines[4] =~ /solver status:\s+(\d+)(\.\d+)?/i) {
			
			my $status = $1;
			$solver_status{$num} = $status;
			
			#update the tallies
			#if already tallied
			if(exists $model_tallies{$status}) {
				
				#update the tally
				$solver_tallies{$status}++;
			
			} else {
				
				#otherwise call it as the first instance of that part
				$solver_tallies{$status} = 1;
				
			}
			
		}
		
		#get the time taken to arrive at the solution
		if($soln_lines[5] =~ /time to solution:\s+(\d+\.\d+)\ssecond/i) {
			
			my $time = $1;
			$soln_times{$num} = $time;
			
			#update min time if nencessary
			if($time < $min_time) {
				
				$min_time = $time;
				
			}
			
			#update max size if necessary
			if($time > $max_time) {
				
				$max_time = $time;
				
			}
			
		}
		
		#element 6 should be blank
		
		#elements 7, 8, and 9 will be simple formatting		
		
		#remaining lines will describe the circuit design
		my $design = "";
		
		#tally of custom parts in the circuit
		my $cust_tally = 0;
		
		for(my $b = 10; $b <= $#soln_lines; $b++) {
			
			#get the solution line
			my $temp = $soln_lines[$b];
			
			#format a bit
			#will be spaces between elements of triads, "|" between triads
			
			#get the part names and part tallies
			my @parts = split /\s+/, $temp;
			
			#update the counts
			for(my $c = 0; $c <= 2; $c++) {
				
				#if already tallied
				if(exists $part_tallies{$parts[$c]}) {
					
					#update the tally
					$part_tallies{$parts[$c]}++;
				
				} else {
					
					#otherwise call it as the first instance of that part
					$part_tallies{$parts[$c]} = 1;
					
				}
				
				#check if the part is a custom one
				if($parts[$c] =~ /^gene_(A|B|C|D|E)$/i) {
					
					#then is a custom part
					$cust_tally++;
					
				}
				
			}
			
			my $design .= "|".$temp;
			
		}
		
		#save the design
		$soln_design{$num} = $design;
		
	} elsif ($soln_lines[0] =~ /No remaining solutions/) {
		
		#line 4 (element 3) will have the solve time
		if($soln_lines[3] =~ /"Time to Determine:\s+(\d+\.\d+) second"/i) {
			
			push @det_times, $1;
			
		}
		
		#get the solver status from element 2
		if($soln_lines[2] =~ /Solver Status:\s+(\d+)(\.\d+)?/i) {
			
			my $status = $1;
			
			#update the tallies
			#if already tallied
			if(exists $model_tallies{$status}) {
				
				#update the tally
				$solver_tallies{$status}++;
			
			} else {
				
				#otherwise call it as the first instance of that part
				$solver_tallies{$status} = 1;
				
			}
			
		}
		
		#get the model status from element 1
		if($soln_lines[1] =~ /model status:\s+(\d+)(\.\d+)?/i) {
			
			my $status = $1;
			
			#update the tallies
			#if already tallied
			if(exists $model_tallies{$status}) {
				
				#update the tally
				$model_tallies{$status}++;
			
			} else {
				
				#otherwise call it as the first instance of that part
				$model_tallies{$status} = 1;
				
			}
			
		}
		
	} else {
		
		#if here then whas a fatal error line
		$fatal_error = "yes";
	
	}
	
}

#report the total number of solutions, remember perl counts from zero
my $total_solns = $#soln_num + 1;

#get the average solution time
my $av_soln_time = 0;

for(my $d = 0; $d <= $#soln_num; $d++) {
	
	$av_soln_time += $soln_times{$soln_num[$d]};
	
}

$av_soln_time = $av_soln_time / $#soln_num;

my $mode_soln_size = 0;

my $biggest_tally = 0;

foreach my $key (keys %size_tallies) {
	
	if($size_tallies{$key} > $biggest_tally) {
		
		$mode_soln_size = $key;
		$biggest_tally = $size_tallies{$key};
		
	}
	
}

#lets now do some basic reporting to make sure 
printf "\n\nNumber of solutions: %s\n\n", $total_solns;
printf "minimum time to solution: %s\n", $min_time; 
printf "maximum time to solution: %s\n", $max_time; 
printf "average solution time: %s\n", $av_soln_time;
printf "minimum size: %s\n", $min_size;
printf "maximum size: %s\n", $max_size;
printf "most common solution size: %s\n", $mode_soln_size;
printf "minimum objective value: %s (soln #%s)\n", $min_obj, $min_obj_soln;
printf "maximum objective value: %s (soln #%s)\n", $max_obj, $max_obj_soln;
printf "\nmodel status tallies: \n";

foreach my $key (keys %model_tallies) {
	
	printf "%s: %s\n", $key, $model_tallies{$key};
	
}

printf "\n";

printf "\nsolver status tallies: \n";

foreach my $key (keys %solver_tallies) {
	
	printf "%s: %s\n", $key, $solver_tallies{$key};
	
}

printf "\n";

printf "fatal error? %s\n\n", $fatal_error;
