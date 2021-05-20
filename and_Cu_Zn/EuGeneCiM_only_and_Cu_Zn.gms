******************************************************************************************************************************
*Code to execute the EUkaryotic GENEtic CIrcuit Design (EUGENECID) and EUkaryotic GENEtic CIrcuit Simulation (EUGENECIM)
*Tools to accompany the publication:
*
*RESEARCH ARTICLE
*TITLE: OPTIMIZATION-BASED EUKARYOTIC GENETIC CIRCUIT DESIGN (EUGENECID) 
*		AND MODELING (EUGENECIM) TOOLS: COMPUTATIONAL APPROACH TO SYNTHETIC BIOLOGY
*AUTHORS: WHEATON L SCHROEDER, ANNA S BABER, RAJIB SAHA
*JOURNAL: ISCIENCE
*YEAR: 2021
*
*ASSOCIATED PROTOCOL PAPER
*TITLE: USING EUGENECID AND EUGENECIM COMPUTATIONAL TOOLS FOR SYNTHETIC BIOLOGY
*AUTHORS: WHEATON L SCHROEDER, ANNA S BABER, RAJIB SAHA
*JOURNAL: STAR PROTOCOLS
*YEAR: 2021
*
*EuGeneCiD/EUGENECIM loosly based on the publication:
*"Course-grained optimization-driven design and piecewise linear modeling of synthetic genetic 
*"circuits" by Ali R. Zomorrodi and Costas D. Maranas 2014
*
*Code Authors: Wheaton L. Schroeder and Anna Baber
*Latest Version: Version 4.0
*Version Date: 05/20/2021  
*
******************************************************************************************************************************

$INLINECOM /*  */
$onlisting
$offdigit
$oneolcom
$onempty

OPTIONS

    decimals = 8
	mip = cplex
	limrow = 1000

;

************************************************* DEFINE SETS ****************************************************************

SETS

	A					set of all molecules includes proteins ligands promotors ribosome binding sites and transcripts
$include "and_Cu_Zn/all_molecules_and_Cu_Zn.txt"

	J(A)				set of protein coding regions (transcripts)
$include "and_Cu_Zn/transcripts_and_Cu_Zn.txt"

	E(A)				set of proteins enzymes or polypeptides in the model
$include "and_Cu_Zn/enzymes_and_Cu_Zn.txt"

	P					set of promotors
$include "and_Cu_Zn/promoters_and_Cu_Zn.txt"

	Ed(E)				set of proteins which are desired to be the response signal to ligand signals
$include "and_Cu_Zn/desired_enzymes_and_Cu_Zn.txt"
 
	L1(A)				set of ligands
$include "and_Cu_Zn/Ligands_and_Cu_Zn.txt"

	Ld(L1)				set of ligands for which a response is desired
$include "and_Cu_Zn/desired_ligands_and_Cu_Zn.txt"

	T					set of terminators
$include "and_Cu_Zn/terminators_and_Cu_Zn.txt"

	tau					set of time values
$include "and_Cu_Zn/time_set_and_Cu_Zn.txt"

;

alias(L2,L1);
alias(A1,A);
alias(J1,J);
alias(Ld1,Ld);
alias(E1,E);
alias(E2,E);

******************************************** DEFINE PARAMETER ****************************************************************

PARAMETERS

	S(P)						strenth of promotor
$include "and_Cu_Zn/promoter_strength_and_Cu_Zn.txt"

	eta(J)						efficiency of ribosome binding sites
$include "and_Cu_Zn/translation_efficiency_and_Cu_Zn.txt"

	Z(P)						normal state of promotor 1 of normally on 0 otherwise
$include "and_Cu_Zn/promoter_normal_state_and_Cu_Zn.txt"

	Zeta(E)						normal state of proteins 1 if normally active 0 otherwise
$include "and_Cu_Zn/protein_normal_state_and_Cu_Zn.txt"

	I(P,A)						promotor molecule interaction -1 if protein j represses promotor p 0 if no effect 1 if activation
$include "and_Cu_Zn/promoter_ligand_interactions_and_Cu_Zn.txt"
	
	B(E,A)						protein molecule interaction -1 if molecule A represses enzyme E 0 if no effect 1 if activation
$include "and_Cu_Zn/protein_ligand_interactions_and_Cu_Zn.txt"

	G(T)						strength of terminators
$include "and_Cu_Zn/terminator_strength_and_Cu_Zn.txt"

	R(E)						enzyme degradation rate
$include "and_Cu_Zn/protein_degradation_and_Cu_Zn.txt"

	Np_max						maximum number of times any promotor can be used in a circuit design
	
	Nj_max						maximum number of times any protein coding region can be used in a circuit design
	
	Nt_max						maximum number of times any terminator can be used in a circuit design
	
	Ncircuit_max				maximum number of promotor protein RBS combinations that can be used in total in a circuit design

	Max_size_allow				maximum circuit size allowed incremented with solutions

	theta(E)					threshold for concentration needed for enzymes to be active 
$include "and_Cu_Zn/expression_threshold_and_Cu_Zn.txt"

	lambda(Ld,Ld1,Ed)			gives the values of W that should exist for desired proteins responses to desired ligands
$include "and_Cu_Zn/logic_table_and_Cu_Zn.txt"

	rho(J,E)					maps transcripts to enzymes
$include "and_Cu_Zn/trans_to_enzyme_and_Cu_Zn.txt"

	sigma(A,A)					determines if enzyme E is the same as enzyme E
$include "and_Cu_Zn/self_and_Cu_Zn.txt"

	temp						temporarily store a number
	
	temp2						temporarily store a number
	
	F(P)						promotor leakiness because usually impossible to completely stop production of an enzyme
$include "and_Cu_Zn/promoter_leakiness_and_Cu_Zn.txt"

	iter						number of iterations for attempting solution
	
	iter_temp					advances iter
	
	iter_max					maximum number of iterations
	
	ECTrans(E,L1,L2) 			enzyme concentration transfered between time points
	
	RNATrans(J,T,L1,L2)			RNA transcript level transfered between time points 
	
	M_design(P,J,T)				Parameter array which stores the M values from EuGeneCiD for EUGENECIM to use so as not to use so many variables
	
	count						just a spare
	
	num_solns					the number of the previous solutions 
	
	H(P,A)						strength of protein and ligand interactions with promoters
$include "and_Cu_Zn/promoter_ligand_strength_and_Cu_Zn.txt"
	
	Q(E,A) 						strength of ligand-protein interactions
$include "and_Cu_Zn/protein_ligand_strength_and_Cu_Zn.txt"

	done						used to exit loop when needed
	
	done_2						used to exit loop when needed
	
	done_3						used to exit loop when needed
	
	max_solns					number of solutions to find
	
	first_time					stores if solving for the first timepoint currently eg where Ms are not fixed
	
	Ed_val(E)					set of ligands for which a response is desired value of 1 if a desired ligand zero otherwise
$include "and_Cu_Zn/desired_enzymes_and_Cu_Zn_vals.txt"

	Ld_val(L1)					set of ligands for which a response is desired value of 1 if a desired ligand zero otherwise
$include "and_Cu_Zn/desired_ligands_and_Cu_Zn_vals.txt"

;

/*right now only look for one solution*/
max_solns = 1000;

/*give initial condition values for enzyme and RNA carry-over*/
/*note: works a bit like an ODE now...*/
ECTrans(E,L1,L2) = 0;
RNATrans(J,T,L1,L2) = 0;

SCALAR epsilon /1E-4/;
SCALAR V /1E4/;

*limit circuit size
Np_max = 3;
Nj_max = 3;
Nt_max = 100;

*limit solution size
Ncircuit_max = 1;
Max_size_allow = 10;

num_solns = 0;

*********************************************** MAJOR UPDATE NOTES ***********************************************************

*now we will have two seperate problems/models
*MODEL 1 EUGENECID (EUkaryotic GENEtic CIrcuit Design tool) will design the circuit with no time consideration
*MODEL 2 EUGENECIM (EUkaryotic GENEtic CIrcuit Simulation tool) will model how this 
*	EUGENECIM will be formulated very similarly to EUGENECID except that the enforcement of the logic table
*	will not be a part of the formulation of this tool and will have time delays between various steps of the 
*	central dogma of biology

******************************************** DEFINE VARIABLES ****************************************************************

VARIABLES

	/*objective variables for CID and CIS formulations respectively*/
	Z_M							objective variable for EUGENECIM
	
	/*these variables are shared by CID and CIS formulations*/
	alpha(P,L1,L2)				alpha determies net effect of all inhibition and activation effects on the given promotor
	gamma(E,L1,L2)				gamma determies net effect of all inhibition and activation effects on the given protein

NONNEGATIVE VARIABLES

	/*these variables are shared by CID and CIS formulations*/
	C(E,L1,L2)					proxy concentration of protein j subject to ligands
	phi(J,T,L1,L2)				RNA level
	xi(P,J,T,L1,L2)				value of deliberate transcription	
	
BINARY VARIABLES
	
	/*these variables are shared by CID and CIS*/
	C_plus(E,L1,L2)				determines if the concentration activation threshold is achieved by comparison to Cjt
	W(E,L1,L2)					determines if the protein is active based on inhibition activation transcription and translation
	alpha_plus(P,L1,L2)			value of 1 if the promotor is active value of zero otherwise
	gamma_plus(E,L1,L2)			value of 1 if the protein is active value of zero otherwise
	omega(E,L1,L2)				value of 1 if enzyme protein or polypeptide E is produced 0 otherwise
	kappa(E,L1,L2)				value of 1 is protein is active (e.g. not inhibited or is activated) 0 otherwise

;

****************************************** HARDCODED SOLUTION FOR EUGENECIM IF NOT RUNNING EUGENECID ***********************************************

*grab a solution from EuGeneCiD to put here, or a repressilator
*set all to a default of zero
M_design(P,J,T) = 0;

*then define the actual design
*P_EXO70B1_11/gene_cI/CaMV25St
M_design("P_EXO70B1_11","gene_cI","CaMV25St") = 1;

*P_FRO2/gene_cI/HSPt
M_design("P_FRO2","gene_cI","HSPt") = 1;

*P_lambda/gene_GFP/CaMV25St
M_design("P_lambda","gene_GFP","CaMV25St") = 1;

******************************************** DEFINE EQUATIONS ****************************************************************

EQUATIONS

	/*simulation equations for EUGENECIM*/
	obj_M							objective function
	prot_conc_M(E,Ld,Ld1)			protein concentration
	trans_level_M(J,T,Ld,Ld1)		find the transcript level for each transcript
	
	produced0_M(E,Ld,Ld1)			determines if a polypeptide is produced
	produced1_M(E,Ld,Ld1)			determines if a polypeptide is produced
	
	act_threshold_0_M(E,Ld,Ld1)		constraint used to trigger the activation threshold
	act_threshold_1_M(E,Ld,Ld1)		constraint used to trigger the activation threshold
	
	findalpha(P,Ld,Ld1)				find value of alpha 
	findalpha_sign0(P,Ld,Ld1)		find sign of alpha
	findalpha_sign1(P,Ld,Ld1)		find sign of alpha
	
	findgamma(E,Ld,Ld1)				Find gamma
	findgamma_sign0(E,Ld,Ld1)		find value of gamma_plus
	findgamma_sign1(E,Ld,Ld1)		find value of gamma_plus
	
	prot_could_act0(E,Ld,Ld1)		determintes if protein could be active based on ligand inputs and other active protiens
	prot_could_act1(E,Ld,Ld1)		determintes if protein could be active based on ligand inputs and other active protiens
	prot_could_act2(E,Ld,Ld1)		determintes if protein could be active based on ligand inputs and other active protiens
	
	prot_is_act0(E,Ld,Ld1)			determines if protein is active based on inhibition activation and concentration
	prot_is_act1(E,Ld,Ld1)			determines if protein is active based on inhibition activation and concentration
	prot_is_act2(E,Ld,Ld1)			determines if protein is active based on inhibition activation and concentration
	
;

*********************************************** EUGENECIM ********************************************************************

*objective function
obj_M..								Z_M =e= sum(Ed, sum(Ld, sum(Ld1, C(Ed,Ld,Ld1))));

*determines total effect of present activator and inhibitors on promotors
findalpha(P,Ld,Ld1)..				alpha(P,Ld,Ld1) =e= Z(P) + sum(E, (W(E,Ld,Ld1) * I(P,E) * H(P,E))) + I(P,Ld) * H(P,Ld) + I(P,Ld1) * H(P,Ld1) - (I(P,Ld1) * H(P,Ld1) * sigma(Ld,Ld1));
findalpha_sign0(P,Ld,Ld1)..			alpha(P,Ld,Ld1) =g= epsilon * alpha_plus(P,Ld,Ld1) - V * (1 - alpha_plus(P,Ld,Ld1));
findalpha_sign1(P,Ld,Ld1)..			alpha(P,Ld,Ld1) =l= V * alpha_plus(P,Ld,Ld1);

*get the value of gamma
findgamma(E,Ld,Ld1)..				gamma(E,Ld,Ld1) =e= Zeta(E) + sum(E1, W(E1,Ld,Ld1) * B(E,E1) * Q(E,E1)) + B(E,Ld) * Q(E,Ld) + B(E,Ld1) * Q(E,Ld1) - B(E,Ld1) * Q(E,Ld1) * sigma(Ld,Ld1);
findgamma_sign0(E,Ld,Ld1)..			gamma(E,Ld,Ld1) =g= epsilon * gamma_plus(E,Ld,Ld1) - V * (1 - gamma_plus(E,Ld,Ld1));
findgamma_sign1(E,Ld,Ld1)..			gamma(E,Ld,Ld1) =l= V * gamma_plus(E,Ld,Ld1);

*determines if the protein is produced
prot_could_act0(E,Ld,Ld1)..			kappa(E,Ld,Ld1) =l= omega(E,Ld,Ld1);
prot_could_act1(E,Ld,Ld1)..			kappa(E,Ld,Ld1) =l= gamma_plus(E,Ld,Ld1);
prot_could_act2(E,Ld,Ld1)..			kappa(E,Ld,Ld1) =g= omega(E,Ld,Ld1) + gamma_plus(E,Ld,Ld1) - 1;

*if protein is both produced and is at the necessary concentration threshold for activation then it is active
prot_is_act0(E,Ld,Ld1)..			W(E,Ld,Ld1) =l= kappa(E,Ld,Ld1);
prot_is_act1(E,Ld,Ld1)..			W(E,Ld,Ld1) =l= C_plus(E,Ld,Ld1);
prot_is_act2(E,Ld,Ld1)..			W(E,Ld,Ld1) =g= kappa(E,Ld,Ld1) + C_plus(E,Ld,Ld1) - 1;

*determine the transcript level
trans_level_M(J,T,Ld,Ld1)..			phi(J,T,Ld,Ld1) =e= sum(P, M_design(P,J,T) * alpha_plus(P,Ld,Ld1) * S(P) + M_design(P,J,T) * F(P));

*determine protein concentration as a function of transcription and translation efficiency
*adds that produced from the previous time points RNA with that already present from previous time point
prot_conc_M(E,Ld,Ld1)..				C(E,Ld,Ld1) =e= sum(J, rho(J,E) * eta(J) * sum(T, RNATrans(J,T,Ld,Ld1)));

*make sure that enzyme is noted as being produced if it does have a non-zero concentration
produced0_M(E,Ld,Ld1)..				omega(E,Ld,Ld1) =l= V * ECTrans(E,Ld,Ld1);
produced1_M(E,Ld,Ld1)..				omega(E,Ld,Ld1) =g= epsilon * ECTrans(E,Ld,Ld1);

*determine if concentration is at the level necessary for activity
act_threshold_0_M(E,Ld,Ld1)..		(theta(E) + epsilon) * C_plus(E,Ld,Ld1) =l= ECTrans(E,Ld,Ld1);
act_threshold_1_M(E,Ld,Ld1)..		ECTrans(E,Ld,Ld1) =l= (V - (theta(E) - epsilon)) * C_plus(E,Ld,Ld1) + (theta(E) - epsilon);

*********************************************** DEFINE MODEL ****************************************************************

MODEL EUGENECIM
/

	obj_M
	
	findalpha
	findalpha_sign0
	findalpha_sign1
	
	produced0_M
	produced1_M
	
	prot_conc_M
	trans_level_M
	
	act_threshold_0_M
	act_threshold_1_M
	
	findgamma
	findgamma_sign0
	findgamma_sign1
	
	prot_could_act0
	prot_could_act1
	prot_could_act2
	
	prot_is_act0
	prot_is_act1
	prot_is_act2
	
/
;

EUGENECIM.holdfixed = 1;
EUGENECIM.optfile = 1;

********************************************** INITIALIZE OUTPUT FILES *******************************************************

/*this file will be essentially a table storing state values for various timepoints*/
/*each subsequent time will be added as a new row*/
FILE STATES / circuit_state_and_Cu_Zn_CiM_only.csv /;
PUT STATES;
STATES.pc = 5;
STATES.lw = 30;
STATES.pw = 32767;

*write document description
PUT "GENETIC CIRCUIT BEHAVIOR FOR EACH SOLUTION AT EACH TIME POINT"//;

*make superheader
PUT "SOLN#","TIME","L1","L2","PROMOTOR BEHAVIOR";

*note: need to skip first value so values align after puttin a superheader text
first_time = 1;

*leave spaces to report alpha, alpha_plus, and by putting the necessary amount of commas
*blank cells for alpha
LOOP(P,

	IF(first_time,
	
		first_time = 0;
		
	ELSE

		PUT " ";
	
	);

);

first_time = 1;

*blank cells for alpha_plus
LOOP(P,

	PUT " ";

);

PUT "TRANSCRIPT BEHAVIOR";

first_time = 1;

*blank cells for phi
LOOP(J,

	IF(first_time,
	
		first_time = 0;
		
	ELSE

		PUT " ";
	
	);

);

*blank cells for RNA carry-over
LOOP(J,

	PUT " ";

);

first_time = 1;

PUT "ENZYME BEHAVIOR";

*blank cells for enzyme carry-over
LOOP(E,

	IF(first_time,
	
		first_time = 0;
		
	ELSE

		PUT " ";
	
	);

);

first_time = 1;

*blank cells for gamma
LOOP(E,

	PUT " ";

);

*blank cells for omega
LOOP(E,

	PUT " ";

);

*blank cells for concentration
LOOP(E,

	PUT " ";

);

*blank cells for Y
LOOP(E,

	PUT " ";

);

*blank cells for kappa
LOOP(E,

	PUT " ";

);

*blank cells for W
LOOP(E,

	PUT " ";

);

first_time = 1;

PUT /;

PUT "SOLN#","TIME","L1","L2";

STATES.pc = 2;

*write headers for alpha
*note that the if else statement exists to fix some formatting issues
LOOP(P,

	IF(first_time,
	
		PUT ",","alpha_",P.tl,",";
		first_time = 0;
	
	ELSE
	
		PUT "alpha_",P.tl,",";
	
	);

);

*write headers for alpha_plus
LOOP(P,

	PUT "alpha_plus_",P.tl,",";

);

*write headers for phi
LOOP(J,

	PUT "phi_",J.tl,",";

);

*write headers for RNA carry-over
LOOP(J,

	PUT "RNA_carry_",J.tl,",";

);

*write headers for enzyme carry-over
LOOP(E,

	PUT "Enzyme_carry_",E.tl,",";

);

*write headers for gamma
LOOP(E,

	PUT "gamma_",E.tl,",";

);

*write headers for omega
LOOP(E,

	PUT "omega_",E.tl,",";

);

*write headers for concentration
LOOP(E,

	PUT "C_",E.tl,",";

);

*write headers for C_plus
LOOP(E,

	PUT "C_plus_",E.tl,",";

);

*write headers for kappa
LOOP(E,

	PUT "kappa_",E.tl,",";

);

*write headers for W
LOOP(E,

	PUT "W_",E.tl,",";

);

PUT /;

PUTCLOSE;

/*this file will be essentially a set of tables storing just concentration data*/
/*each subsequent time will be a new table*/
FILE TABLES / circuit_logic_tables_and_Cu_Zn_CiM_only.txt /;
PUT TABLES;
TABLES.pc = 2;
TABLES.pw = 32767;

PUT 'Input Logic Table'/;
PUT '======================================================================='//;
PUT 'Ligands',system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,'Enzymes'/;

*write the header
PUT "Ligand 1",system.tab,"Ligand 2",system.tab,system.tab,system.tab;

LOOP(Ed,

	PUT Ed.tl;

);

PUT /;

*write the logic table lines
LOOP(Ld,

	LOOP(Ld1,
	
		/*write the ligand conditions*/
		PUT Ld.tl,Ld1.tl;
		
		/*write the desired enzyme response to the system*/
		LOOP(Ed, 
		
			PUT lambda(Ld,Ld1,Ed);
			
		);
		
		/*finish with a newline*/
		PUT /;
	
	);
	
);

PUT //;

PUTCLOSE;

******************** SOLVE EUGENECID/S FOR A GIVEN NUMBER OF SOLUTIONS OR UNTIL INFEASIBLE ***********************************

*stores whether or not we are done looking for new circuit designs
done = 0;

*store whether or not we are incrementing solution size
done_2 = 0;

*store whether or not done solvign EuGeneCiM
done_3 = 0;

/*Inside the loop will be solving the simulation problem EUGENECIM if EUGENECID produced and answer*/
LOOP(tau$(not done_3),

	/*solve EUGENECIM for each time point*/
	SOLVE EUGENECIM USING MIP MAXIMIZING Z_M;
	
	/*check if a solution has been found*/
	/*case 1 solution is found need to write current time point to STATES*/
	IF(((EUGENECIM.modelStat eq 1) OR (EUGENECIM.modelStat eq 2) OR (EUGENECIM.modelStat eq 8) OR (EUGENECIM.modelStat eq 15) OR (EUGENECIM.modelStat eq 16) OR (EUGENECIM.modelStat eq 17)),
	
		/*if here then a solution was found to EUGENECIM*/
		
		/*write the results to output files*/
		
		/*write to states file*/
		STATES.ap = 1;
		PUT STATES;
		STATES.pc = 5;
		STATES.pw = 32767;
		
		LOOP(Ld,
		
			LOOP(Ld1,
			
				PUT "0",tau.tl,Ld.tl,Ld1.tl;
			
				/*write alpha values*/
				LOOP(P,
				
					PUT alpha.l(P,Ld,Ld1);
				
				);
				
				/*write alpha_plus values*/
				LOOP(P,
				
					PUT alpha_plus.l(P,Ld,Ld1);
				
				);

				/*write phi values*/
				LOOP(J,

					PUT sum(T, phi.l(J,T,Ld,Ld1));

				);
				
				/*write RNA carry-over values*/
				LOOP(J,

					PUT sum(T, RNATrans(J,T,Ld,Ld1));

				);

				/*write enzyme carry-over values*/
				LOOP(E,

					PUT ECTrans(E,Ld,Ld1);

				);

				/*write gamma values*/
				LOOP(E,

					PUT gamma.l(E,Ld,Ld1);

				);

				/*write omega values*/
				LOOP(E,

					PUT omega.l(E,Ld,Ld1);

				);

				/*write concentration values*/
				LOOP(E,

					PUT C.l(E,Ld,Ld1);

				);

				/*write C_plus values*/
				LOOP(E,

					PUT C_plus.l(E,Ld,Ld1);

				);

				/*write kappa values*/
				LOOP(E,

					PUT kappa.l(E,Ld,Ld1);

				);

				/*write W values*/
				LOOP(E,

					PUT W.l(E,Ld,Ld1);

				);

				PUT /;
				
			);
			
		);
			
		PUTCLOSE;
			
		/*write the logic table results*/
		TABLES.ap = 1;
		PUT TABLES;
		TABLES.pc = 2;
			
		PUT "TIME POINT: ",tau.tl/;
		PUT "-----------------------------------------------------------------------"//;
		
		PUT 'Output Logic Table at Current Time'/;
		PUT 'Ligands',system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,'Enzymes'/;

		/*write the header*/
		PUT "Ligand 1",system.tab,"Ligand 2",system.tab,system.tab,system.tab;

		LOOP(Ed,

			PUT Ed.tl;

		);

		PUT /;

		/*write the logic table lines*/
		LOOP(Ld,

			LOOP(Ld1,
				
				/*write the ligand conditions*/
				PUT Ld.tl,Ld1.tl;
				
				/*write the desired enzyme response to the system*/
				LOOP(Ed, 
					
					PUT W.l(Ed,Ld,Ld1);
					
				);
					
				/*finish with a newline*/
				PUT /;
				
			);
				
		);

		PUT //;
			
		PUT 'Transcript Level at the Current Time Table'/;
		PUT 'Ligands',system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,'Transcripts'/;

		/*write the header*/
		PUT "Ligand 1",system.tab,"Ligand 2",system.tab,system.tab,system.tab;

		LOOP(J,

			PUT J.tl;

		);

		PUT /;

		/*write the logic table lines*/
		LOOP(Ld,
		
			LOOP(Ld1,
			
				/*write the ligand conditions*/
				PUT Ld.tl,Ld1.tl;
					
				/*write the desired enzyme response to the system*/
				LOOP(J, 
					
					PUT sum(T, RNATrans(J,T,Ld,Ld1));
					
				);
					
				/*finish with a newline*/
				PUT /;
				
			);
			
		);

		PUT //;
			
		PUT 'Concentration at the Current Time Table'/;
		PUT 'Ligands',system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,'Enzymes'/;

		/*write the header*/
		PUT "Ligand 1",system.tab,"Ligand 2",system.tab,system.tab,system.tab;

		LOOP(Ed,

			PUT Ed.tl;

		);

		PUT /;

		/*write the logic table lines*/
		LOOP(Ld,
			
			LOOP(Ld1,
			
				/*write the ligand conditions*/
				PUT Ld.tl,Ld1.tl;
					
				/*write the desired enzyme response to the system*/
				LOOP(Ed, 
					
					PUT ECTrans(Ed,Ld,Ld1);
					
				);
					
				/*finish with a newline*/
				PUT /;
				
			);
				
		);

		PUT //;
			
		PUTCLOSE;
			
		/*use the calculated enzyme production  variable C plus holdover enzyme concentration */
		/*minus the degradation of enzymed based on enzyme concentration that already exists second term*/
		ECTrans(E,Ld,Ld1) = (C.l(E,Ld,Ld1) + ECTrans(E,Ld,Ld1)) * (0.5 ** (1 / ((R(E)/2) + epsilon)));
		
		/*ensure that if degredation is larger than carry over and production that carry over is set to zero*/
		LOOP(E,
		
			LOOP(Ld,
			
				LOOP(Ld1,
				
					/*if less than zero reset to zero*/
					IF((ECTrans(E,Ld,Ld1) < 0),
						
						ECTrans(E,Ld,Ld1) = 0;
						
					);
					
				);
				
			);
				
		);
			
		/*use the calculated mRNA production variable phi plus holdover transcript level*/
		/*after having removed a certain amount of it for */
		RNATrans(J,T,Ld,Ld1) = (phi.l(J,T,Ld,Ld1) + RNATrans(J,T,Ld,Ld1)) * (0.5 ** (1 / (sum(P, M_design(P,J,T) * (G(T)/2)) + epsilon))); 
	
		LOOP(J,
		
			LOOP(T,
		
				LOOP(Ld,
				
					LOOP(Ld1,
					
						/*if less than zero reset to zero*/
						IF((RNATrans(J,T,Ld,Ld1) < 0),
						
							RNATrans(J,T,Ld,Ld1) = 0;
							
						);
						
					);
					
				);
					
			);
			
		);
	
	/*case 2 solution is not found yet probably should be found since the design exists*/
	/*we will treat this like an error by reporting it and end solving*/
	ELSE
	
		/*if here we have a problem in that something which previously had a*/
		/*solution now has none at another timepoint*/
		
		/*if here we are done finding new solutions*/
		done_3 = 1;
		
		TABLES.ap = 1;
		PUT TABLES;
		TABLES.pc = 2;
		TABLES.lw = 4;
		
		PUT //;
		PUT "Solution Number 0 at Time Point ", tau.tl," modeling attempted, was unsuccessful (Fatal Error)"/;
		PUT "model status: ",system.tab,EUGENECIM.modelstat/; 
		PUT "solver status: ",system.tab,EUGENECIM.modelstat;
		
		PUTCLOSE;
		
	);
	
);
	
/*by this point done simulating the genetic circuit behavior*/

/*reset the transfer of enzyme concentration for default initial conditions again*/
ECTrans(E,Ld,Ld1) = 0;

/*reset the transfer of RNA transcription for default initial conditions again*/
RNATrans(J,T,Ld,Ld1) = 0;