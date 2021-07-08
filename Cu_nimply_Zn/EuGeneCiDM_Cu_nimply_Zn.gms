******************************************************************************************************************************
*Code to execute the EUkaryotic GENEtic CIrcuit Design (EUGENECID) and EUkaryotic GENEtic CIrcuit Modeling (EUGENECIM)
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
*EuGeneCiD/EuGeneCiM loosly based on the publication:
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
$include "Cu_nimply_Zn/all_molecules_Cu_nimply_Zn.txt"

	J(A)				set of protein coding regions (transcripts)
$include "Cu_nimply_Zn/transcripts_Cu_nimply_Zn.txt"

	E(A)				set of proteins enzymes or polypeptides in the model
$include "Cu_nimply_Zn/enzymes_Cu_nimply_Zn.txt"

	P					set of promotors
$include "Cu_nimply_Zn/promoters_Cu_nimply_Zn.txt"

	Ed(E)				set of proteins which are desired to be the response signal to ligand signals
$include "Cu_nimply_Zn/desired_enzymes_Cu_nimply_Zn.txt"
 
	L1(A)				set of ligands
$include "Cu_nimply_Zn/Ligands_Cu_nimply_Zn.txt"

	Ld(L1)				set of ligands for which a response is desired
$include "Cu_nimply_Zn/desired_ligands_Cu_nimply_Zn.txt"

	T					set of terminators
$include "Cu_nimply_Zn/terminators_Cu_nimply_Zn.txt"

	tau					set of time values
$include "Cu_nimply_Zn/time_set_Cu_nimply_Zn.txt"

	big_set /0*1000/	set with large number of arbitrary elements

	soln_set(big_set)	set of previous solutions

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
$include "Cu_nimply_Zn/promoter_strength_Cu_nimply_Zn.txt"

	eta(J)						efficiency of ribosome binding sites
$include "Cu_nimply_Zn/translation_efficiency_Cu_nimply_Zn.txt"

	Z(P)						normal state of promotor 1 of normally on 0 otherwise
$include "Cu_nimply_Zn/promoter_normal_state_Cu_nimply_Zn.txt"

	Zeta(E)						normal state of proteins 1 if normally active 0 otherwise
$include "Cu_nimply_Zn/protein_normal_state_Cu_nimply_Zn.txt"

	I(P,A)						promotor molecule interaction -1 if protein j represses promotor p 0 if no effect 1 if activation
$include "Cu_nimply_Zn/promoter_ligand_interactions_Cu_nimply_Zn.txt"
	
	B(E,A)						protein molecule interaction -1 if molecule A represses enzyme E 0 if no effect 1 if activation
$include "Cu_nimply_Zn/protein_ligand_interactions_Cu_nimply_Zn.txt"

	G(T)						strength of terminators
$include "Cu_nimply_Zn/terminator_strength_Cu_nimply_Zn.txt"

	R(E)						enzyme degradation rate
$include "Cu_nimply_Zn/protein_degradation_Cu_nimply_Zn.txt"

	Np_max						maximum number of times any promotor can be used in a circuit design
	
	Nj_max						maximum number of times any protein coding region can be used in a circuit design
	
	Nt_max						maximum number of times any terminator can be used in a circuit design
	
	Ncircuit_max				maximum number of promotor protein RBS combinations that can be used in total in a circuit design

	Max_size_allow				maximum circuit size allowed incremented with solutions

	theta(E)					threshold for concentration needed for enzymes to be active 
$include "Cu_nimply_Zn/expression_threshold_Cu_nimply_Zn.txt"

	lambda(Ld,Ld1,Ed)			gives the values of W that should exist for desired proteins responses to desired ligands
$include "Cu_nimply_Zn/logic_table_Cu_nimply_Zn.txt"

	rho(J,E)					maps transcripts to enzymes
$include "Cu_nimply_Zn/trans_to_enzyme_Cu_nimply_Zn.txt"

	sigma(A,A)					determines if enzyme E is the same as enzyme E
$include "Cu_nimply_Zn/self_Cu_nimply_Zn.txt"

	temp						temporarily store a number
	
	F(P)						promotor leakiness because usually impossible to completely stop production of an enzyme
$include "Cu_nimply_Zn/promoter_leakiness_Cu_nimply_Zn.txt"

	iter						number of iterations for attempting solution
	
	iter_temp					used to advance iter
	
	iter_max					maximum number of iterations
	
	ECTrans(E,L1,L2) 			enzyme concentration transfered between time points
	
	RNATrans(J,T,L1,L2)			RNA transcript level transfered between time points 
	
	pastMs(P,J,T,big_set)		array that stores the past m values for int cut
	
	M_design(P,J,T)				Parameter array which stores the M values from EuGeneCiD for EuGeneCiM to use so as not to use so many variables
	
	count						just a spare
	
	num_solns					the number of the previous solutions 
	
	H(P,A)						strength of protein and ligand interactions with promoters
$include "Cu_nimply_Zn/promoter_ligand_strength_Cu_nimply_Zn.txt"
	
	Q(E,A) 						strength of ligand-protein interactions
$include "Cu_nimply_Zn/protein_ligand_strength_Cu_nimply_Zn.txt"

	done						used to exit loop when needed
	
	done_2						used to exit loop when needed
	
	done_3						used to exit loop when needed
	
	max_solns					number of solutions to find
	
	first_time					stores if solving for the first timepoint currently eg where Ms are not fixed
	
	Ed_val(E)					set of ligands for which a response is desired value of 1 if a desired ligand zero otherwise
$include "Cu_nimply_Zn/desired_enzymes_Cu_nimply_Zn_vals.txt"

	Ld_val(L1)					set of ligands for which a response is desired value of 1 if a desired ligand zero otherwise
$include "Cu_nimply_Zn/desired_ligands_Cu_nimply_Zn_vals.txt"

;

/*right now only look for one solution*/
max_solns = 1000;

/*give initial condition values for enzyme and RNA carry-over*/
/*note: works a bit like an ODE now...*/
ECTrans(E,L1,L2) = 0;
RNATrans(J,T,L1,L2) = 0;

SCALAR epsilon /1E-3/;
SCALAR V /1E3/;

/*right now no solutions exist*/
soln_set(big_set) = no;
pastMs(P,J,T,big_set) = no;

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
*MODEL 2 EUGENECIM (EUkaryotic GENEtic CIrcuit Modeling tool) will model how this 
*	EUGENECIM will be formulated very similarly to EUGENECID except that the enforcement of the logic table
*	will not be a part of the formulation of this tool and will have time delays between various steps of the 
*	central dogma of biology

******************************************** DEFINE VARIABLES ****************************************************************

VARIABLES

	/*objective variables for CID and CIS formulations respectively*/
	Z_D							objective variable for EUGENECID
	Z_M							objective variable for EUGENECIM
	
	/*these variables are shared by CID and CIS formulations*/
	alpha(P,L1,L2)				alpha determies net effect of all inhibition and activation effects on the given promotor
	gamma(E,L1,L2)				gamma determies net effect of all inhibition and activation effects on the given protein

NONNEGATIVE VARIABLES

	/*these variables are shared by CID and CIS formulations*/
	C(E,L1,L2)					proxy concentration of protein j subject to ligands
	phi(J,T,L1,L2)				RNA level
	xi(P,J,T,L1,L2)				value of deliberate transcription	
	
	/*Enzyme attribution to other enzymes partial variables*/
	/*note no longer tracking sign because this might result in sign cancellation*/
	D(E,E1)
	K(E,E1)
	U(E,E1,E1)
	U_prime(E,E1)
	X(E,E1)
	chi(E,E1,P,J,T)
	
	/*Enzyme attribution to other enzymes combined variables*/
	nu(E,E1)
	
BINARY VARIABLES

	/*These variables only used in EUGENECID */
	M(P,J,T)					design variable has value of 1 if protein j expressed from promotor p and Terminator T replaced by parameters in CIS
	
	/*enzyme attribution variables*/
	L(E)						determines if E is encoded in the design
	nu_prime(E,E1)				determine if enzyme E is affected by enzyme E1 in any way shape or form
	beta(E,E1)					determines if E and E1 are encoded in the design
	
	/*these variables are shared by CID and CIS*/
	C_plus(E,L1,L2)				determines if the concentration activation threshold is achieved by comparison to Cjt
	W(E,L1,L2)					determines if the protein is active based on inhibition activation transcription and translation
	alpha_plus(P,L1,L2)			value of 1 if the promotor is active value of zero otherwise
	gamma_plus(E,L1,L2)			value of 1 if the protein is active value of zero otherwise
	omega(E,L1,L2)				value of 1 if enzyme protein or polypeptide E is produced 0 otherwise
	kappa(E,L1,L2)				value of 1 is protein is active (e.g. not inhibited or is activated) 0 otherwise

;

*********************************************** HARDCODED SOLUTIONS FOR DEBUGGING ***********************************************

*hardcode OR solution (used for debugging)
*M.fx(P,J,T) = 0;
*M.fx('P10','gene_GFP','T6') = 1;
*M.fx('P11','gene_GFP','T6') = 1;

*limit solution size
*Ncircuit_max = 2;
*Max_size_allow = 3;

*hardcode AND solution (used for debugging)
*M.fx(P,J,T) = 0;
*M.fx('P10','gene_GFP','T1') = 1;
*M.fx('P11','gene_GFP','T2') = 1;

*limit solution size
*Ncircuit_max = 2;
*Max_size_allow = 3;

*hardcode NOR solution (used for debugging)
*M.fx(P,J,T) = 0;
*M.fx('P10','gene_A','T10') = 1;
*M.fx('P11','gene_A','T10') = 1;
*M.fx('P2','gene_GFP','T10') = 1;

*limit solution size
*Ncircuit_max = 3;
*Max_size_allow = 4;

*hardcode XNOR solution (used for debugging)
*M.fx(P,J,T) = 0;
*M.fx('P10','gene_D','T10') = 1;
*M.fx('P11','gene_D','T10') = 1;
*M.fx('P10','gene_E','T1') = 1;
*M.fx('P11','gene_E','T2') = 1;
*M.fx('P19','gene_GFP','T10') = 1;

*limit solution size
*Ncircuit_max = 5;
*Max_size_allow = 6;

*hardcode NAND solution (used for debugging)
*M.fx(P,J,T) = 0;
*M.fx('P10','gene_D','T3') = 1;
*M.fx('P11','gene_D','T3') = 1;
*M.fx('P19','gene_GFP','T10') = 1;

*limit solution size 
*Ncircuit_max = 3;
*Max_size_allow = 4;

*hardcode NAND solution (used for debugging)
*M.fx(P,J,T) = 0;
*M.fx('P16','gene_GFP','T10') = 1;
*M.fx('P31','gene_GFP','T10') = 1;

*limit solution size 
*Ncircuit_max = 2;
*Max_size_allow = 3;

*hardcode a type of wishful thinking problem to make sure we can catch it
*these three show self-attribution though enzymes and promotors
*this type of wishful thinking should show self-attribution through enzymes
*M.fx(P,J,T) = 0;
*M.fx('P10','gene_F','T10') = 1;
*M.fx('P21','gene_D','T10') = 1;
*M.fx('P24','gene_GFP','T10') = 1;

*limit solution size 
*Ncircuit_max = 3;
*Max_size_allow = 4;

*type of wishful thinking where enzymes regulate each others triads
*M.fx(P,J,T) = 0;
*M.fx('P3','gene_A','T10') = 1;
*M.fx('P2','gene_B','T10') = 1;
*M.fx('P4','gene_GFP','T10') = 1;

*limit solution size 
*Ncircuit_max = 3;
*Max_size_allow = 4;

*another type of wishful thinking working through two enzymes 
*this type of wishful thinking works by one enzyme acting directly on the other
*while the latter works on the triad of the first
*M.fx(P,J,T) = 0;
*M.fx('P9','gene_F','T10') = 1;
*M.fx('P21','gene_D','T9') = 1;
*M.fx('P24','gene_GFP','T10') = 1;

*limit solution size 
*Ncircuit_max = 3;
*Max_size_allow = 4;

*add a buffer to the previous 
*M.fx(P,J,T) = 0;
*M.fx('P9','gene_F','T10') = 1;
*M.fx('P21','gene_A','T10') = 1;
*M.fx('P4','gene_D','T9') = 1;
*M.fx('P24','gene_GFP','T10') = 1;

*limit solution size 
*Ncircuit_max = 4;
*Max_size_allow = 5;

*type of wishful thinking where enzymes regulate each others triads
*M.fx(P,J,T) = 0;
*M.fx('P3','gene_A','T10') = 1;
*M.fx('P2','gene_B','T10') = 1;
*M.fx('P4','gene_GFP','T10') = 1;

*limit solution size 
*Ncircuit_max = 3;
*Max_size_allow = 4;

*Combine 2 types of wishful thinking
*M.fx(P,J,T) = 0;
*M.fx('P9','gene_F','T10') = 1;
*M.fx('P21','gene_D','T9') = 1;
*M.fx('P24','gene_GFP','T10') = 1;
*M.fx('P3','gene_A','T10') = 1;
*M.fx('P2','gene_B','T10') = 1;
*M.fx('P4','gene_GFP','T10') = 1;

*limit solution size 
*Ncircuit_max = 6;
*Max_size_allow = 7;

*type of wishful thinking where enzymes regulate their own triads
*M.fx(P,J,T) = 0;
*M.fx('P4','gene_E','T1') = 1;
*M.fx('P20','gene_GFP','T1') = 1;
*M.fx('P23','gene_E','T1') = 1;

*limit solution size 
*Ncircuit_max = 3;
*Max_size_allow = 4;

*type of wishful thinking where enzymes regulate their own triads
*M.fx(P,J,T) = 0;
*M.fx('P9','gene_D','T1') = 1;
*M.fx('P25','gene_GFP','T1') = 1;
*M.fx('P28','gene_G','T1') = 1;

*limit solution size 
*Ncircuit_max = 3;
*Max_size_allow = 4;


******************************************** DEFINE EQUATIONS ****************************************************************

EQUATIONS

	/*design equations for EuGeneCiD*/
	obj								objective function
	act_threshold_0(E,Ld,Ld1)		constraint used to trigger the activation threshold
	act_threshold_1(E,Ld,Ld1)		constraint used to trigger the activation threshold
	prot_sum(E,Ld,Ld1)				find the sum of transcript levels for each transcript C
	trans_sum(J,T,Ld,Ld1)			find the sum of transcript levels for each transcript
	promotor_limit(P)				limits the number of times a promotor can be used in a given circuit
	trans_limit(J)					limits the number of times a protein coding region can be used in a given circuit
	term_limit(T)					limits the number of times a terminator can be used in a given circuit
	design_limit					limits the number of protein promotor RBS combinations in given circuit

	findalpha(P,Ld,Ld1)				find value of alpha 
	findalpha_sign0(P,Ld,Ld1)		find sign of alpha
	findalpha_sign1(P,Ld,Ld1)		find sign of alpha
	
	transcribed0(P,J,T,Ld,Ld1)		Determines if j is transcribed by p in the current conditions
	transcribed1(P,J,T,Ld,Ld1)		Determines if j is transcribed by p in the current conditions
	transcribed2(P,J,T,Ld,Ld1)		Determines if j is transcribed by p in the current conditions
	
	produced0(E,Ld,Ld1)				determines if a polypeptide is produced
	produced1(E,Ld,Ld1)				determines if a polypeptide is produced
	
	findgamma(E,Ld,Ld1)				Find gamma
	findgamma_sign0(E,Ld,Ld1)		find value of gamma_plus
	findgamma_sign1(E,Ld,Ld1)		find value of gamma_plus
	
	prot_could_act0(E,Ld,Ld1)		determintes if protein could be active based on ligand inputs and other active protiens
	prot_could_act1(E,Ld,Ld1)		determintes if protein could be active based on ligand inputs and other active protiens
	prot_could_act2(E,Ld,Ld1)		determintes if protein could be active based on ligand inputs and other active protiens
	
	prot_is_act0(E,Ld,Ld1)			determines if protein is active based on inhibition activation and concentration
	prot_is_act1(E,Ld,Ld1)			determines if protein is active based on inhibition activation and concentration
	prot_is_act2(E,Ld,Ld1)			determines if protein is active based on inhibition activation and concentration
	
	sat_lambda(Ld,Ld1,Ed)			ensures solution satisfies the logic table
	
	intcut(big_set)					integer cut to prevent repeat solutions
	
	/*necessary for getting rid of the "wishful thinking" problem*/
	/*these determine attribution between enzymes*/
	attribution_0(E)
	attribution_1(E)
	attribution_2(E,E1)
	attribution_3(E,E1)
	attribution_4(E,E1)

	attribution_5(E,E1)
	attribution_6(E,E1)
	attribution_7(E,E1)
	attribution_8(E,E1)
	attribution_9(E,E1,E2)
	attribution_10(E,E1,E2)
	attribution_11(E,E1,E2)
	attribution_12(E,E1)
	
	attribution_13(E,E1)
	attribution_14(E,E1,P,J,T)
	attribution_15(E,E1,P,J,T)
	attribution_16(E,E1,P,J,T)
	attribution_17(E,E1)
	attribution_18(E,E1)
	attribution_19(E,E1)
	
	/*these equations predclude self-attribution*/
	attribution_20(E,E1)
	attribution_21(E,E1)
	
	/*force no meaningless peices*/
	attribution_22(E,E1)
	
	/*add explicit constraints to see if this increases speed*/
	
	explicit_0(E)
	explicit_1(E,Ld,Ld1)
	explicit_2(E,Ld,Ld1)
	explicit_3(E,Ld,Ld1)	
	explicit_4(E,Ld,Ld1)	
	explicit_5(E,Ld,Ld1)	
	
	/*simulation equations for EuGeneCiM*/
	obj_M							objective function
	prot_conc_M(E,Ld,Ld1)			protein concentration
	trans_level_M(J,T,Ld,Ld1)		find the transcript level for each transcript
	
	produced0_M(E,Ld,Ld1)			determines if a polypeptide is produced
	produced1_M(E,Ld,Ld1)			determines if a polypeptide is produced
	
	act_threshold_0_M(E,Ld,Ld1)		constraint used to trigger the activation threshold
	act_threshold_1_M(E,Ld,Ld1)		constraint used to trigger the activation threshold
	
	
;

*********************************************** EUGENECID ********************************************************************

*objective function
obj..								Z_D =e= sum(Ed, sum(Ld, sum(Ld1, C(Ed,Ld,Ld1) * lambda(Ld,Ld1,Ed) - C(Ed,Ld,Ld1) * (1 - lambda(Ld,Ld1,Ed)))));

*limits the size of the designed circuit
promotor_limit(P)..					sum(J, sum(T, M(P,J,T))) =l= Np_max;
trans_limit(J)..					sum(P, sum(T, M(P,J,T))) =l= Nj_max;
term_limit(T)..						sum(P, sum(J, M(P,J,T))) =l= Nt_max;
design_limit..						sum(P, sum(J, sum(T, M(P,J,T)))) =l= Ncircuit_max;

*determines total effect of present activator and inhibitors on promotors
findalpha(P,Ld,Ld1)..				alpha(P,Ld,Ld1) =e= Z(P) + sum(E, (W(E,Ld,Ld1) * I(P,E) * H(P,E))) + I(P,Ld) * H(P,Ld) + I(P,Ld1) * H(P,Ld1) - (I(P,Ld1) * H(P,Ld1) * sigma(Ld,Ld1));
findalpha_sign0(P,Ld,Ld1)..			alpha(P,Ld,Ld1) =g= epsilon * alpha_plus(P,Ld,Ld1) - V * (1 - alpha_plus(P,Ld,Ld1));
findalpha_sign1(P,Ld,Ld1)..			alpha(P,Ld,Ld1) =l= V * alpha_plus(P,Ld,Ld1);

*determine how much transcript J is deliberately transcribed from terminator
transcribed0(P,J,T,Ld,Ld1)..		xi(P,J,T,Ld,Ld1) =l= S(P) * M(P,J,T);
transcribed1(P,J,T,Ld,Ld1)..		xi(P,J,T,Ld,Ld1) =l= S(P) * alpha_plus(P,Ld,Ld1);
transcribed2(P,J,T,Ld,Ld1)..		xi(P,J,T,Ld,Ld1) =g= S(P) * (M(P,J,T) + alpha_plus(P,Ld,Ld1) - 1);

*make sure that enzyme is noted as being produced if it does have a non-zero concentration
produced0(E,Ld,Ld1)..				omega(E,Ld,Ld1) =l= V * C(E,Ld,Ld1);
produced1(E,Ld,Ld1)..				omega(E,Ld,Ld1) =g= epsilon * C(E,Ld,Ld1);

*determine the transcript level
trans_sum(J,T,Ld,Ld1)..				phi(J,T,Ld,Ld1) =e= sum(P, (xi(P,J,T,Ld,Ld1) + M(P,J,T) * F(P)) * ((0.5) ** (1 / (G(T) + epsilon))));

*determine protein concentration as a function of transcription and translation efficiency
*adds that produced from the previous time points RNA with that already present from previous time point
prot_sum(E,Ld,Ld1)..				C(E,Ld,Ld1) =e= sum(J, (rho(J,E) * eta(J) * sum(T, phi(J,T,Ld,Ld1))) * ((0.5) ** (1 / (R(E) + epsilon))));

*determine if concentration is at the level necessary for activity
act_threshold_0(E,Ld,Ld1)..			(theta(E) + epsilon) * C_plus(E,Ld,Ld1) =l= C(E,Ld,Ld1);
act_threshold_1(E,Ld,Ld1)..			C(E,Ld,Ld1) =l= (V - (theta(E) - epsilon)) * C_plus(E,Ld,Ld1) + (theta(E) - epsilon);

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

*need to constrain such that the logic table is satisfied
sat_lambda(Ld,Ld1,Ed)..				W(Ed,Ld,Ld1) =e= lambda(Ld,Ld1,Ed);

*this is the integer cut just on promotors and transcripts since terminators more of a fine tune control not worried about repeat pairs with them
intcut(soln_set)..					sum(P,sum(J,sum(T,pastMs(P,J,T,soln_set)) * sum(T,M(P,J,T)))) =l= sum(P,sum(J,sum(T,pastMs(P,J,T,soln_set)))) - 1;

*these constraints exist to prevent "wishful thinking" circuits
*here we will look at the attribution of enzymes onto other enzymes as this is the overriding problem

*determine if enzyme is encoded
attribution_0(E)..					L(E) =g= epsilon * sum(P, sum(J, sum(T, M(P,J,T) * rho(J,E))));
attribution_1(E)..					L(E) =l= V * sum(P, sum(J, sum(T, M(P,J,T) * rho(J,E))));

*determine if enzyme pairs are encoded
attribution_2(E,E1)..				beta(E,E1) =l= L(E);
attribution_3(E,E1)..				beta(E,E1) =l= L(E1);
attribution_4(E,E1)..				beta(E,E1) =g= L(E) + L(E1) - 1;

*effect of other enzymes directly upon e
attribution_5(E,E1)..				D(E,E1) =e= abs(B(E,E1)) * beta(E,E1);

*effect of other enzymes e1 upon traid of enzyme e
*need to use absolute values
attribution_6(E,E1)..				K(E,E1) =l= V * beta(E,E1);
attribution_7(E,E1)..				K(E,E1) =l= sum(J, sum(P, abs(I(P,E1)) * sum(T, M(P,J,T) * rho(J,E)))) + V * (1 - beta(E,E1));
attribution_8(E,E1)..				K(E,E1) =g= sum(J, sum(P, abs(I(P,E1)) * sum(T, M(P,J,T) * rho(J,E)))) - V * (1 - beta(E,E1));

*effect of other enzymes e2 upon the activity of other enzyme e through the medium of a second enzyme e1
attribution_9(E,E1,E2)..			U(E,E1,E2) =l= V * beta(E,E1);
attribution_10(E,E1,E2)..			U(E,E1,E2) =l= (abs(B(E,E1)) * nu_prime(E1,E2)) + V * (1 - beta(E,E1));
attribution_11(E,E1,E2)..			U(E,E1,E2) =g= (abs(B(E,E1)) * nu_prime(E1,E2)) - V * (1 - beta(E,E1));

*effect of other enzymes e1 upon the activity of other enzymes through the medium of a second enzyme e2
attribution_12(E,E1)..				U_prime(E,E1) =e= sum(E2, U(E,E2,E1) * (1 - sigma(E,E1) * sigma(E1,E2)));

*effect of other enzymes e1 upon the triad of e through the medium of a secon d enzyme e2
attribution_13(E,E1)..				X(E,E1) =e= sum(P, sum(J, sum(T, chi(E,E1,P,J,T))));

attribution_14(E,E1,P,J,T)..		chi(E,E1,P,J,T) =l= V * M(P,J,T);
attribution_15(E,E1,P,J,T)..		chi(E,E1,P,J,T) =l= sum(E2, abs(I(P,E2)) * nu_prime(E2,E1) * rho(J,E)) + V * (1 - M(P,J,T));
attribution_16(E,E1,P,J,T)..		chi(E,E1,P,J,T) =g= sum(E2, abs(I(P,E2)) * nu_prime(E2,E1) * rho(J,E)) - V * (1 - M(P,J,T));

*sum the effects to get attribution of e to e1 attribution 
attribution_17(E,E1)..				nu(E,E1) =e= D(E,E1) + K(E,E1) + U_prime(E,E1) + X(E,E1);
attribution_18(E,E1)..				nu(E,E1) =g= nu_prime(E,E1);
attribution_19(E,E1)..				nu(E,E1) =l= V * nu_prime(E,E1);

*prevent self-attribution and self-feedback cases
attribution_20(E,E1)..				nu_prime(E,E1) =l= 1 - sigma(E,E1);
attribution_21(E,E1)..				nu_prime(E,E1) =g= sigma(E,E1) - 1;

*prevent addition of meaningless pieces
attribution_22(E,E1)..				L(E) =l= sum(Ed, nu_prime(Ed,E)) + Ed_val(E);

explicit_0(E)..						L(E) =g= Ed_val(E);
explicit_1(E,Ld,Ld1)..				W(E,Ld,Ld1) =l= L(E);
explicit_2(E,Ld,Ld1)..				kappa(E,Ld,Ld1) =l= L(E);
explicit_3(E,Ld,Ld1)..				omega(E,Ld,Ld1) =l= L(E);
explicit_4(E,Ld,Ld1)..				C(E,Ld,Ld1) =l= V * L(E);
explicit_5(E,Ld,Ld1)..				C_plus(E,Ld,Ld1) =l= L(E);

*********************************************** EUGENECIM ********************************************************************

*objective function
obj_M..								Z_M =e= sum(Ed, sum(Ld, sum(Ld1, C(Ed,Ld,Ld1))));

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

*********************************************** DEFINE MODELS ****************************************************************

MODEL EUGENECID
/

	obj
	
	promotor_limit
	trans_limit
	term_limit
	design_limit
	
	findalpha
	findalpha_sign0
	findalpha_sign1
	
	transcribed0
	transcribed1
	transcribed2
	
	produced0
	produced1
	
	prot_sum
	
	trans_sum
	
	act_threshold_0
	act_threshold_1
	
	findgamma
	findgamma_sign0
	findgamma_sign1
	
	prot_could_act0
	prot_could_act1
	prot_could_act2
	
	prot_is_act0
	prot_is_act1
	prot_is_act2
	
	sat_lambda
	
	intcut
	
	/*attribution equations to prevent wishful thinking hopefully*/
	/*works by identify attribution between enzymes*/
	attribution_0
	attribution_1
	attribution_2
	attribution_3
	attribution_4
	attribution_5
	attribution_6
	attribution_7
	attribution_8
	attribution_9
	attribution_10
	attribution_11
	attribution_12
	attribution_13
	attribution_14
	attribution_15
	attribution_16
	attribution_17
	attribution_18
	attribution_19
	
	/*force no self attribution*/
	attribution_20
	attribution_21
	
	/*force no meaningless peices*/
	attribution_22

	/*seeing if can increase speed through explicitly stating*/
	/*that which should be implicitly true*/
	explicit_0
	explicit_1
	explicit_2
	explicit_3
	explicit_4
	explicit_5 

/
;

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

EUGENECID.holdfixed = 1;
EUGENECID.optfile = 1;

EUGENECIM.holdfixed = 1;
EUGENECIM.optfile = 1;

********************************************** INITIALIZE OUTPUT FILES *******************************************************

*create a simple results output file for EuGeneCiD
FILE SIMPLECID / a_simple_EUGENECID_report_Cu_nimply_Zn.txt /;
PUT SIMPLECID;
SIMPLECID.pc = 2;
SIMPLECID.lw = 25;

PUT "SIMPLIFIED EUGENECID RESULTS AND STATE FILE"/;

PUTCLOSE;

*create a detailed results output file for EuGeneCiD
FILE DETAILCID / a_detailed_EUGENECID_report_Cu_nimply_Zn.txt /;
PUT DETAILCID;
DETAILCID.pc = 2;
DETAILCID.lw = 25;

PUT "DETAILED EUGENECID RESULTS AND STATE FILE"//;

PUTCLOSE;

/*this file will store just the genetic circuit designs*/
FILE DESIGNS / circuit_designs_file_Cu_nimply_Zn.txt /;
PUT DESIGNS;
DESIGNS.pc = 2;
DESIGNS.lw = 30;

PUT "GENETIC CIRCUIT DESIGNS FOR EACH SOLUTION"//;

PUTCLOSE;

/*this file will be essentially a table storing state values for various timepoints*/
/*each subsequent time will be added as a new row*/
FILE STATES / circuit_state_Cu_nimply_Zn.csv /;
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
FILE TABLES / circuit_logic_tables_Cu_nimply_Zn.txt /;
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

/*for each solution where we aren't done finding solutions yet*/
LOOP(big_set$(not done),

	/*start by solving for a unique design*/
	SOLVE EUGENECID USING MIP MAXIMIZING Z_D;
	
	/*loop through incrementing the size of the allowed circuit in case the limitations are too strict*/
	WHILE(((EUGENECID.modelStat ne 1) AND (EUGENECID.modelStat ne 2) AND (EUGENECID.modelStat ne 8) AND (EUGENECID.modelStat ne 15) AND (EUGENECID.modelStat ne 16) AND (EUGENECID.modelStat ne 17) AND (done_2 ne 1)),
	
		DESIGNS.ap = 1;
		PUT DESIGNS;
		DESIGNS.pc = 2;
		DESIGNS.lw = 30;
		
		PUT ///;
		PUT "No remaining solutions found for circuit size: ",Ncircuit_max/;
		PUT "Model Status: ", EUGENECID.modelStat/;
		PUT "Solver Status: ", EUGENECID.solveStat/;
		PUT "Time to Determine: ",system.tab,EUGENECID.etSolve:0:8," second";
		
		PUT ///;
		
		PUTCLOSE;
	
		/*increment maximum design size*/
		Ncircuit_max = Ncircuit_max + 1;
	
		/*solve with bigger size*/
		SOLVE EUGENECID USING MIP MAXIMIZING Z_D;
		
		/*if at maximum allowed design size call it done regardless*/
		IF((Ncircuit_max ge Max_size_allow),
		
			done_2 = 1;
		
		);
	
	);
	
	/*save the necessary information for the design if a solution was found*/
	/*check if a solution has been found*/
	/*case 1 solution is found and this is the first time point*/
	/*need to write then to both DESIGNS and STATES fix this solution for subsequent timepoints*/
	IF(((EUGENECID.modelStat eq 1) OR (EUGENECID.modelStat eq 2) OR (EUGENECID.modelStat eq 8) OR (EUGENECID.modelStat eq 15) OR (EUGENECID.modelStat eq 16) OR (EUGENECID.modelStat eq 17)),
		
		/*if here then a solution was found*/
		
		/*write to simple results file*/
		SIMPLECID.ap = 1;
		PUT SIMPLECID;
		SIMPLECID.pc = 2;
		SIMPLECID.lw = 30;
		
		PUT //;
		PUT "===================================================================="/;
		PUT "RESULTS FOR SOLUTION NUMBER ",big_set.tl/;
		PUT "===================================================================="//;
		
		PUT "CIRCUIT DESIGN"//;
		PUT "PROMOTOR",system.tab,system.tab,system.tab,system.tab,system.tab,"TRANSCRIPT",system.tab,system.tab,system.tab,system.tab,system.tab,"TERMINATOR"/;
		
		LOOP(P,
		
			LOOP(J,
			
				LOOP(T,
				
					IF((M.l(P,J,T) eq 1),
					
						PUT P.tl, J.tl, T.tl/; 
					
					);
				
				);
				
			);
			
		);
		
		PUT //;
		
		PUT "RESPONDANTS TO LIGANDS"/;
		
		LOOP(P,
		
			LOOP(J,
			
				LOOP(T,
				
					IF((M.l(P,J,T) eq 1),
					
						LOOP(Ld,
					
							IF((I(P,Ld) eq -1),
							
								PUT P.tl," inhibited by ",Ld.tl/;
							
							ELSEIF (I(P,Ld) eq 1),
							
								PUT P.tl," activated by ",Ld.tl/;
							
							);
							
							LOOP(E,
							
								iF((rho(J,E) * B(E,Ld) < 0),
								
									PUT E.tl," inhibited by ",Ld.tl/;
								
								ELSEIF (rho(J,E) * B(E,Ld) > 0),
								
									PUT E.tl," activated by ",Ld.tl/;
								
								);
							
							);
							
						);
					
					);
				
				);
				
			);
			
		);
		
		PUT //;
		
		PUT "ENZYMES ATTRIBUTION TO OTHER ENZYMES"/;
		
		LOOP(E,
		
			LOOP(E1,
			
				IF((nu_prime.l(E,E1) eq 1),
				
					PUT E.tl," is influenced by ",E1.tl/;
				
				);
			
			);
		
		);
		
		PUT //;
		
		/*write the variable values at each ligand combination*/
		LOOP(Ld,

			LOOP(Ld1,
			
				PUT "ligand condition: Ld: ",Ld.tl," Ld1: ",Ld1.tl/;
				PUT "------------------------------------------------------------------------"//;
				
				/*list transcripts produced and their level*/
				PUT "Transcripts Produced"/;
				PUT "transcript",system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,"deliberate",system.tab,system.tab,"leaked",system.tab,system.tab,system.tab,"produced",system.tab,system.tab,"degraded",system.tab,system.tab,"Phi"/;
				
				LOOP(j,
				
					/*check if the transcript is used*/
					IF((sum(P, sum(T, M.l(P,J,T))) ge 1),
					
						/*if used report on transcript levels from various sources*/
						PUT j.tl,(sum(P, sum(T, xi.l(P,J,T,Ld,Ld1)))),system.tab,sum(P, sum(T, M.l(P,j,T) * F(P))),system.tab,(sum(T, sum(P, xi.l(P,J,T,Ld,Ld1) + M.l(P,J,T) * F(P)))),system.tab,((sum(T, sum(P, xi.l(P,J,T,Ld,Ld1) + M.l(P,J,T) * F(P)))) - sum(T,sum(P, (xi.l(P,J,T,Ld,Ld1) + M.l(P,J,T) * F(P)) * ((0.5) ** (1 / (G(T) + epsilon)))))),system.tab,(sum(T,phi.l(J,T,Ld,Ld1)))/;
					
					);
					
				);
				
				PUT //;
				
				/*finally output lists of produced and active proteins*/
				PUT "Proteins Produced:"/;
				
				LOOP(E,
				
					IF((C.l(E,Ld,Ld1) > 0),
					
						PUT E.tl,C.l(E,Ld,Ld1)/;
					
					);
				
				);
				
				PUT //;
				
				PUT "Proteins Expressed:"/;
				
				LOOP(E,
				
					IF((C_plus.l(E,Ld,Ld1) eq 1),
					
						PUT E.tl,C_plus.l(E,Ld,Ld1)/;
					
					);
				
				);
				
				PUT //;
				
				/*finally output lists of produced and active proteins*/
				PUT "Active Proteins:"/;
				
				LOOP(E,
				
					IF((W.l(E,Ld,Ld1) eq 1),
					
						PUT E.tl,W.l(E,Ld,Ld1)/;
					
					);
				
				);
				
				PUT ////;
				
			);
			
		);
		
		PUTCLOSE;
		
		/*write to detailed results file*/
		DETAILCID.ap = 1;
		PUT DETAILCID;
		DETAILCID.pc = 2;
		DETAILCID.lw = 30;
		
		PUT //;
		PUT "===================================================================="/;
		PUT "RESULTS FOR SOLUTION NUMBER ",big_set.tl/;
		PUT "===================================================================="//;
		
		PUT "Designed Genetic Circuit"/;
		PUT "Promotor",system.tab,"Transcript",system.tab,system.tab,system.tab,system.tab,system.tab,"Terminator"/;

		LOOP(P,

			LOOP(J,
			
				LOOP(T,
			
					IF((M.l(P,J,T) eq 1),
				
						PUT P.tl;
					
						DETAILCID.lw = 30;
					
						PUT J.tl;
						PUT T.tl/;
					
						DETAILCID.lw = 12;
						
					);
				);
				
			);
			
		);

		PUT //;
		
		/*since attribution is independent of condition report on the attribution outside the condition loops*/
		
		/*report enzyme pairs encoded which limit enzyme attribution*/ 
		
		PUT "Encoding of E (row): L"/;
		
		LOOP(E,
			
			PUT E.tl,L.l(E)/;
			
		);
		
		PUT //;
		
		PUT "Paired encoding of both E (row) to E1 (col): Beta"/;
		
		PUT system.tab,system.tab,system.tab,system.tab,system.tab;
		
		LOOP(E1,
			
			PUT E1.tl;
			
		);
		
		PUT /;
		
		LOOP(E,
			
			PUT E.tl;
		
			LOOP(E1,
			
				PUT beta.l(E,E1);
			
			);
			
			PUT /;
			
		);
		
		PUT //;
		
		PUT "Partial attribution of E (row) to E1 (col): D"/;
		
		PUT system.tab,system.tab,system.tab,system.tab,system.tab;
		
		LOOP(E1,
			
			PUT E1.tl;
			
		);
		
		PUT /;
		
		LOOP(E,
			
			PUT E.tl;
		
			LOOP(E1,
			
				PUT D.l(E,E1);
			
			);
			
			PUT /;
			
		);
		
		PUT //;
		
		PUT "Partial attribution of E (row) to E1 (col): K"/;
		
		PUT system.tab,system.tab,system.tab,system.tab,system.tab;
		
		LOOP(E1,
			
			PUT E1.tl;
			
		);
		
		PUT /;
		
		LOOP(E,
			
			PUT E.tl;
		
			LOOP(E1,
			
				PUT K.l(E,E1);
			
			);
			
			PUT /;
			
		);
		
		PUT //;
		
		LOOP(E2,
		
			PUT "Partial attribution of E (row) to E1 (col): U, E2: ", E2.tl/;
			
			PUT system.tab,system.tab,system.tab,system.tab,system.tab;
			
			LOOP(E1,
				
				PUT E1.tl;
				
			);
			
			PUT /;
		
			LOOP(E,
				
				PUT E.tl;
			
				LOOP(E1,
				
					PUT U.l(E,E1,E2);
				
				);
				
				PUT /;
				
			);
			
			PUT //;
			
		);
		
		PUT "Partial attribution of E (row) to E1 (col): u_prime"/;
		
		PUT system.tab,system.tab,system.tab,system.tab,system.tab;
		
		LOOP(E1,
			
			PUT E1.tl;
			
		);
		
		PUT /;
		
		LOOP(E,
			
			PUT E.tl;
		
			LOOP(E1,
			
				PUT u_prime.l(E,E1);
			
			);
			
			PUT /;
			
		);
		
		PUT //;
		
		PUT "Partial attribution of E (row) to E1 (col): X"/;
		
		PUT system.tab,system.tab,system.tab,system.tab,system.tab;
		
		LOOP(E1,
			
			PUT E1.tl;
			
		);
		
		PUT /;
		
		LOOP(E,
			
			PUT E.tl;
		
			LOOP(E1,
			
				PUT X.l(E,E1);
			
			);
			
			PUT /;
			
		);
		
		PUT //;
		
		/*report attribution to enzymes*/ 
		PUT "Attribution of E (row) to E1 (col): nu"/;
		
		PUT system.tab,system.tab,system.tab,system.tab,system.tab;
		
		LOOP(E1,
			
			PUT E1.tl;
			
		);
			
		PUT /;
		
		LOOP(E,
			
			PUT E.tl;
		
			LOOP(E1,
			
				PUT nu.l(E,E1);
			
			);
			
			PUT /;
			
		);
		
		PUT //;
		
		/*report attribution to enzymes*/ 
		PUT "Attribution of E (row) to E1 (col): nu_prime"/;
		
		PUT system.tab,system.tab,system.tab,system.tab,system.tab;
		
		LOOP(E1,
			
			PUT E1.tl;
			
		);
		
		PUT /;
		
		LOOP(E,
			
			PUT E.tl;
		
			LOOP(E1,
			
				PUT nu_prime.l(E,E1);
			
			);
			
			PUT /;
			
		);
		
		PUT //;

		/*write the variable values at each ligand combination*/
		LOOP(Ld,

			LOOP(Ld1,
			
				PUT "ligand condition: Ld: ",Ld.tl," Ld1: ",Ld1.tl/;
				PUT "--------------------------------------------------------------------------"/;
				
				PUT "promotor variables:"/;
				PUT "Promotor",system.tab,system.tab,system.tab,"Z(P)",system.tab,system.tab,"sum W*I",system.tab,system.tab,"I(P,L1)",system.tab,system.tab,"I(P,L2)",system.tab,system.tab,"alpha",system.tab,system.tab,"alpha+"/;
				
				LOOP(P,
				
					PUT P.tl,Z(P),(sum(E, (W.l(E,Ld,Ld1) * I(P,E)))),I(P,Ld),I(P,Ld1),alpha.l(P,Ld,Ld1),alpha_plus.l(P,Ld,Ld1)/;

				);
				
				PUT //;
				
				PUT "transcription of transcript J from terminator T:"/;
				PUT "J",system.tab,system.tab,system.tab,system.tab,"T",system.tab,system.tab,system.tab,system.tab,system.tab,"xi",system.tab,system.tab,system.tab,"leaked",system.tab,system.tab,"degraded",system.tab,"total"/;
				
				LOOP(J,
						
					LOOP(T,
						
						PUT J.tl,system.tab,T.tl,sum(P,xi.l(P,J,T,Ld,Ld1)),(sum(P, M.l(P,J,T) * F(P))),(sum(P, xi.l(P,J,T,Ld,Ld1) + M.l(P,J,T) * F(P)) * (1 - (0.5) ** (1 / (G(T) + epsilon)))),phi.l(J,T,Ld,Ld1)/;
					
					);
						
				);	
				
				PUT //;
				
				PUT "transcript variable (phi):"/;
				PUT "transcript",system.tab,system.tab,system.tab;
				
				LOOP(T,
				
					PUT T.tl;
				
				);
				
				PUT /;
				
				LOOP(J,
				
					PUT J.tl;
				
					LOOP(T,
					
						PUT phi.l(J,T,Ld,Ld1);
					
					);
					
					PUT /;
				
				);
				
				PUT //;
				
				PUT "enzyme, protein, and polypeptide variables:"/;
				PUT "enzyme",system.tab,system.tab,system.tab,system.tab,"Zeta",system.tab,system.tab,"sum(W*B)",system.tab,"B(L1)",system.tab,system.tab,"B(L2)",system.tab,system.tab,"gamma",system.tab,system.tab,"gamma+",system.tab,system.tab,"omega",system.tab,system.tab,"kappa",system.tab,system.tab,"C",system.tab,system.tab,system.tab,"C+",system.tab,system.tab,system.tab,"W"/;
				
				LOOP(E,
				
					temp = sum(E1, W.l(E1,Ld,Ld1) * B(E,E1));
					PUT E.tl,Zeta(E),temp,B(E,Ld),B(E,Ld1),gamma.l(E,Ld,Ld1),gamma_plus.l(E,Ld,Ld1),omega.l(E,Ld,Ld1),kappa.l(E,Ld,Ld1),C.l(E,Ld,Ld1),C_plus.l(E,Ld,Ld1),W.l(E,Ld,Ld1)/;
				
				);
				
				PUT //;
				
				/*finally output lists of produced and active proteins*/
				PUT "Proteins Produced:"/;
				
				LOOP(E,
				
					IF((C.l(E,Ld,Ld1) ge 1),
					
						PUT E.tl,C.l(E,Ld,Ld1)/;
					
					);
				
				);
				
				PUT //;
				
				/*finally output lists of produced and active proteins*/
				PUT "Proteins Expressed:"/;
				
				LOOP(E,
				
					IF((C_plus.l(E,Ld,Ld1) eq 1),
					
						PUT E.tl,C_plus.l(E,Ld,Ld1)/;
					
					);
				
				);
				
				PUT //;
				
				/*finally output lists of produced and active proteins*/
				PUT "Active Proteins:"/;
				
				LOOP(E,
				
					IF((W.l(E,Ld,Ld1) eq 1),
					
						PUT E.tl,W.l(E,Ld,Ld1)/;
					
					);
				
				);
				
				PUT ///;
				
			);
			
		);
		
		PUTCLOSE;
		
		/*write the results to an output file*/
		/*note that this first time point will report the design and concentration results*/
		DESIGNS.ap = 1;
		PUT DESIGNS;
		DESIGNS.pc = 2;
		DESIGNS.lw = 30;
		
		PUT //;
		
		PUT "Solution Number: ",big_set.tl/;
		PUT "Objective Value: ",Z_D.l/;
		PUT "Circuit Size: ",Ncircuit_max/;
		PUT "model status: ",system.tab,EUGENECID.modelStat/;
		PUT "solver status: ",system.tab,EUGENECID.solveStat/;
		PUT "time to solution: ",system.tab,EUGENECID.etSolve:0:8," second"//;
	
		PUT "Designed Genetic Circuit:"/;
		PUT "Promoter",system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,"Transcript",system.tab,system.tab,system.tab,system.tab,system.tab,system.tab,"Terminator"/;    
		PUT "--------------------------------------------------------------------------"/;
		
		LOOP(P,

			LOOP(J,
			
				LOOP(T,
				
					IF((M.l(P,J,T) eq 1),
				
						PUT P.tl,system.tab,J.tl,system.tab,T.tl/;
					
					);
				
				);
				
			);
			
		);
		
		PUTCLOSE;
		
		/*write the design logic table results*/
		TABLES.ap = 1;
		PUT TABLES;
		TABLES.pc = 2;
		
		PUT "LOGIC TABLES SOLUTION NUMBER ",big_set.tl/;
		PUT "======================================================================="//;
		
		PUT "DESIGN CONCENTRATION VALUES"/;
		PUT "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"//;
		
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
				
					PUT C.l(Ed,Ld,Ld1);
					
				);
				
				/*finish with a newline*/
				PUT /;
			
			);
			
		);

		PUT //;
		
		PUTCLOSE;
		
		/*since a design was found save the design for use in the simulator*/
		M_design(P,J,T) = M.l(P,J,T);
		
		/*increment the number of solutions*/
		num_solns = num_solns + 1;
	
	ELSE
	
		/*if here we are done looking for new solutions not sure how to handle these types of errors*/
		done = 1;
		
		/*report that designing is finished*/
		DESIGNS.ap = 1;
		PUT DESIGNS;
		DESIGNS.pc = 2;
		DESIGNS.lw = 30;
		
		PUT //;
		
		PUT "Solution Number ", num_solns," Attempted, was unsuccessful"/;
		PUT "model status: ",system.tab,EUGENECID.modelStat/;
		PUT "solver status: ",system.tab,EUGENECID.solveStat;
		
		PUTCLOSE;
		
	);
	
	done_3 = 0;

	/*Inside the loop will be solving the simulation problem EUGENECIM if EUGENECID produced and answer*/
	LOOP(tau$(not done_3),

		/*solve EuGeneCiM for each time point*/
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
				
					PUT big_set.tl,tau.tl,Ld.tl,Ld1.tl;
				
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
			
			DESIGNS.ap = 1;
			PUT DESIGNS;
			DESIGNS.pc = 2;
			DESIGNS.lw = 4;
			
			PUT //;
			PUT "Solution Number ", big_set.tl," at Time Point ", tau.tl," modeling attempted, was unsuccessful (Fatal Error)"/;
			PUT "model status: ",system.tab,EUGENECIM.modelStat/; 
			PUT "solver status: ",system.tab,EUGENECIM.solveStat;
			
			PUTCLOSE;
		
		);
	
	);
	
	/*by this point done simulating the genetic circuit behavior*/
	
	/*save the solution that has been found and save the design to apply the integer cuts*/
	soln_set(big_set) = yes;
	pastMs(P,J,T,soln_set(big_set)) = M.l(P,J,T);
	
	/*reset the transfer of enzyme concentration for default initial conditions again*/
	ECTrans(E,Ld,Ld1) = 0;
	
	/*reset the transfer of RNA transcription for default initial conditions again*/
	RNATrans(J,T,Ld,Ld1) = 0;
	
	/*check to see if we have exceeded our arbitrary solution number limit*/
	IF((num_solns >= max_solns),
		
		done = 1;
		
	);	
	
);