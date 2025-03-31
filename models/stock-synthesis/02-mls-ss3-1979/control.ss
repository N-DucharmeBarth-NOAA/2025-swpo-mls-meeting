#C No comments
#C file created using an r4ss function
#C file write time: 2025-01-09  22:25:06
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
4 # recr_dist_method for parameters
1 # not yet implemented; Future usage:Spawner-Recruitment; 1=global; 2=by area
1 # number of recruitment settlement assignments 
0 # unused option
# for each settlement assignment:
#_GPattern	month	area	age
1	1	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
1 #_Nblock_Patterns
1 #_blocks_per_pattern
#_begin and end years of blocks
1978 1978
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
4 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=Maunder_M;_6=Age-range_Lorenzen
#_ #_Age_natmort_by sex x growthpattern
#_Age_0	Age_1	Age_2	Age_3	Age_4	Age_5	Age_6	Age_7	Age_8	Age_9	Age_10
0.707867	0.4718	0.4007	0.3766	0.3669	0.3629	0.3612	0.3605	0.3602	0.36	0.36	#_natM1
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.25 #_Age(post-settlement)_for_L1;linear growth below this
10 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
2 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
1 #_First_Mature_Age
3 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
    5	     100	    88.3822	    88.3822	 99	0	 -3	0	0	0	0	  0	0	0	#_L_at_Amin_Fem_GP_1 
   50	     600	    210.956	    210.956	 99	0	 -4	0	0	0	0	  0	0	0	#_L_at_Amax_Fem_GP_1 
 0.01	    0.99	     0.8351	     0.8351	 99	0	 -5	0	0	0	0	  0	0	0	#_VonBert_K_Fem_GP_1 
 0.01	      20	     9.9314	     9.9314	 99	0	 -2	0	0	0	0	  0	0	0	#_CV_young_Fem_GP_1  
 0.01	      20	    11.0517	    11.0517	 99	0	 -3	0	0	0	0	  0	0	0	#_CV_old_Fem_GP_1    
   -3	       3	5.39942e-07	5.39942e-07	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_1_Fem_GP_1   
   -3	       5	    3.58378	    3.58378	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_2_Fem_GP_1   
    1	     300	    181.205	    181.205	 99	0	 -3	0	0	0	0	  0	0	0	#_Mat50%_Fem_GP_1    
 -200	       3	  -0.197458	  -0.197458	 99	0	 -3	0	0	0	0	  0	0	0	#_Mat_slope_Fem_GP_1 
   -3	      20	          1	          1	0.5	6	 -3	0	0	0	0	0.5	0	0	#_Eggs_alpha_Fem_GP_1
   -3	       3	          1	          1	0.5	6	 -3	0	0	0	0	0.5	0	0	#_Eggs_beta_Fem_GP_1 
  0.1	      10	          1	          1	  1	6	 -1	0	0	0	0	  0	0	0	#_CohortGrowDev      
1e-06	0.999999	        0.5	        0.5	0.5	0	-99	0	0	0	0	  0	0	0	#_FracFemale_GP_1    
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
1 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
   2	  15	  6	  6	  99	0	 1	0	0	0	0	0	0	0	#_SR_LN(R0)  
 0.2	0.99	0.8	0.8	1000	6	-2	0	0	0	0	0	0	0	#_SR_BH_steep
0.05	 1.9	0.2	0.2	1000	6	-4	0	0	0	0	0	0	0	#_SR_sigmaR  
  -4	   4	  0	  0	  99	0	-1	0	0	0	0	0	0	0	#_SR_regime  
   0	   0	  0	  0	  99	0	-1	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1979 # first year of main recr_devs; early devs can preceed this era
2021 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase
1 # (0/1) to read 13 advanced options
1969 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
2 #_recdev_early_phase
0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1989 #_last_yr_nobias_adj_in_MPD; begin of ramp
2006 #_first_yr_fullbias_adj_in_MPD; begin of plateau
2021 #_last_yr_fullbias_adj_in_MPD
2021 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.2 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
0 #_period of cycles in recruitment (N parms read below)
-5 #min rec_dev
5 #max rec_dev
0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
#Fishing Mortality info
0.2 # F ballpark
-2010 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
5 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms
#_LO	HI	INIT	PRIOR	PRIOR_SD	PR_type	PHASE
0	5	0.153	0.2	99	0	-1	#_1
0	5	0.153	0.2	99	0	-1	#_2
0	5	0.153	0.2	99	0	-1	#_3
0	5	0.153	0.2	99	0	-1	#_4
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
   15	1	0	0	0	1	#_S01_INDEX.1-4
-9999	0	0	0	0	0	#_terminator   
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-25	25	-7.68236	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S01_INDEX.1-4
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	0	0	#_1 F01_LL.JP.1   
24	0	0	0	#_2 F02_LL.JP.2   
24	0	0	0	#_3 F03_LL.JP.3   
24	0	0	0	#_4 F04_LL.JP.4   
 5	0	0	4	#_5 F05_LL.TW.4   
24	0	0	0	#_6 F06_LL.AU.2   
24	0	0	0	#_7 F07_LL.AU.3   
 5	0	0	7	#_8 F08_LL.NZ.3   
24	0	0	0	#_9 F09.REC.AU.3  
 1	0	0	0	#_10 F10_REC.NZ.3 
 5	0	0	1	#_11 F11_LL.ALL.1 
24	0	0	0	#_12 F12_LL.ALL.2 
 5	0	0	3	#_13 F13_LL.ALL.3 
 5	0	0	4	#_14 F14_LL.ALL.4 
24	0	0	0	#_15 S01_INDEX.1-4
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
11	0	0	0	#_1 F01_LL.JP.1   
11	0	0	0	#_2 F02_LL.JP.2   
11	0	0	0	#_3 F03_LL.JP.3   
11	0	0	0	#_4 F04_LL.JP.4   
11	0	0	0	#_5 F05_LL.TW.4   
11	0	0	0	#_6 F06_LL.AU.2   
11	0	0	0	#_7 F07_LL.AU.3   
11	0	0	0	#_8 F08_LL.NZ.3   
11	0	0	0	#_9 F09.REC.AU.3  
11	0	0	0	#_10 F10_REC.NZ.3 
11	0	0	0	#_11 F11_LL.ALL.1 
11	0	0	0	#_12 F12_LL.ALL.2 
11	0	0	0	#_13 F13_LL.ALL.3 
11	0	0	0	#_14 F14_LL.ALL.4 
11	0	0	0	#_15 S01_INDEX.1-4
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_1 
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_2 
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_3 
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_4 
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_5 
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_6 
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_7 
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_8 
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_9 
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_10
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_11
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_12
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_13
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_14
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_15
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_16
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_17
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_18
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_19
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_20
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_21
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_22
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_23
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_24
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_25
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_26
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_27
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_28
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_29
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_30
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_31
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_32
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_33
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_34
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_35
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_36
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_37
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_38
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_39
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_40
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_41
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_42
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_43
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_44
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_45
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_46
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_47
   0	500	  10	  10	99	0	  3	0	0	0	0	0	0	0	#_48
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_49
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_50
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_51
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_52
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_53
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_54
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_55
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_56
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_57
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_58
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_59
 -99	 10	  -1	   1	99	0	-99	0	0	0	0	0	0	0	#_60
  45	325	 100	 100	99	0	  2	0	0	0	0	0	0	0	#_61
  -7	  7	  -7	  -7	99	0	 -3	0	0	0	0	0	0	0	#_62
  -7	  7	   0	   0	99	0	  4	0	0	0	0	0	0	0	#_63
  -7	  7	   0	   0	99	0	  5	0	0	0	0	0	0	0	#_64
-999	  9	-999	-999	99	0	 -4	0	0	0	0	0	0	0	#_65
-999	  9	   0	   0	99	0	  3	0	0	0	0	0	0	0	#_66
#_AgeSelex
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_1 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_2 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_3 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_4 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_5 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_6 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_7 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_8 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_9 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_10
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_11
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_12
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_13
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_14
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_15
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_16
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_17
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_18
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_19
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_20
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_21
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_22
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_23
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_24
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_25
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_26
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_27
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_28
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_29
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_30
#_no timevary selex parameters
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_Factor	Fleet	Value
    1	 1	        0	#_1         
    2	 1	        0	#_2         
    3	 1	        0	#_3         
    4	 1	      0.1	#_4         
    5	 1	        0	#_5         
    6	 1	        0	#_6         
    7	 1	      0.1	#_7         
    1	 2	        0	#_8         
    2	 2	        0	#_9         
    3	 2	        0	#_10        
    4	 2	0.0833333	#_11        
    5	 2	        0	#_12        
    6	 2	        0	#_13        
    7	 2	      0.1	#_14        
    1	 3	        0	#_15        
    2	 3	        0	#_16        
    3	 3	        0	#_17        
    4	 3	      0.1	#_18        
    5	 3	        0	#_19        
    6	 3	        0	#_20        
    7	 3	      0.1	#_21        
    1	 4	        0	#_22        
    2	 4	        0	#_23        
    3	 4	        0	#_24        
    4	 4	      0.1	#_25        
    5	 4	        0	#_26        
    6	 4	        0	#_27        
    7	 4	      0.1	#_28        
    1	 5	        0	#_29        
    2	 5	        0	#_30        
    3	 5	        0	#_31        
    4	 5	      0.1	#_32        
    5	 5	        0	#_33        
    6	 5	        0	#_34        
    7	 5	      0.1	#_35        
    1	 6	        0	#_36        
    2	 6	        0	#_37        
    3	 6	        0	#_38        
    4	 6	0.0555556	#_39        
    5	 6	        0	#_40        
    6	 6	        0	#_41        
    7	 6	      0.1	#_42        
    1	 7	        0	#_43        
    2	 7	        0	#_44        
    3	 7	        0	#_45        
    4	 7	      0.1	#_46        
    5	 7	        0	#_47        
    6	 7	        0	#_48        
    7	 7	      0.1	#_49        
    1	 8	        0	#_50        
    2	 8	        0	#_51        
    3	 8	        0	#_52        
    4	 8	      0.1	#_53        
    5	 8	        0	#_54        
    6	 8	        0	#_55        
    7	 8	      0.1	#_56        
    1	 9	        0	#_57        
    2	 9	        0	#_58        
    3	 9	        0	#_59        
    4	 9	      0.1	#_60        
    5	 9	        0	#_61        
    6	 9	        0	#_62        
    7	 9	0.0833333	#_63        
    1	10	        0	#_64        
    2	10	        0	#_65        
    3	10	        0	#_66        
    4	10	      0.1	#_67        
    5	10	        0	#_68        
    6	10	        0	#_69        
    7	10	0.0434783	#_70        
    1	11	        0	#_71        
    2	11	        0	#_72        
    3	11	        0	#_73        
    4	11	   0.0625	#_74        
    5	11	        0	#_75        
    6	11	        0	#_76        
    7	11	      0.1	#_77        
    1	12	        0	#_78        
    2	12	        0	#_79        
    3	12	        0	#_80        
    4	12	0.0294118	#_81        
    5	12	        0	#_82        
    6	12	        0	#_83        
    7	12	      0.1	#_84        
    1	13	        0	#_85        
    2	13	        0	#_86        
    3	13	        0	#_87        
    4	13	      0.1	#_88        
    5	13	        0	#_89        
    6	13	        0	#_90        
    7	13	      0.1	#_91        
    1	14	        0	#_92        
    2	14	        0	#_93        
    3	14	        0	#_94        
    4	14	      0.1	#_95        
    5	14	        0	#_96        
    6	14	        0	#_97        
    7	14	      0.1	#_98        
    1	15	 0.190777	#_99        
    2	15	        0	#_100       
    3	15	        0	#_101       
    4	15	0.0666667	#_102       
    5	15	        0	#_103       
    6	15	        0	#_104       
    7	15	0.0666667	#_105       
-9999	 0	        0	#_terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 38 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    1	 1	1	0	1	#_1         
    1	 2	1	0	1	#_2         
    1	 3	1	0	1	#_3         
    1	 4	1	0	1	#_4         
    1	 5	1	0	1	#_5         
    1	 6	1	0	1	#_6         
    1	 7	1	0	1	#_7         
    1	 8	1	0	1	#_8         
    1	 9	1	0	1	#_9         
    1	10	1	0	1	#_10        
    1	11	1	0	1	#_11        
    1	12	1	0	1	#_12        
    1	13	1	0	1	#_13        
    1	14	1	0	1	#_14        
    1	15	1	1	1	#_15        
    4	 1	1	1	0	#_16        
    4	 2	1	1	0	#_17        
    4	 3	1	1	0	#_18        
    4	 4	1	0	0	#_19        
    4	 5	1	1	0	#_20        
    4	 6	1	1	0	#_21        
    4	 7	1	1	0	#_22        
    4	 8	1	1	0	#_23        
    4	 9	1	0	0	#_24        
    4	10	1	1	0	#_25        
    4	11	1	1	0	#_26        
    4	12	1	1	0	#_27        
    4	13	1	1	0	#_28        
    4	14	1	1	0	#_29        
    4	15	1	1	0	#_30        
    6	 2	1	1	1	#_31        
    6	10	1	1	1	#_32        
    6	 1	1	1	1	#_33        
    6	 3	1	1	1	#_34        
    6	15	1	1	1	#_35        
    6	 9	1	1	1	#_36        
    6	 6	1	1	1	#_37        
    6	 7	1	1	1	#_38        
-9999	 0	0	0	0	#_terminator
#
0 # 0/1 read specs for more stddev reporting
#
999
