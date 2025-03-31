#V3.30
#C SS3_Control_NP_MAK.xlsx
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
2 # recr_dist_method for parameters
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
1974 1974
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
#_Age_0	Age_1	Age_2	Age_3	Age_4	Age_5	Age_6	Age_7	Age_8	Age_9	Age_10	Age_11	Age_12	Age_13	Age_14	Age_15	Age_16	Age_17	Age_18	Age_19	Age_20	Age_21	Age_22	Age_23	Age_24	Age_25	Age_26	Age_27	Age_28	Age_29	Age_30	Age_31
0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	#_natM1
0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	0.128	#_natM2
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
1 #_First_Mature_Age
2 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
    5	     100	      60	      60	 99	0	 -3	0	0	0	0	  0	0	0	#_L_at_Amin_Fem_GP_1 
   50	     600	   293.1	   293.1	 99	0	 -4	0	0	0	0	  0	0	0	#_L_at_Amax_Fem_GP_1 
 0.01	     0.7	   0.128	   0.128	 99	0	 -5	0	0	0	0	  0	0	0	#_VonBert_K_Fem_GP_1 
 0.01	     0.3	     0.1	     0.1	 99	0	 -2	0	0	0	0	  0	0	0	#_CV_young_Fem_GP_1  
 0.01	     0.3	     0.1	     0.1	 99	0	 -3	0	0	0	0	  0	0	0	#_CV_old_Fem_GP_1    
   -3	       3	4.62e-05	4.62e-05	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_1_Fem_GP_1   
   -3	       5	    2.77	    2.77	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_2_Fem_GP_1   
    1	     300	 233.654	 233.654	 99	0	 -3	0	0	0	0	  0	0	0	#_Mat50%_Fem_GP_1    
 -200	       3	-0.14652	-0.14652	 99	0	 -3	0	0	0	0	  0	0	0	#_Mat_slope_Fem_GP_1 
   -3	      20	       6	       6	0.5	6	 -3	0	0	0	0	0.5	0	0	#_Eggs_alpha_Fem_GP_1
   -3	       3	       0	       0	0.5	6	 -3	0	0	0	0	0.5	0	0	#_Eggs_beta_Fem_GP_1 
    5	     100	      60	      60	 99	0	 -3	0	0	0	0	  0	0	0	#_L_at_Amin_Mal_GP_1 
   50	     600	     232	     232	 99	0	 -4	0	0	0	0	  0	0	0	#_L_at_Amax_Mal_GP_1 
 0.01	     0.7	   0.174	   0.174	 99	0	 -5	0	0	0	0	  0	0	0	#_VonBert_K_Mal_GP_1 
 0.01	     0.3	     0.1	     0.1	 99	0	 -2	0	0	0	0	  0	0	0	#_CV_young_Mal_GP_1  
 0.01	     0.3	     0.1	     0.1	 99	0	 -3	0	0	0	0	  0	0	0	#_CV_old_Mal_GP_1    
   -3	       5	 3.4e-05	 3.4e-05	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_1_Mal_GP_1   
   -3	       5	    2.84	    2.84	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_2_Mal_GP_1   
   -4	       4	       0	       1	 99	0	 -3	0	0	0	0	  0	0	0	#_RecrDist_GP_1      
   -4	       4	       0	       1	 99	0	 -3	0	0	0	0	  0	0	0	#_RecrDist_Area_1    
   -4	       4	       0	       1	 99	0	 -3	0	0	0	0	  0	0	0	#_RecrDist_month_1   
  0.1	      10	       1	       1	  1	6	 -1	0	0	0	0	  0	0	0	#_CohortGrowDev      
1e-06	0.999999	     0.5	     0.5	0.5	0	-99	0	0	0	0	  0	0	0	#_FracFemale_GP_1    
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
   2	  15	7.32641	    7	  99	0	 1	0	0	0	0	0	0	0	#_SR_LN(R0)  
 0.2	0.99	  0.317	0.317	1000	6	-2	0	0	0	0	0	0	0	#_SR_BH_steep
0.05	 1.9	    0.1	  0.1	1000	6	-4	0	0	0	0	0	0	0	#_SR_sigmaR  
  -4	   4	      0	    0	  99	0	-1	0	0	0	0	0	0	0	#_SR_regime  
   0	   0	      0	    0	  99	0	-1	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1975 # first year of main recr_devs; early devs can preceed this era
2016 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase
1 # (0/1) to read 13 advanced options
-5 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
2 #_recdev_early_phase
0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1982 #_last_yr_nobias_adj_in_MPD; begin of ramp
2007.1 #_first_yr_fullbias_adj_in_MPD; begin of plateau
2007.4 #_last_yr_fullbias_adj_in_MPD
2019.9 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.3421 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
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
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0	5	0.0403514	0.2	99	0	1	#_InitF_seas_1_flt_9F9_JPN_SS_II
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
   19	1	0	0	0	1	#_S1_US_SS   
   20	1	0	0	0	1	#_S2_US_DS   
   21	1	0	0	0	1	#_S3_TW_LRG  
   22	1	0	0	0	1	#_S4_JPN_SS  
   23	1	0	0	0	1	#_S5_JPN_RTV 
   24	1	0	0	0	1	#_S6_JPN_OBS 
   25	1	0	0	0	1	#_S7_JPN_GEO 
   26	1	0	0	0	1	#_S8_MEX     
   27	1	0	0	0	1	#_S9_JPN_SS_I
-9999	0	0	0	0	0	#_terminator 
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-25	25	-8.30978	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S1_US_SS(19)   
-25	25	-7.94438	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S2_US_DS(20)   
-25	25	-7.83906	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S3_TW_LRG(21)  
-25	25	-7.90693	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S4_JPN_SS(22)  
-25	25	-7.80228	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S5_JPN_RTV(23) 
-25	25	-7.85765	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S6_JPN_OBS(24) 
-25	25	-7.84295	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S7_JPN_GEO(25) 
-25	25	-8.05204	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S8_MEX(26)     
-25	25	-7.80245	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S9_JPN_SS_I(27)
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	0	 0	#_1 F1_US_CA     
 5	0	0	 1	#_2 F2_US_HI_SS  
24	0	0	 0	#_3 F3_US_HI_DS  
24	0	0	 0	#_4 F4_US_DGN    
 5	0	0	 2	#_5 F5_US_REC    
24	0	0	 0	#_6 F6_TW_LRG_N  
 5	0	0	 6	#_7 F7_TW_LRG_S  
 5	0	0	 6	#_8 F8_TW_SML    
 5	0	0	 6	#_9 F9_JPN_SS_II 
24	0	0	 0	#_10 F10_JPN_DS  
 5	0	0	 9	#_11 F11_JPN_CST 
24	0	0	 0	#_12 F12_JPN_DFN 
 5	0	0	 9	#_13 F13_JPN_OTH 
24	0	0	 0	#_14 F14_MEX_NOR 
24	0	0	 0	#_15 F15_MEX_SOU 
 5	0	0	 1	#_16 F16_WCPFC   
 5	0	0	 1	#_17 F17_IATTC   
 5	0	0	 9	#_18 F18_JPN_SSII
 5	0	0	 2	#_19 S1_US_SS    
 5	0	0	 3	#_20 S2_US_DS    
 5	0	0	 6	#_21 S3_TW_LRG   
 5	0	0	 9	#_22 S4_JPN_SS   
 5	0	0	 9	#_23 S5_JPN_RTV  
 5	0	0	 9	#_24 S6_JPN_OBS  
 5	0	0	 9	#_25 S7_JPN_GEO  
 5	0	0	14	#_26 S8_MEX      
 5	0	0	 6	#_27 S9_JPN_SS_I 
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
11	0	0	0	#_1 F1_US_CA     
11	0	0	0	#_2 F2_US_HI_SS  
11	0	0	0	#_3 F3_US_HI_DS  
11	0	0	0	#_4 F4_US_DGN    
11	0	0	0	#_5 F5_US_REC    
11	0	0	0	#_6 F6_TW_LRG_N  
11	0	0	0	#_7 F7_TW_LRG_S  
11	0	0	0	#_8 F8_TW_SML    
11	0	0	0	#_9 F9_JPN_SS_II 
11	0	0	0	#_10 F10_JPN_DS  
11	0	0	0	#_11 F11_JPN_CST 
11	0	0	0	#_12 F12_JPN_DFN 
11	0	0	0	#_13 F13_JPN_OTH 
11	0	0	0	#_14 F14_MEX_NOR 
11	0	0	0	#_15 F15_MEX_SOU 
11	0	0	0	#_16 F16_WCPFC   
11	0	0	0	#_17 F17_IATTC   
11	0	0	0	#_18 F18_JPN_SSII
11	0	0	0	#_19 S1_US_SS    
11	0	0	0	#_20 S2_US_DS    
11	0	0	0	#_21 S3_TW_LRG   
11	0	0	0	#_22 S4_JPN_SS   
11	0	0	0	#_23 S5_JPN_RTV  
11	0	0	0	#_24 S6_JPN_OBS  
11	0	0	0	#_25 S7_JPN_GEO  
11	0	0	0	#_26 S8_MEX      
11	0	0	0	#_27 S9_JPN_SS_I 
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
  55	297.5	      55	148.87	99	0	 -2	0	0	0	0	  0	0	0	#_SizeSel_P_1_F1_US_CA(1)     
  -9	    4	-7.80488	 -4.56	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_2_F1_US_CA(1)     
  -1	   99	      99	  7.25	99	0	 -3	0	0	0	0	  0	0	0	#_SizeSel_P_3_F1_US_CA(1)     
  -1	   12	  9.8561	  7.61	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_4_F1_US_CA(1)     
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_5_F1_US_CA(1)     
  -9	    9	-3.82876	    -5	99	0	  4	0	0	0	0	  0	0	0	#_SizeSel_P_6_F1_US_CA(1)     
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F2_US_HI_SS(2)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F2_US_HI_SS(2)  
62.5	297.5	 167.798	148.87	99	0	  2	0	0	0	0	  0	0	0	#_SizeSel_P_1_F3_US_HI_DS(3)  
  -6	    4	      -6	 -4.56	99	0	 -3	0	0	0	0	  0	0	0	#_SizeSel_P_2_F3_US_HI_DS(3)  
  -1	    9	 7.39278	  7.25	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_3_F3_US_HI_DS(3)  
  -1	    9	 7.38934	  7.61	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_4_F3_US_HI_DS(3)  
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_5_F3_US_HI_DS(3)  
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_6_F3_US_HI_DS(3)  
62.5	297.5	 96.0996	148.87	99	0	  2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_F4_US_DGN(4)    
  -6	    4	      -6	 -4.56	99	0	 -3	0	0	0	0	  0	0	0	#_SizeSel_P_2_F4_US_DGN(4)    
  -1	    9	 6.89923	  7.25	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_3_F4_US_DGN(4)    
  -1	    9	 7.64571	  7.61	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_4_F4_US_DGN(4)    
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_5_F4_US_DGN(4)    
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_6_F4_US_DGN(4)    
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F5_US_REC(5)    
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F5_US_REC(5)    
62.5	297.5	 140.096	148.87	99	0	  2	0	0	0	0	  0	0	0	#_SizeSel_P_1_F6_TW_LRG_N(6)  
  -6	    4	      -6	 -4.56	99	0	 -3	0	0	0	0	  0	0	0	#_SizeSel_P_2_F6_TW_LRG_N(6)  
  -1	    9	 6.97823	  7.25	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_3_F6_TW_LRG_N(6)  
  -1	    9	 7.29406	  7.61	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_4_F6_TW_LRG_N(6)  
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_5_F6_TW_LRG_N(6)  
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_6_F6_TW_LRG_N(6)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F7_TW_LRG_S(7)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F7_TW_LRG_S(7)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F8_TW_SML(8)    
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F8_TW_SML(8)    
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F9_JPN_SS_II(9) 
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F9_JPN_SS_II(9) 
62.5	297.5	 162.156	148.87	99	0	  2	0	0	0	0	  0	0	0	#_SizeSel_P_1_F10_JPN_DS(10)  
  -6	    4	      -6	 -4.56	99	0	 -3	0	0	0	0	  0	0	0	#_SizeSel_P_2_F10_JPN_DS(10)  
  -1	    9	 6.18822	  7.25	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_3_F10_JPN_DS(10)  
  -1	    9	 7.64275	  7.61	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_4_F10_JPN_DS(10)  
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_5_F10_JPN_DS(10)  
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_6_F10_JPN_DS(10)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F11_JPN_CST(11) 
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F11_JPN_CST(11) 
62.5	297.5	  125.49	148.87	99	0	  2	0	0	0	0	  0	0	0	#_SizeSel_P_1_F12_JPN_DFN(12) 
  -6	    4	      -6	 -4.56	99	0	 -3	0	0	0	0	  0	0	0	#_SizeSel_P_2_F12_JPN_DFN(12) 
  -1	    9	 6.88435	  7.25	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_3_F12_JPN_DFN(12) 
  -1	    9	 6.86101	  7.61	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_4_F12_JPN_DFN(12) 
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_5_F12_JPN_DFN(12) 
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_6_F12_JPN_DFN(12) 
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F13_JPN_OTH(13) 
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F13_JPN_OTH(13) 
62.5	297.5	 100.097	148.87	99	0	  2	0	0	0	0	  0	0	0	#_SizeSel_P_1_F14_MEX_NOR(14) 
  -6	    4	      -6	 -4.56	99	0	 -3	0	0	0	0	  0	0	0	#_SizeSel_P_2_F14_MEX_NOR(14) 
  -1	    9	 8.30474	  7.25	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_3_F14_MEX_NOR(14) 
  -1	    9	 7.53331	  7.61	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_4_F14_MEX_NOR(14) 
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_5_F14_MEX_NOR(14) 
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_6_F14_MEX_NOR(14) 
62.5	297.5	 126.019	148.87	99	0	  2	0	0	0	0	  0	0	0	#_SizeSel_P_1_F15_MEX_SOU(15) 
  -6	    4	      -6	 -4.56	99	0	 -3	0	0	0	0	  0	0	0	#_SizeSel_P_2_F15_MEX_SOU(15) 
  -1	    9	 7.46993	  7.25	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_3_F15_MEX_SOU(15) 
  -1	    9	 6.65529	  7.61	99	0	  3	0	0	0	0	  0	0	0	#_SizeSel_P_4_F15_MEX_SOU(15) 
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_5_F15_MEX_SOU(15) 
-999	    9	    -999	    -5	99	0	 -4	0	0	0	0	  0	0	0	#_SizeSel_P_6_F15_MEX_SOU(15) 
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F16_WCPFC(16)   
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F16_WCPFC(16)   
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F17_IATTC(17)   
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F17_IATTC(17)   
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_F18_JPN_SSII(18)
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_F18_JPN_SSII(18)
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S1_US_SS(19)    
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S1_US_SS(19)    
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S2_US_DS(20)    
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S2_US_DS(20)    
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S3_TW_LRG(21)   
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S3_TW_LRG(21)   
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S4_JPN_SS(22)   
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S4_JPN_SS(22)   
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S5_JPN_RTV(23)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S5_JPN_RTV(23)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S6_JPN_OBS(24)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S6_JPN_OBS(24)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S7_JPN_GEO(25)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S7_JPN_GEO(25)  
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S8_MEX(26)      
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S8_MEX(26)      
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_1_S9_JPN_SS_I(27) 
 -99	   10	      -1	     1	99	0	-99	0	0	0	0	  0	0	0	#_SizeSel_P_2_S9_JPN_SS_I(27) 
#_AgeSelex
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F1_US_CA(1)     
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F1_US_CA(1)     
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F2_US_HI_SS(2)  
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F2_US_HI_SS(2)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F3_US_HI_DS(3)  
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F3_US_HI_DS(3)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F4_US_DGN(4)    
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F4_US_DGN(4)    
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F5_US_REC(5)    
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F5_US_REC(5)    
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F6_TW_LRG_N(6)  
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F6_TW_LRG_N(6)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F7_TW_LRG_S(7)  
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F7_TW_LRG_S(7)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F8_TW_SML(8)    
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F8_TW_SML(8)    
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F9_JPN_SS_II(9) 
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F9_JPN_SS_II(9) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F10_JPN_DS(10)  
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F10_JPN_DS(10)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F11_JPN_CST(11) 
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F11_JPN_CST(11) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F12_JPN_DFN(12) 
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F12_JPN_DFN(12) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F13_JPN_OTH(13) 
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F13_JPN_OTH(13) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F14_MEX_NOR(14) 
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F14_MEX_NOR(14) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F15_MEX_SOU(15) 
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F15_MEX_SOU(15) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F16_WCPFC(16)   
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F16_WCPFC(16)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F17_IATTC(17)   
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F17_IATTC(17)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F18_JPN_SSII(18)
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F18_JPN_SSII(18)
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S1_US_SS(19)    
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S1_US_SS(19)    
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S2_US_DS(20)    
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S2_US_DS(20)    
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S3_TW_LRG(21)   
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S3_TW_LRG(21)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S4_JPN_SS(22)   
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S4_JPN_SS(22)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S5_JPN_RTV(23)  
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S5_JPN_RTV(23)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S6_JPN_OBS(24)  
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S6_JPN_OBS(24)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S7_JPN_GEO(25)  
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S7_JPN_GEO(25)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S8_MEX(26)      
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S8_MEX(26)      
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S9_JPN_SS_I(27) 
10	100	31	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S9_JPN_SS_I(27) 
#_no timevary selex parameters
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_Data_type	Fleet	Value
    5	 1	       0	#_Variance_adjustment_list1 
    6	 1	       0	#_Variance_adjustment_list2 
    5	 2	       0	#_Variance_adjustment_list3 
    6	 2	       0	#_Variance_adjustment_list4 
    5	 3	       0	#_Variance_adjustment_list5 
    6	 3	       0	#_Variance_adjustment_list6 
    5	 4	       0	#_Variance_adjustment_list7 
    6	 4	       0	#_Variance_adjustment_list8 
    5	 5	       0	#_Variance_adjustment_list9 
    6	 5	       0	#_Variance_adjustment_list10
    5	 6	       0	#_Variance_adjustment_list11
    6	 6	       0	#_Variance_adjustment_list12
    5	 7	       0	#_Variance_adjustment_list13
    6	 7	       0	#_Variance_adjustment_list14
    5	 8	       0	#_Variance_adjustment_list15
    6	 8	       0	#_Variance_adjustment_list16
    5	 9	       0	#_Variance_adjustment_list17
    6	 9	       0	#_Variance_adjustment_list18
    5	10	       0	#_Variance_adjustment_list19
    6	10	       0	#_Variance_adjustment_list20
    5	11	       0	#_Variance_adjustment_list21
    6	11	       0	#_Variance_adjustment_list22
    4	12	    0.38	#_Variance_adjustment_list23
    5	12	       0	#_Variance_adjustment_list24
    6	12	       0	#_Variance_adjustment_list25
    5	13	       0	#_Variance_adjustment_list26
    6	13	       0	#_Variance_adjustment_list27
    4	14	     0.4	#_Variance_adjustment_list28
    5	14	       0	#_Variance_adjustment_list29
    6	14	       0	#_Variance_adjustment_list30
    5	15	       0	#_Variance_adjustment_list31
    6	15	       0	#_Variance_adjustment_list32
    5	16	       0	#_Variance_adjustment_list33
    6	16	       0	#_Variance_adjustment_list34
    5	17	       0	#_Variance_adjustment_list35
    6	17	       0	#_Variance_adjustment_list36
    5	18	       0	#_Variance_adjustment_list37
    6	18	       0	#_Variance_adjustment_list38
    1	19	 0.04333	#_Variance_adjustment_list39
    5	19	       0	#_Variance_adjustment_list40
    6	19	       0	#_Variance_adjustment_list41
    1	20	0.010455	#_Variance_adjustment_list42
    5	20	       0	#_Variance_adjustment_list43
    6	20	       0	#_Variance_adjustment_list44
    1	21	     0.2	#_Variance_adjustment_list45
    5	21	       0	#_Variance_adjustment_list46
    6	21	       0	#_Variance_adjustment_list47
    1	22	0.146852	#_Variance_adjustment_list48
    5	22	       0	#_Variance_adjustment_list49
    6	22	       0	#_Variance_adjustment_list50
    1	23	   0.147	#_Variance_adjustment_list51
    5	23	       0	#_Variance_adjustment_list52
    6	23	       0	#_Variance_adjustment_list53
    5	24	       0	#_Variance_adjustment_list54
    6	24	       0	#_Variance_adjustment_list55
    1	25	0.085795	#_Variance_adjustment_list56
    5	25	       0	#_Variance_adjustment_list57
    6	25	       0	#_Variance_adjustment_list58
    1	26	    0.09	#_Variance_adjustment_list59
    5	26	       0	#_Variance_adjustment_list60
    6	26	       0	#_Variance_adjustment_list61
    1	27	     0.1	#_Variance_adjustment_list62
    5	27	       0	#_Variance_adjustment_list63
    6	27	       0	#_Variance_adjustment_list64
-9999	 0	       0	#_terminator                
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 61 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    1	 1	1	0	1	#_Surv_F1_US_CA_Phz1                         
    1	 2	1	0	1	#_Surv_F2_US_HI_SS_Phz1                      
    1	 3	1	0	1	#_Surv_F3_US_HI_DS_Phz1                      
    1	 4	1	0	1	#_Surv_F4_US_DGN_Phz1                        
    1	 5	1	0	1	#_Surv_F5_US_REC_Phz1                        
    1	 6	1	0	1	#_Surv_F6_TW_LRG_N_Phz1                      
    1	 7	1	0	1	#_Surv_F7_TW_LRG_S_Phz1                      
    1	 8	1	0	1	#_Surv_F8_TW_SML_Phz1                        
    1	 9	1	0	1	#_Surv_F9_JPN_SS_II_Phz1                     
    1	10	1	0	1	#_Surv_F10_JPN_DS_Phz1                       
    1	11	1	0	1	#_Surv_F11_JPN_CST_Phz1                      
    1	12	1	0	1	#_Surv_F12_JPN_DFN_Phz1                      
    1	13	1	0	1	#_Surv_F13_JPN_OTH_Phz1                      
    1	14	1	0	1	#_Surv_F14_MEX_NOR_Phz1                      
    1	15	1	0	1	#_Surv_F15_MEX_SOU_Phz1                      
    1	16	1	0	1	#_Surv_F16_WCPFC_Phz1                        
    1	17	1	0	1	#_Surv_F17_IATTC_Phz1                        
    1	18	1	0	1	#_Surv_F18_JPN_SSII_Phz1                     
    1	19	1	1	1	#_Surv_S1_US_SS_Phz1                         
    1	20	1	0	1	#_Surv_S2_US_DS_Phz1                         
    1	21	1	1	1	#_Surv_S3_TW_LRG_Phz1                        
    1	22	1	0	1	#_Surv_S4_JPN_SS_Phz1                        
    1	23	1	1	1	#_Surv_S5_JPN_RTV_Phz1                       
    1	24	1	0	1	#_Surv_S6_JPN_OBS_Phz1                       
    1	25	1	0	1	#_Surv_S7_JPN_GEO_Phz1                       
    1	26	1	1	1	#_Surv_S8_MEX_Phz1                           
    1	27	1	1	0	#_Surv_S9_JPN_SS_I_Phz1                      
    4	 1	1	0	0	#_length_F1_US_CA_sizefreq_method_0_Phz1     
    4	 2	1	1	0	#_length_F2_US_HI_SS_sizefreq_method_0_Phz1  
    4	 3	1	1	0	#_length_F3_US_HI_DS_sizefreq_method_0_Phz1  
    4	 4	1	0	0	#_length_F4_US_DGN_sizefreq_method_0_Phz1    
    4	 5	1	0	0	#_length_F5_US_REC_sizefreq_method_0_Phz1    
    4	 6	1	0	0	#_length_F6_TW_LRG_N_sizefreq_method_0_Phz1  
    4	 7	1	0	0	#_length_F7_TW_LRG_S_sizefreq_method_0_Phz1  
    4	 8	1	0	0	#_length_F8_TW_SML_sizefreq_method_0_Phz1    
    4	 9	1	1	0	#_length_F9_JPN_SS_II_sizefreq_method_0_Phz1 
    4	10	1	1	0	#_length_F10_JPN_DS_sizefreq_method_0_Phz1   
    4	11	1	0	0	#_length_F11_JPN_CST_sizefreq_method_0_Phz1  
    4	12	1	1	0	#_length_F12_JPN_DFN_sizefreq_method_0_Phz1  
    4	13	1	0	0	#_length_F13_JPN_OTH_sizefreq_method_0_Phz1  
    4	14	1	1	0	#_length_F14_MEX_NOR_sizefreq_method_0_Phz1  
    4	15	1	1	0	#_length_F15_MEX_SOU_sizefreq_method_0_Phz1  
    4	16	1	0	0	#_length_F16_WCPFC_sizefreq_method_0_Phz1    
    4	17	1	0	0	#_length_F17_IATTC_sizefreq_method_0_Phz1    
    4	18	1	0	0	#_length_F18_JPN_SSII_sizefreq_method_0_Phz1 
    4	19	1	0	0	#_length_S1_US_SS_sizefreq_method_0_Phz1     
    4	20	1	0	0	#_length_S2_US_DS_sizefreq_method_0_Phz1     
    4	21	1	0	0	#_length_S3_TW_LRG_sizefreq_method_0_Phz1    
    4	22	1	0	0	#_length_S4_JPN_SS_sizefreq_method_0_Phz1    
    4	23	1	0	0	#_length_S5_JPN_RTV_sizefreq_method_0_Phz1   
    4	24	1	0	0	#_length_S6_JPN_OBS_sizefreq_method_0_Phz1   
    4	25	1	0	0	#_length_S7_JPN_GEO_sizefreq_method_0_Phz1   
    4	26	1	0	0	#_length_S8_MEX_sizefreq_method_0_Phz1       
    4	27	1	0	0	#_length_S9_JPN_SS_I_sizefreq_method_0_Phz1  
    6	 4	1	1	2	#_SizeFreq_F4_US_DGN_sizefreq_method_2_Phz1  
    6	 6	1	0	1	#_SizeFreq_F6_TW_LRG_N_sizefreq_method_1_Phz1
    7	 4	1	0	0	#_SizeAge_F4_US_DGN_sizefreq_method_0_Phz1   
    7	 6	1	0	0	#_SizeAge_F6_TW_LRG_N_sizefreq_method_0_Phz1 
    9	 9	1	1	0	#_init_equ_catch_F9_JPN_SS_II_Phz1           
   11	 1	1	0	0	#_parm_prior_F1_US_CA_Phz1                   
   12	 1	1	1	0	#_parm_dev_Phz1                              
-9999	 0	0	0	0	#_terminator                                 
#
0 # 0/1 read specs for more stddev reporting
#
999
