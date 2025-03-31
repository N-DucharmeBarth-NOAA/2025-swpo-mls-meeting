

# Nicholas Ducharme-Barth
# 2025/03/18
# R code to make a baseline stock synthesis model version of the mfcl/1979_20p3 model
# Update to an annual model with length based selectivity
# Keep single sex, and fix initial F (correct the initial F setting)
# Only keep seasonal index for Q4

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(FLR4MFCL)
    library(frqit)
    library(r4ss)

#_____________________________________________________________________________________________________________________________
# define paths
	proj_dir = this.path::this.proj()
	dir_model = paste0(proj_dir,"/models/")
    dir_base_mfcl = paste0(dir_model,"mfcl/1979_20p3/")
    dir_base_stock_synthesis = paste0(dir_model,"stock-synthesis/00-npo-mako-base-file/")
    dir_helper_fns = paste0(proj_dir,"/code/helper-fns/")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(paste0(dir_helper_fns,(list.files(dir_helper_fns))),source)

#_____________________________________________________________________________________________________________________________
# read in baseline mfcl files
    base_frq = parse_frq(paste0(dir_base_mfcl,"mls.frq"))
    base_ini = read.MFCLIni(paste0(dir_base_mfcl,"mls.ini"), nseasons=4)
    base_par = read.MFCLPar(paste0(dir_base_mfcl,"06.par"), first.yr=1979)
    base_rep = read.MFCLRep(paste0(dir_base_mfcl,"plot-06.par.rep"))

#_____________________________________________________________________________________________________________________________
# run new version of stock synthesis
    # file.copy(from=paste0(proj_dir,"/executables/stock-synthesis/3.30.23.1/ss3_win.exe"),to=dir_base_stock_synthesis)
    # run(dir=dir_base_stock_synthesis,exe="ss3_win.exe")

#_____________________________________________________________________________________________________________________________
# read in baseline stock synthesis files
    tmp_starter = SS_readstarter(file=paste0(dir_base_stock_synthesis,"starter.ss_new"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=paste0(dir_base_stock_synthesis,"control.ss_new"),datlist = paste0(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_data = SS_readdat(file=paste0(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_forecast = SS_readforecast(file=paste0(dir_base_stock_synthesis,"forecast.ss_new"),verbose=FALSE)

#_____________________________________________________________________________________________________________________________
# create new directory for stock synthesis mls files
    dir_mls_stock_synthesis_base = paste0(dir_model,"stock-synthesis/25-mls-base-1979-corrected/")
    dir.create(dir_mls_stock_synthesis_base,recursive=TRUE)

#_____________________________________________________________________________________________________________________________
# update starter
    tmp_starter$min_age_summary_bio = 0
    tmp_starter$depl_basis = 1
    tmp_starter$SPR_basis = 4
    SS_writestarter(tmp_starter,dir=dir_mls_stock_synthesis_base,overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update forecast
    tmp_forecast$SPRtarget = 0.3
    tmp_forecast$Btarget = 0.3
    SS_writeforecast(tmp_forecast,dir=dir_mls_stock_synthesis_base,overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update data
    tmp_data$Comments = "#C No comments"
    tmp_data$styr = min(cateffpen(base_frq)$year,na.rm=TRUE)
    tmp_data$endyr = max(cateffpen(base_frq)$year,na.rm=TRUE)
    tmp_data$nseas = 4
    tmp_data$months_per_seas = c(3,3,3,3)
    tmp_data$spawn_month = 1
    tmp_data$spawn_seas = 1 # deprecated
    tmp_data$Nsexes = 1
    tmp_data$Nages = unname(dimensions(base_ini)["agecls"]/4)
    tmp_data$Nareas = 1
    tmp_data$Nfleets = n_fisheries(base_frq)-3
    tmp_data$Nfleet = tmp_data$Nfleets - 1
    tmp_data$Nsurveys = 1

    # define fleet info
    tmp_fleetname = c("F01_LL.JP.1","F02_LL.JP.2","F03_LL.JP.3","F04_LL.JP.4","F05_LL.TW.4","F06_LL.AU.2","F07_LL.AU.3","F08_LL.NZ.3","F09.REC.AU.3","F10_REC.NZ.3","F11_LL.ALL.1","F12_LL.ALL.2","F13_LL.ALL.3","F14_LL.ALL.4","S01_INDEX.1-4")
    tmp_fleetinfo = data.frame(type=c(rep(1,tmp_data$Nfleets-1),3),
                               surveytiming=c(rep(-1,tmp_data$Nfleets-1),1),
                               area=rep(1,tmp_data$Nfleets),
                               units=rep(2,tmp_data$Nfleets),
                               need_catch_mult=rep(0,tmp_data$Nfleets),
                               fleetname=tmp_fleetname)
    tmp_data$fleetinfo = tmp_fleetinfo
    tmp_data$fleetnames = tmp_fleetinfo$fleetname
    tmp_data$surveytiming = tmp_fleetinfo$surveytiming
    tmp_data$units_of_catch = tmp_fleetinfo$units
    tmp_data$areas = tmp_fleetinfo$area

    tmp_data$fleetinfo1 = rbind(tmp_data$surveytiming,tmp_data$areas,c(rep(1,tmp_data$Nfleets-1),3))
    rownames(tmp_data$fleetinfo1) = c("surveytiming","areas","type")
    colnames(tmp_data$fleetinfo1) = tmp_data$fleetnames

    tmp_data$fleetinfo2 = rbind(tmp_data$units_of_catch,rep(0,tmp_data$Nfleets))
    rownames(tmp_data$fleetinfo2) = c("units","need_catch_mult")
    colnames(tmp_data$fleetinfo2) = tmp_data$fleetnames

    # define catch
    # fishery 2 will get the equilibrium catch since it has the greatest catch in 1979
    # equilibrium catch will be the catch averaged across all fisheries for the first 3 years
    # remember catch needs to be divided by 1000 (stock synthesis records catch in 1000s of fish)
    cateffpen_dt = as.data.table(cateffpen(base_frq))
    eq_catch = mean(cateffpen_dt[year %in% c(tmp_data$styr + 0:2),.(catch=sum(catch)/1000),by=year]$catch)
    
    tmp_catch.list = as.list(rep(NA,tmp_data$Nfleets-1))
    for(i in seq_along(tmp_catch.list)){
        tmp_catch.list[[i]] = cateffpen_dt[fishery==i,.(year,month,fishery,catch)] %>%
                         setnames(.,c("month","fishery"),c("seas","fleet")) %>%
                         .[,seas:=c(1,1,1,2,2,2,3,3,3,4,4,4)[seas]] %>%
                         .[,catch:=catch/1000] %>%
                         .[,catch_se:=0.01] %>%
                         .[,.(year,seas,fleet,catch,catch_se)] %>%
                         .[catch>0]
        if(i == 2){
            tmp_catch.list[[i]] = rbind(data.table(year=rep(-999,4),seas=c(1,2,3,4),fleet=rep(i,4),catch=rep(eq_catch/4,4),catch_se=0.01),tmp_catch.list[[i]])
        }
    }
    tmp_catch = rbindlist(tmp_catch.list) %>%
                as.data.frame(.)
    tmp_data$catch = tmp_catch

    # CPUEinfo
    tmp_cpueinfo = data.frame(fleet=1:tmp_data$Nfleets,
                              units=rep(0,tmp_data$Nfleets),
                              errtype=rep(0,tmp_data$Nfleets),
                              SD_report=rep(0,tmp_data$Nfleets))
    rownames(tmp_cpueinfo) = tmp_fleetname
    tmp_data$CPUEinfo = tmp_cpueinfo

    # define CPUE
    # penalty column appears to be cv
    # use fishery 18 (Q4) index
    tmp_cpue = cateffpen_dt[fishery==18,.(cpue=catch/effort,cv=penalty),by=year] %>%
               .[,obs:=cpue/mean(cpue)] %>%
               .[,month:=11] %>%
               .[,se_log:=sqrt(log(1+cv^2))] %>%
               .[,index:=tmp_data$Nfleets] %>%
               .[,.(year,month,index,obs,se_log)] %>%
               as.data.frame(.)
    tmp_data$CPUE = tmp_cpue

    # define population length structure
    mfcl_bin_lower = seq(from=lf_range(base_frq)[3],by=lf_range(base_frq)[4],length.out=lf_range(base_frq)[2])
    mfcl_bin_mid = seq(from=lf_range(base_frq)[3],by=lf_range(base_frq)[4],length.out=lf_range(base_frq)[2]) + 0.5*lf_range(base_frq)[4]
    mfcl_bin_upper = seq(from=lf_range(base_frq)[3],by=lf_range(base_frq)[4],length.out=lf_range(base_frq)[2]) + lf_range(base_frq)[4]

    tmp_data$minimum_size = 5
    tmp_data$maximum_size = max(mfcl_bin_upper)
    tmp_data$lbin_vector_pop = seq(from=tmp_data$minimum_size,to=tmp_data$maximum_size,by=tmp_data$binwidth)
    tmp_data$N_lbinspop = length(tmp_data$lbin_vector_pop)

    # define length info
    tmp_leninfo = data.table(mintailcomp=rep(-1e-04,nrow(tmp_data$fleetinfo))) %>%
                  .[,addtocomp:=1e-04] %>%
                  .[,combine_M_F:=0] %>%
                  .[,CompressBins:=0] %>%
                  .[,CompError:=0] %>%
                  .[,ParmSelect:=0] %>%
                  .[,minsamplesize:=0.001] %>%
                  as.data.frame(.)
    rownames(tmp_leninfo) = tmp_fleetname
    tmp_data$len_info = tmp_leninfo
    tmp_data$N_lbins = length(mfcl_bin_lower)
    tmp_data$lbin_vector = mfcl_bin_lower

    # define length composition data
    tmp_lencomp = as.data.table(lnfrq(base_frq)) %>%
                  na.omit(.) %>%
                  as.data.frame(.)
    tmp_lencomp_a = tmp_lencomp[,c("year","month","fishery")]
    colnames(tmp_lencomp_a) = c("year","month","fleet")
    
    tmp_lencomp_b = tmp_lencomp[,-c(1:4)]
    tmp_lencomp_c = tmp_lencomp_b
    tmp_lencomp_c[tmp_lencomp_c>0] = 0
    colnames(tmp_lencomp_c) = paste0("m",colnames(tmp_lencomp_b))
    colnames(tmp_lencomp_b) = paste0("f",colnames(tmp_lencomp_b))
    
    tmp_lencomp_a$sex = 0
    tmp_lencomp_a$part = 0
    tmp_lencomp_a$Nsamp = rowSums(tmp_lencomp_b)

    # exclude comps from fisheries 15-17
    # rename fishery 18 as 15
    # exclude males

    tmp_lencomp = as.data.table(cbind(tmp_lencomp_a,tmp_lencomp_b))
    tmp_lencomp = tmp_lencomp[!(fleet%in%c(15:17))] %>%
                    .[fleet==18,fleet:=15]

    tmp_data$lencomp = as.data.frame(tmp_lencomp)

    # update generalized size comp
    tmp_data$N_sizefreq_methods_rd = 1
    tmp_data$N_sizefreq_methods = 1
    
    tmp_data$units_per_method = 2
    tmp_data$scale_per_method = 1
    tmp_data$mincomp_per_method = 0.001

    # need to re-aggregate comp data beginning at bin 139 make 4 kg
    tmp_wtcomp_init = as.data.table(wtfrq(base_frq)) %>%
                      melt(.,id.vars=c("year","month","week","fishery")) %>%
                      setnames(.,"variable","bin") %>%
                      .[,bin:=as.numeric(as.character(bin))]
    
    old_wt_bins = seq(from=139,to=249,by=2)
    new_wt_bins = rep(NA,length(old_wt_bins))
    for(i in seq_along(old_wt_bins)){
        if(i %% 2 != 0){
            new_wt_bins[i] = old_wt_bins[i]
        } else {
            new_wt_bins[i] = old_wt_bins[i-1]
        }
    }

    wt_bins_dt = data.table(bin=old_wt_bins,new_bin=new_wt_bins)

    tmp_wtcomp_init_a = tmp_wtcomp_init[bin<139]
    tmp_wtcomp_init_b = tmp_wtcomp_init[bin>=139] %>%
                        merge(.,wt_bins_dt,by="bin") %>%
                        .[,.(value=sum(value)),by=.(year,month,week,fishery,new_bin)] %>%
                        setnames(.,"new_bin","bin")
    tmp_wtcomp = rbind(tmp_wtcomp_init_a,tmp_wtcomp_init_b) %>%
                        dcast(.,year+month+week+fishery~bin) %>%
                  na.omit(.) %>%
                  as.data.frame(.)
    tmp_wtcomp_a = tmp_wtcomp[,c("year","month","fishery")]
    tmp_wtcomp_a = cbind(rep(1,nrow(tmp_wtcomp_a)),tmp_wtcomp_a)
    colnames(tmp_wtcomp_a) = c("method","year","month","fleet")
    
    tmp_wtcomp_b = tmp_wtcomp[,-c(1:4)]
    all_new_wt_bins = as.numeric(colnames(tmp_wtcomp_b))
    tmp_wtcomp_c = tmp_wtcomp_b
    tmp_wtcomp_c[tmp_wtcomp_c>0] = 0
    colnames(tmp_wtcomp_c) = paste0("m",colnames(tmp_wtcomp_b))
    colnames(tmp_wtcomp_b) = paste0("f",colnames(tmp_wtcomp_b))
    
    tmp_wtcomp_a$sex = 0
    tmp_wtcomp_a$part = 0
    tmp_wtcomp_a$Nsamp = rowSums(tmp_wtcomp_b)

    # exclude comps from fisheries 15-17
    # rename fishery 18 as 15
    # exclude males
    tmp_wtcomp = as.data.table(cbind(tmp_wtcomp_a,tmp_wtcomp_b))
    tmp_wtcomp = tmp_wtcomp[!(fleet%in%c(15:17))] %>%
                    .[fleet==18,fleet:=15]
    
    tmp_data$Nobs_per_method = nrow(tmp_wtcomp)
    tmp_data$sizefreq_bins_list = list(all_new_wt_bins)
    tmp_data$nbins_per_method = length(all_new_wt_bins)
    tmp_data$sizefreq_data_list = list(as.data.frame(tmp_wtcomp))

    # misc
    tmp_data$comp_tail_compression = rep(-0.0001,tmp_data$Nfleets)
    tmp_data$add_to_comp = rep(0.0001,tmp_data$Nfleets)
    tmp_data$max_combined_lbin = rep(0,tmp_data$Nfleets)    

    SS_writedat(tmp_data, outfile=paste0(dir_mls_stock_synthesis_base,"data.ss"), overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update control
    tmp_ctl$Comments = "#C No comments"
    tmp_ctl$nseas = tmp_data$nseas
    tmp_ctl$N_areas = tmp_data$N_areas
    tmp_ctl$Nages = tmp_data$Nages
    tmp_ctl$Nsexes = tmp_data$Nsexes
    tmp_ctl$Npopbins = tmp_data$N_lbinspop
    tmp_ctl$Nfleets = tmp_data$Nfleets
    tmp_ctl$fleetnames = tmp_data$fleetnames

    # recruitment timing and distribution
    tmp_ctl$recr_dist_method = 4
    tmp_ctl$Block_Design[[1]] = rep(tmp_data$styr-1,2)

    # natural mortality
    tmp_M_dt = data.table(M = m_at_age(base_rep)) %>%
            .[,age:=seq(from=0.25,by=0.25,length.out=length(m_at_age(base_rep)))] %>%
            .[,int_age:=floor(age)] %>%
            .[,.(M=mean(M)),by=int_age]
    tmp_M = matrix((tmp_M_dt$M)*4,nrow=1)
    rownames(tmp_M) = "natM1"
    colnames(tmp_M) = paste0("Age_",tmp_M_dt$int_age)
    tmp_ctl$natM = as.data.frame(tmp_M)

    # growth
    tmp_laa_dt = as.data.table(mean_laa(base_rep)) %>%
                 .[,.(age,season,value)] %>%
                 .[,dec_age:=as.numeric(age)+(as.numeric(season))/4]
    
    tmp_sdlaa_dt = as.data.table(sd_laa(base_rep)) %>%
                 .[,.(age,season,value)] %>%
                 .[,dec_age:=as.numeric(age)+(as.numeric(season))/4]
    
    
    tmp_ctl$Growth_Age_for_L1 = 0.25
    tmp_ctl$Growth_Age_for_L2 = 10
    tmp_ctl$CV_Growth_Pattern = 2 # make CV a function of length at age; like MFCL

    # fecundity
    tmp_ctl$fecundity_option = 3 # makes fecundity a function of weight

    # add MG params
    # only grab female parameters + last 2 rows
    nrow_mg_params = nrow(tmp_ctl$MG_parms)
    tmp_MG_parms = tmp_ctl$MG_parms[c(grep("_Fem_",rownames(tmp_ctl$MG_parms),fixed=TRUE),nrow_mg_params+(-1:0)),]

    # change growth L1, L2, k, sd1, sd2
    tmp_MG_parms$INIT[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_laa_dt[dec_age==0.25]$value
    tmp_MG_parms$PRIOR[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_laa_dt[dec_age==0.25]$value
    
    tmp_MG_parms$INIT[grep("L_at_Amax_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_laa_dt[dec_age==10]$value
    tmp_MG_parms$PRIOR[grep("L_at_Amax_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_laa_dt[dec_age==10]$value
    
    tmp_MG_parms$INIT[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = growth(base_par)[3,1]*4
    tmp_MG_parms$PRIOR[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = growth(base_par)[3,1]*4
    tmp_MG_parms$HI[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.99
    
    tmp_MG_parms$INIT[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_sdlaa_dt[dec_age==1]$value
    tmp_MG_parms$PRIOR[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_sdlaa_dt[dec_age==1]$value
    tmp_MG_parms$HI[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = 20
    
    tmp_MG_parms$INIT[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_sdlaa_dt[dec_age==10]$value
    tmp_MG_parms$PRIOR[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_sdlaa_dt[dec_age==10]$value
    tmp_MG_parms$HI[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = 20

    # change wtlen1 & wtlen2
    tmp_MG_parms$INIT[grep("Wtlen_1_",rownames(tmp_MG_parms),fixed=TRUE)] = lw_params(base_ini)[1]
    tmp_MG_parms$PRIOR[grep("Wtlen_1_",rownames(tmp_MG_parms),fixed=TRUE)] = lw_params(base_ini)[1]
    
    tmp_MG_parms$INIT[grep("Wtlen_2_",rownames(tmp_MG_parms),fixed=TRUE)] = lw_params(base_ini)[2]
    tmp_MG_parms$PRIOR[grep("Wtlen_2_",rownames(tmp_MG_parms),fixed=TRUE)] = lw_params(base_ini)[2]
    
    # mat 50% & mat slope
    tmp_mat_at_length = as.data.frame(cbind(len = mfcl_bin_mid, mat = mat_at_length(base_par)))

    fit_logistic = nls(mat ~ 1/(1 + exp(mat_slope*(len-mat50))),
                 data = tmp_mat_at_length,
                 start = list(mat_slope = -1, mat50 = 180))
    

    tmp_MG_parms$INIT[grep("Mat50%_Fem_",rownames(tmp_MG_parms),fixed=TRUE)] = summary(fit_logistic)$parameters["mat50","Estimate"]
    tmp_MG_parms$PRIOR[grep("Mat50%_Fem_",rownames(tmp_MG_parms),fixed=TRUE)] = summary(fit_logistic)$parameters["mat50","Estimate"]
    
    tmp_MG_parms$INIT[grep("Mat_slope_",rownames(tmp_MG_parms),fixed=TRUE)] = summary(fit_logistic)$parameters["mat_slope","Estimate"]
    tmp_MG_parms$PRIOR[grep("Mat_slope_",rownames(tmp_MG_parms),fixed=TRUE)] = summary(fit_logistic)$parameters["mat_slope","Estimate"]
    
    # fec a & b; set both to 1
    tmp_MG_parms$INIT[grep("Eggs_alpha_",rownames(tmp_MG_parms),fixed=TRUE)] = 1
    tmp_MG_parms$PRIOR[grep("Eggs_alpha_",rownames(tmp_MG_parms),fixed=TRUE)] = 1
    
    tmp_MG_parms$INIT[grep("Eggs_beta_",rownames(tmp_MG_parms),fixed=TRUE)] = 1
    tmp_MG_parms$PRIOR[grep("Eggs_beta_",rownames(tmp_MG_parms),fixed=TRUE)] = 1
    
    tmp_ctl$MG_parms = tmp_MG_parms

    # update virgin recruitment
    tmp_ctl$SR_parms["SR_LN(R0)","INIT"] = 6
    tmp_ctl$SR_parms["SR_LN(R0)","PRIOR"] = 6

    # update steepness
    tmp_ctl$SR_parms["SR_BH_steep","INIT"] = 0.8
    tmp_ctl$SR_parms["SR_BH_steep","PRIOR"] = 0.8

    # update sigmaR
    tmp_ctl$SR_parms["SR_sigmaR","INIT"] = 0.2
    tmp_ctl$SR_parms["SR_sigmaR","PRIOR"] = 0.2

    # recruitment deviation setup
    tmp_ctl$MainRdevYrFirst = tmp_data$styr
    tmp_ctl$MainRdevYrLast = tmp_data$endyr - 1 # late dev for terminal year
    tmp_ctl$recdev_early_start = tmp_data$styr-10
    tmp_ctl$last_early_yr_nobias_adj = tmp_data$styr + 10 # initial values
    tmp_ctl$first_yr_fullbias_adj = tmp_data$endyr - 16 # initial values
    tmp_ctl$last_yr_fullbias_adj = tmp_data$endyr - 1 # initial values
    tmp_ctl$first_recent_yr_nobias_adj = tmp_data$endyr - 1
    tmp_ctl$max_bias_adj = 0.2

    # initial F
    tmp_initF = matrix(0,nrow=sum(tmp_data$catch$year==-999),ncol=7)
    colnames(tmp_initF) = c("LO","HI","INIT","PRIOR","PRIOR_SD","PR_type","PHASE")
    tmp_initF[,"LO"] = 0
    tmp_initF[,"HI"] = 5
    tmp_initF[,"INIT"] = tail(as.matrix(tmp_ctl$natM)[1,],n=1)/4 * 0.7
    tmp_initF[,"PRIOR"] = 0.2
    tmp_initF[,"PRIOR_SD"] = 99
    tmp_initF[,"PR_type"] = 0
    tmp_initF[,"PHASE"] = -1 # don't estimate

    tmp_ctl$init_F = as.data.frame(tmp_initF)

    # catchability options; define for each survey
    tmp_ctl$Q_options = data.frame(fleet=15,link=1,link_info=0,extra_se=0,biasadj=0,float=1)
    rownames(tmp_ctl$Q_options) = "S01_INDEX.1-4"

    tmp_ctl$Q_parms = tmp_ctl$Q_parms[1,]
    rownames(tmp_ctl$Q_parms) = "LnQ_base_S01_INDEX.1-4"

    # size based selectivity
    tmp_size_selex_types = as.data.frame(matrix(0,nrow=length(tmp_ctl$fleetnames),ncol=4))
    rownames(tmp_size_selex_types) = tmp_ctl$fleetnames
    colnames(tmp_size_selex_types) = colnames(tmp_ctl$size_selex_types)
    tmp_size_selex_types$Pattern = 24
    tmp_size_selex_types$Pattern[c(5,8,11,13,14)] = 5
    tmp_size_selex_types$Pattern[10] = 24 # nz rec fishery
    tmp_size_selex_types$Special[c(5,8,11,13,14)] = c(4,7,1,3,4)
    tmp_ctl$size_selex_types = tmp_size_selex_types

    # age based selectivity
    tmp_age_selex_types = as.data.frame(matrix(0,nrow=length(tmp_ctl$fleetnames),ncol=4))
    rownames(tmp_age_selex_types) = tmp_ctl$fleetnames
    colnames(tmp_age_selex_types) = colnames(tmp_ctl$age_selex_types)
    tmp_age_selex_types$Pattern = 11
    tmp_ctl$age_selex_types = tmp_age_selex_types

    # define selex shapes
    # LO, HI, INIT, PHASE 
    tmp_size_selex_parms = rbind(
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)),
    make_size_selex_par_5(),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)),
    make_size_selex_par_5(),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,7,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,-5),
                          start_logit=c( -999.0,-2,-999.00000,4),
                          end_logit=c( -999.0,9.0,9,-4)),
    make_size_selex_par_5(),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)),
    make_size_selex_par_5(),
    make_size_selex_par_5(),
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-7,7,0,3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999.00000,-4)))
    tmp_ctl$size_selex_parms = as.data.frame(tmp_size_selex_parms)
    # fix selex of certain parameters
    # tmp_ctl$size_selex_parms$PHASE[c(3,8,20,21,29,34,42,56,57)] = -1*tmp_ctl$size_selex_parms$PHASE[c(3,8,20,21,29,34,42,56,57)]
    # tmp_ctl$size_selex_parms$INIT[c(3,8,20,21,29,34,42,56,57)] = c(7,-7,-7,7,-7,-7,-7,-7,7)

    tmp_ctl$age_selex_parms = as.data.frame(rbindlist(lapply(rep(tmp_ctl$Nages,length(tmp_ctl$fleetnames)),make_age_selex_par_11)))

    # variance adjustment list
    tmp_ctl$Variance_adjustment_list = expand.grid(Factor=1:7,Fleet=1:length(tmp_ctl$fleetnames),Value=0)
    tmp_ctl$Variance_adjustment_list$Value[which(tmp_ctl$Variance_adjustment_list$Factor%in%c(4,7))] = 1

    # lambdas
    tmp_lambdas_surv= data.table(like_comp=rep(1,length(tmp_ctl$fleetnames)),fleet=1:length(tmp_ctl$fleetnames)) %>%
                      .[,phase:=1] %>%
                      .[,value:=0] %>%
                      .[,sizefreq_method:=1] %>%
                      .[fleet %in% c(15),value:=1]
    tmp_lambdas_lf = data.table(like_comp=rep(4,length(tmp_ctl$fleetnames)),fleet=1:length(tmp_ctl$fleetnames)) %>%
                      .[,phase:=1] %>%
                      .[,value:=0] %>%
                      .[,sizefreq_method:=0] %>%
                      .[fleet %in% unique(tmp_lencomp$fleet),value:=1]
    tmp_lambdas_gs = data.table(like_comp=rep(6,length(unique(tmp_wtcomp$fleet))),fleet=unique(tmp_wtcomp$fleet)) %>%
                      .[,phase:=1] %>%
                      .[,value:=1] %>%
                      .[,sizefreq_method:=1] 
    tmp_ctl$lambdas = as.data.frame(rbind(tmp_lambdas_surv,tmp_lambdas_lf,tmp_lambdas_gs))
    tmp_ctl$N_lambdas = nrow(tmp_ctl$lambdas)
    SS_writectl(tmp_ctl,paste0(dir_mls_stock_synthesis_base,"control.ss"),overwrite = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# run & summarize
    # run
    file.copy(from=paste0(proj_dir,"/executables/stock-synthesis/3.30.23.1/ss3_win.exe"),to=dir_mls_stock_synthesis_base)
    run(dir=dir_mls_stock_synthesis_base,exe="ss3_win.exe",show_in_console = TRUE)

    # summarize
    summarize_ss_model(dir_mls_stock_synthesis_base)
