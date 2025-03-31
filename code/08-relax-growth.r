

# Nicholas Ducharme-Barth
# 2025/01/22
# Run model with relaxed growth


# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
    library(r4ss)

#_____________________________________________________________________________________________________________________________
# define paths
	proj_dir = this.path::this.proj()
	dir_model = file.path(proj_dir,"models")
    dir_base_stock_synthesis = file.path(dir_model,"stock-synthesis","06-exclude-more-comp")
    dir_helper_fns = file.path(proj_dir,"code","helper-fns")
    executable_dir_stem = file.path(proj_dir,"executables","stock-synthesis")


#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define OS
    linux_ws = is_linux_os() # flag to switch code depending on which workspace is being used

#________________________________________________________________________________________________________________________________________________________________________________________________________
# transfer files

    from_dir = dir_base_stock_synthesis
    to_dir = file.path(dir_model,"stock-synthesis","08-relax-growth")
    dir.create(to_dir,recursive=TRUE)

    ss_transfer(from_dir=from_dir,to_dir=to_dir,executable_dir_stem=executable_dir_stem,linux_ws=is_linux_os())

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in files
    # tmp_starter = SS_readstarter(file=file.path(to_dir,"starter.ss"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=file.path(to_dir,"control.ss"),datlist = file.path(to_dir,"data.ss"))
    tmp_data = SS_readdat(file=file.path(to_dir,"data.ss"))
    # tmp_forecast = SS_readforecast(file=file.path(to_dir,"forecast.ss"),verbose=FALSE)


#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; data
# add new fleet
# switch fleets

    tmp_data$Nfleets = tmp_data$Nfleets + 1
    tmp_data$Nfleet = tmp_data$Nfleet + 1
    tmp_data$Nsurveys = 1

    # define fleet info
    tmp_fleetname = c("F01_LL.JP.1","F02_LL.JP.2","F03_LL.JP.3","F04_LL.JP.4","F05_LL.TW.4","F06_LL.AU.2","F07_LL.AU.3","F08_LL.NZ.3","F09.REC.AU.3","F10_REC.NZ.3","F11_LL.ALL.1","F12_LL.ALL.2","F13_LL.ALL.3","F14_LL.ALL.4","F16_LL.JP.2_early","S01_INDEX.1-4")
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

    # CPUEinfo
    tmp_cpueinfo = data.frame(fleet=1:tmp_data$Nfleets,
                              units=rep(0,tmp_data$Nfleets),
                              errtype=rep(0,tmp_data$Nfleets),
                              SD_report=rep(0,tmp_data$Nfleets))
    rownames(tmp_cpueinfo) = tmp_fleetname
    tmp_data$CPUEinfo = tmp_cpueinfo
    tmp_data$CPUE$index = 16

    # catch
        tmp_catch = as.data.table(tmp_data$catch) %>%
                    .[year<1979&fleet==2,fleet:=15] %>%
                    .[fleet==15&catch>0&seas%in%3:4,catch_se:=0.2] %>%
                    .[fleet==15&catch>0&seas%in%1:2,catch_se:=0.05] %>%
                    .[order(fleet,year,seas)]

        tmp_data$catch = as.data.frame(tmp_catch)
    
    # define length info
    tmp_leninfo = data.table(mintailcomp=rep(1e-02,nrow(tmp_data$fleetinfo))) %>%
                  .[,addtocomp:=1e-04] %>%
                  .[,combine_M_F:=0] %>%
                  .[,CompressBins:=0] %>%
                  .[,CompError:=0] %>%
                  .[,ParmSelect:=0] %>%
                  .[,minsamplesize:=0.001] %>%
                  as.data.frame(.)
    rownames(tmp_leninfo) = tmp_fleetname
    tmp_data$len_info = tmp_leninfo

    tmp_lencomp = as.data.table(tmp_data$lencomp) %>%
                  .[fleet == 15, fleet:=16]
    tmp_data$lencomp = as.data.frame(tmp_lencomp)

    tmp_wtcomp = as.data.table(tmp_data$sizefreq_data_list[[1]]) %>%
                  .[fleet == 15, fleet:=16]
    tmp_data$sizefreq_data_list = list(as.data.frame(tmp_wtcomp))

    # misc
    tmp_data$comp_tail_compression = rep(-0.0001,tmp_data$Nfleets)
    tmp_data$add_to_comp = rep(0.0001,tmp_data$Nfleets)
    tmp_data$max_combined_lbin = rep(0,tmp_data$Nfleets) 


#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; ctl
# F_Method: 4
# add F_4_Fleet_Parms (data.frame): fleet, start_F, first_parm_phase
# F_iter: 4
# mirror selex for new fleet 16
    tmp_ctl$Nfleets = tmp_data$Nfleets
    tmp_ctl$fleetnames = tmp_data$fleetnames

    # fishing mortality
    tmp_ctl$F_Method = 4
    tmp_ctl$F_iter = 4
    tmp_ctl$maxF = 20
    tmp_ctl$F_4_Fleet_Parms = data.frame(fleet=c(1:15),start_F=rep(0.05,15),first_parm_phase=c(rep(99,14),5))

    # catchability options; define for each survey
    tmp_ctl$Q_options = data.frame(fleet=16,link=1,link_info=0,extra_se=1,biasadj=0,float=1)
    rownames(tmp_ctl$Q_options) = "S01_INDEX.1-4"

    tmp_ctl$Q_parms = tmp_ctl$Q_parms[1:2,]
    rownames(tmp_ctl$Q_parms) = c("LnQ","extra se")

    # size based selectivity
    tmp_size_selex_types = as.data.frame(matrix(0,nrow=length(tmp_ctl$fleetnames),ncol=4))
    rownames(tmp_size_selex_types) = tmp_ctl$fleetnames
    colnames(tmp_size_selex_types) = colnames(tmp_ctl$size_selex_types)
    tmp_size_selex_types$Pattern = 24
    tmp_size_selex_types$Pattern[c(5,8,11,12,13,14,15)] = 5
    tmp_size_selex_types$Pattern[10] = 1 # nz rec fishery
    tmp_size_selex_types$Special[c(5,8,11,12,13,14,15)] = c(4,7,1,2,3,4,2)
    tmp_ctl$size_selex_types = tmp_size_selex_types
    tmp_size_selex_parms = tmp_ctl$size_selex_parms

    pointer1 = max(grep("_F14_",rownames(tmp_size_selex_parms),fixed=TRUE))
    pointer2 = min(grep("_S01_",rownames(tmp_size_selex_parms),fixed=TRUE))
    tmp_ctl$size_selex_parms = rbind(tmp_size_selex_parms[1:pointer1,],make_size_selex_par_5(),tmp_size_selex_parms[pointer2:nrow(tmp_size_selex_parms),])
    
    # age based selectivity
    tmp_age_selex_types = as.data.frame(matrix(0,nrow=length(tmp_ctl$fleetnames),ncol=4))
    rownames(tmp_age_selex_types) = tmp_ctl$fleetnames
    colnames(tmp_age_selex_types) = colnames(tmp_ctl$age_selex_types)
    tmp_age_selex_types$Pattern = 11
    tmp_ctl$age_selex_types = tmp_age_selex_types
    tmp_ctl$age_selex_parms = as.data.frame(rbindlist(lapply(rep(tmp_ctl$Nages,length(tmp_ctl$fleetnames)),make_age_selex_par_11)))

    # variance adjustment list
    tmp_ctl$Variance_adjustment_list = expand.grid(Factor=1:7,Fleet=1:length(tmp_ctl$fleetnames),Value=0)
    tmp_ctl$Variance_adjustment_list$Value[which(tmp_ctl$Variance_adjustment_list$Factor%in%c(4,7))] = 1

    # adjust ess for comps
    length_var_adj = 1/c(10,12,10,10,10,18,10,10,10,10,16,34,10,10,15,10)
    weight_var_adj = 1/c(10,10,10,10,10,10,10,10,12,23,10,10,10,10,15,10)

    tmp_ctl$Variance_adjustment_list$Value[which(tmp_ctl$Variance_adjustment_list$Factor%in%c(4))] = length_var_adj
    tmp_ctl$Variance_adjustment_list$Value[which(tmp_ctl$Variance_adjustment_list$Factor%in%c(7))] = weight_var_adj


    # update lambdas
    # lambdas
    tmp_lambdas_surv= data.table(like_comp=rep(1,length(tmp_ctl$fleetnames)),fleet=1:length(tmp_ctl$fleetnames)) %>%
                      .[,phase:=1] %>%
                      .[,value:=0] %>%
                      .[,sizefreq_method:=1] %>%
                      .[fleet %in% c(16),value:=1]
    tmp_ctl$lambdas = as.data.frame(tmp_lambdas_surv)
    tmp_ctl$N_lambdas = nrow(tmp_ctl$lambdas)

    # relax growth
    tmp_ctl$CV_Growth_Pattern = 0 # make CV a function of length at age; 
    tmp_ctl$Growth_Age_for_L1 = 0
    tmp_MG_parms = tmp_ctl$MG_parms
    tmp_MG_parms$INIT[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = 59.9

    tmp_MG_parms$PRIOR[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] =  59.9

    tmp_MG_parms$INIT[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.15

    tmp_MG_parms$PRIOR[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.15

    tmp_MG_parms$HI[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = 20

    tmp_MG_parms$INIT[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.15

    tmp_MG_parms$PRIOR[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.15

    tmp_MG_parms$HI[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = 20
    tmp_ctl$MG_parms = tmp_MG_parms



#________________________________________________________________________________________________________________________________________________________________________________________________________
# write-out files
    # SS_writestarter(tmp_starter,dir=to_dir,overwrite=TRUE)
    # SS_writeforecast(tmp_forecast,dir=to_dir,overwrite=TRUE)
    SS_writedat(tmp_data, outfile=file.path(to_dir,"data.ss"), overwrite=TRUE)
    SS_writectl(tmp_ctl,file.path(to_dir,"control.ss"),overwrite = TRUE)


#________________________________________________________________________________________________________________________________________________________________________________________________________
# run & summarize
    # run
    if(linux_ws){
        run(dir=to_dir,exe="ss3",show_in_console=TRUE,skipfinished=FALSE)
    } else {
        run(dir=to_dir,exe="ss3.exe",show_in_console=TRUE,skipfinished=FALSE)
    }

    # summarize
    summarize_ss_model(to_dir)

    output = SS_output(to_dir)
    SS_plots(output,dir=to_dir)

    # clean directory
    # clean_dir(to_dir)
    # make_ss_output(to_dir)
    
    # make html viewer
    # ss_model_html_viewer(to_dir)


