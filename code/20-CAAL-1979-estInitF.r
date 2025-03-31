

# Nicholas Ducharme-Barth
# 2025/01/21
# Run model excluding bad comp data
    # - drop comp for fisheries 11 (2006-2014); length only no weight
    # - drop comp 12 (mirror to other fleets); length only no weight
    # - remove early fishery 14 composition data (through 2005); length only no weight
    # - drop 2006 - 2014 size comp data for fishery 15; length only (can try dropping all length data later)
    # - remove AU length data (fisheries 6 & 7)
    # - remove NZ wt comp from 1988 on  
    # - remove small sample sizes fleet 1 (length/weight), fleet 8 (length), fleet 13 (length) 

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
    dir_base_stock_synthesis = file.path(dir_model,"stock-synthesis","03-chg-selex-1979")
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
    to_dir = file.path(dir_model,"stock-synthesis","20-CAAL-1979-estInitF")
    dir.create(to_dir,recursive=TRUE)

    ss_transfer(from_dir=from_dir,to_dir=to_dir,executable_dir_stem=executable_dir_stem,linux_ws=is_linux_os())

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in files
    # tmp_starter = SS_readstarter(file=file.path(to_dir,"starter.ss"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=file.path(to_dir,"control.ss"),datlist = file.path(to_dir,"data.ss"))
    tmp_data = SS_readdat(file=file.path(to_dir,"data.ss"))
    # tmp_forecast = SS_readforecast(file=file.path(to_dir,"forecast.ss"),verbose=FALSE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# bring in age data and format it
# convert to eye orbital fork length in cm
    pop_len_bins = tmp_data$lbin_vector_pop
    age_dt = fread(file.path(proj_dir,"data","MLS age data for SPC_04112023.csv")) %>%
            setnames(.,c("Catch date","LJFL (mm)","Decimal age (y)"),c("Date","LJFL","age")) %>%
            .[,year:=sapply(Date,function(x)as.numeric(strsplit(x,"/")[[1]][3]))] %>%
            .[,month:=sapply(Date,function(x)as.numeric(strsplit(x,"/")[[1]][2]))] %>%
            .[COUNTRY=="New Zealand",fleet:=10] %>%
            .[COUNTRY%in%c("New Caledonia","Fiji"),fleet:=12] %>%
            .[COUNTRY=="Australia"&LAT_DEG_MIN<30,fleet:=9] %>%
            .[COUNTRY=="Australia"&LAT_DEG_MIN>=30,fleet:=6] %>%
            .[,sex:=0] %>%
            .[,part:=0] %>%
            .[,ageerr:=1] %>%
            .[,eofl:=0.862069*LJFL/10] %>%
            .[,.(year,month,fleet,sex,part,ageerr,age,eofl)] %>%
            .[,caal_age_bin:=floor(age)] %>%
            .[caal_age_bin>=tmp_data$Nages,caal_age_bin:=tmp_data$Nages] %>%
            .[,Lbin_lo:=cut(eofl,breaks=pop_len_bins,include.lowest=TRUE,right=FALSE,labels=FALSE)] %>%
            .[,Lbin_hi:=Lbin_lo] %>%
            .[,.(year,month,fleet,sex,part,ageerr,caal_age_bin,Lbin_lo,Lbin_hi)] %>%
            setnames(.,"caal_age_bin","age") %>%
            .[,age:=paste0("age",age)] %>%
            .[,age:=factor(age,levels=paste0("age",0:tmp_data$Nages))] %>%
            dcast(.,year+month+fleet+sex+part+ageerr+Lbin_lo+Lbin_hi~age,fun.aggregate=length,drop=c(TRUE,FALSE)) %>%
            .[order(fleet,year,month,Lbin_lo)]

    nsamp_caal = rowSums(as.data.frame(age_dt)[,-c(1:8)])
    caal_data = as.data.frame(age_dt)[,-c(1:8)]
    caal_id_data = as.data.frame(age_dt)[,c(1:8)]
    caal_id_data$Nsamp = nsamp_caal
    caal = cbind(caal_id_data,caal_data)
    

#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; data
# remove comps
    # - drop comp for fisheries 11 (2006-2014); length only no weight
    # - drop comp 12 (mirror to other fleets); length only no weight
    # - remove early fishery 14 composition data (through 2005); length only no weight
    # - drop 2006 - 2014 size comp data for fishery 15; length only (can try dropping all length data later)
    # - remove AU length data (fisheries 6 & 7)
    # - remove small sample sizes fleet 1 (length), fleet 8 (length), fleet 13 (length) 
    tmp_lencomp = as.data.table(tmp_data$lencomp) %>%
                  .[fleet==11 & year%in%2006:2014,year:=abs(year)*-1] %>%
                  .[fleet %in% c(1,6,7,8,12,13), year:=abs(year)*-1] %>%
                  .[fleet==14 & year<=2005, year:=abs(year)*-1] %>%
                  .[fleet==15 & year%in%2006:2014, year:=abs(year)*-1]
    tmp_data$lencomp = as.data.frame(tmp_lencomp)

        # bring in NZ sport comps
        nz_wt = fread(paste0(proj_dir,"/data/nz-all-sport-club-weights.csv"))
        colnames(nz_wt) = c("date","club","species","tagged","weight","boat","locality","year")
        nz_wt = nz_wt %>%
                .[!is.na(weight)&!is.na(year)&tagged=="",.(year,weight)] %>%
                .[year>=1952]

        wtfrq_bins = tmp_data$sizefreq_bins_list[[1]]

        nz_rec_wt_mat = matrix(0,nrow=length(min(nz_wt$year):max(nz_wt$year)),ncol=length(wtfrq_bins))
        colnames(nz_rec_wt_mat) = wtfrq_bins
        rownames(nz_rec_wt_mat) = min(nz_wt$year):max(nz_wt$year)
        nz_nsamp_vec = rep(NA,nrow(nz_rec_wt_mat))

        for(i in 1:nrow(nz_rec_wt_mat)){
            for(j in 1:ncol(nz_rec_wt_mat)){
                if(j<ncol(nz_rec_wt_mat)){
                    tmp_nz = nz_wt[year==as.numeric(rownames(nz_rec_wt_mat)[i])&weight>=as.numeric(colnames(nz_rec_wt_mat)[j])&weight<as.numeric(colnames(nz_rec_wt_mat)[j+1])]
                } else {
                    tmp_nz = nz_wt[year==as.numeric(rownames(nz_rec_wt_mat)[i])&weight>=as.numeric(colnames(nz_rec_wt_mat)[j])]
                }

                nz_rec_wt_mat[i,j] = nrow(tmp_nz)

                if(j==ncol(nz_rec_wt_mat)){
                    nz_nsamp_vec[i] = sum(nz_rec_wt_mat[i,])
                }
                
                rm(list=c("tmp_nz"))
            }
        }
        colnames(nz_rec_wt_mat) = paste0("a",colnames(nz_rec_wt_mat))

        nz_rec_info_mat = matrix(0,nrow=length(min(nz_wt$year):max(nz_wt$year)),ncol=7)
        colnames(nz_rec_info_mat) = c("method","year","month","fleet","sex","part","Nsamp")

        nz_rec_info_mat[,"method"] = 1
        nz_rec_info_mat[,"year"] = min(nz_wt$year):max(nz_wt$year)
        nz_rec_info_mat[,"month"] = 2
        nz_rec_info_mat[,"fleet"] = 10
        nz_rec_info_mat[,"sex"] = 0
        nz_rec_info_mat[,"part"] = 0
        nz_rec_info_mat[,"Nsamp"] = nz_nsamp_vec

        tmp_wtcomp_1952 = cbind(nz_rec_info_mat,nz_rec_wt_mat) %>%
                          as.data.table(.)

        tmp_wtcomp_1979 = as.data.table(tmp_data$sizefreq_data_list[[1]]) %>%
                          .[fleet!=10]

        tmp_wtcomp = rbind(tmp_wtcomp_1952,tmp_wtcomp_1979) %>%
                     .[year>=1979] %>%
                     .[order(method,fleet,year,month)] %>%
                     as.data.frame(.)

        tmp_data$sizefreq_data_list[[1]] = tmp_wtcomp
        tmp_data$Nobs_per_method = nrow(tmp_wtcomp)

    # - remove NZ wt comp from 1988 on  
    # - remove small sample sizes fleet 1 (weight)
    tmp_wtcomp = as.data.table(tmp_data$sizefreq_data_list[[1]]) %>%
                  .[fleet==10 & year>=1988,year:=abs(year)*-1] %>%
                  .[fleet == 1, year:=abs(year)*-1]
    tmp_data$sizefreq_data_list = list(as.data.frame(tmp_wtcomp))

    # remove comps
    # - drop index (15) length comp
    # - drop fishery 2 length comp
    tmp_lencomp = as.data.table(tmp_data$lencomp) %>%
                  .[fleet %in% c(2,15), year:=abs(year)*-1]
    tmp_data$lencomp = as.data.frame(tmp_lencomp)

    # - drop fishery 2 weight comp for 1979
    tmp_wtcomp = as.data.table(tmp_data$sizefreq_data_list[[1]]) %>%
                  .[fleet==2 & year==1979,year:=abs(year)*-1] 
    tmp_data$sizefreq_data_list = list(as.data.frame(tmp_wtcomp))

# add CAAL & ageing error matrix
# integer
tmp_data$N_agebins = 11
# vector of ages
tmp_data$agebin_vector = 0:tmp_data$Nages
# integer
tmp_data$N_ageerror_definitions = 1
# data frame (nrow = 2*N_ageerror_definitions, ncol = N_agebins (mean,sd) for each age)
# colnames 'age0', 'age1', ...
ageerror_mat = matrix(0,nrow=2*tmp_data$N_ageerror_definitions,ncol=tmp_data$N_agebins)
colnames(ageerror_mat) = paste0("age",0:tmp_data$Nages)
ageerror_mat[1,] = (0:tmp_data$Nages) + 0.5
ageerror_mat[2,] = round(seq(from=0.35,to=2,length.out=ncol(ageerror_mat)),digits=2)
tmp_data$ageerror = as.data.frame(ageerror_mat)
# data.frame (nrow=fishery,ncol="mintailcomp","addtocomp","combine_M_F","CompressBins","CompError","ParmSelect","minsamplesize")
tmp_ageinfo = data.table(mintailcomp=rep(-1e-04,nrow(tmp_data$fleetinfo))) %>%
                  .[,addtocomp:=1e-04] %>%
                  .[,combine_M_F:=0] %>%
                  .[,CompressBins:=0] %>%
                  .[,CompError:=0] %>%
                  .[,ParmSelect:=0] %>%
                  .[,minsamplesize:=0.001] %>%
                  as.data.frame(.)
tmp_data$age_info = tmp_ageinfo
# data.frame c("year","month","fleet","sex","part","ageerr","Lbin_lo","Lbin_hi","Nsamp")
tmp_data$agecomp = caal
# integer
tmp_data$Lbin_method = 1

# remove equilibrium catch for non-fleet 2
tmp_catch = as.data.table(tmp_data$catch) %>%
            .[!(fleet!=2&year<0)]
tmp_data$catch = as.data.frame(tmp_catch)


#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; control
     # catchability options; define for each survey
    tmp_ctl$Q_options = data.frame(fleet=15,link=1,link_info=0,extra_se=1,biasadj=0,float=1)
    rownames(tmp_ctl$Q_options) = "S01_INDEX.1-4"

    # extra obs for cpue
    loess_y = tmp_data$CPUE$obs
    loess_x = tmp_data$CPUE$year

    loess_cpue = loess(loess_y~loess_x)
    
    tmp_ctl$Q_parms = rbind(tmp_ctl$Q_parms[1,],tmp_ctl$Q_parms[1,])
    rownames(tmp_ctl$Q_parms) = c("LnQ","extra se")
    tmp_ctl$Q_parms[2,1:3] = c(0,1,mean(abs(loess_cpue$residuals)))

 # size based selectivity
    tmp_size_selex_types = as.data.frame(matrix(0,nrow=length(tmp_ctl$fleetnames),ncol=4))
    rownames(tmp_size_selex_types) = tmp_ctl$fleetnames
    colnames(tmp_size_selex_types) = colnames(tmp_ctl$size_selex_types)
    tmp_size_selex_types$Pattern = 24
    tmp_size_selex_types$Pattern[c(5,8,11,12,13,14)] = 5
    tmp_size_selex_types$Pattern[10] = 1 # nz rec fishery
    tmp_size_selex_types$Special[c(5,8,11,12,13,14)] = c(4,7,1,2,3,4)
    tmp_ctl$size_selex_types = tmp_size_selex_types

    # define selex shapes
    # LO, HI, INIT, PHASE 
    tmp_size_selex_parms = tmp_ctl$size_selex_parms
    pointer1 = max(grep("_F11_",rownames(tmp_size_selex_parms),fixed=TRUE))
    pointer2 = min(grep("_F13_",rownames(tmp_size_selex_parms),fixed=TRUE))
    tmp_ctl$size_selex_parms = rbind(tmp_size_selex_parms[1:pointer1,],make_size_selex_par_5(),tmp_size_selex_parms[pointer2:nrow(tmp_size_selex_parms),])
    
    # update lambdas
    # lambdas
    tmp_lambdas_surv= data.table(like_comp=rep(1,length(tmp_ctl$fleetnames)),fleet=1:length(tmp_ctl$fleetnames)) %>%
                      .[,phase:=1] %>%
                      .[,value:=0] %>%
                      .[,sizefreq_method:=1] %>%
                      .[fleet %in% c(15),value:=1]
    tmp_ctl$lambdas = as.data.frame(tmp_lambdas_surv)
    tmp_ctl$N_lambdas = nrow(tmp_ctl$lambdas)

    # turn on estimation of L1, L2, K
    tmp_ctl$Growth_Age_for_L1 = 0
    tmp_MG_parms = tmp_ctl$MG_parms
    tmp_MG_parms$INIT[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = 59.9
    tmp_MG_parms$PRIOR[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] =  59.9
    
    tmp_MG_parms$PHASE[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = 7
    tmp_MG_parms$PHASE[grep("L_at_Amax_",rownames(tmp_MG_parms),fixed=TRUE)] = 6
    tmp_MG_parms$PHASE[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = 8

    tmp_ctl$MG_parms = tmp_MG_parms
    tmp_initF = matrix(0,nrow=sum(tmp_data$catch$year==-999),ncol=7)
    colnames(tmp_initF) = c("LO","HI","INIT","PRIOR","PRIOR_SD","PR_type","PHASE")
    tmp_initF[,"LO"] = 0
    tmp_initF[,"HI"] = 5
    tmp_initF[,"INIT"] = (tail(as.matrix(tmp_ctl$natM)[1,],n=1)/4 * 1.7) - tail(as.matrix(tmp_ctl$natM)[1,],n=1)/4
    tmp_initF[,"PRIOR"] = 0.2
    tmp_initF[,"PRIOR_SD"] = 99
    tmp_initF[,"PR_type"] = 0
    tmp_initF[,"PHASE"] = 5 # estimate

    tmp_ctl$init_F = as.data.frame(tmp_initF)
    

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


