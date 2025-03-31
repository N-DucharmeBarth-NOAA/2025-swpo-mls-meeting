

# Nicholas Ducharme-Barth
# 2025/01/22
# Run CAAL model


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
    dir_base_stock_synthesis = file.path(dir_model,"stock-synthesis","07-catch-uncertainty")
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
    to_dir = file.path(dir_model,"stock-synthesis","10-CAAL")
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


#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; ctl
# turn on estimation of L1, L2, K
    tmp_ctl$CV_Growth_Pattern = 0 # make CV a function of length at age; 
    tmp_ctl$Growth_Age_for_L1 = 0
    tmp_MG_parms = tmp_ctl$MG_parms
    tmp_MG_parms$INIT[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = 59.9
    tmp_MG_parms$PRIOR[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] =  59.9

    tmp_MG_parms$INIT[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.2
    tmp_MG_parms$PRIOR[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.2
    tmp_MG_parms$HI[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = 20

    tmp_MG_parms$INIT[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.05
    tmp_MG_parms$PRIOR[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.05
    tmp_MG_parms$HI[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = 20
    
    tmp_MG_parms$PHASE[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = 7
    tmp_MG_parms$PHASE[grep("L_at_Amax_",rownames(tmp_MG_parms),fixed=TRUE)] = 6
    tmp_MG_parms$PHASE[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = 8

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


