

# Nicholas Ducharme-Barth
# 2025/01/24
# set up 2 sex model
# also make NZ rec data quarterly


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
    dir_base_stock_synthesis = file.path(dir_model,"stock-synthesis","12-CAAL-old-growth-SD")
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
    to_dir = file.path(dir_model,"stock-synthesis","23-CAAL-2sex")
    dir.create(to_dir,recursive=TRUE)

    ss_transfer(from_dir=from_dir,to_dir=to_dir,executable_dir_stem=executable_dir_stem,linux_ws=is_linux_os())

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in files
    # tmp_starter = SS_readstarter(file=file.path(to_dir,"starter.ss"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=file.path(to_dir,"control.ss"),datlist = file.path(to_dir,"data.ss"))
    tmp_data = SS_readdat(file=file.path(to_dir,"data.ss"))
    # tmp_forecast = SS_readforecast(file=file.path(to_dir,"forecast.ss"),verbose=FALSE)

    # tmp_ctl = SS_readctl(file=file.path(dir_model,"stock-synthesis","00-npo-mako-base-file","control.ss_new"),datlist = file.path(dir_model,"stock-synthesis","00-npo-mako-base-file","data_echo.ss_new"))
    # tmp_data = SS_readdat(file=file.path(dir_model,"stock-synthesis","00-npo-mako-base-file","data_echo.ss_new"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; data
    tmp_data$Nsexes = 2

# duplicate length comp data for males with 0 (sex = 0; aggregated by sex and put with females)
    tmp_lencomp = tmp_data$lencomp
    tmp_lencomp_info = tmp_lencomp[,1:6]
    tmp_lencomp_female = tmp_lencomp[,7:ncol(tmp_lencomp)]
    colnames(tmp_lencomp_female) = gsub("l","f",colnames(tmp_lencomp_female))
    tmp_lencomp_male = as.data.frame(matrix(0,nrow=nrow(tmp_lencomp_female),ncol=ncol(tmp_lencomp_female)))
    colnames(tmp_lencomp_male) = gsub("f","m",colnames(tmp_lencomp_female))

    tmp_data$lencomp = cbind(tmp_lencomp_info,tmp_lencomp_female,tmp_lencomp_male)

# duplicate wt comp
    tmp_wtcomp = tmp_data$sizefreq_data_list[[1]]
    tmp_wtcomp_info = tmp_wtcomp[,1:7]
    tmp_wtcomp_female = tmp_wtcomp[,8:ncol(tmp_wtcomp)]
    colnames(tmp_wtcomp_female) = gsub("a","f",colnames(tmp_wtcomp_female))
    tmp_wtcomp_male = as.data.frame(matrix(0,nrow=nrow(tmp_wtcomp_female),ncol=ncol(tmp_wtcomp_female)))
    colnames(tmp_wtcomp_male) = gsub("f","m",colnames(tmp_wtcomp_female))

    tmp_data$sizefreq_data_list[[1]] = cbind(tmp_wtcomp_info,tmp_wtcomp_female,tmp_wtcomp_male)

# make CAAL sex-specific
    pop_len_bins = tmp_data$lbin_vector_pop
    age_dt = fread(file.path(proj_dir,"data","MLS age data for SPC_04112023.csv")) %>%
            setnames(.,c("Catch date","LJFL (mm)","Decimal age (y)"),c("Date","LJFL","age")) %>%
            .[,year:=sapply(Date,function(x)as.numeric(strsplit(x,"/")[[1]][3]))] %>%
            .[,month:=sapply(Date,function(x)as.numeric(strsplit(x,"/")[[1]][2]))] %>%
            .[COUNTRY=="New Zealand",fleet:=10] %>%
            .[COUNTRY%in%c("New Caledonia","Fiji"),fleet:=12] %>%
            .[COUNTRY=="Australia"&LAT_DEG_MIN<30,fleet:=9] %>%
            .[COUNTRY=="Australia"&LAT_DEG_MIN>=30,fleet:=6] %>%
            .[SEX=="F",sex:=1] %>%
            .[SEX=="M",sex:=2] %>%
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
            .[order(fleet,year,month,sex,Lbin_lo)] %>%
            na.omit(.)

    nsamp_caal = rowSums(as.data.frame(age_dt)[,-c(1:8)])
    caal_data = as.data.frame(age_dt)[,-c(1:8)]
    caal_id_data = as.data.frame(age_dt)[,c(1:8)]
    caal_id_data$Nsamp = nsamp_caal
    
    caal_m = matrix(0,nrow=nrow(caal_data),ncol=ncol(caal_data))
    colnames(caal_m) = colnames(caal_data)
    for(i in 1:nrow(caal_id_data)){
        if(caal_id_data$sex[i]==2){
            caal_m[i,] = as.vector(as.matrix(caal_data[i,]))
            caal_data[i,] = rep(0,ncol(caal_data))
        }
    }
    caal = cbind(caal_id_data,caal_data,caal_m)
    tmp_data$agecomp = caal
    

# add nz rec wt comp in by quarter
# modify for 2 sex
      # bring in NZ sport comps
        nz_wt = fread(paste0(proj_dir,"/data/nz-all-sport-club-weights.csv"))
        colnames(nz_wt) = c("date","club","species","tagged","weight","boat","locality","year")
        nz_wt = nz_wt %>%
                .[!is.na(weight)&!is.na(year)&tagged=="",.(date,year,weight)] %>%
                .[year>=1952&year<1988] %>%
                .[,dd:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][2]))] %>%
                .[,mm:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][1]))] %>%
                .[,yy:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][3]))] %>%
                .[yy>51&yy<100,yy:=yy+1900] %>%
                .[yy<25,yy:=yy+2000] %>%
                .[,year:=yy] %>%
                .[,month:=c(2,2,2,5,5,5,8,8,8,11,11,11)[mm]] %>%
                .[,.(year,month,weight)] %>%
                .[,ts:=paste0(year,"-",month)]


        wtfrq_bins = tmp_data$sizefreq_bins_list[[1]]

        nz_rec_wt_mat_m = nz_rec_wt_mat = matrix(0,nrow=uniqueN(nz_wt$ts),ncol=length(wtfrq_bins))
        colnames(nz_rec_wt_mat) = wtfrq_bins
        rownames(nz_rec_wt_mat) = sort(unique(nz_wt$ts))
        nz_nsamp_vec = rep(NA,nrow(nz_rec_wt_mat))

        for(i in 1:nrow(nz_rec_wt_mat)){
            for(j in 1:ncol(nz_rec_wt_mat)){
                if(j<ncol(nz_rec_wt_mat)){
                    tmp_nz = nz_wt[ts==(rownames(nz_rec_wt_mat)[i])&weight>=as.numeric(colnames(nz_rec_wt_mat)[j])&weight<as.numeric(colnames(nz_rec_wt_mat)[j+1])]
                } else {
                    tmp_nz = nz_wt[ts==(rownames(nz_rec_wt_mat)[i])&weight>=as.numeric(colnames(nz_rec_wt_mat)[j])]
                }

                nz_rec_wt_mat[i,j] = nrow(tmp_nz)

                if(j==ncol(nz_rec_wt_mat)){
                    nz_nsamp_vec[i] = sum(nz_rec_wt_mat[i,])
                }
                
                rm(list=c("tmp_nz"))
            }
        }
        colnames(nz_rec_wt_mat) = paste0("f",colnames(nz_rec_wt_mat))
        colnames(nz_rec_wt_mat_m) = gsub("f","m",colnames(nz_rec_wt_mat))

        nz_rec_info_mat = matrix(0,nrow=uniqueN(nz_wt$ts),ncol=7)
        colnames(nz_rec_info_mat) = c("method","year","month","fleet","sex","part","Nsamp")

        nz_rec_info_mat[,"method"] = 1
        nz_rec_info_mat[,"year"] = sapply(rownames(nz_rec_wt_mat),function(x)as.numeric(strsplit(x,"-")[[1]][1]))
        nz_rec_info_mat[,"month"] = sapply(rownames(nz_rec_wt_mat),function(x)as.numeric(strsplit(x,"-")[[1]][2]))
        nz_rec_info_mat[,"fleet"] = 10
        nz_rec_info_mat[,"sex"] = 0
        nz_rec_info_mat[,"part"] = 0
        nz_rec_info_mat[,"Nsamp"] = nz_nsamp_vec

        tmp_wtcomp_1952 = cbind(nz_rec_info_mat,nz_rec_wt_mat,nz_rec_wt_mat_m) %>%
                          as.data.table(.)

        tmp_wtcomp_1979 = as.data.table(tmp_data$sizefreq_data_list[[1]]) %>%
                          .[fleet!=10]

        tmp_wtcomp = rbind(tmp_wtcomp_1952,tmp_wtcomp_1979) %>%
                     .[order(method,fleet,year,month)] %>%
                     as.data.frame(.)

        tmp_data$sizefreq_data_list[[1]] = tmp_wtcomp
        tmp_data$Nobs_per_method = nrow(tmp_wtcomp)


#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify: ctl
tmp_ctl$Nsexes = 2

# modify natural mortality
# natM (data.frame, nrows= sex, ncol=age)
tmp_ctl$natM = rbind(tmp_ctl$natM,tmp_ctl$natM)

# duplicate growth & wt-len parms for males
# place after female reproductive params but before recr dist
tmp_MG_parms = tmp_ctl$MG_parms
pointer1 = grep("L_at_Amin_Fem_",rownames(tmp_MG_parms),fixed=TRUE)
pointer2 = grep("Wtlen_2_Fem_",rownames(tmp_MG_parms),fixed=TRUE)
pointer3 = grep("Eggs_beta_Fem_",rownames(tmp_MG_parms),fixed=TRUE)

tmp_ctl$MG_parms = rbind(tmp_MG_parms[pointer1:pointer3,],tmp_MG_parms[pointer1:pointer2,],tmp_MG_parms[(pointer3+1):nrow(tmp_MG_parms),])

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


