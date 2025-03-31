

# Nicholas Ducharme-Barth
# 2025/01/23
# put in NZ wt comp by quarter


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
    to_dir = file.path(dir_model,"stock-synthesis","21-CAAL-NZrecwtQtr")
    dir.create(to_dir,recursive=TRUE)

    ss_transfer(from_dir=from_dir,to_dir=to_dir,executable_dir_stem=executable_dir_stem,linux_ws=is_linux_os())

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in files
    # tmp_starter = SS_readstarter(file=file.path(to_dir,"starter.ss"),verbose=FALSE)
    # tmp_ctl = SS_readctl(file=file.path(to_dir,"control.ss"),datlist = file.path(to_dir,"data.ss"))
    tmp_data = SS_readdat(file=file.path(to_dir,"data.ss"))
    # tmp_forecast = SS_readforecast(file=file.path(to_dir,"forecast.ss"),verbose=FALSE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; data
# add nz rec wt comp in by quarter
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

        nz_rec_wt_mat = matrix(0,nrow=uniqueN(nz_wt$ts),ncol=length(wtfrq_bins))
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
        colnames(nz_rec_wt_mat) = paste0("a",colnames(nz_rec_wt_mat))

        nz_rec_info_mat = matrix(0,nrow=uniqueN(nz_wt$ts),ncol=7)
        colnames(nz_rec_info_mat) = c("method","year","month","fleet","sex","part","Nsamp")

        nz_rec_info_mat[,"method"] = 1
        nz_rec_info_mat[,"year"] = sapply(rownames(nz_rec_wt_mat),function(x)as.numeric(strsplit(x,"-")[[1]][1]))
        nz_rec_info_mat[,"month"] = sapply(rownames(nz_rec_wt_mat),function(x)as.numeric(strsplit(x,"-")[[1]][2]))
        nz_rec_info_mat[,"fleet"] = 10
        nz_rec_info_mat[,"sex"] = 0
        nz_rec_info_mat[,"part"] = 0
        nz_rec_info_mat[,"Nsamp"] = nz_nsamp_vec

        tmp_wtcomp_1952 = cbind(nz_rec_info_mat,nz_rec_wt_mat) %>%
                          as.data.table(.)

        tmp_wtcomp_1979 = as.data.table(tmp_data$sizefreq_data_list[[1]]) %>%
                          .[fleet!=10]

        tmp_wtcomp = rbind(tmp_wtcomp_1952,tmp_wtcomp_1979) %>%
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
#________________________________________________________________________________________________________________________________________________________________________________________________________
# write-out files
    # SS_writestarter(tmp_starter,dir=to_dir,overwrite=TRUE)
    # SS_writeforecast(tmp_forecast,dir=to_dir,overwrite=TRUE)
    SS_writedat(tmp_data, outfile=file.path(to_dir,"data.ss"), overwrite=TRUE)
    # SS_writectl(tmp_ctl,file.path(to_dir,"control.ss"),overwrite = TRUE)


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


