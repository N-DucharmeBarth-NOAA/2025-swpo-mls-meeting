

# Nicholas Ducharme-Barth
# 2025/01/10
# Extend model back to 1952
# add catch
# add NZ rec weight comp
# also change parameterization of extra index obs error

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
	dir_model = paste0(proj_dir,"/models/")
    dir_base_mfcl = paste0(dir_model,"mfcl/2019-diagnostic/")
    dir_base_stock_synthesis = paste0(dir_model,"stock-synthesis/03-chg-selex-1979/")
    dir_helper_fns = paste0(proj_dir,"/code/helper-fns/")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(paste0(dir_helper_fns,(list.files(dir_helper_fns))),source)

#_____________________________________________________________________________________________________________________________
# read in baseline mfcl files
    base_frq = parse_frq(paste0(dir_base_mfcl,"cor.au.ll.len.frq"))

#_____________________________________________________________________________________________________________________________
# read in baseline stock synthesis files
    tmp_starter = SS_readstarter(file=paste0(dir_base_stock_synthesis,"starter.ss_new"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=paste0(dir_base_stock_synthesis,"control.ss_new"),datlist = paste0(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_data = SS_readdat(file=paste0(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_forecast = SS_readforecast(file=paste0(dir_base_stock_synthesis,"forecast.ss_new"),verbose=FALSE)

#_____________________________________________________________________________________________________________________________
# create new directory for stock synthesis mls files
    dir_mls_stock_synthesis_base = paste0(dir_model,"stock-synthesis/04-start-1952/")
    dir.create(dir_mls_stock_synthesis_base,recursive=TRUE)

    # write out starter and forecast
    SS_writestarter(tmp_starter,dir=dir_mls_stock_synthesis_base,overwrite=TRUE)
    SS_writeforecast(tmp_forecast,dir=dir_mls_stock_synthesis_base,overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update data
    tmp_data$styr = min(cateffpen(base_frq)$year,na.rm=TRUE)

    # define catch
    # fishery 2 will get the equilibrium catch since it has the greatest catch in 1979
    # equilibrium catch will be the catch averaged across all fisheries for the first 3 years
    # remember catch needs to be divided by 1000 (stock synthesis records catch in 1000s of fish)
    cateffpen_dt = as.data.table(cateffpen(base_frq))
    # eq catch for NZ rec (fishery 10) if needed...
    # eq_catch = mean(cateffpen_dt[year %in% c(tmp_data$styr + 0:2)&fishery==10,.(catch=sum(catch)/1000),by=year]$catch)
    
    tmp_catch.list = as.list(rep(NA,tmp_data$Nfleets-1))
    for(i in seq_along(tmp_catch.list)){
        tmp_catch.list[[i]] = cateffpen_dt[fishery==i,.(year,month,fishery,catch)] %>%
                         setnames(.,c("month","fishery"),c("seas","fleet")) %>%
                         .[,seas:=c(1,1,1,2,2,2,3,3,3,4,4,4)[seas]] %>%
                         .[,catch:=catch/1000] %>%
                         .[,catch_se:=0.01] %>%
                         .[,.(year,seas,fleet,catch,catch_se)] %>%
                         .[catch>0] %>%
                         .[year<1979]
    }
    tmp_catch_1952_dt = rbindlist(tmp_catch.list) 
    tmp_catch_1979_dt = as.data.table(tmp_data$catch) %>%
                       .[year>0]
    tmp_catch = rbind(tmp_catch_1952_dt,tmp_catch_1979_dt) %>%
                .[order(fleet,year,seas)] %>%
                as.data.frame(.)
    tmp_data$catch = tmp_catch

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
                     .[order(method,fleet,year,month)] %>%
                     as.data.frame(.)

        tmp_data$sizefreq_data_list[[1]] = tmp_wtcomp
        tmp_data$Nobs_per_method = nrow(tmp_wtcomp)

    SS_writedat(tmp_data, outfile=paste0(dir_mls_stock_synthesis_base,"data.ss"), overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update control
    # setup default blocks
    tmp_ctl$Block_Design[[1]] = rep(tmp_data$styr-1,2)

    # don't use steepness because starting at unfished
    tmp_ctl$Use_steep_init_equi = 0 

    # recruitment deviation setup
    tmp_ctl$MainRdevYrFirst = tmp_data$styr
    tmp_ctl$MainRdevYrLast = tmp_data$endyr - 1 # late dev for terminal year
    tmp_ctl$recdev_early_start = tmp_data$styr-10
    tmp_ctl$last_early_yr_nobias_adj = tmp_data$styr + 10 # initial values
    tmp_ctl$first_yr_fullbias_adj = tmp_data$endyr - 16 # initial values
    tmp_ctl$last_yr_fullbias_adj = tmp_data$endyr - 1 # initial values
    tmp_ctl$first_recent_yr_nobias_adj = tmp_data$endyr - 1
    tmp_ctl$max_bias_adj = 0.5

    # remove init F
    tmp_ctl$init_F = NULL

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

    SS_writectl(tmp_ctl,paste0(dir_mls_stock_synthesis_base,"control.ss"),overwrite = TRUE)

#_____________________________________________________________________________________________________________________________
# run new version of stock synthesis
    file.copy(from=paste0(proj_dir,"/executables/stock-synthesis/3.30.23.1/ss3_win.exe"),to=dir_mls_stock_synthesis_base)
    run(dir=dir_mls_stock_synthesis_base,exe="ss3_win.exe",show_in_console = TRUE)

    # output = SS_output(dir_mls_stock_synthesis_base)
    # SS_plots(output,dir=dir_mls_stock_synthesis_base)

    # summarize
    summarize_ss_model(dir_mls_stock_synthesis_base)

    # clean directory
    all_files = list.files(to_dir)
    clean_files = all_files[which(!(all_files %in% c(all_files[grep("html_",all_files,fixed=TRUE)],all_files[grep(".csv",all_files,fixed=TRUE)],"ss.par","control.ss","starter.ss_new","forecast.ss_new","control.ss_new","data_echo.ss_new","data.ss_new","ss3.par","warning.sso","executable.txt")))]
    print(paste0(length(clean_files)," files removed taking up ",round(sum(unname(sapply(file.path(to_dir,clean_files),file.size)))/1024^2,digits=2)," MB."))
    unname(sapply(file.path(to_dir,clean_files),unlink))
    print(paste0("Directory size is now ",round(sum(unname(sapply(file.path(to_dir,list.files(to_dir)),file.size)))/1024^2,digits=2)," MB."))


