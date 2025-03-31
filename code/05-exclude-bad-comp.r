

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
    dir_base_stock_synthesis = file.path(dir_model,"stock-synthesis","04-start-1952")
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
    to_dir = file.path(dir_model,"stock-synthesis","05-exclude-bad-comp")
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

    # - remove NZ wt comp from 1988 on  
    # - remove small sample sizes fleet 1 (weight)
    tmp_wtcomp = as.data.table(tmp_data$sizefreq_data_list[[1]]) %>%
                  .[fleet==10 & year>=1988,year:=abs(year)*-1] %>%
                  .[fleet == 1, year:=abs(year)*-1]
    tmp_data$sizefreq_data_list = list(as.data.frame(tmp_wtcomp))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; control
# change selectivity groupings

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


