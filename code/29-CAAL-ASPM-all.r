

# Nicholas Ducharme-Barth
# 2025/04/07
# full ASPM


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
    to_dir = file.path(dir_model,"stock-synthesis","29-CAAL-ASPM-all")
    dir.create(to_dir,recursive=TRUE)

    ss_transfer(from_dir=from_dir,to_dir=to_dir,executable_dir_stem=executable_dir_stem,linux_ws=is_linux_os())

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in files
    # tmp_starter = SS_readstarter(file=file.path(to_dir,"starter.ss"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=file.path(to_dir,"control.ss"),datlist = file.path(to_dir,"data.ss"))
    # tmp_data = SS_readdat(file=file.path(to_dir,"data.ss"))
    # tmp_forecast = SS_readforecast(file=file.path(to_dir,"forecast.ss"),verbose=FALSE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# modify files; ctl
# fix selex
    tmp_ctl$size_selex_parms$PHASE = -1 * abs(tmp_ctl$size_selex_parms$PHASE)
    
# lambda for fleets
    tmp_lambdas_surv = as.data.table(tmp_ctl$lambdas)

    tmp_lambdas_gs1 = as.data.table(expand.grid(like_comp=c(4,5,7),fleet=1:16)) %>%
                      .[order(like_comp,fleet)] %>%
                      .[,phase:=1] %>%
                      .[,value:=0] %>%
                      .[,sizefreq_method:=1]
    
    tmp_lambdas_gs2 = as.data.table(expand.grid(like_comp=c(6),fleet=c(2,3,6,7,9,10,16))) %>%
                      .[order(like_comp,fleet)] %>%
                      .[,phase:=1] %>%
                      .[,value:=0] %>%
                      .[,sizefreq_method:=1]

    tmp_ctl$lambdas = as.data.frame(rbind(tmp_lambdas_surv,tmp_lambdas_gs1,tmp_lambdas_gs2))
    tmp_ctl$N_lambdas = nrow(tmp_ctl$lambdas) 

# turn off growth
    tmp_MG_parms = tmp_ctl$MG_parms
    
    tmp_MG_parms$PHASE[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = -7
    tmp_MG_parms$PHASE[grep("L_at_Amax_",rownames(tmp_MG_parms),fixed=TRUE)] = -6
    tmp_MG_parms$PHASE[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = -8

    tmp_ctl$MG_parms = tmp_MG_parms

#________________________________________________________________________________________________________________________________________________________________________________________________________
# write-out files
    # SS_writestarter(tmp_starter,dir=to_dir,overwrite=TRUE)
    # SS_writeforecast(tmp_forecast,dir=to_dir,overwrite=TRUE)
    # SS_writedat(tmp_data, outfile=file.path(to_dir,"data.ss"), overwrite=TRUE)
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

    # output = SS_output(to_dir)
    # SS_plots(output,dir=to_dir)

    # clean directory
    # clean_dir(to_dir)
    # make_ss_output(to_dir)
    
    # make html viewer
    # ss_model_html_viewer(to_dir)


