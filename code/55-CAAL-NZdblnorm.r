

# Nicholas Ducharme-Barth
# 2025/04/07
# Double-normal selex for NZ recreational fishery; fix terminal selex


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

    selex_vec = c(0.01,0.25,0.5,0.75,0.99)

for(i in 1:length(selex_vec)){
    #________________________________________________________________________________________________________________________________________________________________________________________________________
    # transfer files
        from_dir = dir_base_stock_synthesis

        to_dir = file.path(dir_model,"stock-synthesis",paste0(54+i,"-CAAL-NZdblnorm-",gsub("0.","",as.character(selex_vec[i]),fixed=TRUE)))
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

    tmp_ctl$size_selex_types$Pattern[10] = 24

        new_pars = make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,-10,-3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,10,-5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,boot::logit(selex_vec[i]),-3))

    
    pointer1 = max(grep("F09",rownames(tmp_ctl$size_selex_parms),fixed=TRUE))
    pointer2 = min(grep("F11",rownames(tmp_ctl$size_selex_parms),fixed=TRUE))
    tmp_ctl$size_selex_parms = rbind(tmp_ctl$size_selex_parms[1:pointer1,],
                          new_pars,
                          tmp_ctl$size_selex_parms[pointer2:nrow(tmp_ctl$size_selex_parms),])

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


}

