

# Nicholas Ducharme-Barth
# 2025/01/09
# Update initial MLS synthesis model
# adjust selectivities for
# LL 2
# AU REC 9
# Index 15

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
    dir_base_stock_synthesis = paste0(dir_model,"stock-synthesis/02-mls-ss3-1979/")
    dir_helper_fns = paste0(proj_dir,"/code/helper-fns/")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(paste0(dir_helper_fns,(list.files(dir_helper_fns))),source)

#_____________________________________________________________________________________________________________________________
# read in baseline stock synthesis files
    tmp_starter = SS_readstarter(file=paste0(dir_base_stock_synthesis,"starter.ss_new"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=paste0(dir_base_stock_synthesis,"control.ss_new"),datlist = paste0(dir_base_stock_synthesis,"data_echo.ss_new"))

#_____________________________________________________________________________________________________________________________
# create new directory for stock synthesis mls files
    dir_mls_stock_synthesis_base = paste0(dir_model,"stock-synthesis/03-chg-selex-1979/")
    dir.create(dir_mls_stock_synthesis_base,recursive=TRUE)

#_____________________________________________________________________________________________________________________________
# update starter
    tmp_starter$min_age_summary_bio = 0
    tmp_starter$depl_basis = 2
    tmp_starter$SPR_basis = 2
    SS_writestarter(tmp_starter,dir=dir_mls_stock_synthesis_base,overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update control

    # size based selectivity
    tmp_size_selex_types = as.data.frame(matrix(0,nrow=length(tmp_ctl$fleetnames),ncol=4))
    rownames(tmp_size_selex_types) = tmp_ctl$fleetnames
    colnames(tmp_size_selex_types) = colnames(tmp_ctl$size_selex_types)
    tmp_size_selex_types$Pattern = 24
    tmp_size_selex_types$Pattern[c(5,8,11,13,14)] = 5
    tmp_size_selex_types$Pattern[10] = 1 # nz rec fishery
    tmp_size_selex_types$Special[c(5,8,11,13,14)] = c(4,7,1,3,4)
    tmp_ctl$size_selex_types = tmp_size_selex_types

    # define selex shapes
    # LO, HI, INIT, PHASE 
    tmp_size_selex_parms = rbind(
    # 1 
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,-10,-3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,0,3)),
    # 2
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,0,3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999,-3)),
    # 3
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,-10,-3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,0,3)),
    # 4
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,-10,-3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,0,3)),
    # 5 - m
    make_size_selex_par_5(),
    # 6
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,-10,-3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,0,3)),
    # 7
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,-10,-3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,0,3)),
    # 8 - m
    make_size_selex_par_5(),
    # 9
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,0,3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999,-3)),
    # 10
    make_size_selex_par_1(inflection=c(45,325,100,2),
                          width=c(0,500,10,3)),
    # 11 - m
    make_size_selex_par_5(),
    # 12
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,0,3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999,-3)),
    # 13 - m
    make_size_selex_par_5(),
    # 14 - m
    make_size_selex_par_5(),
    # 15
    make_size_selex_par_24(peak=c(45,325,100,2),
                          top_logit=c(-10,10,0,3),
                          ascend_se=c(-10,10,0,4),
                          descend_se=c(-10,10,0,5),
                          start_logit=c( -999.0,9.0,-999.00000,-4),
                          end_logit=c( -999.0,9.0,-999,-3)))
    tmp_ctl$size_selex_parms = as.data.frame(tmp_size_selex_parms)
    # fix selex of certain parameters
    # tmp_ctl$size_selex_parms$PHASE[c(3,8,20,21,29,34,42,56,57)] = -1*tmp_ctl$size_selex_parms$PHASE[c(3,8,20,21,29,34,42,56,57)]
    # tmp_ctl$size_selex_parms$INIT[c(3,8,20,21,29,34,42,56,57)] = c(7,-7,-7,7,-7,-7,-7,-7,7)

    SS_writectl(tmp_ctl,paste0(dir_mls_stock_synthesis_base,"control.ss"),overwrite = TRUE)

#_____________________________________________________________________________________________________________________________
# run new version of stock synthesis
    file.copy(from=paste0(proj_dir,"/executables/stock-synthesis/3.30.23.1/ss3_win.exe"),to=dir_mls_stock_synthesis_base)
    file.copy(from=paste0(dir_base_stock_synthesis,c("data.ss","forecast.ss")),to=dir_mls_stock_synthesis_base)
    run(dir=dir_mls_stock_synthesis_base,exe="ss3_win.exe",show_in_console = TRUE)

    # output = SS_output(dir_mls_stock_synthesis_base)
    # SS_plots(output,dir=dir_mls_stock_synthesis_base)

    # # directories where models were run need to be defined
    # dir1 = dir_base_stock_synthesis
    # dir2 = dir_mls_stock_synthesis_base

    # # read two models
    # mod1 = SS_output(dir = dir1)
    # mod2 = SS_output(dir = dir2)

    # # create list summarizing model results
    # mod.sum = SSsummarize(list(mod1, mod2))

    # # plot comparisons
    # SSplotComparisons(mod.sum, legendlabels = c("02-mls-ss3-1979", "03-chg-selex-1979"),plotdir=paste0(dir_mls_stock_synthesis_base,"plots/comp/"),plot=FALSE,print=TRUE)
