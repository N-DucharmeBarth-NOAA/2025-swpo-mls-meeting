

# Nicholas Ducharme-Barth
# 2025/03/18
# plots for report

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
library(data.table)
library(magrittr)
library(ggplot2)
library(viridis)
library(FLR4MFCL)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set relative paths
	proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# pres plots
    plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","1979_20p3","01-mls-base-1979","03-chg-selex-1979","04-start-1952","06-exclude-more-comp","07-catch-uncertainty","12-CAAL-old-growth-SD"),
                         c(rep(file.path(proj_dir,"models","mfcl"),3),rep(file.path(proj_dir,"models","stock-synthesis"),6)),
                          categories = c(  "Depletion (D)"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "mfcl-1979-20p3","01-mls-base-1979","03-chg-selex-1979","04-start-1952","06-exclude-more-comp","07-catch-uncertainty","12-CAAL-old-growth-SD"),
                          custom_colors = c("gray80", "gray40", "black",viridis::turbo(6, begin = 0.1, end = 0.8, direction=1)),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.pres_stepwise_6",
                          file_type = "png")

        plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","1979_20p3","01-mls-base-1979","03-chg-selex-1979","04-start-1952","06-exclude-more-comp","07-catch-uncertainty"),
                         c(rep(file.path(proj_dir,"models","mfcl"),3),rep(file.path(proj_dir,"models","stock-synthesis"),5)),
                          categories = c( "Depletion (D)"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "mfcl-1979-20p3","01-mls-base-1979","03-chg-selex-1979","04-start-1952","06-exclude-more-comp","07-catch-uncertainty"),
                          custom_colors = c("gray80", "gray40", "black",viridis::turbo(6, begin = 0.1, end = 0.8, direction=1)[1:5]),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.pres_stepwise_5",
                          file_type = "png")
        plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","1979_20p3","01-mls-base-1979","03-chg-selex-1979","04-start-1952","06-exclude-more-comp"),
                         c(rep(file.path(proj_dir,"models","mfcl"),3),rep(file.path(proj_dir,"models","stock-synthesis"),4)),
                          categories = c( "Depletion (D)"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "mfcl-1979-20p3","01-mls-base-1979","03-chg-selex-1979","04-start-1952","06-exclude-more-comp"),
                          custom_colors = c("gray80", "gray40", "black",viridis::turbo(6, begin = 0.1, end = 0.8, direction=1)[1:4]),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.pres_stepwise_4",
                          file_type = "png")   
        plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","1979_20p3","01-mls-base-1979","03-chg-selex-1979","04-start-1952"),
                         c(rep(file.path(proj_dir,"models","mfcl"),3),rep(file.path(proj_dir,"models","stock-synthesis"),3)),
                          categories = c( "Depletion (D)"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "mfcl-1979-20p3","01-mls-base-1979","03-chg-selex-1979","04-start-1952"),
                          custom_colors = c("gray80", "gray40", "black",viridis::turbo(6, begin = 0.1, end = 0.8, direction=1)[1:3]),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.pres_stepwise_3",
                          file_type = "png")   
        plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","1979_20p3","01-mls-base-1979","03-chg-selex-1979"),
                         c(rep(file.path(proj_dir,"models","mfcl"),3),rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                          categories = c("Depletion (D)"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "mfcl-1979-20p3","01-mls-base-1979","03-chg-selex-1979"),
                          custom_colors = c("gray80", "gray40", "black",viridis::turbo(6, begin = 0.1, end = 0.8, direction=1)[1:2]),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.pres_stepwise_2",
                          file_type = "png") 
        plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","1979_20p3","01-mls-base-1979"),
                         c(rep(file.path(proj_dir,"models","mfcl"),3),rep(file.path(proj_dir,"models","stock-synthesis"),1)),
                          categories = c("Depletion (D)"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "mfcl-1979-20p3","01-mls-base-1979"),
                          custom_colors = c("gray80", "gray40", "black",viridis::turbo(6, begin = 0.1, end = 0.8, direction=1)[1]),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.pres_stepwise_1",
                          file_type = "png") 

    plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","12-CAAL-old-growth-SD"),
                         c(rep(file.path(proj_dir,"models","mfcl"),2),rep(file.path(proj_dir,"models","stock-synthesis"),1)),
                          categories = c( "SSB", "Depletion (D)","Recruitment"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "12-CAAL-old-growth-SD"),
                          custom_colors = c("gray80", "gray40", viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)[2]),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.pres_start1952",
                          file_type = "png")


    plot_model_comparison_selex(c("12-CAAL-old-growth-SD"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),1)),
                          n_col = 4,
                          model_labels = c("12-CAAL-old-growth-SD"),
                          custom_colors = viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)[2],
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "selex_comp.pres",
                          file_type = "png")    

    plot_composition_comparison(c("06-exclude-more-comp"),
                                    c(rep(file.path(proj_dir,"models","stock-synthesis"),1)), 
                                    comp_type = "length", # Can be "length" or "weight"
                                    n_col = 3, 
                                    model_labels = c("06-exclude-more-comp"),
                                    custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)[2]),
                                    save_plot = TRUE, 
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "len_comp.pres_excomp",
                                    file_type = "png") 

    plot_composition_comparison(c("06-exclude-more-comp"),
                                    c(rep(file.path(proj_dir,"models","stock-synthesis"),1)), 
                                    comp_type = "weight", # Can be "length" or "weight"
                                    n_col = 3, 
                                    model_labels = c("06-exclude-more-comp"),
                                    custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)[2]),
                                    save_plot = TRUE, 
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "wt_comp.pres_excomp",
                                    file_type = "png") 
