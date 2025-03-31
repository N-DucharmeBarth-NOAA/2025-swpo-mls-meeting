

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
# extract mfcl in same format as stock synthesis
    # model_id = "1979_20p3"
    # output_dir = file.path(proj_dir,"models","mfcl",model_id)
    # rep_file = file.path(output_dir,"plot-06.par.rep")
    # par_file = file.path(output_dir,"06.par") 

    # new_mfcl = extract_mfcl_quantities(rep_file, par_file, model_id, output_dir, fy = 1979)

    # model_id = "2019-diagnostic"
    # output_dir = file.path(proj_dir,"models","mfcl",model_id)
    # rep_file = file.path(output_dir,"plot-08.par.rep")
    # par_file = file.path(output_dir,"08.par") 

    # diag_2019 = extract_mfcl_quantities(rep_file, par_file, model_id, output_dir, fy = 1952)

    # model_id = "2024-diagnostic"
    # output_dir = file.path(proj_dir,"models","mfcl",model_id)
    # rep_file = file.path(output_dir,"plot-08.par.rep")
    # par_file = file.path(output_dir,"08.par") 

    # diag_2024 = extract_mfcl_quantities(rep_file, par_file, model_id, output_dir, fy = 1979)


    # ss_01 = get_population_numbers(file.path(proj_dir,"models","stock-synthesis","01-mls-base-1979"), "01-mls-base-1979", beginning_of_year = FALSE, exclude_age0 = TRUE)
    # ss_25 = get_population_numbers(file.path(proj_dir,"models","stock-synthesis","25-mls-base-1979-corrected"), "25-mls-base-1979-corrected", beginning_of_year = FALSE, exclude_age0 = TRUE)

    # ss_04 = get_population_numbers(file.path(proj_dir,"models","stock-synthesis","04-start-1952"), "04-start-1952", beginning_of_year = FALSE, exclude_age0 = TRUE)
    
    # Section: pre-meeting
    
    plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","1979_20p3"),
                         c(rep(file.path(proj_dir,"models","mfcl"),3)),
                          categories = c( "SSB", "Depletion (D)","N"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "mfcl-1979-20p3"),
                          custom_colors = c("gray80", "gray40", "black"),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.pre_meeting",
                          file_type = "png")

    # Section: ss3 transition
   
    plot_model_comparison_ts(c("1979_20p3","01-mls-base-1979","25-mls-base-1979-corrected"),
                         c(rep(file.path(proj_dir,"models","mfcl"),1),rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                          categories = c( "SSB", "Depletion (D)","N"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c( "mfcl-1979-20p3", "01-mls-base-1979","25-mls-base-1979-corrected"),
                          custom_colors = c( "black",viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.ss_transition",
                          file_type = "png")

    # Section: early refinements

    plot_model_comparison_selex(c("01-mls-base-1979","02-mls-ss3-1979", "03-chg-selex-1979"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),3)),
                          n_col = 4,
                          model_labels = c("01-mls-base-1979","02-mls-ss3-1979", "03-chg-selex-1979"),
                          custom_colors = c(viridis::turbo(3, begin = 0.1, end = 0.8, direction=1)),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "selex_comp.early",
                          file_type = "png")                      

    plot_index_comparison (c("01-mls-base-1979","02-mls-ss3-1979", "03-chg-selex-1979"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),3)), 
                               n_col = 1, 
                               model_labels = c("01-mls-base-1979","02-mls-ss3-1979", "03-chg-selex-1979"),
                                custom_colors = c(viridis::turbo(3, begin = 0.1, end = 0.8, direction=1)),
                                apply_varadj = FALSE,
                                show_se = FALSE,
                                save_plot = TRUE,
                                save_dir = file.path(proj_dir,"assets","static"),
                                plot_name = "index_comp.early",
                                file_type = "png") 

    plot_composition_comparison(c("02-mls-ss3-1979", "03-chg-selex-1979"),
                                    c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                                    comp_type = "length", # Can be "length" or "weight"
                                    n_col = 4, 
                                    model_labels = c("02-mls-ss3-1979", "03-chg-selex-1979"),
                                    custom_colors = c(viridis::turbo(3, begin = 0.1, end = 0.8, direction=1)[2:3]),
                                    save_plot = TRUE, 
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "len_comp.early",
                                    file_type = "png") 

    plot_composition_comparison(c("02-mls-ss3-1979", "03-chg-selex-1979"),
                                    c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                                    comp_type = "weight", # Can be "length" or "weight"
                                    n_col = 4, 
                                    model_labels = c("02-mls-ss3-1979", "03-chg-selex-1979"),
                                    custom_colors = c(viridis::turbo(3, begin = 0.1, end = 0.8, direction=1)[2:3]),
                                    save_plot = TRUE, 
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "wt_comp.early",
                                    file_type = "png") 

    # Section: start1952
   
    plot_model_comparison_ts(c("2019-diagnostic","2024-diagnostic","04-start-1952"),
                         c(rep(file.path(proj_dir,"models","mfcl"),2),rep(file.path(proj_dir,"models","stock-synthesis"),1)),
                          categories = c( "SSB", "Depletion (D)","Recruitment"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("diagnostic-2019", "diagnostic-2024", "04-start-1952"),
                          custom_colors = c("gray80", "gray40", viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)[2]),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.start1952",
                          file_type = "png")

    plot_index_comparison (c( "03-chg-selex-1979","04-start-1952"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                               n_col = 1, 
                               model_labels = c("03-chg-selex-1979","04-start-1952"),
                                custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                apply_varadj = FALSE,
                                show_se = FALSE,
                                save_plot = TRUE,
                                save_dir = file.path(proj_dir,"assets","static"),
                                plot_name = "index_comp.start1952",
                                file_type = "png") 

    # Section: excomp

    plot_composition_comparison(c("06-exclude-more-comp"),
                                    c(rep(file.path(proj_dir,"models","stock-synthesis"),1)), 
                                    comp_type = "length", # Can be "length" or "weight"
                                    n_col = 4, 
                                    model_labels = c("06-exclude-more-comp"),
                                    custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)[2]),
                                    save_plot = TRUE, 
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "len_comp.excomp",
                                    file_type = "png") 

    plot_composition_comparison(c("06-exclude-more-comp"),
                                    c(rep(file.path(proj_dir,"models","stock-synthesis"),1)), 
                                    comp_type = "weight", # Can be "length" or "weight"
                                    n_col = 4, 
                                    model_labels = c("06-exclude-more-comp"),
                                    custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)[2]),
                                    save_plot = TRUE, 
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "wt_comp.excomp",
                                    file_type = "png") 

    plot_index_comparison (c( "04-start-1952","06-exclude-more-comp"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                               n_col = 1, 
                               model_labels = c("04-start-1952","06-exclude-more-comp"),
                                custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                apply_varadj = TRUE,
                                show_se = TRUE,
                                save_plot = TRUE,
                                save_dir = file.path(proj_dir,"assets","static"),
                                plot_name = "index_comp.excomp",
                                file_type = "png") 

    # catch uncertainty
        plot_catch_comparison(c( "07-catch-uncertainty"),
                               c(rep(file.path(proj_dir,"models","stock-synthesis"),1)),
                               show_fit = TRUE,
                               n_col = 1,
                               fleets = 15,
                               free_y_scale = TRUE,  # Option to have free y scales
                               model_labels = c("07-catch-uncertainty"),
                                custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)[2]),
                               save_plot = TRUE,
                               width = 12,
                               height = 6,
                                save_dir = file.path(proj_dir,"assets","static"),
                                plot_name = "catch_comp.catch_unc",
                                file_type = "png")
        
        plot_model_comparison_ts(c("2019-diagnostic","06-exclude-more-comp","07-catch-uncertainty"),
                         c(rep(file.path(proj_dir,"models","mfcl"),1),rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                          categories = c( "SSB", "Depletion (D)","Recruitment"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c("2019-diagnostic","06-exclude-more-comp","07-catch-uncertainty"),
                          custom_colors = c("gray80", viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.catch_unc",
                          file_type = "png")
        
        plot_index_comparison (c( "06-exclude-more-comp","07-catch-uncertainty"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                               n_col = 1, 
                               model_labels = c("06-exclude-more-comp","07-catch-uncertainty"),
                                custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                apply_varadj = FALSE,
                                show_se = FALSE,
                                save_plot = TRUE,
                                save_dir = file.path(proj_dir,"assets","static"),
                                plot_name = "index_comp.catch_unc",
                                file_type = "png")
        
        # growth

        plot_model_comparison_ts(c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                          categories = c( "SSB", "Depletion (D)","Recruitment"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                          custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.growth",
                          file_type = "png")
        
        plot_composition_comparison(c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                                    comp_type = "length", # Can be "length" or "weight"
                                    n_col = 4, 
                                    model_labels = c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                                    custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                    save_plot = TRUE, 
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "len_comp.growth",
                                    file_type = "png") 

        plot_composition_comparison(c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                                    comp_type = "weight", # Can be "length" or "weight"
                                    n_col = 4, 
                                    model_labels = c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                                    custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                    save_plot = TRUE, 
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "wt_comp.growth",
                                    file_type = "png") 
        
        plot_index_comparison (c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                               n_col = 1, 
                               model_labels = c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                                custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                apply_varadj = TRUE,
                                show_se = TRUE,
                                save_plot = TRUE,
                                save_dir = file.path(proj_dir,"assets","static"),
                                plot_name = "index_comp.growth",
                                file_type = "png")

        plot_model_comparison_bio(c( "07-catch-uncertainty","12-CAAL-old-growth-SD"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                               n_col = 1, 
                               model_labels = c("07-catch-uncertainty","12-CAAL-old-growth-SD"),
                                custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                   show_variability = TRUE, 
                                  save_plot = TRUE,
                                    save_dir = file.path(proj_dir,"assets","static"),
                                    plot_name = "bio_comp.growth",
                                    file_type = "png")

        # plot_model_comparison_ts(c("2024-diagnostic" ,"03-chg-selex-1979","22-1979-estInitF-v2"),
        #                  c(rep(file.path(proj_dir,"models","mfcl"),1),rep(file.path(proj_dir,"models","stock-synthesis"),2)),
        #                   categories = c( "SSB", "Depletion (D)","Recruitment"), 
        #                   show_se = FALSE, n_col = 1,
        #                   model_labels = c("2024-diagnostic" , "03-chg-selex-1979","22-1979-estInitF-v2"),
        #                   custom_colors = c("gray40",viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
        #                   save_plot = TRUE,
        #                   save_dir = file.path(proj_dir,"assets","static"),
        #                   plot_name = "ts_comp.estInitF-v2",
        #                   file_type = "png")


        # sens
        plot_index_comparison (c( "12-CAAL-old-growth-SD","16-CAAL-rm-spike", "17-CAAL-noAU-ASPM"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),3)), 
                               n_col = 1, 
                               model_labels = c( "12-CAAL-old-growth-SD","16-CAAL-rm-spike", "17-CAAL-noAU-ASPM"),
                                custom_colors = c(viridis::turbo(3, begin = 0.1, end = 0.8, direction=1)),
                                apply_varadj = TRUE,
                                show_se = TRUE,
                                save_plot = TRUE,
                                save_dir = file.path(proj_dir,"assets","static"),
                                plot_name = "index_comp.sens_idx",
                                file_type = "png")
        
        plot_model_comparison_ts(c( "12-CAAL-old-growth-SD","23-CAAL-2sex"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                          categories = c( "SSB", "Depletion (D)","Recruitment"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c( "12-CAAL-old-growth-SD","23-CAAL-2sex"),
                          custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.sens_2sex",
                          file_type = "png")
        plot_model_comparison_bio(c( "12-CAAL-old-growth-SD","23-CAAL-2sex"),
                                c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                                    n_col = 1, 
                                    model_labels = c("12-CAAL-old-growth-SD","23-CAAL-2sex"),
                                        custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                        show_variability = TRUE, 
                                        save_plot = TRUE,
                                            save_dir = file.path(proj_dir,"assets","static"),
                                            plot_name = "bio_comp.sens_2sex",
                                            file_type = "png")
        plot_model_comparison_bio(c( "12-CAAL-old-growth-SD","24-CAAL-Richards"),
                                c(rep(file.path(proj_dir,"models","stock-synthesis"),2)), 
                                    n_col = 1, 
                                    model_labels = c("12-CAAL-old-growth-SD","24-CAAL-Richards"),
                                        custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                                        show_variability = TRUE, 
                                        save_plot = TRUE,
                                            save_dir = file.path(proj_dir,"assets","static"),
                                            plot_name = "bio_comp.sens_rich",
                                            file_type = "png")
        plot_model_comparison_ts(c( "12-CAAL-old-growth-SD","21-CAAL-NZrecwtQtr"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                          categories = c( "SSB", "Depletion (D)","Recruitment"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c( "12-CAAL-old-growth-SD","21-CAAL-NZrecwtQtr"),
                          custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.sens_nzrec",
                          file_type = "png")
        
        plot_model_comparison_ts(c( "12-CAAL-old-growth-SD","18-CAAL-cUnc-cv40"),
                         c(rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                          categories = c( "SSB", "Depletion (D)","Recruitment"), 
                          show_se = FALSE, n_col = 1,
                          model_labels = c( "12-CAAL-old-growth-SD","18-CAAL-cUnc-cv40"),
                          custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                          save_plot = TRUE,
                          save_dir = file.path(proj_dir,"assets","static"),
                          plot_name = "ts_comp.sens_catch",
                          file_type = "png")

        plot_catch_comparison(c( "12-CAAL-old-growth-SD","18-CAAL-cUnc-cv40"),
                               c(rep(file.path(proj_dir,"models","stock-synthesis"),2)),
                               show_fit = TRUE,
                               n_col = 1,
                               fleets = 15,
                               free_y_scale = TRUE,  # Option to have free y scales
                               model_labels = c("12-CAAL-old-growth-SD","18-CAAL-cUnc-cv40"),
                                custom_colors = c(viridis::turbo(2, begin = 0.1, end = 0.8, direction=1)),
                               save_plot = TRUE,
                               width = 12,
                               height = 6,
                                save_dir = file.path(proj_dir,"assets","static"),
                                plot_name = "catch_comp.sens_catch",
                                file_type = "png")
