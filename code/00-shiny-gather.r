

# Nicholas Ducharme-Barth
# 2025/01/20
# Create summary files if they do not exist

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(r4ss)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set relative paths
	proj_dir = this.path::this.proj()

    model_stem = file.path(proj_dir,"models","stock-synthesis")
    source_dir_stem = file.path(proj_dir,"code","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(source_dir_stem,(list.files(source_dir_stem))),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# show available directories
    all_dirs = list.files(model_stem,recursive = TRUE)
    all_dirs = all_dirs[grep("executable.txt",all_dirs,fixed=TRUE)]
    all_dirs = gsub("executable.txt","",all_dirs,fixed=TRUE)
    all_dirs = all_dirs[-grep("00-workflow-test",all_dirs,fixe=TRUE)]
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# update summary files
    for(i in seq_along(all_dirs)){
        # check for all needed files
        summary_files= c("comp_len_obs_time.csv","comp_len_exp_time.csv","comp_size_obs_time.csv","comp_size_exp_time.csv","cpue.csv","fleet_summary.csv","lik_tbl.csv","est_par.csv","ssb.csv","f.csv","rec.csv","kobe.csv","dyn_dep.csv","srr.csv","bio.csv","selex_a.csv","selex_l.csv","comp_len.csv","comp_size.csv","summary.csv","catch.csv","bio_age.csv","mature_n.csv","varadj.csv","lambdas.csv","html_parameters.txt","html_parameters_rownames.txt","html_recruitpars.txt","html_recruitpars_rownames.txt","html_estimated_non_dev_parameters.txt","html_estimated_non_dev_parameters_rownames.txt","html_derived_quants.txt","html_derived_quants_rownames.txt")

        if(sum(list.files(file.path(model_stem,all_dirs[i])) %in% summary_files) != length(summary_files)){
            summarize_ss_model(file.path(model_stem,all_dirs[i]))
        }

        # copy to shiny/data dir if they don't already exist there
        dir.create(file.path(proj_dir,"shiny","data",all_dirs[i]),recursive=TRUE)
        if(sum(list.files(file.path(proj_dir,"shiny","data",all_dirs[i])) %in% summary_files) != length(summary_files)){
            file.copy(from=file.path(file.path(model_stem,all_dirs[i],summary_files)),to=file.path(proj_dir,"shiny","data",all_dirs[i],summary_files),overwrite=TRUE)
       }
    }

#_____________________________________________________________________________________________________________________________
# make summary file to seed the shiny app
    all_dirs = list.files(model_stem,recursive = TRUE)
    all_dirs = all_dirs[grep("summary.csv",all_dirs,fixed=TRUE)]
    all_dirs = all_dirs[-grep("fleet_summary.csv",all_dirs,fixed=TRUE)]
    # all_dirs = all_dirs[-grep("surplus-production/runs/",all_dirs,fixed=TRUE)]
    all_dirs = gsub("summary.csv","",all_dirs,fixed=TRUE)

    summary_dt.list = lapply(all_dirs,function(x)as.data.frame(fread(file.path(model_stem,x,"summary.csv"))))

    # round numeric
    for(i in 1:length(summary_dt.list)){
        for(j in 1:ncol(summary_dt.list[[i]])){
            if(class(summary_dt.list[[i]][,j]) == "numeric"){
                # summary_df[,j] = round(summary_df[,j],digits=3)
                summary_dt.list[[i]][,j] = trimws(format(round(summary_dt.list[[i]][,j],digits=3),big.mark=",",scientific=FALSE))
            } else if(class(summary_dt.list[[i]][,j]) == "integer64"){
                summary_dt.list[[i]][,j] = trimws(format(round(as.numeric(summary_dt.list[[i]][,j]),digits=3),big.mark=",",scientific=FALSE))
            }
        }
        summary_dt.list[[i]] = as.data.table(summary_dt.list[[i]])
    }

    summary_dt = rbindlist(summary_dt.list,fill=TRUE)
    fwrite(summary_dt,file=file.path(proj_dir,"data","summary_dt.csv"))
    fwrite(summary_dt,file=file.path(proj_dir,"shiny","data","summary_dt.csv"))
    

