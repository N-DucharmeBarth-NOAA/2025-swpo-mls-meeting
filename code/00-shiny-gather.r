

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

#________________________________________________________________________________________________________________________________________________________________________________________________________
# update summary files
    for(i in seq_along(all_dirs)){
        # check for all needed files
        summary_files= c("lik_tbl.csv","est_par.csv","ssb.csv","f.csv","rec.csv","kobe.csv","dyn_dep.csv","srr.csv","bio.csv","selex_l.csv","comp_len.csv","comp_size.csv","summary.csv","catch.csv","bio_age.csv","mature_n.csv","varadj.csv","lambdas.csv","html_parameters.txt","html_parameters_rownames.txt","html_recruitpars.txt","html_recruitpars_rownames.txt","html_estimated_non_dev_parameters.txt","html_estimated_non_dev_parameters_rownames.txt","html_derived_quants.txt","html_derived_quants_rownames.txt")

        if(sum(list.files(file.path(model_stem,all_dirs[i])) %in% summary_files) != length(summary_files)){
            summarize_ss_model(file.path(model_stem,all_dirs[i]))
        }
    }

    

