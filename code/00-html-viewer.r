

# Nicholas Ducharme-Barth
# 2025/01/17
# Generate and clean-up html files for viewing individual model results

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
    html_location = file.path(proj_dir,"html-dashboard")
    dir.create(html_location,recursive=TRUE)
#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(source_dir_stem,(list.files(source_dir_stem))),source)
#________________________________________________________________________________________________________________________________________________________________________________________________________
# define OS
    linux_ws = is_linux_os() # flag to switch code depending on which workspace is being used
#________________________________________________________________________________________________________________________________________________________________________________________________________
# show available directories
    all_dirs = list.files(model_stem,recursive = TRUE)
    all_dirs = all_dirs[grep("summary.csv",all_dirs,fixed=TRUE)]
    all_dirs = all_dirs[-grep("fleet_summary.csv",all_dirs,fixed=TRUE)]
    all_dirs = gsub("summary.csv","",all_dirs,fixed=TRUE)
    print(all_dirs)
    target_dir = file.path(model_stem,all_dirs)[1] # specify which model you want to make the html for
    
    # make html & plots for desired model run
    ss_model_html_viewer(target_dir,html_dir=html_location,linux_ws=linux_ws, catchasnumbers = TRUE)

    # if viewing on GitHub Codespaces, navigate to the file.
    # Right-click and select preview. This will open the file in a window within Codespaces.

    # clean-up since these plots take up a lot of space
    clean_files = c(list.files(html_location),paste0("plots/",list.files(file.path(html_location,"plots"))))
    print(paste0(length(clean_files)," files removed taking up ",round(sum(unname(sapply(file.path(html_location,clean_files),file.size)))/1024^2,digits=2)," MB."))
    unname(sapply(file.path(html_location,clean_files),unlink,recursive=TRUE))
    print(paste0("Directory size is now ",round(sum(unname(sapply(file.path(html_location,list.files(html_location)),file.size)))/1024^2,digits=2)," MB."))






