

# Nicholas Ducharme-Barth
# 2025/01/20
# script to launch shiny app

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(data.table)
library(markdown)
library(magrittr)
library(ggplot2)
library(viridis)
library(ggthemes)
# library(DT)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set relative paths
	proj_dir = this.path::this.proj()
  model_stem = file.path(proj_dir,"models","stock-synthesis")

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

#_____________________________________________________________________________________________________________________________
# app options
  start_collapsed = FALSE

#_____________________________________________________________________________________________________________________________
# source ui/server
  source(file.path(proj_dir,"code","shiny","css.r"))
  source(file.path(proj_dir,"code","shiny","ui.R"))
  source(file.path(proj_dir,"code","shiny","server.R"))

#_____________________________________________________________________________________________________________________________
# Run the app
  app = shinyApp(ui=ui,server=server)
  runApp(app, port = 8888)
