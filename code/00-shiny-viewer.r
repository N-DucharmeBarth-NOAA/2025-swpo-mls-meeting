

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

#_____________________________________________________________________________________________________________________________
# Run the app
  # app = shinyApp(ui=ui,server=server)
  runApp(file.path(proj_dir,"shiny"), port = 8888)
