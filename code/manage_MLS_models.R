library(r4ss)
library(callr)
library(processx)
library(future)
library(future.callr)

projdir <- here::here()
basedir <- paste0(projdir, "/2025-swpo-mls-meeting/models/stock-synthesis/")
rundirs <- list.dirs(basedir, full.names = TRUE, recursive=FALSE)


for(dirx in rundirs[-1]){
  setwd(dirx)
  if(file.exists("plots")) next
  if(!file.exists("Report.sso")) next
  tmp_out = r4ss::SS_output(dirx, covar = FALSE)
  r4ss::SS_plots(tmp_out)
}

cbind(1:length(rundirs), rundirs)

dirx <- rundirs[15]
setwd(dirx)
tmp_out = r4ss::SS_output(dirx, covar = FALSE)
r4ss::SS_plots(tmp_out)

exepath <- "C:/Users/simon/OneDrive/Consulting/Fisheries_NZ/International/2025_PAW/MLS/2025-swpo-mls-meeting/executables/stock-synthesis/3.30.23.1/ss3_win.exe"

ss_jobs <- list()


max_jobs <- 8  # Adjust based on your CPU capacity

# List to track running processes
running_jobs <- list()

# Helper: Check if model already ran (e.g., check for Report.sso)
model_ran <- function(dir) {
  file.exists(file.path(dir, "CompReport.sso"))
}

# Directories still to run
pending_dirs <- rundirs[!vapply(rundirs, model_ran, logical(1))]

# Launch and monitor
while (length(pending_dirs) > 0 || length(running_jobs) > 0) {
  
  # Clean up any finished jobs
  finished <- sapply(running_jobs, function(p) !p$is_alive())
  if (any(finished)) {
    cat("Cleaning up finished jobs...\n")
    running_jobs <- running_jobs[!finished]
  }
  
  # Start new jobs if under limit
  while (length(running_jobs) < max_jobs && length(pending_dirs) > 0) {
    dir <- pending_dirs[[1]]
    pending_dirs <- pending_dirs[-1]
    
    cat("Starting SS in:", dir, "\n")
    p <- process$new(
      command = exepath,
      args = character(0),
      wd = dir,
      stdout = file.path(dir, "ss_out.txt"),
      stderr = file.path(dir, "ss_err.txt"),
      cleanup = TRUE
    )
    
    running_jobs[[dir]] <- p
  }
  
  cat(length(running_jobs), "job(s) running. Waiting...\n")
  Sys.sleep(5)
}

cat("✅ All selected models are done!\n")

####################################################
library(callr)


library(future.callr)
library(furrr)
plan(callr, workers = 9)  # Change number of parallel jobs here

# Filter dirs needing plots
rd <- rundirs[2:length(rundirs)]
pending_dirs <- rd[!vapply(rd, function(d) file.exists(file.path(d, "plots")), logical(1))]

# The plotting function
plot_model <- function(dir) {
  library(r4ss)
  setwd(dir)
  out <- SS_output(".", covar = FALSE, verbose = FALSE)
  SS_plots(out, verbose = FALSE)
  return(paste("✅ Plotted:", dir))
}

# Launch all plots in parallel futures
future_results <- future_map(pending_dirs, plot_model)

# run12
r12 <- list()
r12$dir <- rundirs[grep("12-CAAL", rundirs)]

## Run a model with selectivity set to 0 for age 0. 
# Load all the data
library(r4ss)
r12$dat = SS_readdat(file = file.path(r12$dir, 'data.ss'))
r12$ctl = SS_readctl(file = file.path(r12$dir, 'control.ss_new'), datlist = r12$dat)
r12$fore = SS_readforecast(file = file.path(r12$dir, 'forecast.ss'))
r12$start = SS_readstarter(file = file.path(r12$dir, 'starter.ss'))

tmp_out = r4ss::SS_output(r12$dir, covar = FALSE)
setwd(r12$dir)
tune_comps(tmp_out, option = 'Francis')

#######################################################

r27 <- r12
r27$dir <- paste0(basedir, "27-CAAL-agesel10")
dir.create(r27$dir)
f_no0 <- c(2,3,5,6,11,14,16)
r27$ctl$age_selex_parms[f_no0 * 2 - 1,"INIT"] <- 2


SS_writedat(r27$dat, outfile = file.path(r27$dir, 'data.ss'), overwrite = T)
SS_writectl(r27$ctl, outfile = file.path(r27$dir, 'control.ss'), overwrite = T)
SS_writeforecast(r27$fore, dir = r27$dir, overwrite = T)
SS_writestarter(r27$start, dir = r27$dir, overwrite = T)

# Run model:
setwd(r27$dir)
tmp_out = r4ss::SS_output(r27$dir, covar = FALSE)
r4ss::SS_plots(tmp_out)

job <- r_bg(
  func = function(ss_exe, dir) {
    setwd(dir)
    system2(ss_exe, stdout = "ss_out.txt", stderr = "ss_err.txt")
  },
  args = list(ss_exe = ss_exe, dir = r27$dir)
)
job$is_alive()
job$get_exit_status()
job$read_output_lines()

#######################################################
## Upwt NZ rec weight data 
# Load all the data
r28 <- r12
r28$dir <- paste0(basedir, "28-CAAL-upwtNZ")
dir.create(r28$dir)
a <- r28$ctl$Variance_adjustment_list
a[a$factor==7 & a$fleet==10,] <- 1
r28$ctl$Variance_adjustment_list <- a

SS_writedat(r28$dat, outfile = file.path(r28$dir, 'data.ss'), overwrite = T)
SS_writectl(r28$ctl, outfile = file.path(r28$dir, 'control.ss'), overwrite = T)
SS_writeforecast(r28$fore, dir = r28$dir, overwrite = T)
SS_writestarter(r28$start, dir = r28$dir, overwrite = T)

# Run model:
job28 <- r_bg(
  func = function(ss_exe, dir) {
    setwd(dir)
    system2(ss_exe, stdout = "ss_out.txt", stderr = "ss_err.txt")
  },
  args = list(ss_exe = ss_exe, dir = r28$dir)
)

setwd(r28$dir)
tmp_out = r4ss::SS_output(r28$dir, covar = FALSE)
r4ss::SS_plots(tmp_out)

#######################################################
## downwt 6.AU weight data
# Load all the data
r29 <- r12
r29$dir <- paste0(basedir, "29-CAAL-dnwt6AU")
dir.create(r29$dir)
a <- r29$ctl$Variance_adjustment_list
a[a$factor==7 & a$fleet==6,"value"] <- 0.001
r29$ctl$Variance_adjustment_list <- a

SS_writedat(r29$dat, outfile = file.path(r29$dir, 'data.ss'), overwrite = T)
SS_writectl(r29$ctl, outfile = file.path(r29$dir, 'control.ss'), overwrite = T)
SS_writeforecast(r29$fore, dir = r29$dir, overwrite = T)
SS_writestarter(r29$start, dir = r29$dir, overwrite = T)

# Run model:
job29 <- r_bg(
  func = function(ss_exe, dir) {
    setwd(dir)
    system2(ss_exe, stdout = "ss_out.txt", stderr = "ss_err.txt")
  },
  args = list(ss_exe = ss_exe, dir = r29$dir)
)

setwd(r29$dir)
tmp_out = r4ss::SS_output(r29$dir, covar = FALSE)
r4ss::SS_plots(tmp_out)

#######################################################
# Profile
  dir_prof = paste0(basedir, '/profiles')
  dir.create(dir_prof, showWarnings = FALSE)
  
  # Write SS files:
  SS_writedat(r12$dat, outfile = file.path(dir_prof, 'data.ss'), overwrite = T)
  SS_writectl(r12$ctl, outfile = file.path(dir_prof, 'control.ss'), overwrite = T)
  SS_writeforecast(r12$fore, dir = dir_prof, overwrite = T)
  SS_writestarter(r12$start, dir = dir_prof, overwrite = T)
  
  # vector of values to profile over
  R0.vec <- seq(5,6,0.1)
  Nprof.R0 <- length(R0.vec)
  #Define directory
  #mydir <- mydir
  
  #Define the starter file
  starter <- SS_readstarter(file.path(dir_prof, "starter.ss"))
  
  #Change control file name in the starter file
  starter$ctlfile <- "control_modified.ss" 
  
  # Make sure the prior likelihood is calculated for non-estimated quantities
  starter$init_values_src <- 0
  starter$prior_like <- 1                                 
  
  SS_writestarter(starter, dir=dir_prof, overwrite=TRUE)
  
  #Run SS_profile command
  ncores <- parallelly::availableCores(omit = 10)
  future::plan(future::multisession, workers = ncores)
  prof.table <- r4ss::profile(
    dir = dir_prof,
    exe=ss_exe,
    oldctlfile="control.ss",
    newctlfile="control_modified.ss",
    string="SR_LN(R0)",
    profilevec=R0.vec,
    conv_criteria = 0.0001,
    extras = "-cbs 400000000 -gbs 2000000000 -nohess")
  
  future::plan(future::sequential)
  



