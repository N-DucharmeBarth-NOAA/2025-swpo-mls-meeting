target_dir = "D:/HOME/SAP/Code/2025-swpo-mls-meeting/models/stock-synthesis/xx-06-simple-3par-noDevs"
tmp_report = SS_output(target_dir,verbose = FALSE,printstats = FALSE)
SS_plots(replist = tmp_report,dir=target_dir,html=FALSE,catchasnumbers = TRUE)
SS_html(tmp_report, filenotes = NULL, plotdir = file.path(target_dir,"plots/"), verbose = TRUE,openfile = FALSE)


n_vec = rep(NA,11)
n_vec[1] = exp(4.56)*1000 # 5.7
for(i in 2:length(n_vec)){
    n_vec[i] = n_vec[i-1]*exp(-0.3)
}
sum(n_vec)

base <- SS_output(dir = "D:/HOME/SAP/Code/2025-swpo-mls-meeting/models/stock-synthesis/xx-01-simple")
alternative <- SS_output(dir = "D:/HOME/SAP/Code/2025-swpo-mls-meeting/models/stock-synthesis/xx-02-simple-prior")
new =  SS_output(dir = "D:/HOME/SAP/Code/2025-swpo-mls-meeting/models/stock-synthesis/xx-04-simple-relax-sigmaR")
no_devs =  SS_output(dir = "D:/HOME/SAP/Code/2025-swpo-mls-meeting/models/stock-synthesis/xx-05-simple-full-prior-noDevs")

# Combine into a summary object
comparison <- SSsummarize(list(alternative,new,no_devs))

# Make comparison plots with custom labels
r4ss::SSplotComparisons(comparison,
                legendlabels = c("xx-02-simple-prior","xx-04-simple-relax-sigmaR","xx-05-simple-full-prior-noDevs"),
                png = TRUE,  # Save as PNG files
                plotdir = file.path("D:/HOME/SAP/Code/2025-swpo-mls-meeting/models/stock-synthesis/xx-05-simple-full-prior-noDevs","comp-plots"))

