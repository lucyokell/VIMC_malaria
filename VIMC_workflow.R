################################################################################
# VIMC orderly workflow
# author: Lydia Haile
# orderly workflow for VIMC model runs. This script: 
# 1) formats input parameters for model runs
# 2) runs models
# 3) processes country outputs
# 4) generates diagnostics
################################################################################

# packages  
library(orderly2)
library(site)
library(data.table)
library(dplyr)
lapply(list.files('functions/', full.names = T), source)

# obtain list of countries to run model for
coverage<- read.csv('src/process_inputs/vimc_inputs/vaccine_coverage/coverage_202310gavi-1_malaria-r3-default.csv')
iso3cs<- unique(coverage$country_code)

dir<- getwd()

# if you have not already, initialize the orderly repository
#orderly2::orderly_init(path = dir)
################################################################################
# 1 prepare and save inputs
# unless inputs change, this only needs to be run once for all countries
# for (iso3c in iso3cs){
# 
#   orderly2::orderly_run(
#     'process_inputs',
#     list(iso3c = iso3c),
#     root = dir)
# }

# PARAMETERS TO CHANGE FOR REPORTS ---------------------------------------------
maps<- make_parameter_maps(
  iso3cs =  iso3cs,                                                                       # Pick 10 countries to begin with
  scenarios= c('malaria-r3-default', 'malaria-rts3-bluesky', 'malaria-rts3-default'),     # if you only want to run reports for certain scenarios. Default is all 7
  population = 100000,                                                                    # population size
  description = 'complete_run',                                                           # reason for model run (change this for every run if you do not want to overwrite outputs)
  parameter_draw = 0,                                                                     # parameter draw to run (0 for central runs)
  burnin= 15,                                                                             # burn-in in years            
  quick_run = FALSE                                                                       # boolean, T or F. If T, makes age groups larger and runs model through 2035.
)

# deduplicate

site_map<- remove_duplicate_reports(report_name = 'process_site',
                                    parameter_map = site_map, day= 20231208)

# check that the preceding report has completed before you launch next report in chronology
site_map<- generate_parameter_map_for_next_report(report_name = 'launch_models',
                                                  parameter_map = maps$site_map, 
                                                  day= 20231208)

# 
# 
# country_map<- remove_duplicate_reports(report_name = 'process_country',
#                                     parameter_map = maps$country_map, day= 20231130)


# check that the preceding report has completed before you launch next report in chronology
# country_map<- generate_parameter_map_for_next_report(report_name = 'process_country',
#                                                   parameter_map = country_map)

site_map<- maps$site_map
sites<- purrr::map(.x = c(1:nrow(site_map)), .f= ~ site_map[.x,])

country_map<- maps$country_map
countries<- purrr::map(.x = 1:nrow(country_map), .f= ~ country_map[.x,])
# 
# # # cluster setup ------------------------------------------------------------
ctx <- context::context_save("ctxs4", sources= 'functions/run_report.R')
config <- didehpc::didehpc_config(
  use_rrq = FALSE,
  cores = 1,
  cluster = "wpia-hn", #"fi--dideclusthn", # , "fi--didemrchnb""fi--didemrchnb"
  template = "AllNodes")

obj <- didehpc::queue_didehpc(ctx, config = config)

# # if you have not already, install packages:  --------------------------------
# pkgs<- c('mrc-ide/orderly2@mrc-4724',
#          'mrc-ide/malariasimulation',
#          'mrc-ide/site_vimc',
#          'data.table',
#          'dplyr',
#          'extrafont',
#          'wesanderson',
#          'ggpubr',
#          'ggplot2',
#          'ggforce',
#           'openxlsx'
#          'mrc-ide/postie@dalys',
#          'countrycode')
# 
# packages you need for postprocessing
pp<- c('mrc-ide/postie@dalys',
       'data.table',
       'dplyr',
       'countrycode',
       'mrc-ide/orderly2',
       'extrafont',
       'ggpubr',
       'wesanderson',
       'scales',
       'openxlsx',
       'ggforce',
       'mrc-ide/scene',
       'ggpubr',
       'ggplot2')

for (pkg in pp){

  obj$install_packages('mrc-ide/malariasimulation')

}


# run report for all sites locally ---------------------------------------------
lapply(
  sites,
  run_report,
  report_name = 'process_site',
  path = dir
)


# 
# # or launch on cluster
pp<- obj$lapply(
  sites,
  run_report,
  report_name = 'process_site',
  path = dir
)

#  # run report for all countries locally  -------------------------------------
# lapply(
#   countries,
#   run_report_country,
#   report_name = 'process_country',
#   path = dir
# )
# 
# # or launch cluster ----------------------------------------------------------
retry<- obj$lapply(
  countries,
  run_report_country,
  report_name = 'country_diagnostics',
  path = dir
)





