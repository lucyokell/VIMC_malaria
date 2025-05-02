################################################################################
## title    VIMC workflow
## author   Lydia Haile
## purpose  run models for VIMC
################################################################################


# setup ------------------------------------------------------------------------
#orderly2::orderly_init()
library(data.table)
library(dplyr)
library(hipercow)
library(vimcmalaria)

coverage<- read.csv('src/process_inputs/vimc_inputs/vaccine_coverage/coverage_202409malaria-1_malaria-r3-r4-default.csv')
iso3cs<- unique(coverage$country_code)
dir<- getwd()

# # generate parameter map for analysis ------------------------------------------
map<- make_parameter_map(iso3cs= c('KEN', 'MRT', 'UGA'),
                         scenarios = c('no-vaccination'),
                          description = 'new_site_files',
                          parameter_draws = 0,
                          quick_run= FALSE)

map$scenario<- 'proxy'
# check function
test<- check_not_a_rerun(report_name = 'process_country', map = map, date = 0)
#map<- check_reports_completed('process_country', map, date_time = 0)
# # STEP 1: run process_inputs report --------------------------------------------
#lapply(iso3cs, function(x) orderly2::orderly_run('process_inputs', parameters = list(iso3c = x)))
#reports<- vimcmalaria::completed_reports('process_country') |> filter(description == 'booster_update')
#
# unique(reports)
# # STEP 2: run process_country for all countries (on cluster)  ------------------
# # if needed, test a report locally before full launch
run_local_reports(map, 'process_country')

# # cluster setup ------
hipercow::hipercow_init(driver = 'windows')
hipercow::hipercow_provision()
hipercow::hipercow_environment_create()
hipercow::hipercow_configuration()

#hipercow::task_log_watch(task)
# # submit groups of jobs by number of cores to submit  ------------------------
#number_order<- c(32, 30, 2, 28, 4, 24, 8, 23, 9, 20, 12, 18, 13, 17, 15, 16,7,3, 1, 5, 6, 7, 10, 11, 14) # intersperse the tasks so the cluster is at optimal usage
lapply(unique(map$site_number), submit_by_core, dt = map, test = FALSE)


for(iso in iso3cs)
task<- hipercow::task_create_expr(
  orderly2::orderly_run(
  "process_country",
  parameters = list(iso3c = 'GHA',
                    description = 'fix_booster_coverage',
                    quick_run = FALSE,
                    scenario = 'malaria-rts3-rts4-default',
                    parameter_draw = 0))
  )
# # launch ethiopia calibrations and save somewhere central ----------------------
# STEP 3: run postprocessing on outputs   --------------------------------------
for(iso in c('CMR', 'COD', 'NGA') ){

message(iso)

#task<- hipercow::task_create_expr(
orderly2::orderly_run(
    "postprocessing",
    parameters = list(
      iso3c = iso,
      description = 'new_site_files',
      quick_run = FALSE
    ))
#)
}


# # STEP 4: run diagnostic reports for outputs  ----------------------------------
for(iso in iso3cs[30:31]){

  message(iso)
#task<- hipercow::task_create_expr(
  orderly2::orderly_run(
    "diagnostics",
    parameters = list(
      iso3c = iso,
      description = 'new_site_files',
      quick_run = FALSE
    ))
#)
}

# compile outputs to shared filepath
compile_diagnostics(descrip = 'new_site_files', date_time = 0)

files<- list.files('montagu/', full.names = TRUE)
files<- files[files %like% 'r3']

for(scen in c('malaria-r3-bluesky', 'malaria-r3-default', 'malaria-r3-r4-bluesky', 'malaria-r3-r4-default')){

  full<- read.csv(paste0('montagu/central-burden-est-', scen, '.csv'))
  message(nrow(full))
  full<- full |>
    filter(country != 'ETH')



  eth<- compile_final_outputs('fix_rtss_booster') |>
    filter(country== 'ETH') |>
    filter(scenario == scen) |>
    select(-run_id) |>
    select(-scenario)


  full<- rbind(full, eth)#

  message(nrow(full))

  write.csv(full, paste0('montagu/update/central-burden-est-', scen, '.csv'), row.names = FALSE)

}




site_output<- lapply(c(1:nrow(completed)), get_site_outputs, map = completed, output_filepath = 'archive/process_country/')








# compile final outputs


task<- hipercow::task_create_expr(vimcmalaria::compile_and_save('booster_update', 'no-vaccination'))
task<- hipercow::task_create_expr(vimcmalaria::compile_and_save('booster_update', 'malaria-r3-r4-default'))
task<- hipercow::task_create_expr(vimcmalaria::compile_and_save('booster_update', 'malaria-rts3-rts4-default'))


hipercow::task_log_watch(task)


