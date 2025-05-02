# process country  -------------------------------------------------------------
# orderly metadata  ----
orderly2::orderly_parameters(iso3c = NULL,
                             scenario = NULL,
                             quick_run = NULL,
                             parameter_draw = NULL,
                             description = NULL)

orderly2::orderly_description('Analyze vaccine impact at the country level')
orderly2::orderly_artefact('Processed output', 'outputs.rds')

# packages and functions ----
library(site)
library(data.table)
library(dplyr)
library(scene)
library(malariasimulation)
library(tidyr)
library(tibble)
library(postie)
library(countrycode)
library(vimcmalaria)

# read in dependencies  ----
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(vimc_input.rds = "vimc_input.rds"))
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(site_file.rds = "site_file.rds"))

vimc_input<- readRDS('vimc_input.rds')
site_data <- readRDS('site_file.rds')

# vimc inputs ----
coverage_data<- vimc_input$coverage_input
le <- vimc_input$le
vimc_pop<- vimc_input$population_input_all_age
pop_single_yr<- vimc_input$population_input_single_yr

# make a map of input parameters for site function
site_df<- remove_zero_eirs(iso3c, site_data)
map<- vimcmalaria::make_analysis_map(site_df, site_data, test = FALSE, run_all = FALSE)


# run analysis function for each site + urban/rural combination ----
cluster_cores <- Sys.getenv("CCP_NUMCPUS")
if (cluster_cores == "") {
  message("running in serial (on a laptop?)")
  output<- lapply(map,
                   analyse_site,
                   site_data= site_data,
                   vimc_input=vimc_input)
} else {
  message(sprintf("running in parallel on %s (on the cluster?)", cluster_cores))
  cl <- parallel::makeCluster(as.integer(cluster_cores))
  invisible(parallel::clusterCall(cl, ".libPaths", .libPaths()))
  parallel::clusterCall(cl, function() {
    message('running')
    library(data.table)
    library(dplyr)
    library(scene)
    library(malariasimulation)
    library(tidyr)
    library(tibble)
    library(postie)
    library(countrycode)
    library(site)
    library(vimcmalaria)
    TRUE
  })
  output<- parallel::clusterApply(cl,
                                  map,
                                  analyse_site,
                                  site_data= site_data,
                                  vimc_input=vimc_input)
  parallel::stopCluster(cl)
}

# reformat outputs into separate data frames
test<- reformat_output(output)
processed_results<- test$processed_full
raw_output<- test$raw_full


# aggregate outputs up to country level
if(scenario == 'no-vaccination'){

  dt<- aggregate_outputs(processed_results, pop_single_yr)

} else{

  dt<- data.table()
}


#save every output to one list
outputs<- list('country_output' = dt,
               'site_output' = output,
               'raw_output' = raw_output)


saveRDS(outputs, 'outputs.rds')
