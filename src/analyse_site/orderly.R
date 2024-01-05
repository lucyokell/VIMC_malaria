# analyse site  ----------------------------------------------------------------
orderly2::orderly_parameters(iso3c = 'MDG', 
                             scenario = 'malaria-rts3-rts4-default',
                             quick_run =TRUE,
                             parameter_draw = 0,
                             description =  'runtime_test')


orderly2::orderly_description('Analyze vaccine impact at the site level')
orderly2::orderly_artefact('Processed output', 'processed_output.rds')

# packages
pkgs <- c( 'site', 'data.table',  'dplyr', 'scene', 'malariasimulation', 'openxlsx', 'ggplot2', 'tidyr', 'tibble', 'postie', 'data.table', 'countrycode')
invisible(lapply(pkgs, library, character.only = TRUE))

# functions
source('site_functions.R')
source('summary_functions.')
# read in dependencies  --------------------------------------------------------
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(vimc_input.rds = "vimc_input.rds"))
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(site_file.rds = "site_file.rds"))

vimc_input<- readRDS('vimc_input.rds') 
site_data <- readRDS('site_file.rds')

# vimc inputs
coverage_data<- vimc_input$coverage_input
le <- vimc_input$le
vimc_pop<- vimc_input$population_input_all_age
pop_single_yr<- vimc_input$population_input_single_yr
more_params<- pull_age_groups_time_horizon(quick_run = quick_run, 
                                           coverage_dt = coverage_data, 
                                           scenario = scenario)

# make a map of site + ur combinations
site_data<- remove_zero_eirs(iso3c, site_data)
site_info<- data.table('site_name' = site_data$sites$name_1, 'ur' = site_data$sites$urban_rural)  


# run analysis function for each site + urban/rural combination
output<- lapply( c(1:length(site_info)),
                  analyse_site,
                  site_df= site_info)
                  



