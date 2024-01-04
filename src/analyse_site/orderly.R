# analyse site  ----------------------------------------------------------------
orderly2::orderly_parameters(iso3c = 'ETH', 
                             site_name = 'Afar',
                             ur = 'rural',
                             scenario = 'malaria-rts3-rts4-default',
                             quick_run =FALSE,
                             parameter_draw = 0,
                             description =  'runtime_test')


orderly2::orderly_description('Analyze vaccine impact at the site level')
orderly2::orderly_artefact('Model input', 'model_input.rds')
orderly2::orderly_artefact('Model output', 'model_output.rds')
orderly2::orderly_artefact('Processed output', 'processed_output.rds')
orderly2::orderly_artefact('Plotting inputs', 'plotting_inputs.rds')

# packages
pkgs <- c( 'site', 'data.table',  'dplyr', 'scene', 'malariasimulation', 'openxlsx', 'ggplot2', 'tidyr', 'tibble', 'postie', 'data.table', 'countrycode')
invisible(lapply(pkgs, library, character.only = TRUE))

# functions
source('site_functions.R')

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


more_params<- pull_age_groups_time_horizon(quick_run = quick_run, coverage_dt = coverage_data, scenario = scenario)

# site data
site <- extract_site(site_file = site_data, site_name = site_name, ur = ur)

# specify vaccine coverage based on forecast  ----------------------------------
site<- expand_intervention_coverage(site, terminal_year = more_params$term_yr)
site<- update_coverage_values(site, coverage_data, scenario_name = scenario)


# add in scenario variable which will be used to implement booster
vaccine_plot_input<- copy(site)
check_eir(site)

# pull parameters for this site ------------------------------------------------
params <- site::site_parameters(
  interventions = site$interventions,
  demography = site$demography,
  vectors = site$vectors,
  seasonality = site$seasonality,
  eir = site$eir$eir[1],
  burnin = more_params$burnin,
  overrides = list(human_population = more_params$pop_val)
)

# set age groups
params$clinical_incidence_rendering_min_ages = more_params$min_ages
params$clinical_incidence_rendering_max_ages = more_params$max_ages
params$severe_incidence_rendering_min_ages = more_params$min_ages
params$severe_incidence_rendering_max_ages = more_params$max_ages
params$age_group_rendering_min_ages = more_params$min_ages
params$age_group_rendering_max_ages = more_params$max_ages

# if this is a stochastic run, set parameter draw ------------------------------
params<- parameterize_stochastic_run(params, parameter_draw)

inputs <- list(
  'param_list' = params,
  'site_name' = site_name,
  'ur' = ur,
  'iso' = iso3c,
  'scenario' = scenario,
  'description' = description,
  'parameter_draw' = parameter_draw
)

message('saving model input')
saveRDS(inputs, 'model_input.rds')

# launch models ----------------------------------------------------------------
model_input <- copy(inputs)

params <- model_input$param_list
params$progress_bar <- TRUE
timesteps <<- model_input$param_list$timesteps

message('running the model')
model <- malariasimulation::run_simulation(timesteps = params$timesteps,
                                           parameters = params)

# add identifying information to output
model <- model |>
  mutate(site_name = site_name,
         urban_rural = ur,
         iso = iso3c,
         description = description, 
         scenario = scenario,
         parameter_draw = parameter_draw,
         population = more_params$pop_val,
         burnin = more_params$burnin)

# save model runs somewhere
message('saving the model')
saveRDS(model, 'model_output.rds')

# process site -----------------------------------------------------------------
# calculate rates
raw_output<- drop_burnin(model, burnin= more_params$burnin* 365)

output <- postie::get_rates(
  raw_output,
  time_divisor = 365,
  baseline_t = 1999,
  age_divisor = 365,
  scaler = 0.215,
  treatment_scaler = 0.517,
)

dt<- vimc_postprocess(output, le,  site_data, vimc_pop, pop_single_yr)


# final formatting  ------------------------------------------------------------
dt<- format_outputs(dt)
saveRDS(dt, 'processed_output.rds')


plotting_input<- pull_plotting_data(scenario)
saveRDS(plotting_input, 'plotting_input.rds')



