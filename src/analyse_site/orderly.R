# analyse site  ----------------------------------------------------------------
orderly2::orderly_parameters(iso3c = 'NGA', 
                             site_name = 'Lagos',
                             ur = 'urban',
                             scenario = 'malaria-rts3-rts4-default',
                             quick_run = TRUE,
                             parameter_draw = 0,
                             description =  'refactor_debug')


orderly2::orderly_description('Analyze vaccine impact at the site level')
orderly2::orderly_artefact('Model input', 'model_input.rds')
orderly2::orderly_artefact('Model output', 'model_output.rds')
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

# site data
site <- extract_site(site_file = site_data, site_name = site_name, ur = ur)
more_params<- pull_age_groups_time_horizon(quick_run = quick_run)

# specify vaccine coverage based on forecast  ----------------------------------
site<- expand_intervention_coverage(site, terminal_year = more_params$term_yr)
site<- update_coverage_values(site, coverage_data, scenario_name = scenario)

# add in scenario variable which will be used to implement booster
saveRDS(site, 'vaccine_plot_input.rds')
check_eir(site)

# pull parameters for this site ------------------------------------------------
params <- site::site_parameters(
  interventions = site$interventions,
  demography = site$demography,
  vectors = site$vectors,
  seasonality = site$seasonality,
  eir = site$eir$eir[1],
  burnin = 15,
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
params<- parameterize_stochastic_run(parameter_draw)

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
         population = pop_val)

# save model runs somewhere
message('saving the model')
saveRDS(model, 'model_output.rds')

# process site -----------------------------------------------------------------
# calculate rates
raw_output<- drop_burnin(model_output, burnin= 15* 365)
saveRDS(raw_output, 'raw_model_output.rds') # save for diagnostics

output <- postie::get_rates(
  raw_output,
  time_divisor = 365,
  baseline_t = 1999,
  age_divisor = 365,
  scaler = 0.215,
  treatment_scaler = 0.517,
)

vimc_postprocess<- function(output, le){
  
  # fill rates out to single year age groups
  output<- output |>
    dplyr::group_by(t) |>
    tidyr::complete(age_lower = c(1:100)) |>
    select(-age_upper) |>
    dplyr::ungroup() |>
    tidyr::fill(clinical, severe, mortality, yld_pp, yll_pp, dalys_pp, .direction = 'down') |>
    select(-prop_n, -n)

  # recalculate YLLs and DALYs based on country-specific life expectancy  --------
  output<- output |>
    select(-yll_pp, -dalys_pp) |>
    rename(year = t)
  
  # merge in inputs for expected remaining years of life (to calculate YLLs)
  le<- expand_life_expectancy(le)
  
  # calculate ylls_pp + dalys per person
  dt<- merge(output, le, by = c('year', 'age_lower'), all.x = TRUE)
  
  dt<- dt |>
    mutate(ylls_pp = mortality * remaining_yrs) |>
    mutate(dalys_pp = ylls_pp + yld_pp) |>
    select(-remaining_yrs)
  
}




# calculate counts  ------------------------------------------------------------

scale_population<- function(site_data, vimc_pop){
  # merge in population from site files (as we only have VIMC inputs for the national level)
  # first, separately sum cases by year
  total_pop<- site_data$population |>
    group_by(year) |>
    summarise(summed_pop = sum(pop))
  
  # pull the population for the site of interest
  pop <- site_data$population |>
    filter(name_1 == site_name & urban_rural == ur) |>
    select(year, pop) |>
    rename(site_file_population = pop)
  
  # merge these two tables together
  pops<- merge(pop, total_pop, by= 'year')
  
  # merge in national population from VIMC (available for entire time period)
  vimc_pop<- vimc_pop |>
    filter(country_code == iso3c,
           year >= 2000)|>
    rename(national_pop = value)|>
    select(year, national_pop)
  
  # merge in vimc population
  pops<- merge(vimc_pop, pops, all.x = T)
  
  # first rescale site file population based on the ratio of (sum of site file pops in country)/ (VIMC country level population)
  # should be more or less the same, but should be done for consistency sake
  
  pops<- pops |>
    mutate(vimc_site_population = (site_file_population * national_pop)/summed_pop)
  
  # calculate population ratio as vimc(site)/ vimc(country)
  pops<- pops |>
    mutate(pop_ratio = vimc_site_population/ national_pop) |>
    tidyr::fill(pop_ratio, .direction = 'down')
  
  # then calculate vimc_site_population by multiplying this ratio to the national population for the final 50 years
  pops<- pops |>
    mutate(vimc_site_population = ifelse(year<= 2050, vimc_site_population, pop_ratio* national_pop))
  
  # subset out site file population for 2000-2100
  site_pop<- pops |>
    select(year, vimc_site_population)
  
  # pull in single year population to calculate proportion_n by age group
  national_pop<- pops |>
    select(year, national_pop)
  
  pop_single_yr<- merge(pop_single_yr, national_pop, by = c('year'))
  pop_single_yr <- pop_single_yr |>
    mutate(prop_n = value/ national_pop) |>
    select(year, age_from, age_to, prop_n) |>
    rename(age_lower = age_from)
  
}


# merge in site population
dt<- merge(dt, site_pop, by= 'year')

# merge in prop_n
dt<- merge(dt, pop_single_yr, by = c('year', 'age_lower'))


# calculate counts for entire time period --------------------------------------
dt<- dt |>
  mutate(
    cases = round(clinical * vimc_site_population * prop_n),
    deaths = round(mortality * vimc_site_population * prop_n),
    dalys = round(dalys_pp * vimc_site_population * prop_n),
    population = round(vimc_site_population * prop_n)) |>
  select(-prop_n)

# final formatting  ------------------------------------------------------------
dt<- format_outputs(dt)
saveRDS(dt, 'processed_output.rds')


if(scenario!="no-vaccination") {
  doses_per_year <- pull_doses_output(raw_output, dt)
  
  saveRDS(doses_per_year, 'doses_per_year.rds')
}


### pull out prevalence
prev <- postie::get_prevalence(raw_output, time_divisor = 365, baseline_t = 1999,
                               age_divisor = 365) 
prev$n_2_10 <- raw_output |>
  mutate(year = floor(timestep/365)) |>
  group_by(year) |>
  summarise(n_2_10 = mean(n_730_3649)) |>
  filter(year<=100) |>
  pull(n_2_10)


saveRDS(prev, 'prevalence_per_year.rds')


