# analyse site  ----------------------------------------------------------------
orderly2::orderly_parameters(iso3c = NULL, 
                             site_name = NULL,
                             ur = NULL,
                             scenario = NULL,
                             quick_run = NULL,
                             parameter_draw = NULL,
                             description =  NULL)


orderly2::orderly_description('Analyze vaccine impact at the site level')

orderly2::orderly_artefact('Model input', 'model_input.rds')
orderly2::orderly_artefact('Model output', 'model_output.rds')
orderly2::orderly_artefact('Vaccine plot input', 'vaccine_plot_input.rds')

# packages
pkgs <- c( 'site', 'data.table',  'dplyr', 'scene', 'malariasimulation', 'openxlsx', 'ggplot2', 'tidyr', 'tibble', 'postie', 'data.table', 'countrycode')
invisible(lapply(pkgs, library, character.only = TRUE))

# functions
lapply(list.files('functions/', full.names = T), source)

# read in dependencies  --------------------------------------------------------
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(vimc_input.rds = "vimc_input.rds"))
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(site_file.rds = "site_file.rds"))

# vimc inputs
vimc_input<- readRDS('vimc_input.rds') 
coverage_data<- vimc_input$coverage_input
le <- vimc_input$le
vimc_pop<- vimc_input$population_input_all_age
pop_single_yr<- vimc_input$population_input_single_yr

# site data
site_data <- readRDS('site_file.rds')
site <- extract_site(site_file = site_data,
                     site_name = site_name,
                     ur = ur)

scen<- scenario # doesn't work when scenario has the same name

if(scen == 'no-vaccination'){
  
  coverage_data<- coverage_data |>           # pull another projection for data table structure
    filter(country_code == iso3c) |>
    filter(scenario == 'malaria-r3-r4-default') |>
    mutate(coverage = 0)  |>
    mutate(scenario = scen,
           set_name = scen) 
  
}else{
  
  coverage_data<- coverage_data |>           # pull another projection for data table structure
    filter(country_code == iso3c) |>
    filter(scenario == scen)
}



# if quick run, set time length to 2035, if not set to 2100
if(quick_run == T){
  term_yr<- 2035
} else{
  term_yr<- 2100
}

# specify vaccine coverage based on forecast  ----------------------------------
site<- expand_intervention_coverage(site, 
                                    terminal_year = term_yr)

site<- update_coverage_values(site, 
                              coverage_data,
                              scenario)

# add in scenario variable which will be used to implement booster
saveRDS(site, 'vaccine_plot_input.rds')
message('formatting')

if(site$eir$eir[[1]] == 0){
  
  stop('Can not model this site beause PfPR EIR is equal to zero. Note this site/ urbanicity combination and exclude from future model runs.')
  
}

# pull parameters for this site ------------------------------------------------
params <- site::site_parameters(
  interventions = site$interventions,
  demography = site$demography,
  vectors = site$vectors,
  seasonality = site$seasonality,
  eir = site$eir$eir[1],
  burnin = 15,
  overrides = list(human_population = 100000)
)


# set age groups  --------------------------------------------------------------
if(quick_run== T) {
  
  year<- 365
  min_ages = c(0:5, 6,15,20) * year
  max_ages = c(1:6, 15,20,200) * year -1
  
}else{
  
  year<- 365
  min_ages = c(seq(0, 19, by= 1), seq(20, 90, by= 10)) * year
  max_ages = c(seq(1, 20, by= 1), seq(30, 100, by= 100)) * year -1
  
}  
params$clinical_incidence_rendering_min_ages = min_ages
params$clinical_incidence_rendering_max_ages = max_ages
params$severe_incidence_rendering_min_ages = min_ages
params$severe_incidence_rendering_max_ages = max_ages
params$age_group_rendering_min_ages = min_ages
params$age_group_rendering_max_ages = max_ages

# if this is a stochastic run, set parameter draw ------------------------------
if (parameter_draw > 0){
  
  params<- params |>
    set_parameter_draw(parameter_draw) |>
    set_equilibrium(init_EIR= params$init_EIR)
  
}

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
         parameter_draw = parameter_draw)


# save model runs somewhere
message('saving the model')
saveRDS(model, 'model_output.rds')

# process site --------------------------------------------------------------
# calculate rates --------------------------------------------------------------
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
le<- le |>
  filter(country_code == iso3c,
         year >= 2000) |>
  dplyr::group_by(year) |>
  tidyr::complete(age_from = c(1:100)) |>
  dplyr::ungroup() |>
  tidyr::fill(dplyr::all_of(names(le)), .direction = "down")

# fill years out (five year age groups)
le<- le |>
  dplyr::group_by(age_from) |>
  tidyr::complete(year = c(2000:2100))|>
  dplyr::ungroup() |>
  tidyr::fill(dplyr::all_of(names(le)), .direction = "down") |>
  rename(age_lower = age_from,
         remaining_yrs = value) |>
  select(year, age_lower, remaining_yrs) 

# calculate ylls_pp + dalys per person
dt<- merge(output, le, by = c('year', 'age_lower'), all.x = TRUE)

dt<- dt |>
  mutate(ylls_pp = mortality * remaining_yrs) |>
  mutate(dalys_pp = ylls_pp + yld_pp) |>
  select(-remaining_yrs)


# calculate counts  ------------------------------------------------------------
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
dt <- dt |>
  mutate(
    disease = 'Malaria',
    country = iso3c,
    country_name = countrycode::countrycode(
      sourcevar = iso3c,
      origin = 'iso3c',
      destination = 'country.name'),
    site_name = site_name,
    urban_rural = ur,
    scenario = scenario,
    description = description
  ) |>
  rename(age = age_lower,
         cohort_size = population) |>
  select(
    disease,
    year,
    age,
    country,
    country_name,
    site_name,
    urban_rural,
    scenario,
    description,
    cohort_size,
    cases,
    dalys,
    deaths,
    clinical,
    mortality,
    dalys_pp
  ) |>
  mutate(
    cases = if_else(is.na(cases), 0, cases),
    deaths = if_else(is.na(deaths), 0, deaths),
    dalys = if_else(is.na(dalys), 0, dalys),
    mortality = if_else(is.na(mortality), 0, mortality),
    clinical = if_else(is.na(clinical), 0, clinical),
    dalys = if_else(is.na(dalys), 0, dalys)
  )


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


