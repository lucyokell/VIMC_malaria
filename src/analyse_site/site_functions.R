# site functions  --------------------------------------------------------------
# utility ----------------------------------------------------------------------
#' Extract a single site-input from a country site file
#' @param site_file  Country site file
#' @param site_name  name of site to extract
#' @param urbanicity urbanicity of site to extract
#' @return Single site
#' @export
extract_site <- function(site_file, site_name, ur){
  
  sites<- data.table(site_file$sites)
  Encoding(sites$name_1) <- "UTF-8"
  
  sites$name_1<- iconv(sites$name_1, from="UTF-8", to="ASCII//TRANSLIT")
  
  index_site <- sites[name_1== site_name & urban_rural== ur]
  
  to_mod <- c("sites", "interventions", "pyrethroid_resistance", "population",
              "vectors", "seasonality", "prevalence", "eir")
  
  site <- site_file
  
  for(level in to_mod){
    mod<- site[[level]]
    Encoding(mod$name_1) <- "UTF-8"
    
    mod$name_1<- iconv(mod$name_1, from="UTF-8", to="ASCII//TRANSLIT")
    
    mc <- intersect(colnames(index_site), colnames(mod))
    site[[level]] <- dplyr::left_join(index_site, mod, by = mc)
  }
  
  return(site)
}



# parameterizing  --------------------------------------------------------------
pull_age_groups_time_horizon<- function(quick_run){
  
  year<- 365
  burnin<- 15
  
  if(quick_run == T){
    
    term_yr<- 2035
    pop_val<- 10000
    
    min_ages = c(0:5, 6,15,20) * year
    max_ages = c(1:6, 15,20,200) * year -1
    
  } else{
    
    term_yr<- 2100
    pop_val<- 100000
    
    min_ages = c(seq(0, 19, by= 1), seq(20, 90, by= 10)) * year
    max_ages = c(seq(1, 20, by= 1), seq(30, 100, by= 100)) * year -1
    
  }
  
  return(list('term_yr' = term_yr, 'pop_val' = pop_val, 'min_ages'= min_ages, 'max_ages' = max_ages, 'burnin' = burnin))
} 

check_eir<- function(site){
  if(site$eir$eir[[1]] == 0){
    
    stop('Can not model this site beause PfPR EIR is equal to zero. Note this site/ urbanicity combination and exclude from future model runs.')
    
  }
}

parameterize_stochastic_run<- function(parameter_draw){
  
  if (parameter_draw > 0){
    
    params<- params |>
      set_parameter_draw(parameter_draw) |>
      set_equilibrium(init_EIR= params$init_EIR)
    
  }
  
  return(params)
}

#' update vaccine coverage based on VIMC inputs from Montagu
#' @param   site             site data file
#' @param   coverage_data    VIMC vaccine forecast for site of interest
#' @returns site file with additional variables 'rtss_coverage', 'rtss_booster_coverage', 'r21_coverage', 'r21_booster_coverage'
update_coverage_values<- function(site, coverage_data, scenario_name){
  
  if(scenario_name == 'no-vaccination'){
    
    coverage_data<- coverage_data |>           # pull another projection for data table structure and fill with zeroes
      filter(country_code == iso3c) |>
      filter(scenario == 'malaria-r3-r4-default') |>
      mutate(coverage = 0)  |>
      mutate(scenario = 'no-vaccination') 
    
  }else{
    
    coverage_data<- coverage_data |>           
      filter(country_code == iso3c) |>
      filter(scenario == scenario_name)
  }
  
  dt<- coverage_data |>
    rename(vaccine_name = vaccine) |>
    data.table()
  
  # add identifying type column for vaccine
  dt[vaccine_name %like% 'RTS', vaccine := 'RTS,S']
  dt[is.na(vaccine), vaccine := 'R21']
  
  vaccine_val<- unique(dt$vaccine)
  
  if (length(vaccine_val) > 1){ stop('Can only implement one type of vaccine at a time. Check vaccine inputs.') }
  
  dt<- dcast(data.table(dt), 
             year + vaccine ~ vaccine_name, 
             value.var= 'coverage')
  
  # if columns for other vaccines or doses are empty, fill them ----------------
  columns_to_check <- c("R3", "R4", "RTS3", "RTS4")
  missing_columns <- setdiff(columns_to_check, names(dt))
  
  dt <- dt |>
    add_column(!!!setNames(rep(0, length(missing_columns)), 
                           missing_columns))
  
  dt <- dt |>
    rename(rtss_coverage = RTS3,
           rtss_booster_coverage = RTS4,
           r21_coverage = R3,
           r21_booster_coverage = R4) 
  
  # transform booster coverage into value per person according to coverage in the preceding year
  for (yr in unique(dt$year)){
    
    dt[year== yr & rtss_coverage!= 0 & rtss_booster_coverage!= 0,
       rtss_booster_coverage := rtss_booster_coverage / dt[year == yr- 1, rtss_coverage]]
    
    dt[year== yr & r21_coverage!= 0 & r21_booster_coverage!= 0,
       r21_booster_coverage := r21_booster_coverage / dt[year == yr- 1, r21_coverage]]
  }
  
  intvns<- data.table(merge(site$interventions, dt, by = 'year', all.x= T))
  
  intvns[is.na(rtss_coverage), "rtss_coverage" := 0]
  intvns[is.na(rtss_booster_coverage), "rtss_booster_coverage" := 0]
  intvns[is.na(r21_coverage), "r21_coverage" := 0]
  intvns[is.na(r21_booster_coverage), "r21_booster_coverage" := 0]  
  intvns[is.na(vaccine), vaccine := vaccine_val]
  
  site$interventions<- intvns 
  
  return(site)
}


#' expand intervention years out to terminal year of forecast using scene package
#' @param   site             site data file
#' @param   terminal_year    terminal year of forecast
#' @returns site file with extrapolated coverage values out to terminal year
expand_intervention_coverage<- function(site, terminal_year){
  
  # first set terminal year to terminal year of forecast
  group_var <- names(site$sites)
  
  first_yr<- max(site$interventions$year)  +1              # first year in site file
  itn_yr<- first_yr- 3                                   # last year to carry over for ITN usage and model input (3 year cycle)
  
  site$interventions <- site$interventions |> 
    scene::expand_interventions(max_year = terminal_year,
                                group_var = group_var)
  
  
  for (yr in c(first_yr:terminal_year)){
    
    comparator<- site$interventions |>
      filter(year == yr - 3)
    
    intvns <-   data.table(site$interventions)
    intvns[year == yr, `:=` (itn_use = comparator$itn_use,
                             itn_input_dist = comparator$itn_input_dist)]
    site$interventions <- intvns
    
  }
  
  
  site$interventions <- site$interventions |>
    scene::fill_extrapolate(group_var = group_var)
  
  return(site)
}

# postprocessing  --------------------------------------------------------------
expand_life_expectancy<- function(le){
  
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
  
  return(le)
}

vimc_postprocess<- function(output, le, site_data, vimc_pop, pop_single_yr){
  
  # fill rates out to single year age groups
  output<- output |>
    dplyr::group_by(t) |>
    tidyr::complete(age_lower = c(1:100)) |>
    select(-age_upper) |>
    dplyr::ungroup() |>
    tidyr::fill(clinical, severe, mortality, yld_pp, yll_pp, dalys_pp, .direction = 'down') |>
    select(-prop_n, -n, -yll_pp, -dalys_pp) |>
    rename(year = t)
  
  # merge in inputs for expected remaining years of life (to calculate YLLs)
  le<- expand_life_expectancy(le)
  
  # calculate ylls_pp + dalys per person
  dt<- merge(output, le, by = c('year', 'age_lower'), all.x = TRUE)
  
  # recalculate YLLs and DALYs based on country-specific life expectancy  ------
  dt<- dt |>
    mutate(ylls_pp = mortality * remaining_yrs) |>
    mutate(dalys_pp = ylls_pp + yld_pp) |>
    select(-remaining_yrs)
  
  # scale site population to be consistent with VIMC population ----------------
  populations<- scale_population(site_data, vimc_pop, pop_single_yr)
  
  # merge in site population + prop_n
  dt<- merge(dt, populations$site_population, by= 'year')
  dt<- merge(dt, populations$age_proportions, by = c('year', 'age_lower'))
  
  # calculate counts for entire time period --------------------------------------
  dt<- dt |>
    mutate(
      cases = round(clinical * vimc_site_population * prop_n),
      deaths = round(mortality * vimc_site_population * prop_n),
      dalys = round(dalys_pp * vimc_site_population * prop_n),
      population = round(vimc_site_population * prop_n)) |>
    select(-prop_n)
  
  return(dt)
  
}

scale_population<- function(site_data, vimc_pop, pop_single_yr){
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
  
  age_proportions<- merge(pop_single_yr, national_pop, by = c('year'))
  age_proportions <- age_proportions |>
    mutate(prop_n = value/ national_pop) |>
    select(year, age_from, age_to, prop_n) |>
    rename(age_lower = age_from)
  
  return(list('site_population' = site_pop, 'age_proportions' = age_proportions))
  
}

format_outputs<- function(dt){
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
  
  return(dt)
}


# save plotting data  ----------------------------------------------------------
pull_plotting_data<- function(scenario){
  
  if(scenario!="no-vaccination") {
    doses_per_year <- pull_doses_output(raw_output, dt)
    
    saveRDS(doses_per_year, 'doses_per_year.rds')
  } else{
    doses_per_year<- 0
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
  # save outputs for plotting
  plotting_inputs<- list('vaccine_plot_input' = vaccine_plot_input,
                         'raw_model_output' = raw_output,
                         'prevalence_per_year' = prev,
                         'doses_per_year' = doses_per_year)
  
  return(plotting_inputs)
}


