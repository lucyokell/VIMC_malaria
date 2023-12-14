# country functions ------------------------------------------------------------
aggregate_outputs<- function(dt, pop){
  dt<- data.table(dt)
  
  dt[, `:=` (
    cases = sum(cases),
    deaths = sum(deaths),
    dalys = sum(dalys)),
    by = c('age', 'year', 'scenario')]
  
  
  # remove cohort size, because for sites with some unmodelled locations, sum of cohort size != national population
  dt<- dt |> 
    select(-cohort_size)
  dt <- unique(dt, by = c('age', 'year', 'scenario'))
  pop<- pop |>
    rename(age = age_from,
           cohort_size = value) |>
    select(year, age, cohort_size)
  dt<- merge(dt, pop, by =c('age', 'year'))
  
  
  # calculate rates --------------------------------------------------------------
  dt[, `:=` (
    clinical = NULL,
    mortality = NULL,
    dalys_pp = NULL,
    site_name = iso3c,
    urban_rural = 'total'
  )]
  
  return(dt)
}

scale_cases<- function(dt, site_data){
  # # scale outputs based on cases from WMR from 2000-2020
  # # first sum cases by year (across all ages) in model output and compare
  pre_scale<- dt |>
    group_by(year) |>
    summarise(cases = sum(cases))
  
  #average site file cases across last three years
  site_file_cases<- data.table::data.table(site_data$cases_deaths[, c('year', 'wmr_cases')])
  site_file_cases<- site_file_cases[year >= 2018]
  average_value<- mean(site_file_cases$wmr_cases)
  
  # calculate ratio in comparison to year 2020 cases in output
  output_cases<- pre_scale |>
    filter(year == 2020) |>
    select(cases)
  
  ratio<- average_value/output_cases$cases
  
  # test this scaling with an extra column and review
  dt<- dt |>
    mutate(pre_scaled_cases = cases)
  
  dt<- dt |>
    mutate(cases = cases * ratio)
  
  dt<- dt|>
    mutate(clinical= cases/cohort_size,
           mortality = deaths/ cohort_size,
           dalys_pp = dalys/ cohort_size) |>
    select(-site_name, -urban_rural)
  
  return(dt)
}

remove_zero_eirs<- function(iso3c, sites, eirs){
  
  eirs<- data.table::data.table(eirs)
  eirs<- eirs[spp == 'pf' & eir == 0]
  remove<- eirs[, c('name_1', 'urban_rural')]
  
  sites<- data.table::data.table(sites)
  if(nrow(remove) > 0){
    for (i in 1:nrow(remove)){
      
      message(paste0('removing site ', i))
      
      sites<- sites[!(name_1== remove[i, name_1] & urban_rural == remove[i, urban_rural])]
      
    }
    
  } else{
    message('No zero eir sites to remove')
  }
  return(sites)
}

