site_counts<- rbindlist(lapply(1:nrow(full_map), pull_site_numbers, map= full_map))


#' Pull the number of sites to model for country job
#' @param index observation index for parameter map
#' @param map   parameter map
#' @export
pull_site_numbers<- function(index, map){

  map<- map[index,]
  iso3c<- map$iso3c
  scenario<- map$scenario

  site<- readRDS(paste0('src/process_inputs/site_files/', iso3c, '_new_eir.rds'))

  if(scenario == 'no-vaccination'){
    site_number<- nrow(site$sites)

  } else{

    site$prevalence <- site$prevalence |>
      mutate(run_model = ifelse(pfpr > 0.10, TRUE, FALSE))|>
      mutate(run_model = ifelse(name_1 == 'Toliary', TRUE, run_model))|>   # Madagascar exception
      mutate(run_model = ifelse(name_1 %like% 'Gambela', TRUE, run_model))|> # Ethiopia exception
      mutate(run_model = ifelse(name_1 == 'South Darfur', TRUE, run_model)) |> # Sudan exceptions
      mutate(run_model = ifelse(name_1 == 'West Kurdufan', TRUE, run_model))

    site_number<- nrow(site$prevalence |> filter(run_model == TRUE, year == 2019))
  }

  return(data.table('iso3c' = iso3c, 'scenario' = scenario, 'site_number' = site_number))
}
