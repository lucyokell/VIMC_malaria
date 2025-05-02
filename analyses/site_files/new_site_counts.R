# generate site counts
coverage<- read.csv('src/process_inputs/vimc_inputs/vaccine_coverage/coverage_202409malaria-1_malaria-r3-r4-default.csv')
iso3cs<- unique(coverage$country_code)


grid <- expand.grid(
  "iso3c" = as.character(iso3cs),
  "scenario" = c("malaria-r3-r4-default", "no-vaccination", "malaria-rts3-rts4-default")
)

#' Pull the moderate to high transmission sites in a site file
#' @param iso3c   iso3c
#' @return sites component of site file with only sites that will be modelled
#' @export
pull_site_counts<- function(iso, scen){

  site_data<- site::fetch_site(iso3c = iso)

  site_data$prevalence<- site_data$prevalence |>
    dplyr::filter(year == 2024) |>
    dplyr::mutate(run_model = ifelse(pfpr > 0.10, TRUE, FALSE)) |>
    mutate(run_model = ifelse(name_1 == 'Toliary', TRUE, run_model)) |> #hardcoded exceptions
    mutate(run_model = ifelse(name_1 %like% 'Gambela', TRUE, run_model)) |>
    mutate(run_model = ifelse(name_1 %like% 'Bay', TRUE, run_model)) |>
    mutate(run_model = ifelse(name_1 %like% 'Nouakchott', TRUE, run_model)) |>
    mutate(run_model = ifelse(name_1 %like% 'Bolama', TRUE, run_model)) |>
    mutate(run_model = ifelse(name_1 == 'South Darfur', TRUE, run_model)) |>
    mutate(run_model = ifelse(name_1 == 'West Kurdufan', TRUE, run_model))

  if(scen == 'no-vaccination'){site_data$prevalence$run_model<- TRUE}


return(data.frame('iso3c'= iso, 'scenario'= scen, 'site_count'= nrow(site_data$prevalence |> filter(run_model == TRUE))))

}
testing <- lapply(1:nrow(grid), function(x) {
  message(x)  # Print the index to track progress
  sub <- grid[x, ]  # Access the x-th row of the data frame
  counts <- pull_site_counts(iso= as.character(sub$iso3c), scen =sub$scenario)  # Use the row's values for the columns
  return(counts)
})


counts<- rbindlist(testing)
counts<- counts |>
  rename(site_number = site_count)
saveRDS(counts, 'site_counts.rds')
