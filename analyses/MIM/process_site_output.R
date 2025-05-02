# postprocess outputs at the site level for MIM poster
output_filepath<- 'J:/VIMC_malaria/archive/process_country/'

# process inputs for each country modelled  (must only be run once)  ----------
coverage<- read.csv('src/process_inputs/vimc_inputs/vaccine_coverage/coverage_202310gavi-7_malaria-rts3-rts4-bluesky.csv')
iso3cs<- unique(coverage$country_code)

source('J:/VIMC_malaria/src/postprocess_bulk/postprocessing_functions.R')
#source('diagnostic_report_functions.R')


# packages
library(site)
library(data.table)
library(dplyr)
library(malariasimulation)
library(openxlsx)
library(ggplot2)
library(tidyr)
library(tibble)
library(postie)
library(data.table)
library(countrycode)
library(scene)
library(scales)

postprocess_site_output<- function(iso3c){
  message(iso3c)
  # pull metadata for completed outputs for this country
  completed<- pull_most_recent_output({{iso3c}})
  completed<- completed |>
    filter(parameter_draw <= 50)

  # vimc inputs ----
  inputs<- completed_reports('process_inputs') |>
    filter(iso3c == {{iso3c}}) |>
    arrange(desc(date_time))
  inputs<- unique(inputs, by = 'iso3c')

  site_data<- readRDS(paste0('J:/VIMC_malaria/archive/process_inputs/', inputs$directory_name, '/merged_site_file.rds'))

  # pull prevalence data at 2019
  prev<- site_data$prevalence |>
    filter(year == 2019) |>
    rename(site_name = name_1,
           pfpr_2019 = pfpr) |>
    select(-country, -iso3c, -year)


  # pull model outputs for all baseline scenarios (as a separate input into intervention processing)
  outputs<- rbindlist(lapply(c(1:nrow(completed)), get_site_output, map = completed))

  no_vaccine<- outputs |>
    filter(scenario == 'no-vaccination') |>
    rename(cases_novax = cases,
           severe_novax = severe,
           dalys_novax = dalys,
           ylls_novax = ylls,
           deaths_novax = deaths,
           mortality_novax = mortality,
           clinical_novax = clinical,
           dalys_pp_novax = dalys_pp) |>
    select(-scenario, -cohort_size)

  intvn<- outputs |>
    filter(scenario != 'no-vaccination')

  intvn<- merge(intvn, no_vaccine, by = c('disease', 'year', 'age', 'country', 'country_name', 'site_name', 'urban_rural', 'parameter_draw'))

  dose_output<- rbindlist(lapply(c(1:nrow(completed)), get_site_dose_output, map = completed), fill = TRUE)
  dose_output<- dose_output |>
      filter(year %in% c(2000:2100)) |>
    select(-cohort_size, -doses_model, -doses, -n_model)

  intvn<- merge(intvn, dose_output, by = c('year', 'site_name', 'urban_rural', 'scenario', 'parameter_draw'))
  intvn<- merge(intvn, prev, by = c('site_name', 'urban_rural'))

  #merge in cohort size at corresponding age to calculate doses delivered

  four_dose<- intvn |> filter(scenario %like% 'rts4')
  three_dose<- intvn |> filter(!scenario %like% 'rts4')



saveRDS(intvn, paste0('outputs/site_outputs/', iso3c, '_site_output.rds'))
}


completed_reports<- function(report_name){


  meta <- orderly2::orderly_metadata_extract(name = report_name, extract = c('time', 'parameters'),  options = orderly2::orderly_search_options(allow_remote = TRUE))

  meta<- meta|>
    mutate(directory_name = id) |>
    tidyr::separate(col = id, into = c('date', 'time'), sep = '-')|>
    mutate(date= as.numeric(date)) |>
    mutate(date_time = as.numeric(paste0(date, time)))

  meta<- data.table(meta)
  meta<- meta[, index:= c(1:.N) ]


  unique(lapply(meta$parameters, names))
  nms <- names(meta$parameters[[1]])
  pars <- do.call("data.frame", setNames(lapply(nms, function(nm) sapply(meta$parameters, function(x) x[[nm]])), nms))
  pars<- data.table(pars)
  pars<- pars[, index:= c(1:.N)]

  meta<- meta |>
    select(directory_name, index, date_time)
  map<- merge(pars, meta, by = 'index')
  map<- map |>
    select(-index)

  return(map)
}

pull_most_recent_output<- function(iso3c){
  completed<- completed_reports('process_country') |>
    filter(iso3c == {{iso3c}},
           pfpr10 == TRUE,
           description == 'test_round2_changes',
           quick_run == FALSE) |>
    arrange(desc(date_time)) |>
    dplyr::distinct(iso3c, scenario, quick_run, parameter_draw, description, pfpr10, .keep_all = TRUE) |>
    arrange(iso3c, scenario, parameter_draw)


  return(completed)
}


get_site_output<- function(index, map){

  metadata<- map[ index,]

  directory<- metadata$directory_name
  draw<- metadata$parameter_draw

  message(directory)

  output<- readRDS(paste0(output_filepath, directory, '/outputs.rds'))             # get output file
  sites<- rbindlist(lapply(output$site_output, function(x) return(x$processed_output))) #pull out processed site_level output
  sites<- sites |>
    mutate(parameter_draw = draw)
  return(sites)
}

get_site_dose_output<- function(index, map){

  metadata<- map[ index,]

  directory<- metadata$directory_name
  draw<- metadata$parameter_draw

  message(directory)

  output<- readRDS(paste0(output_filepath, directory, '/outputs.rds'))             # get output file

  pull_doses<- function(x){

    dose_output<- x$doses
    site_output<- x$processed_output
    dose_output<-  dose_output |>
      mutate(site_name = unique(site_output$site_name),
             urban_rural= unique(site_output$urban_rural),
             scenario = unique(site_output$scenario))

    return(dose_output)

  }
  sites<- rbindlist(lapply(output$site_output, pull_doses))       # pull out processed site_level output
  sites<- sites |>
    mutate(parameter_draw = draw)
  return(sites)
}





lapply(iso3cs, postprocess_site_output)
