# postprocess in bulk
orderly2::orderly_parameters(iso3c = 'NGA',
                             description =  'new_site_files',
                             quick_run = FALSE)


# packages  --------------------------------------------------------------------
library(site)
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
library(wesanderson)
library(vimcmalaria)

# read in dependencies  --------------------------------------------------------
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(vimc_input.rds = "vimc_input.rds"))
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(site_file.rds = "site_file.rds"))

vimc_input<- readRDS('vimc_input.rds')
site_data <- readRDS('site_file.rds')
coverage_data<- vimc_input$coverage_input
le <- vimc_input$le
vimc_pop<- vimc_input$population_input_all_age
pop_single_yr<- vimc_input$population_input_single_yr
pop_data<- vimc_input$population_input_all_age

# if this is ivory coast, pull WMR cases and deaths from a flat file for scaling portion
# this will be fixed in upcoming versions of the site file
wmr<- read.xlsx('CIV_cases_deaths.xlsx')|> select(year, cases_mean) |> rename(wmr_cases = cases_mean)
wmr<- list('cases_deaths' = wmr )

# pull metadata for completed outputs for this country
completed<- completed_reports('process_country') |>
  dplyr::filter(iso3c == {{iso3c}},
                description == {{description}},
                quick_run == {{quick_run}}) |>
  dplyr::arrange(dplyr::desc(date_time)) |>
  dplyr::distinct(iso3c, scenario, quick_run, parameter_draw, description, .keep_all = TRUE) |>
  dplyr::arrange(iso3c, scenario, parameter_draw)

max_draw<- max(completed[scenario == 'no-vaccination', parameter_draw])

# vimc inputs
inputs<- completed_reports('process_inputs') |>
  filter(iso3c == {{iso3c}}) |>
  arrange(desc(date_time))
inputs<- unique(inputs, by = 'iso3c')

# pull model outputs for scaling to WMR cases: no-vaccination draw 0
scaling_filepath<- completed |> filter(parameter_draw == 0,
                                       scenario == 'no-vaccination')
scaling<- readRDS(paste0('J:/VIMC/VIMC_malaria/archive/process_country/', scaling_filepath$directory_name, '/outputs.rds'))
scaling<- scaling$country_output


# postprocess by draw ----------------------------------------------------------
final_postprocessing<- function(draw){

  message(paste0('postprocessing draw ', draw))

  bl_filepaths<- completed |> filter(parameter_draw == draw,
                                     scenario == 'no-vaccination')
  intvn_filepaths<- completed |> filter(parameter_draw == draw)

  # pull model outputs for all scenarios
  intvn<- rbindlist(lapply(c(1:nrow(intvn_filepaths)), get_site_output, map = intvn_filepaths, output_filepath = 'J:/VIMC/VIMC_malaria/archive/process_country/' ))
  
  # pull model outputs for all baseline scenarios (as a separate input into intervention processing)
  bl<- rbindlist(lapply(c(1:nrow(bl_filepaths)), get_site_output, map = bl_filepaths, output_filepath = 'J:/VIMC/VIMC_malaria/archive/process_country/'))

  # commenting out as now modelling introduction in all sites regardless of transmission intensity
  message('adding low transmission sites')
  
  low<- pull_low_transmission_sites(iso3c, site_data, bl, threshold = 0.10, threshold = .10) # pull sites below PFPR threshold
  intvn<- subset_high_transmission_sites(intvn, threshold = .10, site_data) #subset to sites meeting PFPR threshold
  print(nrow(low))
  intvn<- append_low_transmission_sites(low_transmission = low, intvn)


  message('aggregating')
  dt<- aggregate_outputs(intvn, pop_single_yr)

  message('scaling cases')
  output<- scale_cases(output, scaling_data= scaling,  site_data = site_data)


  # scale cases based on difference between site file PAR and VIMC PAR
  message('scaling PAR')
  processed_output<- scale_par(output, iso3c= {{iso3c}})

  message('formatting')
  # format and save (
  processed_output<- processed_output |>
    mutate(cases = round(cases),
           deaths = round(deaths),
           yll= round(ylls),
           dalys = round(dalys),
           cohort_size = cohort_size) |>
    rename(run_id = parameter_draw) |>
    select(-model_proportion_risk,
           -proportion_risk,
           -scaling_ratio,
           -prop_severe,
           -prop_deaths,
           -dalys_pp,
           -ylls,
           -life_expectancy) |>
    select(run_id,
           scenario,
           disease,
           year,
           age,
           country,
           country_name,
           cohort_size,
           cases,
           dalys,
           deaths,
           yll) |>
    mutate(pfpr_threshold = .10)



  return(processed_output)

}

# save -------------------------------------------------------------------------
max_draw<- 50
outputs<- lapply(c(0:max_draw), final_postprocessing)
saveRDS(outputs, 'final_output.rds')



dose_postprocessing<- function(draw){
  intvn_filepaths<- completed |> filter(parameter_draw == draw)
  
  # pull model outputs for all scenarios
  raw<- bind_rows(lapply(c(1:nrow(intvn_filepaths)), get_raw_output, map = intvn_filepaths, output_filepath = 'J:/VIMC/VIMC_malaria/archive/process_country/' ))
  raw<- raw |>
    mutate(site = paste0(site_name, '_', urban_rural, '_', scenario),
           site_ur = paste0(site_name, '_', urban_rural))
  
  ids<- unique(raw$site)
  
  case_output<- rbindlist(lapply(ids, site_postprocessing, dt = raw)) 
             
  case_output<- case_output |>
             rename(year = t,
                    age = age_lower) |>
             mutate(cases = round(clinical * n),
                    deaths = round(mortality * n)) |>
             group_by(site, year, site_ur) |>
             summarise(cases = sum(cases),
                       deaths= sum(deaths),
                       .groups = 'keep') 
           
           #calculate cases averted by year and site
           novax<- case_output|> 
             filter(site %like%'no-vaccination') |>
             rename(cases_novax = cases,
                    deaths_novax = deaths) |>
             ungroup() |>
             select(-site)
         
           vax<- case_output |> filter(!site %like% 'no-vaccination')
           
           averted<- merge(vax, novax, by = c('year', 'site_ur')) |>
             mutate(cases_averted = cases_novax- cases,
                    deaths_averted = deaths_novax - deaths)
         
             
           doses<- raw |>
             select(timestep, n_pev_epi_dose_1, n_pev_epi_dose_2, n_pev_epi_dose_3, n_pev_epi_booster_1, scenario, site) |>
             mutate(year = as.integer(timestep/365)) |>
             group_by(site, year, scenario) |>
             summarise(n_pev_epi_dose_1 = sum(n_pev_epi_dose_1),
                       n_pev_epi_dose_2 = sum(n_pev_epi_dose_2),
                       n_pev_epi_dose_3 = sum(n_pev_epi_dose_3),
                       n_pev_epi_booster_1 = sum(n_pev_epi_booster_1),
                       .groups = 'keep') |>
               mutate(fvp = n_pev_epi_dose_3) |> # to align with Nora and HIllary's paper
               mutate(doses_total = n_pev_epi_dose_1 + n_pev_epi_dose_2 +n_pev_epi_dose_3 + n_pev_epi_booster_1) |>
             select(-n_pev_epi_dose_1, -n_pev_epi_dose_2, -n_pev_epi_dose_3, -n_pev_epi_booster_1) |>
             mutate(year = year +1998) #there was misalignment here
             
           dose_output<- merge(averted, doses, by = c('year', 'site')) 
           dose_output <- dose_output |>
             mutate(parameter_draw = draw)
         
  return(dose_output)
  
}

dose_output<- lapply(c(0:max_draw), dose_postprocessing)

saveRDS(dose_output, 'dose_output.rds')

message('done with postprocessing')



