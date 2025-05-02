### calibrate Ethiopia  
##  purpose calibrate Ethiopia sites to lower baseline transmission intensity


library(vimcmalaria)
library(data.table)
library(dplyr)
library(cali)
scenario = 'no-vaccination'
parameter_draw = 0
gfa = FALSE
quick_run = FALSE
iso3c = 'ETH'
calibrate_eth_sites<- function(site_name){
  
  site_data<- readRDS('src/process_inputs/site_files/ETH_new_eir.rds')
  site_data$interventions$irs_cov = 0
  
  completed<- vimcmalaria::completed_reports('process_inputs') |>
    filter(iso3c== 'ETH')
  vimc_input<- readRDS(paste0('archive/process_inputs/', completed$directory_name, '/vimc_input.rds'))  
  
  
  # vimc inputs ----
  coverage_data<- vimc_input$coverage_input
  le <- vimc_input$le
  vimc_pop<- vimc_input$population_input_all_age
  pop_single_yr<- vimc_input$population_input_single_yr
  
  # make a map of input parameters for site function
  site_df<- remove_zero_eirs(iso3c, site_data)
  map<- vimcmalaria::make_analysis_map(site_df, site_data, test = FALSE)
  
  
  
  message('parameterizing')
  # site data
  site <- vimcmalaria::extract_site(site_file = site_data,
                       site_name = site_name,
                       ur = 'both')
  
  
  run_params<- vimcmalaria::pull_age_groups_time_horizon(quick_run)
  
  # specify vaccine coverage based on forecast  ----------------------------------
  site<- vimcmalaria::expand_intervention_coverage(site,
                                      terminal_year = run_params$term_yr)
  site<- vimcmalaria::update_coverage_values(site,
                                iso3c = iso3c,
                                coverage_data,
                                scenario_name = scenario)
  
  # check the site has a non-zero EIR
  vimcmalaria::check_eir(site)
  
  # pull parameters for this site ------------------------------------------------
  params <- site::site_parameters(
    interventions = site$interventions,
    demography = site$demography,
    vectors = site$vectors,
    seasonality = site$seasonality,
    eir = site$eir$eir[1],
    burnin = run_params$burnin,
    overrides = list(human_population = run_params$pop_val)
  )
  
  
  # set age groups
  params$clinical_incidence_rendering_min_ages = run_params$min_ages
  params$clinical_incidence_rendering_max_ages = run_params$max_ages
  params$severe_incidence_rendering_min_ages = run_params$min_ages
  params$severe_incidence_rendering_max_ages = run_params$max_ages
  params$age_group_rendering_min_ages = run_params$min_ages
  params$age_group_rendering_max_ages = run_params$max_ages
  
  
  cal_params<- vimcmalaria::recalibrate(params,
                           site_name= site_name,
                           site_dt = site_data)
  
  
  saveRDS(cal_params, paste0('analyses/ethiopia/calibrations/calibrated_site', site_name, '.rds'))
  
}

