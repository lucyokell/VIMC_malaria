#' update vaccine coverage based on VIMC inputs from Montagu
#'
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
  
  site$interventions<- intvns # keep in mind these booster values will not be fed into the model
  
  if(scenario == 'malaria-rts3-rts4-bluesky'){
    
    site$interventions$scenario_type<- 'bluesky'
  }
  
  if(scenario_name == 'malaria-r3-r4-default' | scenario_name == 'malaria-rts3-rts4-default'){
    
    site$interventions$scenario_type<- 'routine'
    
  }

  if (scenario == 'malaria-rts3-default' | scenario == 'malaria-r3-default' | scenario == 'malaria-rts3-bluesky'){
    
    site$interventions$scenario_type<- 'single_dose'
  }
  
  return(site)
}

