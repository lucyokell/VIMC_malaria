# summary checks
# purpose: write up quick checks of model outputs based on key issues + plot 
################################################################################


library(data.table)
library(vimcmalaria)
library(dplyr)
library(ggpubr)
library(ggforce)

coverage<- read.csv('src/process_inputs/vimc_inputs/vaccine_coverage/coverage_202409malaria-1_malaria-r3-r4-default.csv')
iso3cs<- unique(coverage$country_code)
source('src/diagnostics/diagnostic_report_functions.R')
outputs<- rbindlist(lapply(iso3cs , vimcmalaria::pull_postprocessed_output, description = 'fix_booster_coverage', quick_run= FALSE))

retrieve_final<- function(directory){
  
  output<- rbindlist(readRDS(paste0('archive/postprocessing/', directory, '/final_output.rds')))

  
  return(output)
}



retrieve_dose<- function(directory){
  
  output<- readRDS(paste0('archive/postprocessing/', directory, '/dose_output.rds'))
                   
  return(output)
}


final<- rbindlist(lapply(outputs$directory_name, retrieve_final))
dose<-  lapply(outputs$directory_name, retrieve_dose))
doses<- rbindlist(dose, fill= TRUE)


plot_incidence<- function(country_nm){


  message(country_nm)
test<- final |> filter(country_name == country_nm) |>
  group_by(scenario, year, run_id) |>
  summarise(cases = sum(cases),
            deaths= sum(deaths),
          dalys= sum(dalys),
          cohort_size= sum(cohort_size),
        .groups = 'keep') |>
  data.table()
  
  
  p1<- ggplot() +
    #geom_line(data = test[run_id == 0], mapping = aes(x= year, y= cases/cohort_size, color =scenario))  +
    geom_line(data = test[run_id!= 0], mapping = aes(x= year, y= cases/cohort_size, color = scenario, group = run_id), alpha = 0.02)  +
    stat_summary(data = test[run_id!= 0], mapping = aes(x= year, y= cases/cohort_size, color = scenario), fun.y = median, geom = "line") +
    plotting_theme +
    labs(title= 'Incidence over time',
          subtitle= country_nm)
  
  test<- final |> filter(country_name == country_nm) |>
  group_by(scenario, age, run_id) |>
  summarise(cases = sum(cases),
            deaths= sum(deaths),
          dalys= sum(dalys),
          cohort_size= sum(cohort_size),
        .groups = 'keep') |>
    data.table()
  
  
p2<- ggplot() +
  #geom_line(data = test[run_id == 0], mapping = aes(x= age, y= cases/cohort_size, color =scenario))  +
    geom_line(data = test[run_id!= 0], mapping = aes(x= age, y= cases/cohort_size, color = scenario, group = run_id), alpha = 0.02)  +
    stat_summary(data = test[run_id!= 0], mapping = aes(x= age, y= cases/cohort_size, color = scenario), fun.y = median, geom = "line") +
  plotting_theme +
  labs(title= 'Incidence over age',
        subtitle= country_nm)
  
  
  plots<- ggarrange(p1, p2, ncol =2 )

  print(plots)
}


plot_deaths<- function(country_nm){


  message(country_nm)
test<- final |> filter(country_name == country_nm) |>
  group_by(scenario, year, run_id) |>
  summarise(cases = sum(cases),
            deaths= sum(deaths),
          dalys= sum(dalys),
          cohort_size= sum(cohort_size),
        .groups = 'keep') |>
  data.table()
  
  
  p1<- ggplot() +
    #geom_line(data = test[run_id == 0], mapping = aes(x= year, y= deaths/cohort_size, color =scenario))  +
    geom_line(data = test[run_id!= 0], mapping = aes(x= year, y= deaths/cohort_size, color = scenario, group = run_id), alpha = 0.02)  +
    stat_summary(data = test[run_id!= 0], mapping = aes(x= year, y= deaths/cohort_size, color = scenario), fun.y = median, geom = "line") +
    #facet_wrap(~scenario) +
    plotting_theme +
    labs(title= 'Mortality over time',
          subtitle= country_nm)
  
  test<- final |> filter(country_name == country_nm) |>
  group_by(scenario, age, run_id) |>
  summarise(cases = sum(cases),
            deaths= sum(deaths),
          dalys= sum(dalys),
          cohort_size= sum(cohort_size),
        .groups = 'keep') |>
    data.table()
  
  
p2<- ggplot() +
 # geom_line(data = test[run_id == 0], mapping = aes(x= age, y= deaths/cohort_size, color =scenario))  +
    geom_line(data = test[run_id!= 0], mapping = aes(x= age, y= deaths/cohort_size, color = scenario, group = run_id), alpha = 0.02)  +
    stat_summary(data = test[run_id!= 0], mapping = aes(x= age, y= deaths/cohort_size, color = scenario), fun.y = median, geom = "line") +
  plotting_theme +
  labs(title= 'Mortality over age',
        subtitle= country_nm)
  
  
  plots<- ggarrange(p1, p2, ncol =2 )

  print(plots)
}

# run a few checks
# plot WMR cases over incidence over time
# pull WMR cases
pull_wmr<-function (iso3c){
  
  site<- readRDS(paste0('src/process_inputs/site_files/', iso3c, '_new_eir.rds'))
  wmr<- site$cases_deaths
  
  return(wmr)

}

wmr_cases<- rbindlist(lapply(iso3cs, pull_wmr))


wmr_diagnostic<- function(){
  
  outputs<- final[scenario == 'no-vaccination'] |>
    rename(iso3c = country) |>
    group_by(iso3c, year, run_id) |>
    summarise(cases= sum(cases), .groups = 'keep') |>
    data.table()
  
  p1<- ggplot()+
    geom_line(data = outputs[run_id == 0], mapping= aes(x= year, y= cases))+
    geom_line(data = wmr_cases, mapping  = aes(x= year, y= wmr_cases, color = 'red'))+
    geom_vline(xintercept = 2018, linetype= 'dotted')+
    geom_vline(xintercept = 2020, linetype= 'dotted')+
    facet_wrap_paginate(~iso3c, nrow = 4, ncol = 4, scales = 'free', page = 1) +
    # plotting_theme + 
    theme(legend.position = 'none')+
    labs(title = 'World Malaria Report cases vs. modelled outputs')
  
  p2<- ggplot()+
    geom_line(data = outputs, mapping= aes(x= year, y= cases))+
    geom_line(data = wmr_cases, mapping  = aes(x= year, y= wmr_cases, color = 'red'))+
    geom_vline(xintercept = 2018, linetype= 'dotted')+
    geom_vline(xintercept = 2020, linetype= 'dotted')+
    facet_wrap_paginate(~iso3c, nrow = 4, ncol = 4, scales = 'free', page = 2) +
    # plotting_theme + 
    theme(legend.position = 'none')+
    labs(title = 'World Malaria Report cases vs. modelled outputs')
  
  print(p1)
  print(p2)
}



pdf('wmr_diagnostic_new_runs.pdf', width = 12, height = 12)
wmr_diagnostic()
dev.off()



pdf('incidence_diagnostics_medians.pdf', width = 12, height = 12)
lapply(unique(final$country_name), plot_incidence)
dev.off()



pdf('mortality_diagnostic_medians.pdf', width = 12, height = 12)
lapply(unique(final$country_name), plot_deaths)
dev.off()



# sensitivity analysis



for(scen in unique(final$scenario)){

  message(scen)
  subset<- final |> filter(scenario == scen) |>
    filter(run_id== 0) |>
    select(-scenario,
            -run_id)

  write.csv(subset, paste0('montagu/central-burden-est-', scen, '.csv'), row.names = FALSE)

}
