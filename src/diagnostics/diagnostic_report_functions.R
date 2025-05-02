# Diagnostic report functions
# load plotting theme
library(ggplot2)



plotting_theme<- theme_bw(base_size = 14) +
  theme( legend.position = 'bottom',
         strip.text.x = element_text(size = rel(0.8)),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
         #text = element_text(family= 'Calibri'),
         axis.ticks.y= element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())


compile_diagnostics<- function(description, date_time){

  completed<- completed_reports('diagnostics')
  completed<- completed |>
    filter(description == {{description}})|>
    filter(date_time >= {{date_time}}) |>
    arrange(desc(date_time)) |>
    dplyr::distinct(iso3c, description, .keep_all = TRUE) |>
    arrange(iso3c, description)

  copy_report<- function(index, map){

    message(index)
    map<- map[ index,]
    directory_name<- map$directory_name
    iso3c<- map$iso3c
    file.copy(from= paste0('J:/VIMC_malaria/archive/diagnostics/', directory_name, '/diagnostic_report_', iso3c, '.html'),
              to= paste0('J:/VIMC_malaria/diagnostics/', iso3c, '.html'))
  }

  lapply(c(1:nrow(completed)), copy_report, map = completed)
}


prepare_averted_outputs<- function(output){
    output<- data.table(output)

  # by year
  dt<- output |>
    filter(!is.na(scenario)) |>
    group_by(scenario, parameter_draw, year, country_name) |>
    summarise(cases= sum(cases),
              deaths = sum(deaths),
              cohort_size = sum(cohort_size),
              dalys = sum(dalys), .groups = 'keep')  |>
    mutate(cases = ifelse(is.na(cases), 0, cases),
           deaths = ifelse(is.na(deaths), 0, deaths),
           dalys = ifelse(is.na(dalys), 0, dalys))
  novax<- dt |>
    filter(scenario == 'no-vaccination') |>
    rename(cases_novax = cases,
           deaths_novax = deaths,
           dalys_novax = dalys) |>
    ungroup() |>
    select(-scenario)
  dt<- merge(dt, novax, by = c('year', 'parameter_draw'))
  dt<- dt |>
    mutate(cases_averted = cases_novax - cases,
           deaths_averted = deaths_novax - deaths,
           dalys_averted = dalys_novax - dalys) |>
    filter(scenario != 'no-vaccination')
  by_year<- copy(dt)

  # by age
  dt<- output |>
    filter(!is.na(scenario)) |>
    group_by(scenario, parameter_draw, age, country_name) |>
    summarise(cases= sum(cases),
              deaths = sum(deaths),
              cohort_size = sum(cohort_size),
              dalys = sum(dalys), .groups = 'keep')  |>
    mutate(cases = ifelse(is.na(cases), 0, cases),
           deaths = ifelse(is.na(deaths), 0, deaths),
           dalys = ifelse(is.na(dalys), 0, dalys))
  novax<- dt |>
    filter(scenario == 'no-vaccination') |>
    rename(cases_novax = cases,
           deaths_novax = deaths,
           dalys_novax = dalys) |>
    ungroup() |>
    select(-scenario)
  dt<- merge(dt, novax, by = c('age', 'parameter_draw'))
  dt<- dt |>
    mutate(cases_averted = cases_novax - cases,
           deaths_averted = deaths_novax - deaths,
           dalys_averted = dalys_novax - dalys) |>
    filter(scenario != 'no-vaccination')

  by_age<- copy(dt)
  return(list('by_year' = by_year, 'by_age' = by_age))

}
prepare_burden_outputs<- function(output){
  dt<- output |>
    group_by(scenario, year, country_name, parameter_draw) |>
    select(country_name, scenario, year, cases, deaths, dalys, cohort_size, parameter_draw) |>
    summarise(cases = sum(cases),
              deaths= sum(deaths),
              dalys = sum(dalys),
              cohort_size = sum(cohort_size),
              .groups= 'keep') |>
    mutate(cases = ifelse(is.na(cases), 0, cases),
           deaths = ifelse(is.na(deaths), 0, deaths),
           dalys = ifelse(is.na(dalys), 0, dalys))

  by_year<- copy(dt)

  dt<- output |>
    group_by(scenario, age, country_name, parameter_draw) |>
    select(country_name, scenario, age, cases, deaths, dalys, cohort_size, parameter_draw) |>
    summarise(cases = sum(cases),
              deaths= sum(deaths),
              dalys = sum(dalys),
              cohort_size = sum(cohort_size),
              .groups= 'keep') |>
    mutate(cases = ifelse(is.na(cases), 0, cases),
           deaths = ifelse(is.na(deaths), 0, deaths),
           dalys = ifelse(is.na(dalys), 0, dalys))


  by_age<- copy(dt)
  return(list('by_year' = by_year, 'by_age' = by_age))
}

incidence_over_time<- function(output, scenario, site_data, coverage_data){

  intro_yr<- min(coverage_data[scenario == {{scenario}}& coverage> 0, year])
  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  stochastic_run<- output |>
    filter(parameter_draw != 0)

  median_run<- output |>
    group_by(year, scenario) |>
    summarise(cases = mean(cases),
              deaths = mean(deaths),
              dalys = mean(dalys),
              .groups = 'keep')

  p<-   ggplot()+
    geom_line(data = median_run, mapping = aes(x= year, y= cases, color =scenario))  +
    geom_line(data = stochastic_run, mapping = aes(x= year, y= cases, color = scenario, group = parameter_draw), alpha = 0.1)  +
    geom_line(data= site_data$cases_deaths, mapping = aes(x= year, y= wmr_cases), color= 'darkgreen')+
    geom_vline(xintercept= intro_yr, linetype= "dotted") +
    labs(x= 'Time (in years)', y= 'Cases',
         title= paste0('Cases over time'),
         subtitle = {{scenario}},
         color= 'Scenario', fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Darjeeling1', n= 2)) +
    scale_fill_manual(values= wes_palette('Darjeeling1', n= 2))  +
    plotting_theme

  return(p)
}
incidence_over_age<- function(output, scenario, coverage_data){
  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  stochastic_run<- output |>
    filter(parameter_draw != 0)
  median_run<- output |>
    group_by(age, scenario) |>
    summarise(cases = mean(cases),
              deaths = mean(deaths),
              dalys = mean(dalys),
              .groups = 'keep')
  p<-   ggplot()+
    geom_line(data = median_run, mapping = aes(x= age, y= cases, color =scenario))  +
    geom_line(data = stochastic_run, mapping = aes(x= age, y= cases, color = scenario,
                                                   group = parameter_draw), alpha = 0.05)  +
    #facet_wrap_paginate(~scenario) +
    labs(x= 'Time (in years)', y= 'Cases',
         title= paste0('Cases over age'),
         subtitle = unique(median_run$country_name),
         color= 'Scenario', fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Darjeeling1', n= 2)) +
    scale_fill_manual(values= wes_palette('Darjeeling1', n= 2))  +
    plotting_theme

  return(p)
}
incidence_rate_over_time<- function(output, scenario, coverage_data){

  intro_yr<- min(coverage_data[scenario == {{scenario}} & coverage> 0, year])
  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  stochastic_run<- output |>
    filter(parameter_draw != 0)
  median_run<- output |>
    group_by(year, scenario) |>
    summarise(cases = mean(cases),
              deaths = mean(deaths),
              dalys = mean(dalys),
              cohort_size = mean(cohort_size),
              .groups = 'keep')

  p<-   ggplot()+
    geom_line(data = median_run, mapping = aes(x= year, y= cases/cohort_size, color =scenario))  +
    geom_line(data = stochastic_run, mapping = aes(x= year, y= cases/cohort_size, color = scenario, group = parameter_draw), alpha = 0.05)  +
    geom_vline(xintercept= intro_yr, linetype= "dotted") +
    labs(x= 'Time (in years)', y= 'Incidence',
         title= paste0('Incidence over time'),
         subtitle = {{scenario}},
         color= 'Scenario', fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Darjeeling1', n= 2)) +
    scale_fill_manual(values= wes_palette('Darjeeling1', n= 2))  +
    plotting_theme
  return(p)
}
incidence_rate_over_age<- function(output, scenario, coverage_data){
  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  stochastic_run<- output |>
    filter(parameter_draw != 0)
  median_run<- output |>
    group_by(age, scenario) |>
    summarise(cases = mean(cases),
              deaths = mean(deaths),
              dalys = mean(dalys),
              cohort_size= mean(cohort_size),
              .groups = 'keep')
  p<-   ggplot()+
    geom_line(data = median_run, mapping = aes(x= age, y= cases/cohort_size, color =scenario))  +
    geom_line(data = stochastic_run, mapping = aes(x= age, y= cases/cohort_size, color = scenario, group = parameter_draw), alpha = 0.05)  +
    #facet_wrap_paginate(~scenario) +
    labs(x= 'Time (in years)', y= 'Incidence',
         title= paste0('Incidence over age'),
         subtitle = unique(median_run$country_name),
         color= 'Scenario', fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Darjeeling1', n= 2)) +
    scale_fill_manual(values= wes_palette('Darjeeling1', n= 2))  +
    plotting_theme

  return(p)
}
plot_intvn_coverage<- function(site_data){
  p<- plot_interventions_combined(
    interventions = site_data$interventions,
    population = site_data$population$population_total,
    group_var = c("country", "name_1"),
    include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov",'irs_cov', "pmc_cov", "rtss_cov"),
    labels = c("ITN usage", "ITN model input", "Treatment","SMC", "IRS", "PMC", "RTSS")
  )

  return(p)
}
plot_vaccine_coverage<- function(coverage_data,country_name){

  coverage_data<- coverage_data

  p<- ggplot(data= coverage_data, mapping= aes(x= year, y= coverage, color= vaccine))+
    geom_line()+
    facet_wrap(~scenario, ncol= 4) +
    plotting_theme +
    labs(y= 'Coverage value', x= 'Year', title= 'Vaccine coverage over time (VIMC input)',
         subtitle = paste0('Country: ', unique(coverage_data$country)))


  print(p)
}
scenario_comparison_plot<- function(output, country_name){

  output<- data.table(output)
  dt<- output |>
    group_by(scenario, parameter_draw) |>
    summarise(cases = sum(cases),
              cases_novax= sum(cases_novax),
              .groups= 'keep') |>
    mutate(prop_cases_averted = 1- cases/ cases_novax) |>
    mutate(scenario = factor(scenario, levels =c("malaria-rts3-default",
                                                 "malaria-rts3-bluesky",
                                                 "malaria-rts3-rts4-default",
                                                 "malaria-rts3-rts4-bluesky",
                                                 "no-vaccination")))

  dt<- data.table(dt)
  p<- ggplot(data = dt, mapping = aes(x= scenario, y= prop_cases_averted, fill = scenario)) +
    geom_col(position= 'dodge') +
    #facet_wrap(~parameter_draw, scales = 'free') +
    labs(title= 'Proportion of cases averted over 100-year period by age and scenario',
         subtitle = unique(dt$country_name))+
    plotting_theme

  return(p)
}
deaths_over_time<- function(output, scenario, coverage_data){

  intro_yr<- min(coverage_data[scenario == {{scenario}}& coverage> 0, year])
  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  stochastic_run<- output |>
    filter(parameter_draw != 0)
  median_run<- output |>
    group_by(year, scenario) |>
    summarise(cases = mean(cases),
              deaths = mean(deaths),
              dalys = mean(dalys),
              .groups = 'keep')

  p<-   ggplot()+
    geom_line(data = median_run, mapping = aes(x= year, y= deaths, color =scenario))  +
    geom_line(data = stochastic_run, mapping = aes(x= year, y= deaths, color = scenario, group = parameter_draw), alpha = 0.05)  +
    geom_vline(xintercept= intro_yr, linetype= "dotted") +
    labs(x= 'Time (in years)', y= 'Deaths',
         title= paste0('Deaths over time'),
         subtitle = {{scenario}},
         color= 'Scenario', fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Darjeeling1', n= 2)) +
    scale_fill_manual(values= wes_palette('Darjeeling1', n= 2))  +
    plotting_theme

  return(p)
}
deaths_over_age<- function(output, scenario, coverage_data){
  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  stochastic_run<- output |>
    filter(parameter_draw != 0)
  median_run<- output |>
    group_by(age, scenario) |>
    summarise(cases = mean(cases),
              deaths = mean(deaths),
              dalys = mean(dalys),
              .groups = 'keep')

  p<-   ggplot()+
    geom_line(data = median_run, mapping = aes(x= age, y= deaths, color =scenario))  +
    geom_line(data = stochastic_run, mapping = aes(x= age, y= deaths, color = scenario, group = parameter_draw), alpha = 0.05)  +
    #facet_wrap_paginate(~scenario) +
    labs(x= 'Time (in years)', y= 'Deaths',
         title= paste0('Deaths over age'),
         subtitle = unique(median_run$country_name),
         color= 'Scenario', fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Darjeeling1', n= 2)) +
    scale_fill_manual(values= wes_palette('Darjeeling1', n= 2))  +
    plotting_theme

  return(p)
}
dalys_over_time<- function(output, scenario, coverage_data){

  intro_yr<- min(coverage_data[scenario == {{scenario}} & coverage> 0, year])
  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  stochastic_run<- output |>
    filter(parameter_draw != 0)
  median_run<- output |>
    group_by(year, scenario) |>
    summarise(cases = mean(cases),
              deaths = mean(deaths),
              dalys = mean(dalys),
              .groups = 'keep')

  p<-   ggplot()+
    geom_line(data = median_run, mapping = aes(x= year, y= dalys, color =scenario))  +
    geom_line(data = stochastic_run, mapping = aes(x= year, y= dalys, color = scenario, group = parameter_draw), alpha = 0.05)  +
    geom_vline(xintercept= intro_yr, linetype= "dotted") +
    labs(x= 'Time (in years)', y= 'DALYs',
         title= paste0('DALYs over time'),
         subtitle = {{scenario}},
         color= 'Scenario', fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Darjeeling1', n= 2)) +
    scale_fill_manual(values= wes_palette('Darjeeling1', n= 2))  +
    plotting_theme


  return(p)
}
dalys_over_age<- function(output, scenario, coverage_data){
  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  stochastic_run<- output |>
    filter(parameter_draw != 0)
  median_run<- output |>
    group_by(age, scenario) |>
    summarise(cases = mean(cases),
              deaths = mean(deaths),
              dalys = mean(dalys),
              .groups = 'keep')

  p<-   ggplot()+
    geom_line(data = median_run, mapping = aes(x= age, y= dalys, color =scenario))  +
    geom_line(data = stochastic_run, mapping = aes(x= age, y= dalys, color = scenario, group = parameter_draw), alpha = 0.05)  +
    #facet_wrap_paginate(~scenario) +
    labs(x= 'Time (in years)', y= 'DALYs',
         title= paste0('DALYs over age'),
         subtitle = unique(median_run$country_name),
         color= 'Scenario', fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Darjeeling1', n= 2)) +
    scale_fill_manual(values= wes_palette('Darjeeling1', n= 2))  +
    plotting_theme

  return(p)
}
cases_deaths_diagnostic<- function(site){

  dt<- data.table(site$cases_deaths)

  site_name<- site$sites$name_1
  ur<- site$sites$urban_rural

  lower <-
    melt(dt,
         measure.vars = c('wmr_cases_l', 'wmr_deaths_l', 'wmr_incidence_l', 'wmr_mortality_l'))  |>
    select(year, variable, value) |>
    rename(lower = value)

  upper <-
    melt(dt,
         measure.vars = c('wmr_cases_u', 'wmr_deaths_u', 'wmr_incidence_u', 'wmr_mortality_u')) |>
    select(year, variable, value) |>
    rename(upper = value)

  median <-
    melt(dt,
         measure.vars = c('wmr_cases', 'wmr_deaths', 'wmr_incidence', 'wmr_mortality')) |>
    select(year, variable, value)

  lower[variable %like% 'cases', variable := 'cases']
  lower[variable %like% 'deaths', variable := 'deaths']
  lower[variable %like% 'mortality', variable := 'mortality']
  lower[variable %like% 'incidence', variable := 'incidence']

  upper[variable %like% 'cases', variable := 'cases']
  upper[variable %like% 'deaths', variable := 'deaths']
  upper[variable %like% 'mortality', variable := 'mortality']
  upper[variable %like% 'incidence', variable := 'incidence']

  median[variable %like% 'cases', variable := 'cases']
  median[variable %like% 'deaths', variable := 'deaths']
  median[variable %like% 'mortality', variable := 'mortality']
  median[variable %like% 'incidence', variable := 'incidence']

  plot <- Reduce(merge, list(lower, median, upper))

  p<- ggplot(data = plot) +
    geom_line(mapping = aes(x = year, y = value, color = variable)) +
    geom_ribbon(
      mapping = aes(
        x = year,
        ymin = lower,
        ymax = upper,
        fill = variable
      ),
      alpha = 0.2
    ) +
    facet_wrap( ~ variable, ncol = 4, scales = 'free') +
    scale_y_continuous(labels = comma) +
    scale_color_manual(values = wes_palette('Darjeeling1', n = 4)) +
    scale_fill_manual(values = wes_palette('Darjeeling1', n = 4))  +
    plotting_theme +
    labs(
      y = 'Value',
      x = 'Year',
      title = 'Cases, deaths, mortality, and incidence over time')

  return(p)
}
cases_averted_over_time<- function(output, scenario, country_name){

  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  ggplot(data = output, mapping = aes(x= year, y= cases_averted), alpha = 0.05) +
    geom_line(aes( color = scenario, group= parameter_draw)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='Cases averted over time',
         subtitle = country_name) +
    plotting_theme
}
cases_averted_over_age<- function(output, scenario, country_name){

  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  ggplot(data = output, mapping = aes(x= age, y= cases_averted), alpha = 0.05) +
    geom_line(aes( color = scenario, group= parameter_draw)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='Cases averted over age',
         subtitle = country_name) +
    plotting_theme
}
deaths_averted_over_time<- function(output, scenario, country_name){

  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  ggplot(data = output, mapping = aes(x= year, y= deaths_averted), alpha = 0.05) +
    geom_line(aes( color = scenario, group= parameter_draw)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='Cases averted over time',
         subtitle = country_name) +
    plotting_theme
}
deaths_averted_over_age<- function(output, scenario, country_name){

  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')


  ggplot(data = output, mapping = aes(x= age, y= deaths_averted), alpha = 0.05) +
    geom_line(aes( color = scenario, group= parameter_draw)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='Deaths averted over age',
         subtitle = country_name) +
    plotting_theme
}
dalys_averted_over_time<- function(output, scenario, country_name){

  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')

  ggplot(data = output, mapping = aes(x= year, y= dalys_averted), alpha = 0.05) +
    geom_line(aes( color = scenario, group= parameter_draw)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='DALYs averted over time',
         subtitle = country_name) +
    plotting_theme
}
dalys_averted_over_age<- function(output, scenario, country_name){

  output<- output |>
    filter(scenario == {{scenario}} | scenario == 'no-vaccination')


  ggplot(data = output, mapping = aes(x= age, y= dalys_averted), alpha = 0.05) +
    geom_line(aes( color = scenario, group= parameter_draw)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='DALYs averted over age',
         subtitle = country_name) +
    plotting_theme
}


prepare_dose_output<- function(dose_output, site_data, coverage_data){
  intro_yr<- min(coverage_data[scenario == 'malaria-rts3-rts4-default' & coverage> 0, year])
  
  # first aggregate cases + deaths averted by site_ur over entire simulation period
  site_doses<- dose_output |>
    group_by(site_ur, scenario) |>
    #filter(parameter_draw == 0)|>
    filter(year %in% c(intro_yr:(intro_yr+15))) |>
    summarise(cases_averted = sum(cases_averted),
              deaths_averted = sum(deaths_averted), 
              doses_total = sum(doses_total),
              fvp = sum(fvp),
              .groups = 'keep')
  
  # then merge on site data based on site name and ur
  sites<- site_data$prevalence |> 
    filter(year == 2019) |>
    mutate(site_ur = paste0(name_1, '_', urban_rural)) |>
    select(site_ur, pfpr)
  site_doses<- merge(site_doses, sites, by = 'site_ur')
  
  
  # then sum all doses + cases + deaths averted up to country level by year
  country_doses<- dose_output |>
    group_by(year, scenario) |>
    #filter(parameter_draw== 0) |>
    summarise(cases_averted = sum(cases_averted),
              deaths_averted = sum(deaths_averted), 
              doses_total = sum(doses_total),
              fvp = sum(fvp),
              .groups = 'keep') |>
    filter(doses_total != 0)
  
  return(list('site_doses' = site_doses, 'country_doses' = country_doses)) 
}





plot_case_site_comparison<- function(site_doses){
  p<- ggplot(data = site_doses, mapping = aes(x= pfpr, y = (cases_averted/fvp) *100000, color = scenario, shape = site_ur)) +
    geom_point()+
    plotting_theme +
    theme(legend.position = 'none')+
    labs(title = 'Unscaled cases averted per 100,000 FVP, 15 yrs after introduction',
         x= 'PFPR(2-10) in 2019',
         y = 'Cases averted per 100k FVPs')
  
  ggplotly(p)
  
}

plot_death_site_comparison<-  function(site_doses){
  p<- ggplot(data = site_doses, mapping = aes(x= pfpr, y = deaths_averted/fvp *100000, color = scenario, shape = site_ur)) +
    geom_point()+
    plotting_theme +
    theme(legend.position = 'none')+
    labs(title = 'Unscaled deaths averted per 100,000 FVP, 2000-2100',
         x= 'PFPR(2-10) in 2019',
         y = 'Deaths averted per 100k FVPs')
  
  
  ggplotly(p)
  
  
}


cases_averted_per_dose<- function(country_doses, scenario, country_name){
  plot<- country_doses |> filter(scenario == {{scenario}})
  
  ggplot(data = plot, mapping = aes(x= year, y= cases_averted * 100000/doses_total), alpha = 0.05) +
    geom_line(aes( color = scenario)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='Unscaled cases averted per 100k doses delivered over time',
         subtitle = country_name) +
    plotting_theme
  
}


cases_averted_per_fvp<- function(country_doses, scenario, country_name){
  
  
  plot<- country_doses |> filter(scenario == {{scenario}})
  
  ggplot(data = plot, mapping = aes(x= year, y= cases_averted* 100000/fvp), alpha = 0.05) +
    geom_line(aes( color = scenario)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='Unscaled cases averted per 100k FVP over time',
         subtitle = country_name) +
    plotting_theme
  
}


deaths_averted_per_dose<- function(country_doses, scenario, country_name){
  plot<- country_doses |> filter(scenario == {{scenario}})
  
  ggplot(data = plot, mapping = aes(x= year, y= deaths_averted* 100000/doses_total), alpha = 0.05) +
    geom_line(aes( color = scenario)) +
    stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='Unscaled deaths averted per 100k doses delivered over time',
         subtitle = country_name) +
    plotting_theme
  
}



deaths_averted_per_fvp<- function(country_doses, scenario, country_name){
  
  plot<- country_doses |> filter(scenario == {{scenario}})
  ggplot(data = plot, mapping = aes(x= year, y= deaths_averted* 100000/fvp), alpha = 0.05) +
    geom_line(aes( color = scenario)) +
    #stat_summary(fun.y = mean, geom = "line", color = "black") +
    geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
    facet_wrap(~scenario)+
    labs(title ='Unscaled deaths averted per 100k FVP over time',
         subtitle = country_name) +
    plotting_theme
  
}