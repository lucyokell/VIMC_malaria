# site file diagnostics --------------------------------------------------------

library(scales)
library(ggpubr)
library(wesanderson)
library(tidyr)


# plot cases + deaths

# cases and deaths  ------------------------------------------------------------
cases_deaths_diagnostic<- function(site){
  
  dt<- data.table::data.table(site$cases_deaths)
  
  site_name<- site$sites$name_1
  ur<- site$sites$urban_rural
  
  lower <-
    melt(dt,
         measure.vars = c('wmr_cases_l', 'wmr_deaths_l', 'incidence_l', 'mortality_l'))  |>
    select(year, variable, value) |>
    rename(lower = value)
  
  upper <-
    melt(dt,
         measure.vars = c('wmr_cases_u', 'wmr_deaths_u', 'incidence_u', 'mortality_u')) |>
    select(year, variable, value) |>
    rename(upper = value)
  
  median <-
    melt(dt,
         measure.vars = c('wmr_cases', 'wmr_deaths', 'incidence', 'mortality')) |>
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
    facet_wrap( ~ variable, scales = 'free', ncol = 4) +
    scale_y_continuous(labels = comma) +
    scale_color_manual(values = wes_palette('Royal2', n = 4)) +
    scale_fill_manual(values = wes_palette('Royal2', n = 4))  +
    plotting_theme +
    labs(
      y = 'Value',
      x = 'Year',
      title = 'Cases, deaths, mortality, and incidence over time',
      subtitle = paste('Site: ', site_name, ' ', ur)
    )
  
  return(p)
}


# prevalence -------------------------------------------------------------------
prevalence_diagnostic<- function(site){
  
  dt<- melt(data.table::data.table(site$prevalence), measure.vars= c('pfpr', 'pvpr'))
  
  site_name<- site$sites$name_1
  ur<- site$sites$urban_rural
  
  p<- ggplot(dt, mapping = aes(x= year, y= value, color= variable))+
    geom_line()+
    plotting_theme+
    labs(x= 'Year', y= 'Prevalence value', color= 'Prevalence type',
         title= 'Parasite prevalence over time',
         subtitle = paste0('Site: ', site_name, ' ', ur))
  
  print(p)
}



# population  ------------------------------------------------------------------

population_diagnostic<- function(site){
  
  dt<- melt(data.table::data.table(site$population), measure.vars = c('par_pf', 'par_pv'))
  
  p1<- ggplot(data=dt)+
    geom_col(position = 'stack', mapping= aes(x= year, y= value, fill= variable)) +
    scale_fill_manual(values = wes_palette('Royal2', n = 2))  +
    geom_line(mapping= aes(x= year, y=par)) +
    plotting_theme+
    labs(title= 'Population at risk over time', 
         subtitle= paste0('Site: ', site_name, ' ', ur),
         y= 'Population', x= 'Year',
         fill= 'Population at risk') +
    scale_y_continuous(labels= comma)
  
  
  dt<- melt(data.table::data.table(site$population), measure.vars = c('pop', 'par'))
  
  p2<- ggplot(data=dt)+
    geom_line(mapping= aes(x= year, y= value, color= variable)) +
    plotting_theme+
    labs(title= 'Population over time', 
         subtitle= paste0('Site: ', site_name, ' ', ur),
         y= 'Population', x= 'Year') +
    scale_y_continuous(labels= comma)
  
  plots<- ggarrange(p1, p2, ncol= 2)
  
  return(plots)         
}
