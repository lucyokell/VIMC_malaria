# pull outputs for slides

library(vimcmalaria)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(plotly)
library(htmlwidgets)

cols<- c('#CACF85', '#8CBA80', '#658E9C', '#4D5382', '#514663')
cols2<- c('#658E9C', '#4D5382')

plotting_theme<- theme_bw(base_size = 14) +
  theme( legend.position = 'bottom',
         strip.text.x = element_text(size = rel(0.8)),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
         #text = element_text(family= 'Arial'),
         axis.ticks.y= element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#  pull new VIMC runs

outputs<- compile_final_outputs('new_site_files') 
outputs<- outputs |>
  mutate(run = 'new') |>
  group_by(country_name, year, scenario, run) |>
  summarise(cases = sum(cases),
            deaths= sum(deaths),
            cohort_size = sum(cohort_size),
          .groups= 'keep')


previous_r21<- read.csv('J:/september_runs/VIMC_malaria/montagu/central-burden-est-malaria-r3-r4-default.csv') |> mutate(scenario = 'malaria-r3-r4-default')
previous_rtss<- read.csv('J:/september_runs/VIMC_malaria/montagu/central-burden-est-malaria-rts3-rts4-default.csv')|> mutate(scenario = 'malaria-rts3-rts4-default')
previous_bl<- read.csv('J:/september_runs/VIMC_malaria/montagu/central-burden-est-no-vaccination.csv') |> mutate(scenario = 'no-vaccination')

prev<- bind_rows(previous_r21, previous_rtss, previous_bl) |>
  mutate(run = 'previous_site_file') |>
group_by(country_name, year, scenario, run) |>
  summarise(cases = sum(cases),
            deaths= sum(deaths),
            cohort_size = sum(cohort_size),
          .groups = 'keep')

merged<- merge(prev, outputs, by = c('scenario', 'year', 'country_name'))

merged<- rbind(prev, outputs)

pdf('analyses/site_files/site_file_comparison_plots.pdf')
for (country in unique(merged$country_name)){
  message(country)

  subset<- merged |> filter(country_name == country)

  p1<- ggplot(subset)+
    geom_line(aes(x=year, y= cases, color= run)) +
    facet_wrap(~scenario)+
    labs(title= 'Annual incidence, new site files vs. old',
          subtitle= country,
         x= 'Year',
        y= 'Cases') + 
    plotting_theme +
    scale_color_manual(values= cols) 


  p2<- ggplot(subset)+
    geom_line(aes(x=year, y= deaths, color= run)) +
    facet_wrap(~scenario)+
    labs(title= 'Annual mortality, new site files vs. old',
          subtitle= country,
         x= 'Year',
        y= 'Deaths') + 
    plotting_theme +
    scale_color_manual(values= cols[3:4]) 


  
  print(p1)
  print(p2)
}

dev.off()


# now do a general percent change plot (all-age, all years)
# do this for
prev<- bind_rows(previous_r21, previous_rtss, previous_bl) |>
group_by(country_name, year, scenario) |>
  summarise(cases_pre = sum(cases),
            deaths_pre= sum(deaths),
          .groups = 'keep')

merged<- merge(prev, outputs, by = c('scenario', 'year', 'country_name'))


merged<- merged |>
  group_by(scenario, country_name) |>
  summarise(cases_pre= sum(cases_pre),
            deaths_pre = sum(deaths_pre),
            cases= sum(cases),
            deaths= sum(deaths),
          .groups = 'keep') |>
  mutate(pct_change_cases= (cases-cases_pre)/cases_pre,
         pct_change_deaths = (deaths-deaths_pre)/ deaths_pre)


p<- ggplot(merged)+
  geom_point(aes(x= pct_change_cases*100, y= pct_change_deaths*100, color= country_name)) +
  labs(title= 'Percent change in cases and deaths between site files, 2000-2100')
# pre-post diagnostic
ggplotly()

saveWidget(ggplotly(p), file = "analyses/site_files/pre_post_diagnostic.html")
