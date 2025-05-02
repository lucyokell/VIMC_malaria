#   MIM Poster Plots  ----------------------------------------------------------


setwd('J:/VIMC_malaria/MIM')
source('J:/VIMC_malaria/src/diagnostics/diagnostic_report_functions.R')

windowsFonts(myFont = windowsFont("Arial"))
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)
library(scene)
library(scales)
library(extrafont)
loadfonts(device = "pdf", quiet = T)
inputs<- readRDS('diagnostic_inputs.rds')

cols<- c('#CACF85', '#8CBA80', '#658E9C', '#4D5382', '#514663')
cols2<- c('#658E9C', '#4D5382')

plotting_theme<- theme_bw(base_size = 18) +
  theme( legend.position = 'bottom',
         strip.text.x = element_text(size = rel(0.8)),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
         #text = element_text(family= 'Arial'),
         axis.ticks.y= element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())


vetting<- inputs$vetting
vetting<- vetting |>
  filter(!scenario %like% 'default')
output<- vetting

site_output<- inputs$site_output
averted_output<- inputs$averted_output
burden_output<- inputs$burden_output
coverage_data<- inputs$coverage_data |> filter (!scenario %like% 'default')
site_data<- inputs$site_data
population <- inputs$population
population_by_age<- inputs$population_by_age
country_name<- inputs$country_name

pfprs<- readRDS('pfpr10plus_admins.rds') |> filter(iso3c == 'COD')


plot_intvn_coverage(site_data)
cases_deaths_diagnostic(site_data)

# incidence over time and age --------------------------------------------------
intro_yr<- min(coverage_data[coverage> 0, year])

output<- burden_output$by_year |> filter(!scenario %like% 'default')

stochastic_run<- output |>
  filter(parameter_draw != 0)

mean_run<- stochastic_run |>
  group_by(year, scenario) |>
  summarise(cases = mean(cases),
            deaths = mean(deaths),
            dalys = mean(dalys),
            cohort_size = mean(cohort_size),
            .groups = 'keep')

p<-   ggplot()+
  geom_line(data = mean_run, mapping = aes(x= year, y= cases, color =scenario))  +
  #geom_line(data = stochastic_run, mapping = aes(x= year, y= cases, color = scenario, group = parameter_draw), alpha = 0.1)  +
  geom_line(data= site_data$cases_deaths, mapping = aes(x= year, y= wmr_cases), color= 'darkgreen')+
  geom_vline(xintercept= intro_yr, linetype= "dotted") +
  labs(x= 'Time (in years)', y= 'Cases',
       title= paste0('Cases over time'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values= cols) +
  scale_fill_manual(values= cols)  +
  plotting_theme

p1<- ggplot()+
  geom_line(data = mean_run, mapping = aes(x= year, y= cases/cohort_size, color =scenario))  +
  #geom_line(data = stochastic_run, mapping = aes(x= year, y= cases, color = scenario, group = parameter_draw), alpha = 0.1)  +
  #geom_line(data= site_data$cases_deaths, mapping = aes(x= year, y= wmr_cases), color= 'darkgreen')+
  geom_vline(xintercept= intro_yr, linetype= "dotted") +
  labs(x= 'Time (in years)', y= 'Incidence rate',
       title= paste0('Incidence rate over time'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values=cols) +
  scale_fill_manual(values= cols)  +
  plotting_theme


p2<- ggplot()+
  geom_line(data = mean_run, mapping = aes(x= year, y= deaths/cohort_size, color =scenario))  +
  #geom_line(data = stochastic_run, mapping = aes(x= year, y= cases, color = scenario, group = parameter_draw), alpha = 0.1)  +
  #geom_line(data= site_data$cases_deaths, mapping = aes(x= year, y= wmr_cases), color= 'darkgreen')+
  geom_vline(xintercept= intro_yr, linetype= "dotted") +
  labs(x= 'Time (in years)', y= 'Mortality',
       title= paste0('Mortality rate over time'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values= cols) +
  scale_fill_manual(values= cols)  +
  plotting_theme

ggarrange(p, p1, p2, ncol= 3,  widths = c(1, 1, 2))

output<- burden_output$by_age |> filter(!scenario %like% 'default') |> filter(!scenario == 'malaria-rts3-bluesky')


output<- vetting |> filter(!scenario %like% 'default') |> filter(!scenario == 'malaria-rts3-bluesky')|> filter(year >= 2038 & year <2089)

stochastic_run<- output |>
  filter(parameter_draw != 0)

mean_run<- stochastic_run |>
  group_by(scenario, age, run_id) |>
  summarise(cases = sum(cases),
            deaths = sum(deaths),
            dalys = sum(dalys),
            cohort_size = sum(cohort_size),
            .groups = 'keep')

mean_run<- stochastic_run |>
  group_by(scenario, age) |>
  summarise(cases = mean(cases),
            deaths = mean(deaths),
            dalys = mean(dalys),
            cohort_size = mean(cohort_size),
            .groups = 'keep')

p5<- ggplot()+
  geom_line(data = mean_run, mapping = aes(x= age, y= cases/cohort_size, color =scenario))  +
  labs(x= 'Age (in years)', y= 'Incidence rate',
       #title= paste0('Incidence rate over age'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values= wes_palette('Darjeeling1', n= 5)) +
  scale_fill_manual(values= wes_palette('Darjeeling1', n= 5))  +
  plotting_theme +
  theme(legend.position = 'none')


p6<- ggplot()+
  geom_line(data = mean_run, mapping = aes(x= age, y= deaths/cohort_size, color =scenario))  +
  labs(x= 'Age (in years)', y= 'Mortality rate',
       #title= paste0('Mortality rate over age'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values= wes_palette('Darjeeling1', n= 5)) +
  scale_fill_manual(values= wes_palette('Darjeeling1', n= 5))  +
  plotting_theme +
  theme(legend.position = 'bottom')

ggarrange(p5, p6, ncol = 2, widths = c(1, 2))



# incidence over time and age
output<- vetting |> filter(!scenario %like% 'default')

output<- output |>
  mutate(age_group = ifelse(age <= 5, '0_5', NA)) |>
  mutate(age_group = ifelse(age > 5 & age <=15, '5_15', age_group)) |>
  mutate(age_group = ifelse(age > 15, '15+', age_group))

stochastic_run<- output |>
  filter(parameter_draw != 0)

mean_run<- stochastic_run |>
  group_by(year, scenario, age_group, run_id) |>
  summarise(cases = sum(cases),
            deaths = sum(deaths),
            dalys = sum(dalys),
            cohort_size = sum(cohort_size),
            .groups = 'keep')

mean_run<- mean_run |>
  group_by(year, scenario, age_group) |>
  summarise(cases = mean(cases),
            deaths = mean(deaths),
            dalys = mean(dalys),
            cohort_size = mean(cohort_size),
            .groups = 'keep')

mean_run$age_group<- factor(mean_run$age_group, levels = c('0_5', '5_15', '15+'))
mean_run<- mean_run |>
  filter(scenario %in% c('no-vaccination', 'malaria-rts3-rts4-bluesky'))

mean_run<- data.table(mean_run)
mean_run[scenario %like% 'bluesky', scenario := 'vaccine']



p<-   ggplot()+
  geom_smooth(data = mean_run, mapping = aes(x= year, y= cases/cohort_size, color =scenario, fill= scenario))  +
  #geom_line(data = stochastic_run, mapping = aes(x= year, y= cases, color = scenario, group = parameter_draw), alpha = 0.1)  +
  geom_vline(xintercept= intro_yr, linetype= "dotted") +
  facet_wrap(~age_group, scales = 'free')+
  labs(x= 'Year', y= 'Incidence rate',
       #title= paste0('Incidence rate over time'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values= cols) +
  scale_fill_manual(values= cols)  +
  plotting_theme


p2<-   ggplot()+
  geom_smooth(data = mean_run, mapping = aes(x= year, y= deaths/cohort_size, color =scenario, fill= scenario))  +
  #geom_line(data = stochastic_run, mapping = aes(x= year, y= cases, color = scenario, group = parameter_draw), alpha = 0.1)  +
  geom_vline(xintercept= intro_yr, linetype= "dotted") +
  facet_wrap(~age_group, scales = 'free')+
  labs(x= 'Year', y= 'Mortality rate',
       #title= paste0('Mortality rate over time'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values= cols2) +
  scale_fill_manual(values= cols2)  +
  plotting_theme


ggarrange(p, p2)


# cases averted
output<- vetting |> filter(!scenario %like% 'default')

output<- output |>
  mutate(age_group = ifelse(age <= 5, '0_5', NA)) |>
  mutate(age_group = ifelse(age >5 & age <=15, '5_15', age_group)) |>
  mutate(age_group = ifelse(age > 15, '15+', age_group))

stochastic_run<- output |>
  filter(parameter_draw != 0)

mean_run<- stochastic_run |>
  group_by(year, scenario, age_group, run_id) |>
  summarise(cases = sum(cases),
            deaths = sum(deaths),
            dalys = sum(dalys),
            cohort_size = sum(cohort_size),
            .groups = 'keep')

novax<- mean_run |>
  filter(scenario == 'no-vaccination') |>
  rename(cases_novax = cases,
         deaths_novax = deaths,
         dalys_novax = dalys) |>
  ungroup() |>
  select(-cohort_size, -scenario)

mean_run<- mean_run |>
  filter(scenario!= 'no-vaccination')


merged<- merge(mean_run, novax, by = c('year', 'age_group', 'run_id'))

merged<- merged |>
  group_by(year, scenario, age_group) |>
  mutate(cases_averted = cases_novax- cases,
         deaths_averted = deaths_novax - deaths,
         dalys_averted = dalys_novax- dalys) |>
  summarise(cases = mean(cases),
            deaths = mean(deaths),
            dalys = mean(dalys),
            cohort_size = mean(cohort_size),
            .groups = 'keep')

p<-   ggplot()+
  geom_smooth(data = merged |> filter(!scenario == 'malaria-rts3-bluesky'), mapping = aes(x= year, y= cases_averted/cohort_size, color =scenario, fill= scenario))  +
  #geom_line(data = stochastic_run, mapping = aes(x= year, y= cases, color = scenario, group = parameter_draw), alpha = 0.1)  +
  geom_vline(xintercept= intro_yr, linetype= "dotted") +
  facet_wrap(~age_group, scales = 'free')+
  labs(x= 'Time (in years)', y= 'Incidence rate',
       #title= paste0('Incidence rate over time'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values= cols) +
  scale_fill_manual(values= cols)  +
  plotting_theme


p2<-   ggplot()+
  geom_smooth(data = merged |> filter(!scenario == 'malaria-rts3-bluesky'), mapping = aes(x= year, y= deaths_averted/cohort_size, color =scenario, fill= scenario))  +
  #geom_line(data = stochastic_run, mapping = aes(x= year, y= cases, color = scenario, group = parameter_draw), alpha = 0.1)  +
  geom_vline(xintercept= intro_yr, linetype= "dotted") +
  facet_wrap(~age_group, scales = 'free')+
  labs(x= 'Time (in years)', y= 'Mortality rate',
       #title= paste0('Mortality rate over time'),
       color= 'Scenario', fill= 'Scenario') +
  scale_color_manual(values= cols2) +
  scale_fill_manual(values= cols2)  +
  plotting_theme
