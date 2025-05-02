
library(ggplot2)
library(dplyr)
plotting_theme<- theme_bw(base_size = 14) +
  theme( legend.position = 'bottom',
         strip.text.x = element_text(size = rel(0.8)),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
         #text = element_text(family= 'Calibri'),
         axis.ticks.y= element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())


no_vax<- read.csv('montagu/central-burden-est-no-vaccination.csv')
full<- rbind(no_vax, rtss)

full<- full |>
  filter(country== 'ETH') |>
  group_by(scenario, year) |>
  summarise(cases = sum(cases),
            .groups = 'keep')



source('src/diagnostics/diagnostic_report_functions.R')
rtss_3<- read.csv('montagu/central-burden-est-malaria-rts3-bluesky.csv')
rtss_3<- rtss_3 |> mutate(scenario = 'rtss_3_dose_bluesky') |>
  rename(cases_rtss_3 = cases)


rtss_4<- read.csv('montagu/central-burden-est-malaria-rts3-rts4-bluesky.csv')
rtss_4<- rtss_4 |> mutate(scenario = 'rtss_4_dose_bluesky') |>
  rename(cases_rtss_4 = cases)


full<- merge(rtss_3, rtss_4, by = c('year', 'age', 'country'))
full<- full |>
  filter(country== 'ETH') |>
  group_by(year) |>
  summarise(cases_rtss_3 = sum(cases_rtss_3),
            cases_rtss_4 = sum(cases_rtss_4),
            .groups = 'keep')

ggplot(data = full, mapping = aes(x= year, y=  cases_rtss_3 - cases_rtss_4)) +
  geom_line(color = 'darkgreen', size= 0.75) +
  plotting_theme +
  labs(title= 'Cases averted by calendar year, malaria-rts3-rts4-bluesky')+
  theme(text = element_text(family = 'Calibri'))


ggplot(data = full, mapping = aes(x= year, y=  cases_rtss_4 - cases_rtss_3)) +
  geom_line() 


rtss_3<- read.csv('montagu/central-burden-est-malaria-rts3-bluesky.csv')
rtss_3<- rtss_3 |> mutate(scenario = 'rtss_3_dose_bluesky')

rtss_4<- read.csv('montagu/central-burden-est-malaria-rts3-rts4-bluesky.csv')
rtss_4<- rtss_4 |> mutate(scenario = 'rtss_4_dose_bluesky') 
  
# plot cases
outputs<- compile_final_outputs('fix_rtss_booster')
#outputs_2<- compile_final_outputs('recalibrated_irs')

rtss_3<- outputs |> filter(scenario == 'no-vaccination') |>
  rename(cases_rtss_3 = cases)
rtss_4<- outputs |> filter(scenario == 'malaria-r3-r4-bluesky') |>
  rename(cases_rtss_4 = cases)

full<- merge(rtss_3, rtss_4, by = c('year', 'age', 'country', 'run_id'))

full<- full |>
  #filter(country== 'ETH') |>
  #filter(run_id != 0) |>
  group_by(year, country, run_id) |>
  summarise(cases_rtss_3 = sum(cases_rtss_3),
            cases_rtss_4 = sum(cases_rtss_4),
            .groups = 'keep') |>
  mutate(averted = cases_rtss_3- cases_rtss_4) 


ggplot() +
  geom_line(data = full, mapping = aes(x= year, y=  averted, group = country), color = 'red', alpha = 0.5) +
  #geom_line(data = mean_run, mapping = aes(x= year, y= averted), color = 'black') +
  plotting_theme +
  labs(title= 'Cases averted by calendar year, malaria-rts3-rts4-bluesky')+
  theme(text = element_text(family = 'Calibri'))



full<- rbind(rtss_3, rtss_4)

full<- full |>
  filter(country== 'ETH') |>
  group_by(scenario, year) |>
  summarise(cases = sum(cases),
            .groups = 'keep')

ggplot(data = full, mapping = aes(x= year, y=  cases, color= scenario)) +
  geom_line() 
