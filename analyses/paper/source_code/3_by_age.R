# by age
library(vimcmalaria)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(ggpubr)
library(data.table)
library(stringr)



  cols<- c('#CE5137', '#1B4F72', '#45B39D', '#BA4A00')

  plotting_theme<- theme_bw(base_size = 14) +
    theme( legend.position = 'bottom',
           strip.text.x = element_text(size = rel(0.8)),
           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
           text = element_text(family= 'TT Arial'),
           axis.ticks.y= element_blank(),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# FIGURE 2: cases + deaths averted by age
vax<- readRDS('analyses/paper/inputs/vax.rds')

vax_age<- copy(vax)
vax_age$age<- cut(vax_age$age, breaks = c(0, 5, 15, 101), include.lowest = TRUE)
table_2<- vax_age |>
  filter(scenario== 'Routine') |>
  group_by(scenario, age, year, run_id) |>
  summarise(cases_averted = sum(cases_averted),
            deaths_averted = sum(deaths_averted),
            dalys_averted = sum(dalys_averted),
            cohort_size = sum(cohort_size),
            .groups = 'keep') |>
  data.table()  |>
    ungroup() |>
    group_by(scenario, age, year) |>
    summarise(
      cases_mean = mean(cases_averted, na.rm = TRUE),
      cases_lower = quantile(cases_averted, 0.025, na.rm = TRUE),
      cases_upper = quantile(cases_averted, 0.975, na.rm = TRUE),
      deaths_mean = mean(deaths_averted, na.rm = TRUE),
      deaths_lower = quantile(deaths_averted, 0.025, na.rm = TRUE),
      deaths_upper = quantile(deaths_averted, 0.975, na.rm = TRUE),
      cohort_size = mean(cohort_size, na.rm = TRUE),
      .groups= 'keep'
    ) |>
    unique() |>
  as.data.table()

table_2<- data.table(table_2)

  table_2[age== '[0,5]', age:= 'under 5']
  table_2[age== '(5,15]', age:= '5 to 15']
  table_2[age== '(15,101]', age:= 'over 15']

p3<- ggplot(table_2) +
  geom_line(mapping = aes(x= year, y= cases_mean, color = age)) +
  geom_ribbon(aes(x= year, ymin = cases_lower, ymax = cases_upper, fill= age) , alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  xlim(2020, 2100) +
  scale_y_continuous(labels = scales::comma) +
  plotting_theme+
    theme(legend.position = 'none') +
  labs(title = '',
       x= '',
       y= 'Cases averted')

p4<- ggplot(table_2) +
  geom_line(mapping = aes(x= year, y= deaths_mean, color = age)) +
  geom_ribbon(aes(x= year, ymin = deaths_lower, ymax = deaths_upper, fill= age) , alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  xlim(2020, 2100)+
  scale_y_continuous(labels = scales::comma) +
  plotting_theme+
  labs(title = '',
       x= 'Year',
       y= 'Deaths averted')

ggarrange(p3, p4, nrow = 2, heights= c(1, 1.25))
# try same plot but as a stacked bar chart
# should calculate mean + cis after all calculations but update this later
stacked<- table_2 |>
  mutate(year_cut = cut(year,
                        breaks = c(1999, 2005, 2010, 2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060, 2065, 2070, 2075, 2080, 2085, 2090, 2095, 2101), ordered_result = TRUE))
stacked<- stacked |>
  mutate(year_cut = as.character(year_cut)) |>
  mutate(year_cut = str_replace(year_cut, ',', '-')) |>
 # mutate(year_cut = str_replace(year_cut, "(", ""))
  #mutate(year_cut = str_replace(year_cut, ']', ''))
group_by(year_cut,age, scenario) |>
  summarise(cases_mean = sum(cases_mean),
            cases_lower = sum(cases_lower),
            cases_upper = sum(cases_upper),
            deaths_mean = sum(deaths_mean),
            deaths_lower = sum(deaths_lower),
            deaths_upper = sum(deaths_upper),
            .groups = 'keep')

p3<- ggplot() +
  geom_col(stacked, position = 'stack', mapping = aes(x= year_cut, y= cases_mean, fill = age)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::comma) +
  #xlim((1999,2005], (2095,2101]) +
  plotting_theme+
  theme(legend.position = 'none') +
  labs(title = '',
       x= '',
       y= 'Cases averted')

p4<- ggplot() +
  geom_col(stacked, position = 'stack', mapping = aes(x= year_cut, y= deaths_mean, fill = age)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::comma) +
  plotting_theme+
  theme(legend.position = 'none') +
  labs(title = '',
       x= '',
       y= 'Deaths averted')

ggarrange(p3, p4, ncol = 2)




p5<- ggplot(table_2) +
  geom_line(mapping = aes(x= year, y= cohort_size, color = age)) +
  #geom_ribbon(aes(x= year, ymin = deaths_lower, ymax = deaths_upper, fill= age) , alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::comma) +
  plotting_theme+
    theme(legend.position = 'none')+
  labs(title = '',
       x= 'Year',
       y= 'Cohort size')

props<- table_2 |>
  group_by(year) |>
  summarise(total_pop = sum(cohort_size))
table_2<- merge(table_2, props, by = 'year')
table_2<- table_2 |>
  mutate(prop_population = cohort_size/total_pop)

p6<- ggplot(table_2) +
  geom_line(mapping = aes(x= year, y= prop_population, color = age)) +
  #geom_ribbon(aes(x= year, ymin = deaths_lower, ymax = deaths_upper, fill= age) , alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::comma) +
  plotting_theme+
  labs(title = '',
       x= 'Year',
       y= 'Proportion of total population')

ggarrange(p5, p6, ncol = 2)

# calculate figures summarized by age
table_2<- vax_age |>
  filter(scenario== 'Routine') |>
  group_by(scenario, age, run_id) |>
  summarise(cases_averted = sum(cases_averted),
            deaths_averted = sum(deaths_averted),
            dalys_averted = sum(dalys_averted),
            .groups = 'keep') |>
  data.table()  |>
    ungroup() |>
    group_by(scenario, age) |>
    summarise(
      cases_mean = mean(cases_averted, na.rm = TRUE),
      cases_lower = quantile(cases_averted, 0.025, na.rm = TRUE),
      cases_upper = quantile(cases_averted, 0.975, na.rm = TRUE),
      deaths_mean = mean(deaths_averted, na.rm = TRUE),
      deaths_lower = quantile(deaths_averted, 0.025, na.rm = TRUE),
      deaths_upper = quantile(deaths_averted, 0.975, na.rm = TRUE),
      .groups= 'keep'
    ) |>
    unique()


  get_values(table_2 |> filter(age== '[0, 5]') |> pull(cases_averted))

# pull cases averted summarized over age
table_2<- vax_age |>
  filter(scenario== 'Routine') |>
  group_by(scenario, age) |>
  summarise(cases_averted = sum(cases_averted),
            deaths_averted = sum(deaths_averted),
            dalys_averted = sum(dalys_averted),
            .groups = 'keep') |>
  data.table()



##Figure 3: Age pattern pre- and post- vaccination
age_pattern <- total |>
  group_by(scenario, age) |>
  summarise(
    cases = sum(cases),
    deaths = sum(deaths),
    .groups = "keep"
  )


p1<- ggplot(data = age_pattern, aes(x = age, y = cases, color = scenario)) +
  geom_line()+
    scale_y_continuous(labels = scales::comma) +
      scale_color_manual(values = cols)+
#labs(title = 'Incident cases by age and vaccine scenario') +
plotting_theme

p2<- ggplot(data = age_pattern, aes(x = age, y = deaths, color = scenario)) +
  geom_line()+
    scale_y_continuous(labels = scales::comma) +
      scale_color_manual(values = cols)+
#labs(title = 'Deaths by age and vaccine scenario') +
plotting_theme +
  theme( legend.position = 'none')



ggarrange(p1, p2, ncol = 2)
