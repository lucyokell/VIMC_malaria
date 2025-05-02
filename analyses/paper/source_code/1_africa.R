# Africa overall
library(vimcmalaria)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(ggpubr)
library(data.table)

cols<- c('#CE5137', '#1B4F72', '#45B39D', '#BA4A00')

plotting_theme<- theme_bw(base_size = 14) +
    theme( legend.position = 'bottom',
           strip.text.x = element_text(size = rel(0.8)),
           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
           text = element_text(family= 'TT Arial'),
           axis.ticks.y= element_blank(),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank())

get_values <- function(col) {
  mean <- mean(col, na.rm = TRUE)
  lower <- quantile(col, 0.025, na.rm = TRUE)
  upper <- quantile(col, 0.975, na.rm = TRUE)

  return(paste0(mean, "  (", lower, " - ", upper, ")"))
}
extract_site_data <- function(directory) {
  site_file <- readRDS(paste0("archive/process_inputs/", directory, "/site_file.rds"))

  # extract prevalence data
  prev <- site_file$prevalence |>
    filter(year %in% c(2000, 2019)) |>
    mutate(id = paste0(name_1, "_", urban_rural)) |>
    data.table()

  prevs <- dcast(prev, id ~ year, value.var = "pfpr") |>
    rename(
      pfpr_2000 = `2000`,
      pfpr_2019 = `2019`
    )


  # extract intervention coverage data for last available year
  intvns <- site_file$interventions |>
    filter(year == 2022) |>
    mutate(id = paste0(name_1, "_", urban_rural)) |>
    select(id, itn_use, tx_cov, smc_cov, pmc_cov, irs_cov)

  site_info <- merge(prevs, intvns, by = "id")

  return(site_info)
}

# metrics of interest:
completed <- completed_reports("process_country") |>
  dplyr::filter(
    description == "new_site_files",
    quick_run == FALSE
  ) |>
  dplyr::arrange(dplyr::desc(date_time)) |>
  dplyr::distinct(iso3c, scenario, quick_run, parameter_draw, description, .keep_all = TRUE) |>
  dplyr::arrange(iso3c, scenario, parameter_draw)



# pull country level output -------------- --------------------------------------unique
default<- compile_final_outputs('new_site_files')
default <- default |>
  filter(scenario == "proxy" | scenario == "no-vaccination") |>
  data.table()


default[scenario %like% 'proxy', scenario := 'Routine']
default[scenario == 'no-vaccination', scenario := 'Control']
total<- copy(default)


# generate a table of annual cases and deaths for both scenarios (for supplement)
supplement<- total |>
  filter(run_id < 51) |>
  group_by(scenario, run_id, year) |>
  summarise(cases = sum(cases, na.rm= TRUE),
            deaths = sum(deaths, na.rm = TRUE),
            .groups = 'keep') |>
  ungroup() |>
  group_by(scenario, year) |>
  summarise(cases = median(cases, na.rm = TRUE),
            deaths= median(deaths, na.rm = TRUE),
            .groups = 'keep') |>
  data.table()
supplement2<- dcast(supplement, year ~ scenario, value.var = c('cases', 'deaths'))
write.csv(supplement2, 'analyses/paper/tables/scenario_outputs_supplement.csv')

# RESULTS SECTION
  # calculate cases averted by year and age
  novax<- total |> filter(scenario == 'Control') |>
    rename(cases_novax = cases,
           deaths_novax = deaths,
           dalys_novax = dalys) |>
    select(year, age, run_id, country, cases_novax, deaths_novax, dalys_novax)

  vax<- total |> filter(scenario == 'Routine')
  vax<- merge(vax, novax, by = c('country', 'age', 'year', 'run_id'))

  vax<- vax |>
    filter(run_id < 51) |>
    #filter(age< 6) |>
    #filter(year < 2051) |>
    mutate(cases_averted = cases_novax- cases,
           deaths_averted = deaths_novax- deaths,
           dalys_averted = dalys_novax - dalys)


  # overall impact---------------------------------------------------------------------------------------
  # baseline scenario

  bl<- vax |>
    group_by(scenario,run_id) |>
    summarise(cases_baseline = sum(cases_novax), # pull ivory coast site file
              deaths_baseline = sum(deaths_novax),
              .groups = 'keep')

#baseline cases and deaths
  get_values(bl |> filter(scenario == 'Routine') |> pull(cases_baseline))
  get_values(bl |> filter(scenario == 'Routine') |> pull(deaths_baseline))


  africa<- vax |>
    group_by(scenario, run_id) |>
    summarise(cases_averted = sum(cases_averted),
              deaths_averted = sum(deaths_averted),
              dalys_averted = sum(dalys_averted),
              cases_novax = sum(cases_novax),
              deaths_novax = sum(deaths_novax),
              dalys_novax = sum(dalys_novax),
              .groups = 'keep') |>
    ungroup() |>
    group_by(scenario) |>
    filter(scenario!= 'Control') |>
    mutate(percent_cases_averted = cases_averted/cases_novax,
           percent_deaths_averted = deaths_averted/ deaths_novax,
           percent_dalys_averted = dalys_averted/ dalys_novax) |>
    data.table()
  print(africa)




# cases and deaths averted 
get_values(africa |> filter(scenario == 'Routine') |> pull(cases_averted))
get_values(africa |> filter(scenario == 'Routine') |> pull(deaths_averted))

# these numbers present X% of cases and deaths
get_values(africa |> filter(scenario == 'Routine') |> pull(percent_cases_averted))
get_values(africa |> filter(scenario == 'Routine') |> pull(percent_deaths_averted))


total_deaths<- vax |>
  group_by(run_id, year, scenario) |>
  summarise(deaths = sum(deaths), .groups= 'keep')

get_values(total_deaths |> filter(scenario == 'Routine', year == 2100) |> pull(deaths))

# FIGURE 1: line plot of cases + deaths averted over time
fig_1<- vax |>
  group_by(scenario,year, run_id) |>
  summarise(cases_averted = sum(cases_averted),
            deaths_averted = sum(deaths_averted),
            .groups = 'keep') |>
  ungroup() |>
  group_by(scenario, run_id) |>
  mutate(cumulative_cases = cumsum(cases_averted),
         cumulative_deaths = cumsum(deaths_averted)) |>
  data.table() |>
  ungroup() |>
  group_by(scenario, year) |>
  mutate(cases_mean = mean(cases_averted, na.rm= TRUE),
        cases_lower = quantile(cases_averted, 0.025, na.rm= TRUE),
      cases_upper =  quantile(cases_averted, 0.975, na.rm= TRUE),
      cases_cumulative_mean = mean(cumulative_cases, na.rm= TRUE),
      cases_cumulative_lower = quantile(cumulative_cases, 0.025, na.rm= TRUE),
      cases_cumulative_upper =  quantile(cumulative_cases, 0.975, na.rm= TRUE),
      deaths_mean = mean(deaths_averted, na.rm= TRUE),
      deaths_lower = quantile(deaths_averted, 0.025, na.rm= TRUE),
    deaths_upper =  quantile(deaths_averted, 0.975, na.rm= TRUE),
    deaths_cumulative_mean = mean(cumulative_deaths, na.rm= TRUE),
    deaths_cumulative_lower = quantile(cumulative_deaths, 0.025, na.rm= TRUE),
    deaths_cumulative_upper =  quantile(cumulative_deaths, 0.975, na.rm= TRUE))


p1<- ggplot(fig_1) +
  geom_line(aes(x = year, y = cases_mean/1000000, color = scenario)) +
  geom_ribbon(aes(x= year, ymin= cases_lower/1000000, ymax= cases_upper/1000000, fill= scenario), alpha = 0.25) +
  scale_color_manual(values = '#CE5137') +
  scale_fill_manual(values = '#CE5137') +
  plotting_theme +
  theme(legend.position = 'none') +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "",
    y = "Cases averted (millions)",
    x = ""
  )

p2<- ggplot(fig_1) +
  geom_line(aes(x = year, y = deaths_mean/1000, color = scenario)) +
  geom_ribbon(aes(x= year, ymin= deaths_lower/1000, ymax= deaths_upper/1000, fill= scenario), alpha = 0.25)+
  scale_color_manual(values = '#1B4F72') +
  scale_fill_manual(values = '#1B4F72') +
  plotting_theme +
    theme(legend.position = 'none') +
  scale_y_continuous(labels = scales::comma) +
  labs(
    y = "Deaths averted (thousands)",
    x = "Year"
  )

ggarrange(p1, p2, nrow= 2)


p3<- ggplot(fig_1) +
  geom_line(aes(x = year, y = cases_cumulative_mean/1000000, color = scenario)) +
  geom_ribbon(aes(x= year, ymin= cases_cumulative_lower/1000000, ymax= cases_cumulative_upper/1000000, fill= scenario), alpha = 0.25) +
  scale_color_manual(values = '#CE5137') +
  scale_fill_manual(values = '#CE5137') +
  plotting_theme +
  theme(legend.position = 'none') +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "",
    y = "Cumulative cases averted (millions)",
    x = ""
  )

p4<- ggplot(fig_1) +
  geom_line(aes(x = year, y = deaths_cumulative_mean/1000, color = scenario)) +
  geom_ribbon(aes(x= year, ymin= deaths_cumulative_lower/1000, ymax= deaths_cumulative_upper/1000, fill= scenario), alpha = 0.25)+
  scale_color_manual(values = '#1B4F72') +
  scale_fill_manual(values = '#1B4F72') +
  plotting_theme +
  theme(legend.position = 'none') +
  scale_y_continuous(labels = scales::comma) +
  labs(
    y = "Cumulative deaths averted (thousands)",
    x = "Year")


ggarrange(p3, p4, nrow= 2)
ggarrange(p1, p3, p2, p4, nrow= 2, ncol= 2)
