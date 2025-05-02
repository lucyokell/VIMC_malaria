
# TITLE:   Create vaccine scenario
# PURPOSE: Devise vaccine scenario as a proxy for VIMC coverage given forecast sensitivites
#          Generates vaccine scenario by country based on DTP coverage
# AUTHOR   Lydia Haile
#################################################################################################################

#packages and objects
library(data.table)
library(ggplot2)
library(dplyr)

cols<- c('#CE5137', '#1B4F72', '#45B39D', '#BA4A00')

# data
population<- read.csv('src/process_inputs/vimc_inputs/demography/202409malaria-1_dds-202208_tot_pop_both.csv') # total population
pop_age<- read.csv('src/process_inputs/vimc_inputs/demography/202409malaria-1_dds-202208_int_pop_both.csv') #single year population
dtp<- read.csv('who_dtp_coverage.csv')
vimc<- read.csv('src/process_inputs/vimc_inputs/vaccine_coverage/coverage_202409malaria-1_malaria-r3-r4-default.csv') # vimc vaccine scenario

# data on vaccine introduction year by country
inputs<- data.table('iso3c' = c('AGO', 'BDI', 'BEN', 'BFA', 'CAF', 'CIV', 'CMR', 'COD', 'COG', 'ETH', 'GHA', 'GIN', 'GNB', 'KEN', 'LBR', 'MDG', 'MLI', 'MOZ', 'MRT', 'MWI', 'NER', 'NGA', 'SDN', 'SLE', 'SOM', 'SSD', 'TCD', 'TGO', 'TZA', 'UGA', 'ZMB'),
                    'intro_yr'= c(2026,  2025,   2024,  2024, 2024,  2024,   2024,  2024,  2026,   2025, 2019, 2025,   2026,  2019,  2024,  2026,  2025, 2024,  2026,   2019, 2024,   2024,  2024,  2024,   2026, 2024, 2024, 2026,   2026,  2025,  2026))


#pull out the number of children from 0-1 yrs old from 2023-2030--------------------------------------------------
# pop_age<- pop_age |>
  #filter(year %in% c(2026:2030)) |>
  #filter(age_from < 2)

# sum the total population of children under 2 from 2026-2030 (the comparator period for gavi estimate)
total_par<- pop_age |>
  group_by(age_from) |>
  summarise(value = sum(value))

# estimate par based on moderate-to-high transmission cutoff (will be smaller)--------------------------------------
# read in population data from site files
read_pop<- function(file_name){

  message(file_name)
  site_data<- readRDS(file_name)

  site_data$prevalence<- site_data$prevalence |>
    dplyr::filter(year == 2019) |>
    mutate(run_model = ifelse(pfpr > 0.10, TRUE, FALSE))

  # make exceptions for Madagascar, Ethiopia, and Sudan
  # hardcode for time's sake but operationalize later
  if(unique(site_data$country) == 'MDG'){

    site_data$prevalence <- site_data$prevalence |>
      mutate(run_model = ifelse(name_1 == 'Toliary', TRUE, run_model))
  }

  if(unique(site_data$country) == 'ETH'){

    site_data$prevalence <- site_data$prevalence |>
      mutate(run_model = ifelse(name_1 %like% 'Gambela', TRUE, run_model))


  }

  if(unique(site_data$country) == 'SDN'){

    site_data$prevalence <- site_data$prevalence |>
      mutate(run_model = ifelse(name_1 == 'South Darfur', TRUE, run_model)) |>
      mutate(run_model = ifelse(name_1 == 'West Kurdufan', TRUE, run_model))
  }

  prevalence<- site_data$prevalence |>
    select(name_1, urban_rural, iso3c, run_model)

  #now pull annual population and merge this on
  pop<- site_data$population


  pop<- merge(pop, prevalence, by = c('name_1', 'urban_rural', 'iso3c'))

  return(pop)
}

files<- list.files('src/process_inputs/site_files/', full.names = T)
pops<- rbindlist(lapply(files, read_pop))

#pull out identifying info for countries, clean names for merging
ids<- pops |> select(country, iso3c) |> unique()
ids <- ids |>
  mutate(country = ifelse(country== 'Republic of the Congo', 'Congo', country)) |>
  mutate(country = ifelse(country == 'Tanzania', 'United Republic of Tanzania', country))

# subset to moderate to high transmission areas
pop_ar<- pops |>
  filter(run_model== TRUE)

# calculate target population by admin 1 unit through 2100
# this assumes that the national age distribution is the same as the subnational age distribution for any given year
#pull out the number of children from 0-1 yrs old from 2023-2030
pop_under_1<- pop_age |>
  filter(age_from < 1) |>
  #filter(year %in% c(2026:2030)) |>
  group_by(country, year, country_code) |>
  rename(iso3c= country_code) |>
  summarise(value_1 = sum(value), .groups = 'keep')

pop_at_2<- pop_age |>
  filter(age_from == 1) |>
  #filter(year %in% c(2026:2030)) |>
  group_by(country, year, country_code) |>
  rename(iso3c= country_code) |>
  summarise(value_2 = sum(value), .groups = 'keep')

# pull total population
total_pop<- population |>
  group_by(country, year, country_code) |>
  rename(iso3c= country_code) |>
  #filter(year %in% c(2026:2030)) |>
  summarise(total= sum(value), .groups = 'keep')

pop_prop<- merge(pop_under_1, total_pop, by = c('iso3c', 'year'))
pop_prop<- merge(pop_at_2, pop_prop, by = c('iso3c', 'year')) |>
  mutate(prop_1 = value_1/total,
         prop_2 = value_2/total ) |>
  select(country, year, iso3c, prop_1, prop_2)

# merge onto site specific data set
pops<- merge(pop_ar, pop_prop, by = c('iso3c', 'year')) |>
  mutate(pop_1 = pop * prop_1,
        pop_2 = pop * prop_2)

summary<- pops|>
  summarise(pop_1 = sum(pop_1),
            pop_2 = sum(pop_2))


#clean dtp covearge data set  ----------------------------------------------------------------------------
unique(dtp$DIM_TIME)
dtp<- data.table(dtp)
dtp<- dtp |>
  filter(DIM_TIME == 2022 ) |>
  rename(country = GEO_NAME_SHORT,
          year = DIM_TIME,
        coverage = RATE_PER_100_N) |>
  select(year, country, coverage)
dtp<- merge(dtp, ids, by = 'country')


# generate vaccine scenario for each country  --------------------------------------------------------------------------
devise_vaccine_scenario <- function(iso) {
  message(iso)
  # pull the year of vaccine introduction
  intro <- inputs |>
    filter(iso3c == iso) |>
    pull(intro_yr)

  if (intro < 2022) { intro <- 2022 } # for MVIP countries, vaccine coverage pre-2022 is already in the malariasimulation site files and this is for RTS,S

  dtp_coverage <- dtp |>
    filter(iso3c == iso) |>
    pull(coverage)

  # interpolate between 65% and 100% of this value from 2022 to 2030
  # Define the years and the corresponding values
  years <- intro:2100
  value_intro <- (.4 * dtp_coverage) / 100 # Value in intro year
  value_2030 <-  (.85* dtp_coverage) / 100 # Value in 2022
  value_2060<-  0.95 #assume every country reaches 95% threshold by 2060
#
#   # make sure the maximum is no more than 95%
  if (value_2030 > 0.90) {value_2030 <- 0.90 } 

  # Perform linear extrapolation
  values <- approx(x = c(intro, 2030, 2060, 2100), y = c(value_intro, value_2030, value_2060, 0.95), xout = years)$y

  # Combine years and values into a data frame for easier viewing
  result <- data.table(year = years, coverage = values, iso3c = iso)

  if (nrow(result[year == 2022]) != 1) {
    append <- data.frame(
      year = c(2022:intro - 1),
      coverage = rep(0, length(c(2022:intro - 1))),
      iso3c = rep(iso, length(c(2022:intro - 1)))
    )

    result <- data.table(rbind(append, result))
  }


  # add in booster coverage, which should be 80% of total coverage in the preceding year
  result <- result |>
    mutate(booster = .80) |>
    mutate(booster = ifelse(year <= intro, 0, booster)) |> # if it is the year of introduction, no booster
    data.table()

  return(result)
}

# run for all countries
proxy_scenario<- rbindlist(lapply(ids$iso3c, devise_vaccine_scenario))

calculate_children_protected<- function(coverage){

  doses<- merge(pops, coverage, by = c('iso3c', 'year'), all.x = TRUE)
  # calculate doses at the admin- 1 level
  # assuming that coverage is uniform in all moderate to high transmission areas
  #first calculate the number of children protected with the primary series (first 3 doses)
  # number of children protected with the booster dose corresponds to 
  # booster coverage * children_protected_primary 

  doses<- doses |>
    group_by(year) |>
    mutate(coverage= ifelse(is.na(coverage), 0, coverage),
           booster= ifelse(is.na(booster), 0, booster)) |>
    mutate(children_protected_primary = pop_1 * coverage) |>
    mutate(children_fully_protected = children_protected_primary * booster) |>
    summarise(children_protected = sum(children_fully_protected))


  return(doses)
}
calculate_children_protected_summary<- function(coverage){

  doses<- merge(pops, coverage, by = c('iso3c', 'year'), all.x = TRUE)
  # calculate doses at the admin- 1 level
  # assuming that coverage is uniform in all moderate to high transmission areas
  doses<- doses |>
    filter(year %in%  c(2026:2030)) |>
      mutate(coverage= ifelse(is.na(coverage), 0, coverage),
    booster= ifelse(is.na(booster), 0, booster)) |>
mutate(children_protected_primary = pop_1 * coverage) |>
mutate(children_fully_protected = children_protected_primary * booster) |>
summarise(children_protected = sum(children_fully_protected))

  return(doses$children_protected)
}

#pull in VIMC scenario for comparison
vimc_cov <- vimc |>
  #filter(age_first== 0) |>
  rename(iso3c = country_code) |>
  mutate(scenario = 'vimc') |>
 mutate(coverage = coverage / proportion_risk) |>  # convert coverage to per population instead of per population at risk
  mutate(coverage= ifelse(is.na(coverage), 0, coverage)) |>
  select(year,  coverage, iso3c, age_first, scenario)

initial<- vimc_cov |> filter(age_first == 0)   |> select(-age_first) 
boosted<- vimc_cov |> filter(age_first == 1) |> rename(booster = coverage) |> select(-age_first)
vimc_cov<- merge(initial, boosted, by = c('year', 'iso3c', 'scenario'))

cov_scenario<- calculate_children_protected(proxy_scenario)
vimc_scenario<- calculate_children_protected(vimc_cov)

cov_scenario<- calculate_children_protected_summary(proxy_scenario)
vimc_scenario<- calculate_children_protected_summary(vimc_cov)

proxy_scenario<- proxy_scenario |>
  mutate(scenario= 'proxy') |>
  select(year, coverage, booster, iso3c, scenario)

cov<- data.table(rbind(vimc_cov, proxy_scenario))


compare_coverage<- function(iso){

  message(iso)
  p<- ggplot(cov[iso3c == iso], mapping = aes(x= year, y= coverage, color= scenario)) +
    geom_line()+
    geom_vline(mapping= aes(xintercept= 2022), linetype= 'dotted')+
      geom_vline(mapping= aes(xintercept= 2030), linetype= 'dotted')+
    scale_color_manual(values = cols) +
    labs(title = 'Comparison of coverage scenarios',
          subtitle = iso)


  print(p)

}

pdf('comparison_scenario_final2.pdf')

lapply(unique(cov$iso3c), compare_coverage)

dev.off()


saveRDS(proxy_scenario, 'proxy_scenario.rds')
