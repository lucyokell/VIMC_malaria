# process country --------------------------------------------------------------
orderly2::orderly_parameters(iso3c = NULL,
                             description = NULL,
                             population = NULL,
                             burnin = NULL,
                             parameter_draw = NULL,
                             scenario = NULL,
                             quick_run = NULL)


library(postie)
library(dplyr)
library(data.table)

source('remove_zero_eirs.R')

orderly2::orderly_description('Process model outputs')
orderly2::orderly_artefact('Processed output', 'country_output.rds')


# read in model outputs for all sites in country
orderly2::orderly_dependency("process_inputs",
                             "latest(parameter:iso3c == this:iso3c)",
                             c(site_file.rds = "site_file.rds"))


site_data <- readRDS('site_file.rds')


sites<- site_data$sites

sites<- remove_zero_eirs(iso3c, sites, site_data$eir)

Encoding(sites$name_1) <- "UTF-8"
sites$name_1<- iconv(sites$name_1, from="UTF-8", to="ASCII//TRANSLIT")

output<- data.table()
doses<-data.table()

for (i in 1:nrow(sites)) {
    
    site<- sites[i,]
    site_name<- site$name_1
    ur<- site$urban_rural


    message(i)
    metadata<-orderly2::orderly_dependency("process_site", quote(latest(parameter:iso3c == this:iso3c &&
                                                                   parameter:description == this:description &&
                                                                   parameter:population == this:population &&
                                                                   parameter:scenario == this:scenario &&
                                                                   parameter:burnin == this:burnin &&
                                                                   parameter:site_name == environment:site_name &&
                                                                   parameter:ur == environment:ur &&
                                                                   parameter:parameter_draw == this:parameter_draw &&
                                                                   parameter:quick_run == this:quick_run)),
                                           c('processed_output_${site_name}_${ur}.rds' = "processed_output.rds"))
    dt<- readRDS(metadata$files$here)
    output<- rbind(output, dt, fill = T)
    
    
    scenario<-dt$scenario[1]
    
    if(scenario!="no-vaccination") {
      metadata<-orderly2::orderly_dependency("process_site", quote(latest(parameter:iso3c == this:iso3c &&
                                                                            parameter:description == this:description &&
                                                                            parameter:population == this:population &&
                                                                            parameter:scenario == this:scenario &&
                                                                            parameter:burnin == this:burnin &&
                                                                            parameter:site_name == environment:site_name &&
                                                                            parameter:ur == environment:ur &&
                                                                            parameter:parameter_draw == this:parameter_draw &&
                                                                            parameter:quick_run == this:quick_run)),
                                             c('doses_per_year_${site_name}_${ur}.rds' = "doses_per_year.rds"))
      
      doses_site<- readRDS(metadata$files$here)
      doses<- rbind(doses, doses_site, fill = T)
      
    }
}


orderly2::orderly_dependency("process_inputs",
                             "latest(parameter:iso3c == this:iso3c)",
                             c(population_input_single_yr.rds = "population_input_single_yr.rds"))

pop<- readRDS("population_input_single_yr.rds")

# sum cases up to country level ------------------------------------------------
dt<- copy(output)
dt<- data.table(dt)

dt[, `:=` (
  cases = sum(cases),
  deaths = sum(deaths),
  dalys = sum(dalys)),
by = c('age', 'year', 'scenario')]


# remove cohort size, because for sites with some unmodelled locations, sum of cohort size != national population

dt<- dt |> 
  select(-cohort_size)

dt <- unique(dt, by = c('age', 'year', 'scenario'))

pop<- pop |>
  rename(age = age_from,
         cohort_size = value) |>
  select(year, age, cohort_size)


dt<- merge(dt, pop, by =c('age', 'year'))


# calculate rates --------------------------------------------------------------
dt[, `:=` (
  clinical = NULL,
  mortality = NULL,
  dalys_pp = NULL,
  site_name = iso3c,
  urban_rural = 'total'
)]


# # scale outputs based on cases from WMR from 2000-2020
# # first sum cases by year (across all ages) in model output and compare
pre_scale<- dt |>
  group_by(year) |>
  summarise(cases = sum(cases))

#average site file cases across last three years
site_file_cases<- data.table::data.table(site_data$cases_deaths[, c('year', 'wmr_cases')])
site_file_cases<- site_file_cases[year >= 2018]
average_value<- mean(site_file_cases$wmr_cases)

# calculate ratio in comparison to year 2020 cases in output
output_cases<- pre_scale |>
  filter(year == 2020) |>
  select(cases)


ratio<- average_value/output_cases$cases

# test this scaling with an extra column and review
dt<- dt |>
  mutate(pre_scaled_cases = cases)

dt<- dt |>
  mutate(cases = cases * ratio)

dt<- dt|>
  mutate(clinical= cases/cohort_size,
         mortality = deaths/ cohort_size,
         dalys_pp = dalys/ cohort_size) |>
  select(-site_name, -urban_rural)

if(scenario!="no-vaccination") {
  doses_per_year <- doses |>
    dplyr::group_by(year) |>
    summarise(doses=sum(doses)) |>
    mutate(scenario=scenario)
  
  saveRDS(doses_per_year, 'country_doses_per_year.rds')
}

# save outputs  ----------------------------------------------------------------
saveRDS(dt, 'country_output.rds')
message('done')

