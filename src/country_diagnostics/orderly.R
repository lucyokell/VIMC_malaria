# produce country diagnostic report ------------------------------------------------------------------
orderly2::orderly_parameters(iso3c = NULL,
                             description = NULL,
                             population = NULL,
                             burnin = NULL,
                             parameter_draw = NULL,
                             scenario = NULL,
                             quick_run= NULL)


orderly2::orderly_description('Produce country diagnostic report for site')


print("Country ISO = ")
print(iso3c)
print("Scenario = ")
print(scenario)

library(dplyr)
library(data.table)
lapply(list.files('functions/', full.names = T), source)

# pull processed output
intvn_metadata<- orderly2::orderly_dependency("process_country",
                                              "latest(parameter:iso3c == this:iso3c 
                                                          && parameter:population == this:population
                                                          && parameter:description == this:description
                                                          && parameter:scenario == this:scenario
                                                          && parameter:parameter_draw == this:parameter_draw
                                                          && parameter:quick_run == this:quick_run)",
                                              c(intvn_output.rds = "country_output.rds"))

# pull baseline output
baseline_scenario_name<- 'no-vaccination'
baseline_metadata<- orderly2::orderly_dependency("process_country",
                                                 "latest(parameter:iso3c == this:iso3c
                                                          && parameter:population == this:population
                                                          && parameter:description == this:description
                                                          && parameter:scenario == environment:baseline_scenario_name
                                                          && parameter:parameter_draw == this:parameter_draw
                                                          && parameter:quick_run == this:quick_run)",
                                                 c(baseline_output.rds = "country_output.rds"))


# pull doses output
doses_metadata<- orderly2::orderly_dependency("process_country",
                                              "latest(parameter:iso3c == this:iso3c 
                                                          && parameter:population == this:population
                                                          && parameter:description == this:description
                                                          && parameter:scenario == this:scenario
                                                          && parameter:parameter_draw == this:parameter_draw
                                                          && parameter:quick_run == this:quick_run)",
                                              c(doses_output.rds = "country_doses_per_year.rds"))

# VIMC inputs
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(site_file.rds = "vimc_input.rds"))

# read in VIMC inputs ----------------------------------------------------------
site_data<- readRDS('site_file.rds')
vimc_input<- readRDS('vimc_input.rds')
coverage_data<- vimc_input$coverage_input
mort<- vimc_input$mort_rate_input
population_data<- vimc_input$population_input_all_age

# read in model inputs ---------------------------------------------------------
intvn_output<- readRDS('intvn_output.rds')
baseline_output<- readRDS('baseline_output.rds')
doses_output <- readRDS('doses_output.rds')



format_for_country_report<- function(){
  # bind intervention and baseline outputs together
  processed_output<- rbind(baseline_output, intvn_output, fill = T)
  processed_output<- data.table(processed_output)
  processed_output<- processed_output[scenario!= TRUE]
  
  # make summary output for all ages
  agg_output<- processed_output |>
    group_by(year, scenario) |>
    summarise(cases = sum(cases),
              dalys = sum(dalys),
              deaths = sum(deaths),
              pre_scaled_cases = sum(pre_scaled_cases),
              cohort_size = sum(cohort_size)) |>
    mutate(mortality = deaths/cohort_size,
           clinical = cases/ cohort_size,
           dalys_pp = dalys/ cohort_size) 
  
  
  # subset VIMC inputs to country of interest ------------------------------------
  intvn_scenario<- scenario
  coverage_data<- coverage_data |>
    filter(scenario == intvn_scenario,
           country_code == iso3c)
  
  population_data<- population_data |>
    filter(country_code == iso3c)
  message('read inputs successfully')
  
  # pull some numbers out for the report
  # Total cases in the country in WMR 2020 
  cases_wmr_2020 <- site_data$cases_deaths$wmr_cases[which(site_data$cases_deaths$year==2020)]
  
  # Total modelled cases in 2020 were
  cases_mod_2020 <- baseline_output |>
    filter(year==2020) |>
    summarise(cases = sum(cases)) |>
    pull(cases)
  
  # Up to 2100, xxxx % of cases were averted by vaccine
  percent_cases_averted_2100 <- 
    round(100 * (1 - sum(intvn_output$cases) / sum(baseline_output$cases)), 2)
  
  # MAP PfPR 2-10 2023 was... TO DO EXTRACT PREV IN PROCESS SITE AND COUNTRY TO POP WEIGHT.
  #site_data$prevalence$pfpr
  
  # Cases averted per 100,000 vaccinated were:
  outcomes_averted <- round(pull_outcomes_averted_per_100k_vacc(intvn_output , baseline_output , doses_output ))
  
  return(list('outcomes_averted' = outcomes_averted, 'agg_output' = agg_output, ))
}


# render report
rmarkdown::render(input= 'diagnostic_report_country.Rmd',
                  output_file = 'country_diagnostic_report',
                  output_format = 'html_document',
                  params= list('iso3c' = iso3c,
                               'description'= description,
                               'population' = population,
                               'quick_run' = quick_run,
                               'burnin' = burnin,
                               'parameter_draw' = parameter_draw,
                               'processed_output' = processed_output,
                               'site_data' = site_data,
                               'coverage_data' = coverage_data,
                               'population_data' = population_data,
                               'outcomes_averted' = outcomes_averted,
                               'agg_output' = agg_output,
                               'mort' = mort))

