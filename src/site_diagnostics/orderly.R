# produce diagnostic report ------------------------------------------------------------------
orderly2::orderly_parameters(iso3c = NULL,
                             description = NULL,
                             site_name = NULL,
                             ur = NULL,
                             parameter_draw = NULL,
                             scenario = NULL,
                             quick_run = NULL)


orderly2::orderly_description('Produce diagnostic report for site')

library(dplyr)
library(data.table)

lapply(list.files('functions/', full.names = T), source)

# pull processed output
intvn_metadata<- orderly2::orderly_dependency("process_site",
                                              "latest(parameter:iso3c == this:iso3c 
                                                          && parameter:site_name == this:site_name 
                                                          && parameter:ur == this:ur 
                                                          && parameter:population == this:population
                                                          && parameter:description == this:description
                                                          && parameter:scenario == this:scenario
                                                          && parameter:parameter_draw == this:parameter_draw
                                                          && parameter:quick_run == this:quick_run)",
                                              c(intvn_output.rds = "processed_output.rds"))

# pull baseline output
baseline_scenario_name<- 'no-vaccination'
baseline_metadata<- orderly2::orderly_dependency("process_site",
                                              "latest(parameter:iso3c == this:iso3c
                                                          && parameter:site_name == this:site_name
                                                          && parameter:ur == this:ur
                                                          && parameter:population == this:population
                                                          && parameter:description == this:description
                                                          && parameter:scenario == environment:baseline_scenario_name
                                                          && parameter:parameter_draw == this:parameter_draw
                                                          && parameter:quick_run == this:quick_run)",
                                              c(baseline_output.rds = "processed_output.rds"))
# raw model output for intervention
raw_output<- orderly2::orderly_dependency("process_site",
                                          "latest(parameter:iso3c == this:iso3c 
                                                          && parameter:site_name == this:site_name 
                                                          && parameter:ur == this:ur 
                                                          && parameter:population == this:population
                                                          && parameter:description == this:description
                                                          && parameter:scenario == this:scenario
                                                          && parameter:parameter_draw == this:parameter_draw
                                                          && parameter:quick_run == this:quick_run)",
                                          c(raw_model_output.rds = "raw_model_output.rds"))
# site file inputs
vaccine_coverage_input<- orderly2::orderly_dependency("set_parameters",
                                                      "latest(parameter:iso3c == this:iso3c 
                                                          && parameter:site_name == this:site_name 
                                                          && parameter:ur == this:ur 
                                                          && parameter:population == this:population
                                                          && parameter:description == this:description
                                                          && parameter:scenario == this:scenario
                                                          && parameter:parameter_draw == this:parameter_draw
                                                          && parameter:quick_run == this:quick_run)",
                                                      c(vaccine_plot_input.rds = "vaccine_plot_input.rds"))
# VIMC inputs
orderly2::orderly_dependency("process_inputs",
                                        "latest(parameter:iso3c == this:iso3c)",
                                        c(site_file.rds = "site_file.rds"))

site_data<- readRDS('site_file.rds')
coverage_data<- readRDS('coverage_input.rds')
mort<- readRDS('mort_rate_input.rds')

# read in inputs to pass into report as parameters
model_input<- readRDS('vaccine_plot_input.rds')
raw_output<- readRDS('raw_model_output.rds')
intvn_output<- readRDS('intvn_output.rds')
baseline_output<- readRDS('baseline_output.rds')

# bind intervention and baseline outputs together
processed_output<- rbind(baseline_output, intvn_output, fill = T)
processed_output<- data.table(processed_output)
processed_output<- processed_output[scenario!= TRUE]

message('read inputs successfully')

site_data<-extract_site(site_file = site_data,
                                 site_name = site_name,
                                 ur = ur)

# make summary output for all ages
agg_output<- processed_output |>
  group_by(year, scenario) |>
  summarise(cases = sum(cases),
            dalys = sum(dalys),
            deaths = sum(deaths),
            cohort_size = sum(cohort_size)) |>
  mutate(mortality = deaths/cohort_size,
         clinical = cases/ cohort_size,
         dalys_pp = dalys/ cohort_size) 


key_outcomes<-pull_outcomes_averted_per_100k_vacc(raw_output, processed_output)

# summarize outputs by age as wee
# render report
rmarkdown::render(input= 'diagnostic_report.Rmd',
                  output_file = 'site_diagnostic_report',
                  output_format = 'html_document',
                  params= list('iso3c' = iso3c,
                               'description'= description,
                               'site_name' = site_name,
                               'ur' = ur,
                               'population' = population,
                               'quick_run' = quick_run,
                               'burnin' = burnin,
                               'parameter_draw' = parameter_draw,
                               'model_input' = model_input,
                               'raw_output' = raw_output,
                               'processed_output' = processed_output,
                               'agg_output' = agg_output,
                               'site_data' = site_data,
                               'coverage_data' = coverage_data,
                               'key_outcomes' = key_outcomes,
                               'mort' = mort))

