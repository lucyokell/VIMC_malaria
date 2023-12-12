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

source('diagnostic_functions.R')

# pull processed output
intvn_metadata<- orderly2::orderly_dependency("process_site",
                                              "latest(parameter:iso3c == this:iso3c 
                                               && parameter:site_name == this:site_name 
                                               && parameter:ur == this:ur 
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
                                              && parameter:description == this:description
                                              && parameter:scenario == environment:baseline_scenario_name
                                              && parameter:parameter_draw == this:parameter_draw
                                              && parameter:quick_run == this:quick_run)",
                                              c(baseline_output.rds = "processed_output.rds"))
# plotting inputs for intervention scenario
plotting_input<- orderly2::orderly_dependency("process_site",
                                          "latest(parameter:iso3c == this:iso3c 
                                          && parameter:site_name == this:site_name 
                                          && parameter:ur == this:ur 
                                          && parameter:description == this:description
                                          && parameter:scenario == this:scenario
                                          && parameter:parameter_draw == this:parameter_draw
                                          && parameter:quick_run == this:quick_run)",
                                          c(plotting_input.rds = "plotting_input.rds"))

# VIMC inputs
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)",c(site_file.rds = "site_file.rds"))
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)",c(vimc_input.rds = "vimc_input.rds"))

site_data<- readRDS('site_file.rds')
vimc_input<- readRDS('vimc_input.rds')
plotting_input<- readRDS('plotting_input.rds')
intvn_output<- readRDS('intvn_output.rds')
baseline_output<- readRDS('baseline_output.rds')

# render report
rmarkdown::render(input= 'diagnostic_report.Rmd',
                  output_file = 'site_diagnostic_report',
                  output_format = 'html_document',
                  params= list('iso3c' = iso3c,
                               'description'= description,
                               'site_name' = site_name,
                               'ur' = ur,
                               'population' = report_input$population,
                               'quick_run' = quick_run,
                               'burnin' = report_input$burnin,
                               'parameter_draw' = parameter_draw,
                               'model_input' = report_input$model_input,
                               'raw_output' = report_input$raw_output,
                               'processed_output' = report_input$processed_output,
                               'agg_output' = report_input$agg_output,
                               'site_data' = report_input$site_data,
                               'coverage_data' = report_input$coverage_data,
                               'key_outcomes' = report_input$key_outcomes,
                               'mort' = report_input$mort))

