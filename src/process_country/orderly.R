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

source('country_functions.R')

orderly2::orderly_description('Process model outputs')
orderly2::orderly_artefact('Processed output', 'country_output.rds')

# read in model outputs for all sites in country
orderly2::orderly_dependency("process_inputs", "latest(parameter:iso3c == this:iso3c)", c(site_file.rds = "site_file.rds"))
orderly2::orderly_dependency("process_inputs","latest(parameter:iso3c == this:iso3c)", c(vimc_input.rds = "vimc_input.rds"))

site_data <- readRDS('site_file.rds')
vimc_input<- readRDS('vimc_input.rds')
pop<- vimc_input$population_input_single_yr

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
      metadata<-orderly2::orderly_dependency("analyse_site", quote(latest(parameter:iso3c == this:iso3c &&
                                                                            parameter:description == this:description &&
                                                                            parameter:population == this:population &&
                                                                            parameter:scenario == this:scenario &&
                                                                            parameter:burnin == this:burnin &&
                                                                            parameter:site_name == environment:site_name &&
                                                                            parameter:ur == environment:ur &&
                                                                            parameter:parameter_draw == this:parameter_draw &&
                                                                            parameter:quick_run == this:quick_run)),
                                             c('plotting_input_${site_name}_${ur}.rds' = "plotting_input.rds"))
      
      input<- readRDS(metadata$files$here)
      doses_site<- input$doses_per_year
      
      doses<- rbind(doses, doses_site, fill = T)
      
    }
}

# sum cases up to country level ------------------------------------------------

dt<- aggregate_outputs(output, pop)
dt<- scale_cases(dt, site_data)


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

