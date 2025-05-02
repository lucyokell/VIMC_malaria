# title   ethiopia recalibration
# author  Lydia Haile
# purpose recalibrate sites in Ethiopia to match MAP prevalence after turning IRS off

library(vimcmalaria)
library(cali)
library(malariasimulation)
library(dplyr)
library(data.table)
# vimc inputs ----
site_data<- readRDS('src/process_inputs/site_files/ETH_new_eir.rds')


# vimc inputs
inputs<- completed_reports('process_inputs') |>
  filter(iso3c == 'ETH') |>
  arrange(desc(date_time))
inputs<- unique(inputs, by = 'iso3c')
vimc_inputs<- readRDS(paste0('archive/process_inputs/', inputs$directory_name, '/vimc_input.rds'))
coverage_input<- vimc_inputs$coverage_input

# turn IRS coverage off
site_data$interventions$irs_cov = 0

# write recalibration function for each site of ethiopia


recalibrate<- function(site){
  
  
  params<- vimcmalaria::pull_input_params(site_name = site,
                                          ur = 'both',
                                          site_data= site_data,
                                          iso3c= 'ETH',
                                          coverage_data = coverage_input,
                                          scenario = 'no-vaccination',
                                          gfa= FALSE,
                                          parameter_draw = 0,
                                          quick_run = FALSE)
  
  
  
  summary_mean_pfpr_2_10 <- function (x) {
    message('calibrating')
    x<- data.table(x)
    # Calculate the PfPR2-10:
    prev_2_10 <- mean(x[timestep %in% c((10*365):(11*365))]$n_detect_pcr_730_3649/x[timestep %in% c((10*365):(11* 365))]$n_730_3649) # average over a year 
    
    # Return the calculated PfPR2-10:
    return(prev_2_10)
  }
  
  
  summary_mean_pfpr_2_10(output)
  # pull target pfpr from 2010 for corresponding site
  target_pfpr <- site_data$prevalence |> filter(year == 2010, name_1 == site) |> pull(pfpr)
  
  print(paste0('target pfpr ', target_pfpr ))
  
  # Add a parameter to the parameter list specifying the number of timesteps 
  simparams<- params$param_list
  simparams$timesteps <- 12 * 365
  
  # Establish a tolerance value:
  pfpr_tolerance <- 0.01
  
  # Set upper and lower EIR bounds for the calibrate function to check
  lower_EIR <- 0.01; upper_EIR <- 60
  
  output<- run_simulation(timesteps = simparams$timesteps, parameters = simparams)
  # Run the calibrate() function:
  cali_EIR <- calibrate(target = target_pfpr,
                        summary_function = summary_mean_pfpr_2_10,
                        parameters = simparams,
                        tolerance = pfpr_tolerance, 
                        low = lower_EIR, high = upper_EIR)
  
  print(paste0('calibrated EIR for site ', site, ' :', cali_EIR))
  simparams<- set_equilibrium(simparams, init_EIR = cali_EIR)
}


recalibrate('Afar')





cols <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7","#F0E442", "#0072B2", "#D55E00")
# pull map prevalence for ETH
map_prev<- site_data$prevalence
map_prev<- map_prev |>
  mutate(site_ur = paste0(name_1, '_', urban_rural))

# plot map prevalence by site
ggplot()+
  geom_line(data= map_prev, mapping = aes(x= year, y = pfpr), color = '56B4E9') +
  facet_wrap(~site_ur) +
  plotting_theme

plot<- raw |> filter(scenario == 'no-vaccination') |>
  mutate(year = as.integer(timestep/365)+ 1999) |>
  group_by(year, site_ur) |>
  summarise(pfpr = mean(n_detect_pcr_730_3649 /n_730_3649),
            .groups = 'keep')

# plot modelled prevalence by site
ggplot()+
  geom_line(data= plot, mapping = aes(x= year, y = pfpr), color = 'darkblue') +
  geom_line(data= map_prev, mapping = aes(x= year, y = pfpr), color = 'darkred') +
  facet_wrap(~site_ur) +
  #plotting_theme +
  scale_fill_manual(values = cols) +
  labs(title = 'MAP prevalence vs modelled prevalence, pre-IRS fix')


plot<- raw_new |> filter(scenario == 'no-vaccination') |>
  mutate(year = as.integer(timestep/365)+ 1999) |>
  group_by(year, site_ur) |>
  summarise(pfpr = mean(n_detect_pcr_730_3649 /n_730_3649),
            .groups = 'keep')

# plot modelled prevalence by site
ggplot()+
  geom_line(data= plot, mapping = aes(x= year, y = pfpr), color = 'darkblue') +
  geom_line(data= map_prev, mapping = aes(x= year, y = pfpr), color = 'darkred') +
  facet_wrap(~site_ur) +
  plotting_theme +
  scale_fill_manual(values = cols) +
  labs(title = 'MAP prevalence vs modelled prevalence, post-IRS fix')
