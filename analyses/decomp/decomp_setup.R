##  title quick decomp
##  author Lydia Haile
##  purpose run an ad-hoc decomp analysis to assess relative contribution of the following to VIMC outputs:
##  -fixed demography
##  - historical intercention coverage
##  - both or none


# packages  --------------------------------------------------------------------
library(data.table)
library(vimcmalaria)
library(ggplot2)
library(dplyr)
library(site)
library(malariasimulation)
# vimc inputs
site_file<- readRDS('src/process_inputs/site_files/COD_new_eir.rds')
inputs<- completed_reports('process_inputs') |>
  filter(iso3c == 'COD') |>
  arrange(desc(date_time))
inputs<- unique(inputs, by = 'iso3c')
vimc_input<- readRDS(paste0('archive/process_inputs/', inputs$directory_name, '/vimc_input.rds'))
site<- extract_site(site_file,
                    site_name = 'Ituri',
                    ur = 'both')

# plot prevalence to get sense of pfpr towards 2019/2020
ggplot2::ggplot(data= site$prevalence, mapping = aes(x= year, y= pfpr))+
  geom_point() +
  facet_wrap(~name_1)

# cluster setup ----------------------------------------------------------------
hipercow::hipercow_init(driver = 'windows')
hipercow::hipercow_provision()
hipercow::hipercow_environment_create(sources= 'analyses/decomp/run_and_save.R')
hipercow::hipercow_configuration()

#launch a model as baseline ----------------------------------------------------
mod<- pull_input_params(site_name = 'Ituri',
                                ur = 'both',
                                site_data = site_file,
                                coverage_data = vimc_input$coverage_input,
                                scenario = 'no-vaccination',
                                iso3c = 'COD',
                                parameter_draw = 0,
                                gfa = FALSE,
                                quick_run = FALSE)
test2<- hipercow::task_create_expr(run_and_save(model_input = mod, description = 'no-vaccine'))


# launch a model with fixed demography  ----------------------------------------
site_file_fixed<- copy(site_file)
fixed_dem<- site_file$demography |> filter(year == 2024)
fixed_dem <- fixed_dem |> select(-year) |> rename(new_mort = mortality_rate)
site_file_fixed$demography<- merge(site_file_fixed$demography, fixed_dem, by = c('iso3c', 'age_upper', 'country'))
site_file_fixed$demography <- site_file_fixed$demography |>
  mutate(mortality_rate = ifelse(year > 2024, new_mort, mortality_rate)) |> 
  select(-new_mort)
site_file_fixed$demography<- site_file_fixed$demography[order(site_file_fixed$demography$year, site_file_fixed$demography$age_upper), ]

mod_dem<- pull_input_params(site_name = 'Ituri',
                        ur = 'both',
                        site_data = site_file_fixed,
                        coverage_data = vimc_input$coverage_input,
                        scenario = 'malaria-r3-r4-bluesky',
                        iso3c = 'COD',
                        parameter_draw = 0,
                        gfa = FALSE,
                        quick_run = FALSE)

test<- hipercow::task_create_expr(run_and_save(model_input = mod_dem, description = 'fixed_demography'))
hipercow::task_log_watch(test)


# launch model with historical interevention coverage turned off  --------------
site_file_intvn<- copy(site_file)
intvns<- site_file_intvn$interventions
intvns<- intvns |>
  mutate(pmc_cov = 0,
         smc_cov = 0,
         tx_cov = 0,
         itn_use = 0,
         itn_input_dist = 0,
         irs_cov = 0)

site_file_intvn$interventions<- intvns
mod_intvn<- pull_input_params(site_name = 'Ituri',
                            ur = 'both',
                            site_data = site_file_intvn,
                            coverage_data = vimc_input$coverage_input,
                            scenario = 'malaria-r3-r4-bluesky',
                            iso3c = 'COD',
                            parameter_draw = 0,
                            gfa = FALSE,
                            quick_run = FALSE)


# Prepare a summary function that returns the endpoint PfPR2-10 from each simulation output: 
library(cali)
summary_mean_pfpr_2_10 <- function (x) {
  message('calibrating')
  x<- data.table(x)
  # Calculate the PfPR2-10:
  prev_2_10 <- mean(x[timestep %in% c((19*365):(20*365))]$n_detect_730_3649/x[timestep %in% c((19*365):(20* 365))]$n_730_3649) # average over a year 
  
  # Return the calculated PfPR2-10:
  return(prev_2_10)
}

# Establish a target PfPR2-10 value:
target_pfpr <- 0.5869696

# Add a parameter to the parameter list specifying the number of timesteps to
# simulate over. Note, increasing the number of steps gives the simulation longer
# to stablise/equilibrate, but will increase the runtime for calibrate().  
simparams<- copy(mod_intvn$param_list)
simparams$timesteps <- 20 * 365

# Establish a tolerance value:
pfpr_tolerance <- 0.01

# Set upper and lower EIR bounds for the calibrate function to check (remembering EIR is
# the variable that is used to tune to the target PfPR):
lower_EIR <- 15; upper_EIR <- 100

# Run the calibrate() function:
cali_EIR <- calibrate(target = target_pfpr,
                      summary_function = summary_mean_pfpr_2_10,
                      parameters = simparams,
                      tolerance = pfpr_tolerance, 
                      low = lower_EIR, high = upper_EIR)


mod_intvn$param_list <- set_equilibrium(mod_intvn$param_list, init_EIR = 65.79351)
test<- hipercow::task_create_expr(run_and_save(model_input = mod_intvn, description = 'no_historical_interventions'))

# Use the match_EIR_to_PfPR() function to return the EIR predicted to be required under the
# malariasimulation method:
malsim_EIR <- match_EIR_to_PfPR(x = target_pfpr)

# View the dataframe containing the EIR and matching PfPR2-10 values:
malSim_P2E


# run postprocessing
dem<- vimcmalaria::process_output(fixed_demography,
                              ur = 'both',
                              site_name = 'Ituri',
                              description = 'fixed_demography',
                              site_data = site_file,
                              vimc_input = vimc_input,
                              scenario = 'malaria-r3-r4-bluesky',
                              iso3c = 'COD',
                              parameter_draw = 0,
                              gfa = FALSE,
                              quick_run = FALSE)

bl<- vimcmalaria::process_output(baseline,
                                  ur = 'both',
                                  site_name = 'Ituri',
                                  description = 'baseline',
                                  site_data = site_file,
                                  vimc_input = vimc_input,
                                  scenario = 'malaria-r3-r4-bluesky',
                                  iso3c = 'COD',
                                  parameter_draw = 0,
                                  gfa = FALSE,
                                  quick_run = FALSE)

no_intvn<-  vimcmalaria::process_output(no_historical_interventions,
                                     ur = 'both',
                                     site_name = 'Ituri',
                                     description = 'no_historical_intvns',
                                     site_data = site_file,
                                     vimc_input = vimc_input,
                                     scenario = 'malaria-r3-r4-bluesky',
                                     iso3c = 'COD',
                                     parameter_draw = 0,
                                     gfa = FALSE,
                                     quick_run = FALSE)

no_vax<-  vimcmalaria::process_output(`no-vaccine`,
                                        ur = 'both',
                                        site_name = 'Ituri',
                                        description = 'no-vaccine',
                                        site_data = site_file,
                                        vimc_input = vimc_input,
                                        scenario = 'malaria-r3-r4-bluesky',
                                        iso3c = 'COD',
                                        parameter_draw = 0,
                                        gfa = FALSE,
                                        quick_run = FALSE)
comparison<- rbind(dem$processed_output, bl$processed_output) 
comparison<- rbind(comparison, no_intvn$processed_output)
comparison<- rbind(comparison, no_vax$processed_output)

comparison<- comparison |>
  group_by(description, year) |>
  summarise(cases = sum(cases),
            cohort_size = sum(cohort_size), 
            .groups = 'keep',
            )


ggplot(data = comparison, mapping = aes(x = year, y= cases, color = description)) +
  geom_line() +
  #plotting_theme +
  labs(title = 'Cases in Ituri, COD, vaccine scenario over time')


ggplot(data = comparison, mapping = aes(x = year, y= cases/cohort_size, color = description)) +
  geom_line() +
  #plotting_theme +
  labs(title = 'Incidence in Ituri, COD, vaccine scenario over time')



# Calculate Pf PR 2-10
cols  <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
baseline$PfPR2_10 <- baseline$n_detect_730_3649/baseline$n_730_3649

# Plot Pf PR 2-10
plot(x = baseline$timestep/365, y = baseline$PfPR2_10, type = "l",
     col = cols[7], ylim = c(0,1), lwd = 2,
     ylab = expression(paste(italic(Pf),"PR"[2-10])), xlab = "Years",
     xaxs = "i", yaxs = "i")
grid(lty = 2, col = "grey80", lwd = 0.5)
