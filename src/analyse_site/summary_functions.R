#summary functions

pull_input_params<- function(site_name, ur){
  
  message('parameterizing')
  # site data
  site <- extract_site(site_file = site_data, 
                       site_name = site_name,
                       ur = ur)
  
  # specify vaccine coverage based on forecast  ----------------------------------
  site<- expand_intervention_coverage(site, 
                                      terminal_year = more_params$term_yr)
  site<- update_coverage_values(site, 
                                coverage_data, 
                                scenario_name = scenario)
  
  
  # add in scenario variable which will be used to implement booster
  vaccine_plot_input<- copy(site)
  check_eir(site)
  
  # pull parameters for this site ------------------------------------------------
  params <- site::site_parameters(
    interventions = site$interventions,
    demography = site$demography,
    vectors = site$vectors,
    seasonality = site$seasonality,
    eir = site$eir$eir[1],
    burnin = more_params$burnin,
    overrides = list(human_population = more_params$pop_val)
  )
  
  # set age groups
  params$clinical_incidence_rendering_min_ages = more_params$min_ages
  params$clinical_incidence_rendering_max_ages = more_params$max_ages
  params$severe_incidence_rendering_min_ages = more_params$min_ages
  params$severe_incidence_rendering_max_ages = more_params$max_ages
  params$age_group_rendering_min_ages = more_params$min_ages
  params$age_group_rendering_max_ages = more_params$max_ages
  
  # if this is a stochastic run, set parameter draw ------------------------------
  params<- parameterize_stochastic_run(params, parameter_draw)
  
  inputs <- list(
    'param_list' = params,
    'site_name' = site_name,
    'ur' = ur,
    'iso' = iso3c,
    'scenario' = scenario,
    'description' = description,
    'parameter_draw' = parameter_draw
  )
  
  return(inputs)
  
}

run_model<- function(model_input){
  message('running the model')
  
  params <- model_input$param_list
  params$progress_bar <- TRUE
  timesteps <<- model_input$param_list$timesteps
  
  model <- malariasimulation::run_simulation(timesteps = params$timesteps,
                                             parameters = params)
  
  # add identifying information to output
  model <- model |>
    mutate(site_name = model_input$site_name,
           urban_rural = model_input$ur,
           iso = model_input$iso3c,
           description = model_input$description, 
           scenario = model_input$scenario,
           parameter_draw = model_input$parameter_draw,
           population = more_params$pop_val,
           burnin = more_params$burnin)
  
  # save model runs somewhere
  message('saving the model')
  return(model)
}

process_output<- function(model){
  
  message('postprocessing')
  
  # calculate rates
  raw_output<- drop_burnin(model, burnin= more_params$burnin* 365)
  
  output <- postie::get_rates(
    raw_output,
    time_divisor = 365,
    baseline_t = 1999,
    age_divisor = 365,
    scaler = 0.215,
    treatment_scaler = 0.517,
  )
  
  dt<- vimc_postprocess(output, le,  site_data, vimc_pop, pop_single_yr)
  
  
  # final formatting  ------------------------------------------------------------
  output<- format_outputs(dt)
  
  return(output)
}


analyse_site<- function(index, site_df){
  
  site_name<- site_df[index, site_name]
  ur<- site_df[index, ur]
  
  model_input<- pull_input_params(site_name, ur)
  model<- run_model(model_input)
  output<- process_output(model)
  
  return(output)
}