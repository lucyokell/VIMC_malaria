make_parameter_maps <- function(iso3cs,
                                population,
                                scenarios = c(
                                  'no-vaccination',
                                  'malaria-r3-default',
                                  'malaria-r3-r4-default',
                                  'malaria-rts3-bluesky',
                                  'malaria-rts3-default',
                                  'malaria-rts3-rts4-bluesky',
                                  'malaria-rts3-rts4-default'
                                ),
                                description,
                                parameter_draw,
                                quick_run,
                                burnin,
                                metadata) {
  
  # Function to read sites for a given iso3c
  get_sites <- function(iso3c) {
    sites <- data.table(readRDS(paste0('src/process_inputs/site_files/', iso3c, '.rds'))$sites)
    remove_zero_eirs(iso3c, sites)
  }
  
  # Create a data.table of iso3c, site_name, and ur
  map <- rbindlist(lapply(iso3cs, function(iso3c) {
    sites <- get_sites(iso3c)
    data.table('iso3c' = iso3c, 'site_name' = sites$name_1, 'ur' = sites$urban_rural)
  }))
  
  # Expand the grid to include input parameters
  map[, c('population', 'description', 'parameter_draw', 'quick_run', 'burnin') := 
        .(population, description, parameter_draw, quick_run, burnin)]
  
  # Create a data.table of iso3cs
  country_map <- data.table('iso3c' = iso3cs)
  
  # Expand the grid to include input parameters
  country_map[, c('population', 'description', 'parameter_draw', 'quick_run', 'burnin') := 
                .(population, description, parameter_draw, quick_run, burnin)]
  
  # Create a full map data.table
  full_map <- rbindlist(lapply(scenarios, function(scen) {
    map[, scenario := scen]
  }))
  
  # Create a full country map data.table
  full_country_map <- rbindlist(lapply(scenarios, function(scen) {
    country_map[, scenario := scen]
  }))
  
  # Convert site_name to ASCII
  Encoding(full_map$site_name) <- "UTF-8"
  full_map$site_name <- iconv(full_map$site_name, from="UTF-8", to="ASCII//TRANSLIT")
  
  return(list('site_map' = full_map, 'country_map' = full_country_map))
}
