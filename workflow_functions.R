# workflow functions  ----------------------------------------------------------

# mapping ----------------------------------------------------------------------
make_parameter_maps<- function(iso3cs,
                               population,
                               scenarios =  c('no-vaccination', 'malaria-r3-default', 'malaria-r3-r4-default', 'malaria-rts3-bluesky', 'malaria-rts3-default', 'malaria-rts3-rts4-bluesky', 'malaria-rts3-rts4-default'),
                               description,
                               parameter_draw,
                               quick_run,
                               metadata){
  
  # could be cleaner but will do for now
  
  map<- data.table()
  for (iso3c in iso3cs){
    
    sites<- data.table(readRDS(paste0('src/process_inputs/site_files/', iso3c, '.rds'))$sites)  # sites for country of interest
    sites<- remove_zero_eirs(iso3c, sites)
    
    site_num<- data.table('iso3c' = iso3c, 'site_name' = sites$name_1, 'ur' = sites$urban_rural)  
    
    map<- rbind(site_num, map)
  }
  
  # expand grid out to include input parameters-
  map<- map |>
    mutate(description = description,
           parameter_draw = parameter_draw,
           quick_run = quick_run)
  
  
  
  full_map<- data.table()
  for (scen in scenarios){
    subset<- map|>
      mutate(scenario = scen)
    
    full_map<- rbind(subset, full_map)
  }
  
  
  country_map<- data.table('iso3c' = iso3cs)
  
  # expand grid out to include input parameters-
  country_map<- country_map |>
    mutate(description = description,
           parameter_draw = parameter_draw,
           quick_run = quick_run)
  
  
  full_country_map<- data.table()
  for (scen in scenarios){
    subset<- country_map|>
      mutate(scenario = scen)
    
    full_country_map<- rbind(subset, full_country_map)
  }
  
  Encoding(full_map$site_name) <- "UTF-8"
  full_map$site_name<- iconv(full_map$site_name, from="UTF-8", to="ASCII//TRANSLIT")
  
  
  return(list('site_map' = full_map, 'country_map'= full_country_map))
}
# metadata functions  ----------------------------------------------------------
remove_duplicate_reports<- function(report_name, parameter_map, day= NULL){
  # check if you have run this report before; if so remove from list of parameters to run
  # note you may want to rerun a report with the same parameters if you have changed source code; in that case do not use this function
  
  meta <- orderly2::orderly_metadata_extract(name = report_name, extract = c('time', 'parameters'), options = orderly2::orderly_search_options(allow_remote = TRUE))
  
  meta<- meta|>
    tidyr::separate(col = id, into = c('date', 'other'), sep = '-')|>
    mutate(date= as.numeric(date))
  
  if(day){
    
    meta<- meta |>
      filter(date >= day)
  }
  
  unique(lapply(meta$parameters, names))
  
  nms <- names(meta$parameters[[1]])
  pars <- do.call("data.frame", setNames(lapply(nms, function(nm) sapply(meta$parameters, function(x) x[[nm]])), nms))
  
  # observations that are in the parameter map, but have not already been run (in report metadata)
  # these are the observations you will actually want to run
  diff<- setdiff(parameter_map, pars)
  
  
  return(diff)
}



generate_parameter_map_for_next_report<- function(report_name, parameter_map, day= NULL){
  
  # before you run a report, check that the preceding reports have already been run
  
  meta <- orderly2::orderly_metadata_extract(name = report_name, extract = c('time', 'parameters'),  options = orderly2::orderly_search_options(allow_remote = TRUE))
  
  meta<- meta|>
    tidyr::separate(col = id, into = c('date', 'other'), sep = '-')|>
    mutate(date= as.numeric(date))
  
  if(day){
    
    meta<- meta |>
      filter(date >= day)
  }
  
  unique(lapply(meta$parameters, names))
  nms <- names(meta$parameters[[1]])
  pars <- do.call("data.frame", setNames(lapply(nms, function(nm) sapply(meta$parameters, function(x) x[[nm]])), nms))
  
  # pars contains completed reports
  # these are the observations you will actually want to run
  can_run<- intersect(parameter_map, pars)
  
  return(can_run)
  
}



# check that the sites you seek to launch have non-zero EIRs
remove_zero_eirs<- function(iso3c, sites){
  eirs<- data.table(readRDS(paste0('src/process_inputs/site_files/', iso3c, '.rds'))$eir)  # sites for country of interest
  eirs<- eirs[spp == 'pf' & eir == 0]
  remove<- eirs[, c('name_1', 'urban_rural')]
  
  if(nrow(remove) > 0){
    for (i in 1:nrow(remove)){
      
      message(paste0('removing site ', i))
      
      sites<- sites[!(name_1== remove[i, name_1] & urban_rural == remove[i, urban_rural])]
      
    }
    
  } else{
    message('No zero eir sites to remove')
  }
  return(sites)
}


# launch reports  --------------------------------------------------------------
run_site_report<- function(sites, 
                           report_name,
                           path){
  
  orderly2::orderly_run(report_name,
                        list(
                          iso3c = sites$iso3c,
                          site_name = sites$site_name,
                          ur = sites$ur,
                          description = sites$description,
                          scenario = sites$scenario,
                          parameter_draw = sites$parameter_draw,
                          quick_run = sites$quick_run),
                        root = path)
  
  message('report complete')
  
}



run_country_report<- function(countries, 
                              report_name,
                              path){
  
  orderly2::orderly_run(report_name,
                        list(
                          iso3c = countries$iso3c,
                          description = countries$description,
                          scenario = countries$scenario,
                          parameter_draw = countries$parameter_draw,
                          quick_run = countries$quick_run),
                        root = path)
  
  message('report complete')
  
}
