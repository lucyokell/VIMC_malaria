#############
# COUNTRY LEVEL
#############

path <-"N:/Lydia/VIMC_malaria/archive/process_country"
filenames <-list.files(path)
library(purrr)
library(dplyr)
library(ggplot2)

date = as.numeric(substring(filenames,1,8))
filenames <- filenames[which(date>20231130)]

list.files(paste0(path,"/",filenames[1]))

# some functions
get_country <- function(filename) {
  site_data <-readRDS(paste0(path,"/",filename,"/site_file.rds"))
  return(data.frame(iso = site_data$sites$iso3c[1], country=site_data$sites$country[1]))
}


get_scenario <- function(filename) {
  country_output <- readRDS(paste0(path,"/",filename,"/country_output.rds"))
  return(data.frame(scenario = country_output$scenario[1], 
                    description = country_output$description[1]))
}

#### compile country data
cres <-data.frame(filenames = filenames, date = date[which(date>20231130)],
                  date_time = as.numeric(paste0(substring(filenames,1,8), substring(filenames, 10,13)))
)
cres <- cbind(cres, bind_rows(map(filenames, get_country)))
cres <- cbind(cres, bind_rows(map(filenames, get_scenario)))

dim(cres)
cres <- cres |>
  arrange(desc(date_time)) |>
  dplyr::distinct(iso, country, scenario, description, .keep_all = TRUE) |>
  arrange(country, scenario)


############################
#### Country impact_metrics
get_tot_impacts <- function(filename) {
  country_output <- readRDS(paste0(path,"/",filename,"/country_output.rds"))
  return(data.frame(cases=sum(country_output$cases),
                    deaths=sum(country_output$deaths),
                    dalys= sum(country_output$dalys),
                    cases_model = sum(country_output$clinical*country_output$cohort_size)))
  
}

cres <- cbind(cres, bind_rows(map(cres$filenames, get_tot_impacts)))
cres <- cres |> 
  group_by(iso) |>
  mutate(cases_averted = cases[scenario=='no-vaccination'] - cases,
         prop_cases_averted =  1- cases / cases[scenario=='no-vaccination'],
         deaths_averted = deaths[scenario=='no-vaccination'] - deaths,
         prop_deaths_averted = 1 - deaths / deaths[scenario=='no-vaccination'])

cres <-cres |>
  arrange(prop_cases_averted) 

## get clinical and mortality outputs

## arrange in expected order of impact and re order factor scenario
cres$scenario <- factor(cres$scenario, levels=c("malaria-rts3-default",
                                                "malaria-rts3-bluesky",
                                                "malaria-rts3-rts4-default",
                                                "malaria-rts3-rts4-bluesky",
                                                "malaria-r3-default"      ,
                                                "malaria-r3-r4-default" ))

ggplot(cres |> filter(scenario!='no-vaccination'), aes(y=scenario, x=prop_cases_averted)) +
  geom_bar(stat = "identity") +
  facet_wrap(~iso) + #, scales = "free") +
  theme_bw()


dt<- cres |> filter(scenario!= 'no-vaccination')
check_bluesky_routine_order<- function(dt){

  full_rerun<- data.table()
  for (cty in (unique(dt$iso))){
    subset<- dt |>
      filter(iso == cty) |>
      data.table()
    
    if(subset[scenario== 'malaria-rts3-rts4-bluesky', prop_cases_averted]
       < subset[scenario== 'malaria-rts3-rts4-default', prop_cases_averted]){
      
      message(paste0('malaria-rts3-rts4-bluesky is less than malaria-rts3-rts4-default for ', cty))
      
      rerun<- data.table('scenarios' = c('malaria-rts3-rts4-bluesky', 'malaria-rts3-rts4-default'), 'iso3c' = cty)
      full_rerun<- rbind(rerun, full_rerun)
    }
    
    if(subset[scenario== 'malaria-rts3-bluesky', prop_cases_averted]
       < subset[scenario== 'malaria-rts3-default', prop_cases_averted]){
      
      message(paste0('malaria-rts3-bluesky is less than malaria-rts3-default for ', cty))
      
      rerun<- data.table('scenarios' = c('malaria-rts3-bluesky', 'malaria-rts3-default'), 'iso3c' = cty)
      full_rerun<- rbind(rerun, full_rerun)
      
    }
  }
  
  return(full_rerun)
}
check_three_four_dose_order<- function(dt){
  
  
  full_rerun<- data.table()
  for (cty in (unique(dt$iso))){
    
    subset<- dt |>
      filter(iso == cty) |>
      data.table()
    
    if(subset[scenario== 'malaria-rts3-rts4-bluesky', prop_cases_averted]
       < subset[scenario== 'malaria-rts3-bluesky', prop_cases_averted]){
      
      message(paste0('malaria-rts3-rts4-bluesky is less than malaria-rts3-bluesky for ', cty))
      
      rerun<- data.table('scenarios' = c('malaria-rts3-rts4-bluesky', 'malaria-rts3-bluesky'), 'iso3c' = cty)
      full_rerun<- rbind(rerun, full_rerun)
    }
    
    if(subset[scenario== 'malaria-r3-r4-default', prop_cases_averted]
       < subset[scenario== 'malaria-r3-default', prop_cases_averted]){
      
      message(paste0('malaria-r3-r4-default is less than malaria-r3-default for ', cty))
      
      rerun<- data.table('scenarios' = c('malaria-r3-r4-default', 'malaria-r3-default'), 'iso3c' = cty)
      full_rerun<- rbind(rerun, full_rerun)
    }
    
    if(subset[scenario== 'malaria-rts3-rts4-default', prop_cases_averted]
       < subset[scenario== 'malaria-rts3-default', prop_cases_averted]){
      
      message(paste0('malaria-rts3-rts4-default is less than malaria-rts3-default for ', cty))
      
      rerun<- data.table('scenarios' = c('malaria-rts3-rts4-default', 'malaria-rts3-default'), 'iso3c' = cty)
      full_rerun<- rbind(rerun, full_rerun)
    }
  }
  return(full_rerun)
}
rerun<- check_bluesky_routine_order(dt)
rerun2<- check_three_four_dose_order(dt)
reruns<- rbind(rerun, rerun2)
reruns<- unique(reruns)
