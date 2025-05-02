# TITLE    MVIP mortality plots
# PURPOSE  pull out mortality reduction estimates by pfpr / admin 1 unit to compare to MVIP all-cause mortality reduction
# AUTHOR   Lydia Haile
##########################################################################################################################

setwd('J:/september_runs/VIMC_malaria')

library(dplyr)
library(data.table)
library(ggplot2)
library(vimcmalaria)
cols<- c('#CE5137', '#1B4F72', '#45B39D', '#BA4A00')

completed <- completed_reports("process_country") |>
  dplyr::filter(
    description == "booster_update",
    scenario %in% c('malaria-r3-r4-bluesky', 'no-vaccination'),
    quick_run == FALSE,
    parameter_draw == 0
  ) |>
  dplyr::arrange(dplyr::desc(date_time)) |>
  dplyr::distinct(iso3c, scenario, quick_run, parameter_draw, description, .keep_all = TRUE) |>
  dplyr::arrange(iso3c, scenario, parameter_draw)


#' Pull site level processed output based on metadata input
#' @param index           observation in metadata df
#' @param map             metadata df
#' @param output_filepath filepath where outputs live
#' @export
get_site_outputs<- function(index, map, output_filepath){

  metadata<- map[ index,]

  directory<- metadata$directory_name
  iso3c<- metadata$iso3c
  scenario<- metadata$scenario

  message(directory)

  output<- readRDS(paste0(output_filepath, directory, '/outputs.rds'))                  # get output file
  sites<- data.table::rbindlist(lapply(output$site_output, function(x) return(x$processed_output))) #pull out processed site_level output


  saveRDS(sites, paste0('J:/september_runs/VIMC_malaria/analyses/draft/site_output/', iso3c, '_', scenario, '_site_output.rds'))

  return(sites)
}
extract_site_data <- function(directory) {
  site_file <- readRDS(paste0("archive/process_inputs/", directory, "/site_file.rds"))

  # extract prevalece data
  prev <- site_file$prevalence |>
    filter(year %in% c(2000, 2019)) |>
    mutate(id = paste0(name_1, "_", urban_rural, '_', iso3c)) |>
    data.table()

  prevs <- dcast(prev, id ~ year, value.var = "pfpr") |>
    rename(
      pfpr_2000 = `2000`,
      pfpr_2019 = `2019`
    )


  # extract intervention coverage data for last available year
  intvns <- site_file$interventions |>
    filter(year == 2022) |>
    mutate(id = paste0(name_1, "_", urban_rural, '_', iso3c)) |>
    select(id, itn_use, tx_cov, smc_cov, pmc_cov, irs_cov)

  site_info <- merge(prevs, intvns, by = "id") |>
    select(pfpr_2000, pfpr_2019, id)

  return(site_info)
}

# pull VIMC mortality results
# admin level site results
#site_output<- lapply(c(1:nrow(completed)), get_site_outputs, map = completed, output_filepath = 'archive/process_country/')
files<- list.files('analyses/draft/site_output/', full.names = T)
files<- files[files %like% 'no-vaccination']
output<- rbindlist(lapply(files, readRDS))
output<- output |>
  mutate(id = paste0(site_name, '_', urban_rural, '_', country))

# Extract site file prevalence data for 2000 and 2019 and merge onto site level data-set
completed<- completed_reports('process_inputs') |>
  dplyr::arrange(dplyr::desc(date_time)) |>
  dplyr::distinct(iso3c, .keep_all = TRUE) 


prevs<- rbindlist(lapply(completed$directory_name, extract_site_data))
output_prev<- merge(output, prevs, by = 'id', allow.cartesian = TRUE)

pull_mortality_metrics<- function(id_name){
  message(id_name)
# summarize under 5 mortality 
test<- output_prev |> 
  filter(id == id_name,
         age < 5,
         year < 2051) 

pfpr_2000<- unique(test$pfpr_2000)
pfpr_2019<- unique(test$pfpr_2019)
iso3c<- unique(test$country)
  
test <- test |>
  group_by(scenario, year, id, pfpr_2019) |>
  summarise(cohort_size= sum(cohort_size),
            deaths= sum(deaths),
            .groups = 'keep') |>
  mutate(mortality= (deaths/cohort_size) *1000) |>
  filter(year == 2019) |>   #pull mortality from 2030 (at this point the vaccine has been introduced in all locations)
  as.data.table() 
  
  return(test)
}


mort_outputs<- rbindlist(lapply(unique(output_prev$id), pull_mortality_metrics))


ggplot(mort_outputs, mapping = aes(x= pfpr_2019, y= mortality))+
geom_point()+
scale_color_manual(values= cols) +
  theme_minimal() +
labs(title = 'U5 mortality rate per 1000 vs. PFPR in 2019',
    y= 'U5 mortality rate per 1000',
  x= 'PFPR in 2019 (MAP)')


# should run same plots:
# summarized across parameter draws
# for bluesky scenario so scaleup is the same across all settings (coverage is already > 85% in MVIP countries in 2030)

#look at cause fraction data from GBD in 2019
# Ghana: malaria cause fraction is 16% in U5s
# Malawi: CF is 12%
# Kenya introduced the vaccine in the lake-endemic region
# in this region CF is between 11-16%

# if we assume a 33% reduction in malaria specific mortality, and mortality accounts for
# between 12- 16% of all-cause mortality in these countries, 
# what reduction in all-cause mortality do we estimate?

# all cause mortality = sum(cause specific mortality * cause fraction)
# reduction in all-cause mortality due to malaria = malaria-specific mortality * malaria CF
# .33 * .16 = 0.0528, which is about 50% of the 13% reduction in all-cause mortality observed in Malaria Vaccine Implementation Program
