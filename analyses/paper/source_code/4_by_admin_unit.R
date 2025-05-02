
library(vimcmalaria)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(ggpubr)
library(data.table)


  cols<- c('#CE5137', '#1B4F72', '#45B39D', '#BA4A00')

  plotting_theme<- theme_bw(base_size = 14) +
    theme( legend.position = 'bottom',
           strip.text.x = element_text(size = rel(0.8)),
           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
           text = element_text(family= 'TT Arial'),
           axis.ticks.y= element_blank(),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank())

completed <- completed_reports("process_country") |>
  dplyr::filter(
    description == "new_site_files",
    scenario== 'no-vaccination',
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

  message(directory)

  output<- readRDS(paste0(output_filepath, directory, '/outputs.rds'))                  # get output file
  sites<- data.table::rbindlist(lapply(output$site_output, function(x) return(x$processed_output))) #pull out processed site_level output


  saveRDS(sites, paste0('J:/VIMC/VIMC_malaria/analyses/paper/inputs/site_output/', iso3c, '_site_output_novax.rds'))

  return(sites)
}

lapply(c(1:nrow(completed)), get_site_outputs, map= completed, output_filepath= 'J:/VIMC/VIMC_malaria/archive/process_country/'))


# admin level site results
#site_output<- lapply(c(1:nrow(completed)), get_site_outputs, map = completed, output_filepath = 'archive/process_country/')
files<- list.files('analyses/paper/inputs/site_output/', full.names = TRUE)
outputs<- rbindlist(lapply(files, readRDS))
outputs<- outputs |>
  mutate(id = paste0(site_name, '_', urban_rural))


#or try country level results and see
# incidence rate before vaccine is introduced versus 15 years after
inc_22<- outputs |>
  group_by(scenario, year, id) |>
  filter(year == 2022) |>
  mutate(cases_2022 = cases,
         cohort_22 = cohort_size)

inc_37<- outputs |>
  group_by(scenario, year, id) |>
  filter(year == 2037) |>
  rename(cases_2037 = cases,
         cohort_37 = cohort_size) |>
  ungroup() |>
  select(cases_2037, cohort_37, scenario, id)

# inc_52<- vax |>
#   group_by(scenario, year, id) |>
#   filter(year == 2043) |>
#   rename(cases_52 = cases,
#          cohort_52 = cohort_size) |>
#   ungroup() |>
#   select(cases_52, cohort_52, scenario, id)

# inc_67<- vax |>
#   group_by(scenario, year, id) |>
#   filter(year == 2067) |>
#   rename(cases_67 = cases,
#          cohort_67 = cohort_size) |>
#   ungroup() |>
#   select(cases_67, cohort_67, scenario, id)

inc<- merge(inc_22, inc_37, by = c('id', 'scenario'))
# inc<- merge(inc, inc_52, by = c('id', 'scenario'))
#inc<- merge(inc, inc_67, by = c('id', 'scenario'))

inc<- inc |>
  group_by(scenario, id) |>
  summarise(cases_2022 = sum(cases_2022),
            cases_2037 = sum(cases_2037),
            cohort_22 = sum(cohort_22),
            cohort_37 = sum(cohort_37),
            #cases_52 = sum(cases_52),
            #cohort_52 = sum(cohort_52),
            # cases_67 = sum(cases_67),
            # cohort_67 = sum(cohort_67),
            .groups = 'keep')

#inc<- inc |> filter(scenario != 'Control')
inc<- data.table(inc)
inc<- inc[ scenario == 'malaria-r3-r4-default', scenario:= 'Routine']
inc<- inc[scenario == 'no-vaccination', scenario:= 'Control']

ggplot(data = inc, mapping = aes(x= cases_2022/cohort_22, y= cases_2037/cohort_37, color= scenario)) +
  geom_point()+
  #facet_wrap(~scenario) +
  geom_abline(slope = 1, intercept = 0) +
  plotting_theme +
  xlim(0, 1.2) +
  ylim(0, 1.2)+
  theme(legend.position = 'bottom') +
  scale_color_manual(values = cols)+
    labs(color= 'Scenario',
         x= 'Incidence rate before vaccine introduction',
         y= 'Incidence rate 15 years after vaccine introduction')


vax<- readRDS('analyses/draft/vax.rds')


#or try country level results and see
# incidence rate before vaccine is introduced versus 15 years after
inc_22<- vax |>
  group_by(scenario, year, country) |>
  filter(year == 2022) |>
  mutate(cases_2022 = cases,
         cohort_22 = cohort_size)

inc_37<- vax |>
  group_by(scenario, year, country) |>
  filter(year == 2037) |>
  rename(cases_2037 = cases,
         cohort_37 = cohort_size) |>
  ungroup() |>
  select(cases_2037, cohort_37, scenario, country)



inc<- merge(inc_22, inc_37, by = c('country', 'scenario'))

inc<- inc |>
  group_by(scenario, country) |>
  summarise(cases_2022 = sum(cases_2022),
            cases_2037 = sum(cases_2037),
            cohort_22 = sum(cohort_22),
            cohort_37 = sum(cohort_37),

            .groups = 'keep')

#inc<- inc |> filter(scenario != 'Control')
ggplot(data = inc, mapping = aes(x= cases_2022/cohort_22, y= cases_2037/cohort_37, color= scenario)) +
  geom_point()+
  #facet_wrap(~scenario) +
  geom_abline(slope = 1, intercept = 0) +
  plotting_theme +
  theme(legend.position = 'none') +
  scale_color_manual(values = cols)+
    labs(
         x= 'Incidence rate before vaccine introduction',
         y= 'Incidence rate 15 years after vaccine introduction')

