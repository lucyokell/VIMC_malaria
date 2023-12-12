# process inputs for all sites -------------------------------------------------
orderly2::orderly_parameters(iso3c = NULL)

pkgs <- c( 'data.table', 'countrycode',  'dplyr')
invisible(lapply(pkgs, library, character.only = TRUE))

# outputs from this report  ----------------------------------------------------
orderly2::orderly_artefact('Site file input', 'site_file.rds')
orderly2::orderly_artefact('VIMC input', 'vimc_input.rds')


# pull coverage data -----------------------------------------------------------
coverage_files<- list.files('vimc_inputs/vaccine_coverage/', full.names = T)
coverage_dt<- rbindlist(lapply(coverage_files, read.csv)) 

coverage_dt <- coverage_dt |>
  filter(country_code == iso3c)


# pull population data (single year) -------------------------------------------
demog_single_yr<- read.csv('vimc_inputs/demography/202310gavi-1_dds-202208_int_pop_both.csv') |>
  filter(country_code == iso3c) |>
  filter(year >= 2000)


# pull population data (all ages) ----------------------------------------------
demog_all_ages<- read.csv('vimc_inputs/demography/202310gavi-1_dds-202208_tot_pop_both.csv') |>
  filter(country_code == iso3c) |>
  filter(year >= 2000)


# pull life expectancy data ----------------------------------------------------
le<- read.csv('vimc_inputs/demography/202310gavi-1_dds-202208_life_ex_both.csv') |>
  filter(country_code == iso3c) |>
  filter(year >= 2000)

# pull mortality rate data  ----------------------------------------------------
mort<- read.csv('vimc_inputs/demography/202310gavi-1_dds-202208_mort_rate_both.csv') |>
  filter(country_code == iso3c) |>
  filter(year >= 2000)


# make a list of all dataframes needed for one country to reduce dependency calls
vimc_input<- list('coverage_input' = coverage_dt,
                  'population_input_single_yr' = demog_single_yr,
                  'population_input_all_age' = demog_all_ages,
                  'le_input' = le,
                  'mort_rate_input' = mort)

saveRDS(vimc_input, 'vimc_input.rds')

# pull site data  --------------------------------------------------------------
site_data <- readRDS(paste0('site_files/', iso3c, '.rds'))
saveRDS(site_data, 'site_file.rds')

