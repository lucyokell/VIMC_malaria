# by country
library(vimcmalaria)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(ggpubr)
library(data.table)
library(purrr)
library(sf)

cols<- c('#CE5137', '#1B4F72', '#45B39D', '#BA4A00')

saveRDS(vax, 'analyses/paper/inputs/vax.rds')
vax<- readRDS('analyses/paper/inputs/vax.rds')

table_1 <- vax |>
  group_by(scenario, country, country_name, run_id) |>
  summarise(
    cases_averted = sum(cases_averted),
    deaths_averted = sum(deaths_averted),
    dalys_averted = sum(dalys_averted),
    .groups = "keep"
  ) |>
  data.table() |>
  ungroup() |>
  group_by(scenario, country) |>
  mutate(
    cases_mean = mean(cases_averted, na.rm = TRUE),
    cases_lower = quantile(cases_averted, 0.025, na.rm = TRUE),
    cases_upper = quantile(cases_averted, 0.975, na.rm = TRUE),
    deaths_mean = mean(deaths_averted, na.rm = TRUE),
    deaths_lower = quantile(deaths_averted, 0.025, na.rm = TRUE),
    deaths_upper = quantile(deaths_averted, 0.975, na.rm = TRUE)
  )

table_1 <- arrange(table_1, -cases_averted)
table_1 <- table_1 |>
  filter(!scenario == "Control")

# plug text
get_values(table_1 |> filter(scenario == 'Routine', country_name == 'Nigeria') |> pull(cases_averted))
get_values(table_1 |> filter(scenario == 'Routine', country_name %like% 'Kinshasa') |> pull(cases_averted))
get_values(table_1 |> filter(scenario == 'Routine', country_name == 'Nigeria') |> pull(deaths_averted))
get_values(table_1 |> filter(scenario == 'Routine', country_name %like% 'Kinshasa') |> pull(deaths_averted))

p1 <- ggplot(data = table_1, aes(x = reorder(country, -cases_mean), group = scenario)) +
  geom_col(position = "dodge", width = 0.8, mapping = aes(y = cases_mean, color = scenario, fill = scenario)) +
  geom_errorbar(aes(ymin = cases_lower, ymax = cases_upper), position = position_dodge(0.8)) +
  plotting_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Cases averted",
    x = "Country",
    y = "Cases averted"
  )


p2 <- ggplot(data = table_1, aes(x = reorder(country, -deaths_mean), group = scenario)) +
  geom_col(position = "dodge", width = 0.8, mapping = aes(y = deaths_mean, color = scenario, fill = scenario)) +
  geom_errorbar(aes(ymin = deaths_lower, ymax = deaths_upper), position = position_dodge(0.8)) +
  plotting_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Deaths averted",
    x = "Country",
    y = "Deaths averted"
  )


ggarrange(p1, p2, ncol = 2)


#  same table, but per FVP over the simulation period
# pull counts of fvps and total doses delivered by country and scenario
total_doses<- compile_dose_outputs('new_site_files', iso = 'all') |> data.table() |> filter(scenario %in% c('proxy', 'no-vaccination'))
total_doses[scenario %like% 'no-vaccination', scenario := 'Control']



total_doses<- total_doses |>
  filter(year <= 2037)|>
  group_by(iso3c, scenario, run_id) |>
  summarise(cases_averted = sum(cases_averted),
            deaths_averted = sum(deaths_averted),
            fvp = sum(fvp),
            .groups = 'keep') |>
  mutate(cases_averted_per_fvp= cases_averted*100000/ fvp,
         deaths_averted_per_fvp = deaths_averted*100000/ fvp) |>
  select(iso3c, scenario, run_id, cases_averted_per_fvp, deaths_averted_per_fvp) |>
  rename(country = iso3c) |>
  ungroup() |>
  group_by(scenario, country) |>
  summarise(
    cases_mean = mean(cases_averted_per_fvp, na.rm = TRUE),
    cases_lower = quantile(cases_averted_per_fvp, 0.025, na.rm = TRUE),
    cases_upper = quantile(cases_averted_per_fvp, 0.975, na.rm = TRUE),
    deaths_mean = mean(deaths_averted_per_fvp, na.rm = TRUE),
    deaths_lower = quantile(deaths_averted_per_fvp, 0.025, na.rm = TRUE),
    deaths_upper = quantile(deaths_averted_per_fvp, 0.975, na.rm = TRUE),
    .groups= 'keep'
  ) |>
  unique()


# plot cases averted by country
p3<- ggplot(data = total_doses, aes(x = reorder(country, -cases_mean), group = scenario)) +
  geom_col(position = "dodge", width = 0.8, mapping = aes(y = cases_mean, color = scenario, fill = scenario)) +
  geom_errorbar(aes(ymin = cases_lower, ymax = cases_upper), position = position_dodge(0.8)) +
  plotting_theme +
  scale_fill_manual(values = cols)+
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = 'none') +
  labs(title = 'Cases averted per 100,000 fully vaccinated people',
       x= 'Country',
       y= 'Cases averted per 100,000 fully vaccinated people')


p4<- ggplot(data = total_doses, aes(x = reorder(country, -deaths_mean), group = scenario)) +
  geom_col(position = "dodge", width = 0.8, mapping = aes(y = deaths_mean, color = scenario, fill = scenario)) +
  geom_errorbar(aes(ymin = deaths_lower, ymax = deaths_upper), position = position_dodge(0.8)) +
  plotting_theme +
  scale_fill_manual(values = cols)+
    theme(legend.position = 'none') +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Deaths averted per 100,000 fully vaccinated people',
       x= 'Country',
       y= 'Deaths averted per 100,000 fully vaccinated people')

ggarrange(p1, p2, p3, p4, ncol =2, nrow= 2)


# results section results per dose (average and IQR)
proxt<- get_values(total_doses |> filter(scenario == 'proxy') |> pull(cases_mean))



# generate maps
Africa <- c("DZA","AGO","BEN","BWA","BFA","BDI","CMR","CPV","CAF","TCD","COM",
            "COD","DJI","EGY","GNQ","ERI","ETH","GAB","GMB","GHA","GIN","GNB",
            "CIV","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS","MYT",
            "MAR","MOZ","NAM","NER","NGA","COG","REU","RWA","SHN","STP","SEN",
            "SYC","SLE","SOM","ZAF","SSD","SDN","SWZ","TZA","TGO","TUN","UGA",
            "ESH","ZMB","ZWE")

# save copies of gadm SpatVector as .rds files in data folder
import_gadm <- function(ISO, level){
  geodata::gadm(country = ISO, level = level, path = "J:/september_runs/VIMC_malaria/analyses/gadm", version = "4.0")
}

map2(Africa, 0, import_gadm)    # admin0
map2(Africa, 1, import_gadm)    # admin1

# create a list of all countries
countries_africa <- lapply(c(Africa),
                           function(x){
                             list.files(path = "J:/september_runs/VIMC_malaria/analyses/gadm",
                                        pattern = paste0("*", x , "_0_pk.rds"),
                                        full.names = TRUE)
                           }) |>
  unlist()


admin1s_africa <- lapply(c(Africa),
                         function(x){
                           list.files(path = "J:/september_runs/VIMC_malaria/analyses/gadm",
                                      pattern = paste0("*", x , "_1_pk.rds"),
                                      full.names = TRUE)
                         }) |>
  unlist()
# unpack gadms
unpack_gadm <- function(file){

  object <- readRDS(file) # read in object
  object <- terra::vect(object) # unpack SpatVector
  st_as_sf(object) # transform to sf object

}

countries_africa <- map_dfr(countries_africa, unpack_gadm) # loop over each country
admin1s_africa <- map_dfr(admin1s_africa, unpack_gadm) # loop over each admin1
st_crs(countries_africa); st_crs(countries_africa) # view CRS


# select out African countries with malaria (i.e. countries in foresite package)
countries <- countries_africa[!(countries_africa$ID_0 %in% c("DZA", "CPV", "EGY", "LSO", "LBY", "MUS", "MYT", "MAR", "REU", "SHN", "STP", "SYC", "TUN", "ESH")), ]
admin1s <- admin1s_africa[!(admin1s_africa$ID_0 %in% c("DZA", "CPV", "EGY", "LSO", "LBY", "MUS", "MYT", "MAR", "REU", "SHN", "STP", "SYC", "TUN", "ESH")), ]

mapdat <- total_doses |>
  left_join(countries_africa |> dplyr::select(ID_0, geometry),
            by = c("country" = "ID_0")) |>
  st_as_sf()


# make admin-1 boolean variable
read_pop<- function(file_name){

  message(file_name)
  site_data<- readRDS(file_name)

  site_data$prevalence<- site_data$prevalence |>
    dplyr::filter(year == 2019) |>
    mutate(run_model = ifelse(pfpr > 0.10, TRUE, FALSE))

  # make exceptions for Madagascar, Ethiopia, and Sudan
  # hardcode for time's sake but operationalize later
  if(unique(site_data$country) == 'MDG'){

    site_data$prevalence <- site_data$prevalence |>
      mutate(run_model = ifelse(name_1 == 'Toliary', TRUE, run_model))
  }

  if(unique(site_data$country) == 'ETH'){

    site_data$prevalence <- site_data$prevalence |>
      mutate(run_model = ifelse(name_1 %like% 'Gambela', TRUE, run_model))


  }

  if(unique(site_data$country) == 'SDN'){

    site_data$prevalence <- site_data$prevalence |>
      mutate(run_model = ifelse(name_1 == 'South Darfur', TRUE, run_model)) |>
      mutate(run_model = ifelse(name_1 == 'West Kurdufan', TRUE, run_model))
  }

  prevalence<- site_data$prevalence |>
    select(name_1, urban_rural, iso3c, run_model)

  #now pull annual population and merge this on
  pop<- site_data$population


  pop<- merge(pop, prevalence, by = c('name_1', 'urban_rural', 'iso3c'))

  return(pop)
}

files<- list.files('src/process_inputs/site_files/', full.names = T)
pops<- rbindlist(lapply(files, read_pop))
pops<- pops |>
  unique(by= c('name_1', 'urban_rural'))|>
  data.table()
pops<- pops[, run_model := as.character(run_model)]
pops<- pops[run_model == TRUE, run_model := 'Introduced']
pops<- pops[run_model == FALSE, run_model := 'Not introduced']

pop_reduced<- remove_urbanicity(pops) |>
  mutate(id = paste0(name_1, '_', iso3c))

# output this as a dataset to add to the appendix
pop_format<- pop_reduced |>
  select(name_1, country, run_model) |>
  rename(admin_unit = name_1,
         vaccine_introduced= run_model)
write.csv(pop_format, 'admin_units_modelled.csv')

admin1s<- admin1s |>
  mutate(id= paste0(NAME_1, '_', ID_0))
#' site file is split by admin1 and urbanicity, for the sake of mapping need to reduce down to admin1
#' @param admin_data
#'
#' @returns admin1 data reduced down to urbanicity
#' @export
#'
#' @examples
remove_urbanicity<- function(admin_data){

  reduced<- data.frame()
  for(nm in unique(admin_data$name_1)){

    message(nm)
    subset<- admin_data |> filter(name_1 == nm)

    if(nrow(subset) > 1 & length(unique(subset$run_model))== 1){

      subset<- subset[1]
    }
    if(nrow(subset) > 1 & length(unique(subset$run_model))== 2){

      subset<- subset[run_model == 'Introduced']
    }
    subset<- subset |> select(-urban_rural)
    reduced<- rbind(subset, reduced)
  }

  return(reduced)
}


mapdat_admin<- pop_reduced |>
  left_join(admin1s |> dplyr::select(id, geometry),
            by = 'id')|>
  st_as_sf()

# Country borders --------------------------------------------------------------
# define borders between country polygons
borders <- st_boundary(countries) # will find all outside and inside borders
border_admins <- admin1s[unique(unlist(st_touches(borders, admin1s))),]

borders_inner <- countries |> # will find inner borders
  rmapshaper::ms_innerlines() |>
  # as_tibble() |>
  st_as_sf()

admin1s_intersect <- st_intersects(admin1s, borders_inner)

index <- tibble(index = as.numeric(admin1s_intersect %>% lengths > 0))
border_admins <- bind_cols(admin1s, index) |>
  mutate(nrow = row_number()) |>
  # remove admins not touching the inner borders
  filter(index == 1) |>
  # remove admin1s missed by the function (Madagascar, South Africa, Nigeria)
  filter(!(nrow %in% c(482, 490, 305, 307, 424)))


#  make a plot of all of the admin-1 units modelled for study (methods) -----
p <- ggplot() +
  geom_sf(data = countries_africa) +
  geom_sf(data = countries, fill = "cornsilk2", color = "cornsilk3") +
  #geom_sf(data = border_admins, fill = "tomato") +
  geom_sf(data = borders_inner, color = "black", size = 1) +
  geom_sf(data= mapdat_admin, aes(fill = run_model), show.legend = TRUE, linewidth = 0.005) +
  scale_fill_manual(values = cols) +
  labs(fill = "Vaccine introduced") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())



p1 <- ggplot() +
#geom_sf(data = countries, fill = "cornsilk2", color = "cornsilk3") +
geom_sf(data= mapdat, aes(fill = cases_mean), color = "lightgrey", show.legend = TRUE, linewidth = 0.005) +
scale_fill_distiller(palette = "Oranges", direction = 1, labels= scales::comma) + # "RdYlBu"
#scale_fill_continuous(labels = scales::comma) +
# scale_fill_gradient(low = "#ffcccc", high = "#990000") +
labs(fill = "Cases averted per 100,000 FVP") +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.ticks = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_blank(),
panel.border = element_blank(),
plot.background = element_blank(),
panel.background = element_blank())


p2 <- ggplot() +
#geom_sf(data = countries, fill = "cornsilk2", color = "cornsilk3") +
geom_sf(data= mapdat, aes(fill = deaths_mean), color = "lightgrey", show.legend = TRUE, linewidth = 0.005) +
scale_fill_distiller(palette = "Oranges", direction = 1, labels= scales::comma) + # "RdYlBu"
#scale_fill_continuous(labels = scales::comma) +
# scale_fill_gradient(low = "#ffcccc", high = "#990000") +
labs(fill = "Deaths averted per 100,000 FVP") +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.ticks = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_blank(),
panel.border = element_blank(),
plot.background = element_blank(),
panel.background = element_blank())



table_summary <- table_1 |>
  group_by(scenario, country, country_name) |>
  summarise(
    cases_averted = mean(cases_averted, na.rm = TRUE),
    deaths_averted = mean(deaths_averted, na.rm= TRUE),
    dalys_averted = mean(dalys_averted, na.rm = TRUE),
    .groups = "keep"
  )
mapdat <- table_summary |>
  left_join(countries_africa |> dplyr::select(ID_0, geometry),
            by = c("country" = "ID_0")) |>
  st_as_sf()


ggpubr::ggarrange(p1, p2, nrow= 2)

p3 <- ggplot() +
#geom_sf(data = countries, fill = "cornsilk2", color = "cornsilk3") +
geom_sf(data= mapdat, aes(fill = cases_averted), color = "lightgrey", show.legend = TRUE, linewidth = 0.005) +
scale_fill_distiller(palette = "PuOr", direction = -1, labels= scales::comma) + # "RdYlBu"
#scale_fill_continuous(labels = scales::comma) +
# scale_fill_gradient(low = "#ffcccc", high = "#990000") +
labs(fill = "Cases averted") +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.ticks = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_blank(),
panel.border = element_blank(),
plot.background = element_blank(),
panel.background = element_blank())


p4 <- ggplot() +
#geom_sf(data = countries, fill = "cornsilk2", color = "cornsilk3") +
geom_sf(data= mapdat, aes(fill = deaths_averted), color = "lightgrey", show.legend = TRUE, linewidth = 0.005) +
scale_fill_distiller(palette = "PuOr", direction = -1, labels= scales::comma) + # "RdYlBu"
#scale_fill_continuous(labels = scales::comma) +
# scale_fill_gradient(low = "#ffcccc", high = "#990000") +
labs(fill = "Deaths averted") +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.ticks = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_blank(),
panel.border = element_blank(),
plot.background = element_blank(),
panel.background = element_blank())



