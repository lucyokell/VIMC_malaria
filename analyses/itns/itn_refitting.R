##  update ITN usage in site files


library(netz)

test_site<- readRDS('src/process_inputs/site_files/new_site_files/AGO_new_eir.rds')


plot_itn_use<- function(site_name, test_site){

  message(site_name)

  intvns<- test_site$interventions |> filter(name_1 == {{site_name}})
  intvns<- intvns |>
    mutate(t = (year - 1999) * 365)

  year <- 365
  timesteps <- 23 * year
  population <-5000

  distribution <- intvns$itn_input_dist
  distribution_timesteps <- year * (intvns$year - 1999)
  net_hl <- 5 * year



  pu <- population_usage(
    timesteps = timesteps,
    distribution = distribution,
    distribution_timesteps = distribution_timesteps,
    half_life = net_hl)


  pd <- data.frame(t = 1:timesteps, usage = pu)

  p<- ggplot() +
    geom_line(pd, mapping = aes(x = t, y = usage), color = 'blue') +
    ylim(0, 1) +
    theme_bw() +
    geom_line(intvns, mapping = aes(x= t, y = itn_input_dist)) +
    geom_line(intvns, mapping = aes(x= t, y = itn_use), color = 'green') +
    ggtitle(unique(intvns$name_1))+
    labs(subtitle = paste0(unique(intvns$country), ', green is original ITN use, blue is refit ITN use, black is ITN input distribution'))

  print(p)
}


plot_itn_use('Moxico', test_site)
site_names<- unique(test_site$sites$name_1)


