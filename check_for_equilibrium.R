# check when outputs equilibriate by scenario
# may help reduce run time for scenarios
plotting_theme<- theme_bw(base_size = 12) +
  theme( legend.position = 'bottom',
         strip.text.x = element_text(size = rel(0.8)),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1), 
         text = element_text(family= 'Arial Narrow'),
         axis.ticks.y= element_blank(), 
         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

library(ggplot2)
library(wesanderson)


path <-"N:/Lydia/VIMC_malaria/archive/process_site"



# some functions
get_country <- function(filename) {
  test <- try(readRDS(paste0(path,"/",filename,"/processed_output.rds")), silent=TRUE)
  
  if (class(test) %in% 'try-error') {
    next
  } else {
    site_output <- readRDS(paste0(path, "/", filename, "/processed_output.rds"))
  }
  
  return(data.frame(scenario = site_output$scenario[1], 
                    description = site_output$description[1],
                    country = site_output$country[1],
                    site_name = site_output$site_name[1],
                    urban_rural = site_output$urban_rural[1]))
}




#### compile country data
cres <-data.frame(filenames = filenames, date = date[which(date>20231130)],
                  date_time = as.numeric(paste0(substring(filenames,1,8), substring(filenames, 10,13)))
)
cres <- cbind(cres, bind_rows(map(filenames, get_country)))

dim(cres)
cres <- cres |>
  arrange(desc(date_time)) |>
  dplyr::distinct(country, site_name, urban_rural, scenario, description, .keep_all = TRUE) |>
  arrange(country, scenario)



get_site_output <- function(filename) {
  site_output <- readRDS(paste0(path,"/",filename,"/processed_output.rds"))
  return(site_output)
  
}

outputs<- rbindlist(lapply(filenames[1:500], get_site_output))

outputs<- data.table(outputs)
outputs<- outputs[urban_rural== 'urban']

# plot incidence by site
pdf('incidence_runtime_plots2.pdf')

for (site in unique(outputs$site_name)){
  
  message(site)
  output<- outputs[site_name == site]
  
  
  p<- ggplot(data= output[age < 25], mapping = aes(x= year, y= clinical, color= scenario, fill= scenario))+
    geom_line()  +
    facet_wrap(~age, scales= 'free') +
  labs(x= 'Time (in years)', 
         y= 'Incidence rate', 
         title= paste0('Incidence rate over time: ', unique(output$site_name), ', urban'),
         color= 'Scenario', 
         fill= 'Scenario') +
    scale_color_manual(values= wes_palette('Royal2' )) +
    scale_fill_manual(values= wes_palette('Royal2'))  +
    plotting_theme
  
  print(p)
  
  
}

dev.off()
