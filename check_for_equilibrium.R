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
filenames <-list.files(path)

get_site_output <- function(filename) {
  message(filename)
  site_output <- readRDS(paste0(path,"/",filename,"/processed_output.rds"))
  return(site_output)
  
}

check_if_output_exists<- function(filename){
  
  filepath<- paste0(path, '/', filename, '/processed_output.rds')
  
  if(file.exists(filepath)){
    message('file exists')
  }else{
    
    message('no file')
    return(filename)
  }
}

filenames_to_avoid<- lapply(filenames, check_if_output_exists)

outputs<- rbindlist(lapply(filenames[1:8000], get_site_output))
orig<- copy(outputs)

outputs<- data.table(orig)
outputs<- outputs[urban_rural== 'urban']

# plot incidence by site
pdf('incidence_runtime_plots_urban_sites2.pdf', height= 12, width = 12)

for (site in unique(outputs$site_name)){
  
  message(site)
  output<- outputs[site_name == site]
  
  
  p<- ggplot(data= output[age < 25], mapping = aes(x= year, y= clinical, color= scenario, fill= scenario))+
    geom_line()  +
    facet_wrap(~age, scales= 'free') +
  labs(x= 'Time (in years)', 
         y= 'Incidence rate', 
         title= paste0('Incidence rate over time: ', unique(output$site_name)),
         color= 'Scenario', 
         fill= 'Scenario') +
    #scale_color_manual(values= wes_palette('Royal2' )) +
    #scale_fill_manual(values= wes_palette('Royal2'))  +
    plotting_theme
  
  print(p)
  
  
}

dev.off()
