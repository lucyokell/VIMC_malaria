##  run site-level summary plots
library(site)
library(scene)
library(data.table)
library(dplyr)

# intervention coverage --------------------------------------------------------
site_files<- list.files('src/process_inputs/site_files/new_site_files/', full.names = TRUE)

pdf('MIM/intervention_coverage.pdf', width= 12)
for(site in site_files){
  
  message(site)
  site_file<- readRDS(site)
  
  
  p<- plot_intvn_coverage(site_file)
  print(p)
}
dev.off()



# incidence over time and age --------------------------------------------------
plot_incidence<- function(iso3c){
  
  output<- readRDS(paste0('outputs/site_outputs/', iso3c, '_site_output_draw_0.rds'))
  
  sites<- data.table(unique(output[, c('site_name', 'urban_rural')]))
  
  for(index in c(1:nrow(sites))){
    
    site_name<- sites[index, site_name]
    ur<- sites[index, urban_rural]
    
    plot<- output[site_name == {{site_name}} & urban_rural == {{ur}}]
    plot<- plot |>
      group_by(year, scenario, parameter_draw) |>
      summarise(cases = sum(cases),
                cases_novax = sum(cases_novax),
                doses = doses,
                .groups = 'keep')
    
    ggplot(data = plot, mapping = aes(x= year, y= (cases_novax - cases)/doses, alpha = 0.05)) +
      geom_line(aes( color = scenario)) +
      geom_abline(slope = 0 , intercept = 0, linetype = 'dashed') +
      facet_wrap(~scenario)+
      labs(title ='Cases averted over age',
           subtitle = paste0(site_name, ' ', ur)) +
      plotting_theme
    
  }

}

# incidence over age  ----------------------------------------------------------