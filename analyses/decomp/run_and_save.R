run_and_save<- function(model_input, description){

  
  library(data.table)
  library(vimcmalaria)
  library(ggplot2)
  library(dplyr)
  library(site)
  
  message('modelling')
  
model<- vimcmalaria::run_model(model_input) 
  message('saving')
  
  saveRDS(model, paste0('J:/VIMC_malaria/analyses/decomp/', description, '.rds'))
}