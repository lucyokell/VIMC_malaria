# launch orderly reports (locally or on the cluster)
run_site_report<- function(sites, 
                      report_name,
                      path){
  
 orderly2::orderly_run(report_name,
                        list(
                          iso3c = sites$iso3c,
                          site_name = sites$site_name,
                          ur = sites$ur,
                          description = sites$description,
                          scenario = sites$scenario,
                          parameter_draw = sites$parameter_draw,
                          quick_run = sites$quick_run),
                          root = path)
  
 message('report complete')
  
}



run_country_report<- function(countries, 
                      report_name,
                      path){
  
  orderly2::orderly_run(report_name,
                        list(
                          iso3c = countries$iso3c,
                          description = countries$description,
                          scenario = countries$scenario,
                          parameter_draw = countries$parameter_draw,
                          quick_run = countries$quick_run),
                        root = path)
  
  message('report complete')
  
}
