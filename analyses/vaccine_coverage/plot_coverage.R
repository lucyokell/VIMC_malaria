### plot coverage inputs

library(data.table)
library(tidyr)
library(ggplot2)

r21<- read.csv('src/process_inputs/vimc_inputs/sept_inputs/coverage_202409malaria-1_malaria-r3-r4-default.csv')
rtss<- read.csv('src/process_inputs/vimc_inputs/sept_inputs/coverage_202409malaria-1_malaria-rts3-rts4-default.csv')

data<- rbind(r21, rtss)
ctys<- unique(r21$country)

plot_coverage<- function(cty){


  subset<- data |> filter(country == cty) |>
    mutate(coverage = coverage / proportion_risk)  # convert coverage to per population instead of per population at risk


  p<- ggplot(data= subset, mapping = aes(x= year, y= coverage, color = vaccine))+
    geom_line()+
    facet_wrap(~scenario) +
    labs(title = 'Vaccine coverage, Montagu input 202409malaria-1',
          subtitle = cty) +
    theme_classic()

  print(p)
}

pdf('coverage_inputs_montagu_per_pop.pdf')
lapply(ctys, plot_coverage)

dev.off()
