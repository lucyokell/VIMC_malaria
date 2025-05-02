## VIMC malaria modelling
This repository contains code used to estimate the impact of malaria vaccines on cases, deaths, and DALYs from 2000-2100 for the Vaccine Impact Modelling Consortium (VIMC). More information on the VIMC can be found [here](https://www.vaccineimpact.org/). This workflow is written in [orderly2](https://mrc-ide.github.io/orderly2/), a package designed at ICL to facilitate reproducible analysis. Source code is found in the [/src](https://github.com/mrc-ide/VIMC_malaria/tree/main/src) folder, with the name of each folder corresponding to the name of each report and the corresponding script containing code for execution. The following reports + actions are run chronologically to produce final estimates, with more detail below. 

The helper functions developed for this workflow can be found in the [vimcmalaria](https://github.com/mrc-ide/vimcmalaria) package. This package is a holding place for helper functions used in this workflow, and therefore may not be useful when applied to other contexts. More documentation + testing for this package is forthcoming.

## Directory
```
.
├── analyses                             # Ad-hoc analyses for project-specific work
|   ├── MIM                              # Generate plots for Pan-African Malaria Conference poster (MIM)
|   ├── decomp                           # Run a generic vaccine impact model and compare how impact metrics vary if non-vaccine malaria interventions are included (was not used)
|   ├── ethiopia                         # recalibrate Ethiopia admin-1 units to lower baseline transmission intensity (no longer used)
|   ├── gf_data_request                  # Data request from BMGF regarding subnational tailoring (see gf_data_request branch for more code)
|   ├── itns                             # refit ITN usage in consolidated site files using netz package (diagnostic)
|   ├── mvip                             # Pulled mortality metrics from previous model runs to compare to impact metrics from Malaria Vaccine Implementation Programme
|   ├── paper                            # Generate tables and figures for manuscript
|   ├── site_files                       # compare outputs run with new vs. old site files
|   ├── vaccine_scenario                 # Generate proxy vaccine coverage scenario for manuscript based on DTP3 coverage, vaccine type, and year of introduction
├── src                                  # Source code for VIMC malaria main workflow
|   ├── process_inputs                   # Format VIMC demography, vaccine coverage, and malariasimulation site file inputs for analysis
|   |   ├── process_inputs.R                   
|   ├── process_country                  # Paramaterize, run models, and do initial postprocessing
|   |   ├── process_country.R                   
|   ├── postprocessing                   # Do additional postprocessing, including scaling national outputs to World Malaria Report data and aggregation
|   |   ├── postprocessing.R                   
|   ├── diagnostics                      # Generate an Rmarkdown as a diagnostic report of vaccine impact results
|   |   ├── diagnostics.R                   
├── VIMC_malaria.Rproj                   # R.Studio project file
├── orderly_config.yml                   # Configuration for orderly workflow. Do not modify
├── pkgdepends.txt                       # Package dependencies for cluster installation
├── postie_0.1.2.txt                     # Postie package version for VIMC workflow (because DALYs are recalculated based on VIMC inputs). Might be a separate branch, should be updated in the future
├── workflow.R                           # Script to initialize orderly workflow, set up cluster environment with hipercow package, and run analysis. Start here
└── README.md                            # Project overview

```

###  Quick Start
In order to run this workflow, run the ["workflow.R"](https://github.com/mrc-ide/VIMC_malaria/blob/main/workflow.R) script. 
- `Process_inputs` formats the input files needed for this workflow.
- `Process_country` parameterizes, models, and postprocesses outputs used to estimate vaccine impact at the admin1 level.
- `Postprocessing` conducts basic postprocessing, including scaling and aggregation.
- `Diagnostics` generates a basic report with plots for vetting.

#### Save input files
VIMC model inputs are saved locally and not tracked on this repository due to large size and data privacy contraints. VIMC inputs and site files should be saved under `/src/process_inputs/vimc_inputs` and `/src/process_inputs/site_files`, respectively. Contact Lydia for access to these files.

#### Run process inputs report for each country of interest
This only needs to be run one time for each country, if VIMC inputs and site files do not change. This report pulls in VIMC inputs and site files, reformatting for later use. 

####  Change input parameters
The following parameters must be changed for each run:
- `iso3c`: country/countries you would like to run models for
- `parameter_draw`: parameter draw value for model run. For median parameter values, set to 0.
- `description`: description of the reason for a certain model run. Make sure to change this value for each run, unless you seek to overwrite pre-existing outputs.
- `quick_run`: Boolean, for whether you would like to run a test model or the full simulation. A test model produces outputs with wider age bands and a shorter time horizon (2000-2035), in order to optimize model run time. Preferable to set to `TRUE` when testing models locally, debugging changes, etc.
- `scenario`: vaccine scenario you would like to run models for. Options:
    * `no-vaccination`: No vaccines implemented. 
    * `r3-r4-default`: Full series of R21 (with booster), based on GAVI forecasts.
    * `rts3-rts4-default`: Full series of RTS,S (with booster), based on GAVI forecasts.
    * `proxy`: proxy routine coverage scenario based on public GAVI estimates of protected children, routine DTP3 coverage, and vaccine choice by country (R21 vs. RTS,S). 
      
The `make_parameter_map` function will create input parameter data frames (at the site and country level) for all of the sites in the 31 VIMC-modelled countries, as well as each VIMC vaccination scenario. 

#### Run process_country  for all countries, parameter draws, and scenarios of interest
This reports must be run in chronological order for all of the sites in a country for the vaccine scenario of interest. The code to run these reports can be found [here](https://github.com/mrc-ide/VIMC_malaria/blob/main/workflow.R). 
Note that given extended run time, you will likely prefer to launch models on the cluster. Ensure orderly2, malariasimulation, dplyr, and data.table are installed in your cluster environment before launching models or they will fail. 

#### Other notes
This workflow can be quite space-intensive (particularly when running models for all sites in the 31 VIMC input countries). It is worthwhile to monitor the size of this repository and regularly clean out non-final reports. 

---

## Methods documentation
###  Format inputs
Process external inputs for country of interest, including site file, demographic, and vaccine coverage data [("process_inputs")](https://github.com/mrc-ide/VIMC_malaria/blob/main/src/process_inputs/process_inputs.R). Reads in and saves the following inputs by country:
- Site file (characterizing the level and pattern of transmission in an admin 1 unit, in addition to intervention coverage and vector species)
- Life expectacy by age (from VIMC)
- Mortality rate (from VIMC)
- Population size (from VIMC)

These reports should only be run once, then rerun if any inputs change.

### Parameterize model
Parameterize model. Models are run with single year age groups from 0 to 20, with 10-year age groups from 20 through 100. Models were run from 2000-2100, with a 15-year burn-in period.

The `process_country` script pulls in the corresponding site file for an admin 1 unit, which characterizes the pattern on malaria transmission in this area. The interventions component of the site file is modified based on the VIMC scenario of interest, specifying the coverage and booster coverage for each year of the simulation (for either RTS,S or R21). These values are added to the site file via the following columns: **rtss_cov**, **r21_cov**. We additionally added the columns **vaccine** to specify whether the vaccine scenario is RTS,S or R21. 

R21 vaccine profile parameters are obtained from ["The Public Health Impact and Cost-Effectiveness of the R21/Matrix-M Malaria Vaccine: A Mathematical Modelling Study"]([https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4597985](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(23)00816-2/fulltext)), in which a semi-mechanistic vaccine efficacy model was fit to phase 2b R21 trial data on children observed over 12-18 months of follow-up in Nanoro, Burkina Faso.

RTS,S vaccine profile parameters are obtained from [Public health impact and cost-effectiveness of the RTS,S/AS01 malaria vaccine: a systematic comparison of predictions from four mathematical models](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(23)00816-2/fulltext). 

We additionally carry over intervention coverage from the last observed year (typically 2023) out to 2100, assuming constant values for the remainder of the simulation period. Note that insecticide-treated net (ITN) usage follows a 3-year cyclical pattern based on administrated and time-based waning of net efficacy-- the pattern of the last 3 year cycle observed is carried out for the remainder of the simulation period, to capture this temporal trend. 

We utilize the [site package](https://mrc-ide.github.io/site/) to translate site file inputs into malariasimulation parameters. The site package contains data on parasite prevalence, demography, coverage of non-vaccine malaria interventions, seasonality of transmission, and the relative abundance of different mosquito species for each admin-1 unit in Sub-Saharan Africa. For access to these site files, please contact Pete Winskill. 

### Launch model
Run malariasimulation model. If you would like to obtain an estimate of model run time, test this report locally before launching on the cluster.

### Postprocess outputs 
Outputs are processed in the following steps.

* `get_rates()` from the [postie](https://github.com/mrc-ide/postie) package is used to estimate incidence, mortality, YLD, YLL, and DALY rates from model outputs. Note that these rates are based on the population size used to run the model, not the real-world population size of the site modelled. DWs for malaria are sourced from the Global Burden of Disease study. For the purposes of VIMC modelling, we ignore malaria GBD disability weights for comorbid conditions such as anemia and motor impairiment. For documentation on the DW values used, see source code [here](https://github.com/mrc-ide/postie/blob/dalys/R/epi.R#L36-L85).
* Because the VIMC utilizes country-specific life-expectancy, YLLs and DALY rates are recalculated based on these inputs.
* Rates are multiplied by site population to estimate cases, deaths, and DALYs from 2000-2100.
* We scale site file populations such that the sum of site file populations is equivalent to the national VIMC population.

Aggregate outputs up to country level [("process_country")](https://github.com/mrc-ide/VIMC_malaria/blob/main/src/process_country/orderly.R) via simple summation.

### Visualize results
Produce diagnostic report. Note that for this report to run properly, you must have processed outputs (from "process_country") for the no-vaccination scenario in addition to the intervention scenario you specify. You cannot run a diagnostic report for the no-vaccination scenario alone.
