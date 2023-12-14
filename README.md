# VIMC Malaria modelling
This repository contains code used to estimate the impact of malaria vaccines on cases, deaths, and DALYs from 2000-2100 for the Vaccine Impact Modelling Consortium (VIMC). More information on the VIMC can be found [here](https://www.vaccineimpact.org/). This workflow is written in [orderly2](https://mrc-ide.github.io/orderly2/), a package designed to facilitate reproducible analysis. Source code is found in the [/src](https://github.com/mrc-ide/VIMC_malaria/tree/main/src) folder, with the name of each folder corresponding to the name of each report. The following reports + actions are run chronologically to produce final estimates, with more detail below. 

##  Quick Start
In order to run this workflow, run the ["VIMC_workflow.R"](https://github.com/mrc-ide/VIMC_malaria/blob/main/VIMC_workflow.R) script. Carry out the following steps:

### Initialize orderly repository
Install orderly2, if you have not already, with the command `install_github('mrc-ide\orderly2')`.
You may need to initialize this repository by setting your working directory to the repository and running the command `orderly2::orderly_init()`. This will install supplemental folders in the repository that are not tracked via Git but are essential to run orderly reports.

###  Install packages 
This workflow uses a modified version of the site package [(link here)](https://github.com/mrc-ide/site_vimc), which is used to translate site files into malariasimulation model paramters. The key change is that this version allows the user to parameterize the R21 vaccine, and also implements flat booster coverage of 90% for the routine scenarios and 100% for the blue sky scenarios. The modified version of this package should be installed prior to running this workflow using the following command:
`install_github('mrc-ide/site_vimc')`. 

This workflow also uses the postie package to postprocess outputs, with DALY functionality for postie is currently on the [DALYs](https://github.com/mrc-ide/postie/tree/dalys/R) repository branch; you will need to install this version of the postie package to run this code using the function call `install_github('mrc-ide/postie@dalys')`. 

### Save input files
VIMC model inputs are saved locally and not tracked on this repository due to large size and privacy issues. VIMC inputs and site files should be saved under `/src/process_inputs/vimc_inputs` and `/src/process_inputs/site_files`, respectively. Contact Lydia for access to these files.

### Run process inputs report
This only needs to be run one time, if VIMC inputs and site files do not change.

###  Change input parameters
The following parameters must be changed for each run:
* **iso3c:** country/countries you would like to run models for
- **draw:** draw value for model run. For median parameter values, set to 0.
- **description:** description of the reason for a certain model run. Make sure to change this value for each run, unless you seek to overwrite pre-existing outputs.
- **quick_run**: Boolean, for whether you would like to run a test model or the full simulation. A test model produces outputs with wider age bands and a shorter time horizon (2000-2035), in order to optimize model run time. Preferable to set to TRUE when testing models locally, debugging changes, etc.

The `make_parameter_maps` function will create input parameter data frames (at the site and country level) for all of the sites in the 31 VIMC-modelled countries, as well as each VIMC vaccination scenario. If you would like to change these parameters specifically, or run models on only a subset of sites or scenarios, they are as follows:
- **site_name:** name of site to run
- **ur:** urban/rural split for site of interest.
* **scenario:** scenario you would like to run models for. Options:
    * 'malaria-no-vaccination': No vaccines implemented
    * 'malaria-r3-default': Initial series of R21, based on GAVI forecasts.
    * 'malaria-r3-r4-default': Full series of R21 (with booster), based on GAVI forecasts.
    * 'malaria-rts3-bluesky': 90% coverage of initial series of RTS,S for entire modelling period
    * 'malaria-rts3-default': Initial series of RTS,S, based on GAVI forecasts.
    * 'malaria-rts3-rts4-bluesky': 90% of full series of RTS,S (90% coverage for entire modelling time period)
 
### Run analyse_site for all sites
This report must be run for all of the sites in a country for the vaccine scenario of interest.
Note that given long run time, you will likely prefer to launch models on the cluster. Ensure packages are installed in your cluster environment before launching models or they will fail.

### Run process_country and country_diagnostics for all countries
Run these reports using the code linked [here](https://github.com/mrc-ide/VIMC_malaria/blob/main/VIMC_workflow.R#L93-L104). 

### Other notes
This workflow can be quite space-intensive (particularly when running models for all sites in the 31 VIMC input countries). It is worthwhile to monitor the size of this repository and regularly clean out non-final reports. 

# Methods documentation
##  Process inputs
Process external inputs for country of interest, including site file, demographic, and vaccine coverage data [("process_inputs")](https://github.com/mrc-ide/VIMC_malaria/blob/main/src/process_inputs/orderly.R). Reads in and saves the following inputs by country:
- Site file (characterizing the level and pattern of transmission in an admin 1 unit, in addition to intervention coverage and vector species)
- Life expectacy by age (from VIMC)
- Mortality rate (from VIMC)
- Population size (from VIMC)

These reports should only be run once, then rerun if inputs change.

## Analyse site
### Parameterize
Models are run with single year age groups from 0 to 100, from 2000-2100. 
This analysis pulls in the corresponding site file for an admin 1 unit, which characterizes the pattern on malaria transmission in this area. The interventions component of the site file is modified based on the VIMC scenario of interest, specifying the coverage and booster coverage for each year of the simulation (for either RTS,S or R21). These values are added to the site file via the following columns: **rtss_coverage**, **rtss_booster_coverage**, **r21_coverage**, and **r21_booster_coverage**. We additionally added the columns **vaccine** to specify whether the vaccine scenario is RTS,S or R21; and **scenario_type** to specify whether a routine or optimal (blue-sky) scenario is being implemented. 

R21 vaccine profile parameters are obtained from the pre-print ["The Public Health Impact and Cost-Effectiveness of the R21/Matrix-M Malaria Vaccine: A Mathematical Modelling Study"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4597985), in which a semi-mechanistic vaccine efficacy model was fit to phase 2b R21 trial data on children observed over 12-18 months of follow-up in Nanoro, Burkina Faso. The median vaccine efficacy parameters used are found in the table below.

![image](https://github.com/mrc-ide/VIMC_malaria/assets/55333260/f5935495-0bf0-48f1-a68c-d5962c2fae7b)

We additionally carry over intervention coverage from the last observed year (typically 2023) out to 2100, assuming constant values for the remainder of the simulation period. Note that insecticide-treated net (ITN) usage follows a 3-year cyclical pattern based on administrated and time-based waning of net efficacy-- the pattern of the last 3 year cycle observed is carried out for the remainder of the simulation period, to capture this temporal trend. 

We utilize a modified version of the site package to translate site file inputs into malariasimulation parameters. The package is found [here](https://github.com/mrc-ide/site_vimc). The `add_interventions()` function in the site file was modified to impose varying booster coverage depending on scenario and pull in median R21 vaccine profile parameters in addition to RTS,S. 

### Postprocessing
Outputs are processed in the following steps.

* `get_rates()` from the [postie](https://github.com/mrc-ide/postie) package is used to estimate incidence, mortality, YLD, YLL, and DALY rates from model outputs. Note that these rates are based on the population size used to run the model, not the real-world population size of the site modelled. DWs for malaria are sourced from the Global Burden of Disease study. For the purposes of VIMC modelling, we ignore malaria GBD disability weights for comorbid conditions such as anemia and motor impairiment. For documentation on the DW values used, see source code [here](https://github.com/mrc-ide/postie/blob/dalys/R/epi.R#L36-L85).
* Because the VIMC utilizes country-specific life-expectancy, YLLs and DALY rates are recalculated based on these inputs.
* Rates are multiplied by site population to estimate cases, deaths, and DALYs from 2000-2050.
* We scale site file populations such that the sum of site file populations is equivalent to the national VIMC population.
* From 2050 onward, site population is estimated as (national population in year) * [(scaled site population in 2050)/ (national population in 2050)]. We assume that the proportional breakdown of population by site is fixed from 2050-2100, in the absence of other data.

## Produce site diagnostics
Produce diagnostic report (at the site level) [("site_diagnostics")](https://github.com/mrc-ide/VIMC_malaria/blob/main/src/process_site/orderly.R). Note that for this report to run properly, you must have processed outputs (from "analyse_site") for the no-vaccination scenario in addition to the intervention scenario you specify. You cannot run a diagnostic report for the no-vaccination scenario alone.

## Process country
Aggregate outputs up to country level [("process_country")](https://github.com/mrc-ide/VIMC_malaria/blob/main/src/process_country/orderly.R) via simple summation. We additionally scale cases based on the average number of cases in the World Malaria report for the country of interest for the years 2018-2020.

## Produce country diagnostics
Produce similar diagnostic report for country outputs [("country_diagnostics")](https://github.com/mrc-ide/VIMC_malaria/blob/main/src/country_diagnostics/orderly.R). Note that for this report to run properly, you must have processed outputs (from "process_country") for the no-vaccination scenario in addition to the intervention scenario you specify. You cannot run a diagnostic report for the no-vaccination scenario alone.
