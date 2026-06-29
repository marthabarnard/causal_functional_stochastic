# Code and data for Causal Inference for Functional Treatments with Stochastic Policies
By Martha Barnard, Jared D. Huling, and Julian Wolfson (https://arxiv.org/abs/2606.27518)

## Code
 R version 4.4.0 and 4.5.0 was used for all analysis

### Replicate results
simulations_figures.qmd, nhanes_toy_figures.qmd, nhanes_real_data_analysis.qmd replicate/recreate paper figures using saved results

### Re-run analysis
The folders simulations/ and nhanes_analysis/ contain the code to re-run all analyses from the paper. We recommend the use of a computing cluster for all code, however a computer cluster is **necessary** to run the simulation code. In these folders:
* simulations/
  * Files starting with `input_` have functions for generating the data
  * Files starting with `mid_` have functions for running our proposed method
  * Files starting with `run_` are the main code files for running the simulations
  * Files starting with `init_` are sample bash scripts for running the `run_` files on a Slurm managed cluster
* nhanes_analysis/
  * Files starting with `process_` go through the data cleaning process
  * `real_data_helper_funcs.R` contains functions to run our proposed method on the NHANES data
  * Files starting with `run_` are the main code files that run our proposed method on the NHANES data (takes ~1-3 hours for each file on a single machine)


## Data
All data files are linked in the `nhanes_analysis/process_` files however the core data files are:
* Acceleromter data: https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/PAXMIN_G.xpt
* Acceleromter data: https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/PAXMIN_H.xpt
* Step count data: https://physionet.org/files/minute-level-step-count-nhanes/1.0.1/csv/nhanes_1440_oaksteps.csv.xz
* Mortality data: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/
* All other covariate data can be found at: https://wwwn.cdc.gov/nchs/nhanes/
