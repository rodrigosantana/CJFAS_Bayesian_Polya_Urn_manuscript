# Counting the Uncountable: a Bayesian Polya Urn model for monitoring artisanal fisheries

Repository containing **data, scripts, and figures** that support the manuscript:

The project focuses on evaluating and illustrating the use of a Bayesian **Polya Urn** model for **finite population sampling inference** in artisanal fisheries monitoring.

## Overview

This repository brings together:

- a virtual reference population (synthetic census) of fishers;
- R scripts for descriptive analysis and simulation-based calibration;
- intermediate/final `.Rda` objects used to reproduce tables and outputs;
- `.png` figures used for inspection and communication of results.

In practice, the workflow implements a two-stage sampling design (clusters and units within clusters), generates Polya posterior draws, and computes performance metrics such as accuracy, interval coverage, and variability (CV/RCV).

## Repository structure

```text
.
├── README.md
├── script_data_analysis _CJFAS.R
├── script_calibracao_CJFAS.R
├── dkft_CensoVirtual_perCapita.Rda
├── PopVir_list.Rda
├── dd.Rda
├── Rplot_PopPairs.png
├── Rplot_cpue_G.png
└── Rplot_PolyaPost.png
```

### Main files

- `script_data_analysis _CJFAS.R`  
  Main exploratory-analysis script and preparation of outputs for the manuscript. It loads the virtual census, summarizes effort/catch by week and cluster, organizes tables, and runs a sampling example for Polya-posterior-based estimation.

- `script_calibracao_CJFAS.R`  
  Simulation-based calibration script (repeated loops) to evaluate estimator performance in finite populations under two-stage sampling, including statistics such as mean, median, quantiles, CV/RCV, interval width, and coverage.

- `dkft_CensoVirtual_perCapita.Rda`  
  Core data object containing the virtual population/census used by the scripts.

- `PopVir_list.Rda` and `dd.Rda`  
  Auxiliary R objects that support analyses and consistency checks.

- `Rplot_PopPairs.png`, `Rplot_cpue_G.png`, `Rplot_PolyaPost.png`  
  Figures generated in the analysis context (relationships among variables, CPUE behavior, and posterior-related outputs).

## Methodological workflow (summary)

1. **Define the virtual finite population** with cluster/community structure and weekly observations.
2. **Run two-stage sampling**:
   - stage 1: select clusters with probability linked to the number of possible within-cluster sampling combinations;
   - stage 2: sample fishers within each selected cluster.
3. **Apply Bayesian Polya Urn inference** to expand sampled information and obtain posterior distributions for monthly population totals (effort and catch categories).
4. **Assess estimator performance** via repeated simulation (calibration), focusing on:
   - accuracy of central estimates;
   - interval width and coverage (percentile and/or HDI intervals);
   - stability through classical and robust variability measures.

## Requirements

- **R** (recommended: R >= 4.2)
- R packages used in the scripts:
  - `polyapost`
  - `HDInterval`

> Note: scripts assume `.Rda` files are available in the current working directory.

## Reproducibility

From the project directory, run in R:

```r
# 1) Descriptive analysis and table/object preparation
source("script_data_analysis _CJFAS.R")

# 2) Simulation-based calibration
source("script_calibracao_CJFAS.R")
```

### Execution notes

- The calibration script can be computationally intensive because it uses many simulation loops (`nloop` and `sim`).
- For quick exploratory runs, temporarily reduce these parameters.
- Keep `set.seed()` when strict reproducibility of random draws is required.

## Scope and status

This repository stores computational support material for the manuscript (analysis and calibration). It is not structured as an R package; instead, it prioritizes transparency and traceability of the study workflow via scripts.

## Authors

- **Paul Gerhard Kinas**
- **Rodrigo Sant'Ana**

## Citation

If this material supports academic or technical work, please cite the associated manuscript:

*Counting the Uncountable: a Bayesian Polya Urn model for monitoring artisanal fisheries*.

## License

 
                  Creative Commons License 4.0
                       (CC BY-NC-SA 4.0)
 
  This is a humam-readable summary of (and not a substitute for) the
  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 
  You are free to:
 
  Share - copy and redistribute the material in any medium or format.
 
  The licensor cannot revoke these freedoms as long as you follow the
  license terms.
 
  Under the following terms:
 
  Attribution - You must give appropriate credit, provide a link to
  license, and indicate if changes were made. You may do so in any
  reasonable manner, but not in any way that suggests the licensor
  endorses you or your use.
 
  NonCommercial - You may not use the material for commercial
  purposes.
 
  ShareAlike - If you remix, transform, or build upon the material,
  you must distributive your contributions under the same license
  as the  original.
 
  No additional restrictions — You may not apply legal terms or
  technological measures that legally restrict others from doing
  anything the license permits.
 
