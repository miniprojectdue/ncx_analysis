# ncx_analysis

Sex-specific dofetilide-induced repolarisation instability in a human ventricular model of type 2 diabetes (T2DM), using the [ToR-ORd / ToRORd-Land framework](https://github.com/MaxxHolmes/Sex_Specific_Human_Electromechanics/tree/main) with NCX-focused follow-up analyses.

## Overview

This repository contains the MATLAB code, scripts, and output files used for investigating whether compounded Na+/Ca2+ exchange (NCX) upregulation in sex-parameterised T2DM ventricular cell models increases sensitivity to dofetilide-induced repolarisation instability, and whether this effect is attenuated by resetting NCX to control level.

The analysis uses:
- sex- and cell-type-specific ToRORd-Land baselines from [Holmes et al.](https://www.biorxiv.org/content/10.1101/2025.03.14.643310v1.full.pdf)
- a compact T2DM remodelling layer
- dofetilide concentration-response simulations
- instability and strict EAD threshold analyses
- NCX-reset follow-up simulations

## Repository structure

- `Sex_Specific_Human_Electromechanics/`  
  Public Holmes sex-specific [ToR-ORd / ToRORd-Land framework](https://github.com/MaxxHolmes/Sex_Specific_Human_Electromechanics/tree/main) model files used as the simulation framework.

- `functions/`  
  Helper MATLAB functions for pacing, biomarker extraction, threshold detection, and analysis.

- `scripts/`  
  MATLAB scripts for running baseline simulations, dofetilide threshold searches, follow-up analyses, and figure generation.

- `results/`  
  Generated output files, summary CSVs, and exported figures.

- `setup_paths.m`  
  MATLAB path setup script for this project.

## Main study outputs

The main outputs reported in the project are
- baseline AP/CaT summaries across sex, disease state, and cell type
- strict EAD threshold summary (`C_EAD`)
- broader repolarisation-instability threshold summary (`C_instab`)
- NCX-reset follow-up results in endocardial T2DM models
- final figure showing representative action-potential traces at 1000 nM dofetilide

## How to run

1. Clone this repository.
2. Open MATLAB in the project root.
3. Run

```matlab
setup_paths
