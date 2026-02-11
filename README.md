# HPV MSD Power-Law Fit Analysis

This repository contains the MATLAB workflow used to compute ensemble mean squared displacement (MSD) curves from TrackMate spot trajectories and to quantify segmented power-law diffusion behavior during early Human Papillomavirus (HPV) trafficking.

The analysis reproduces the MSD plots and fitted diffusion exponents (α) reported in the associated manuscript.

---

## Requirements

- MATLAB (tested on R2023 or newer)
- msdanalyzer package installed and added to the MATLAB path  
  https://tinevez.github.io/msdanalyzer/

---

## Repository structure

scripts/  
 run_msd_analysis.m  Main analysis script  

data/trackmate_csv/  
 TrackMate trajectory input files (.csv)  

figures/  
 Final exported figures  

---

## Data placement

The analysis script reads files from the current working directory.

Before running, either:

• place the TrackMate input files in the same directory as the script, **or**  
• set MATLAB’s working directory to the folder containing the input files.

By default, the provided repository stores inputs in:

data/trackmate_csv/

If using this structure, change into that folder before running the script.

---

## How to run

1. Open MATLAB  
2. Navigate to the directory containing the input files and script  
3. Run:

run_msd_analysis

The script will compute the ensemble MSD, perform power-law fitting, and export results.

---

## Outputs generated

Running the script produces:

- msd_analysis_powerlaw_fits.png  
  Ensemble MSD curve with power-law fits

- msd_analysis_powerlaw_fit_params.csv  
  Estimated diffusion exponents (α) for each fitted segment

---

## Customization

User-adjustable parameters are defined at the top of the script and may be modified to control:

- input filenames  
- time rounding and frame interval  
- maximum lag limits  
- minimum contributing tracks per lag   
- figure appearance (axes limits, ticks, styles, formatting)

These options allow adaptation of the analysis without modifying the core computation.

---

## Citation

For full experimental design, biological context, and methodological details, please refer to:

Love M, Dang RC, Xie J, Zhang P (2026)  
The Rab5 effector Rabankyrin-5 mediates endosomal fusion and trafficking of Human Papillomavirus during early entry.  
