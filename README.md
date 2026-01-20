## Mathematica and R code for mechano-electro-osmotic (MEO) spheroid model: Figure 5B and Figure S5B

This repository contains code used to reproduce results associated with the accompanying manuscript.

### What this code reproduces

**Figure 5B (Mathematica)**  
The Mathematica notebook reproduces the results shown in Figure 5B, including:
- Cell area (radial profile)
- Cell aspect ratio (radial profile)

For completeness, the notebook also computes and plots the solid stress field versus radial distance, although this quantity is not displayed in Figure 5B.

**Figure S5B (R)**  
The R script reproduces Figure S5B by computing and plotting the tension-responsive transcriptional program as a function of spheroid radial position (Smart-seq3 dataset), including:
- Single-cell values (background points)
- Binned means (10 radial bins) with error bars
- A fitted trend curve from a negative binomial model

### Contents

**Mathematica**
- `MEO_spheroid_model_Fig5B.nb`  
  Main Mathematica notebook (parameters, equations, numerical solution, and plotting for Figure 5B and the solid stress field).

**R**
- `figS5B_tension_responsive.R`  
  R script for Figure S5B (tension-responsive program radial trend).

### Software requirements

**Mathematica**
- Wolfram Mathematica 13.1 or later (the notebook was prepared in Mathematica 13.1).

**R**
- R (>= 4.2 recommended)
- Required R packages: `dplyr`, `tibble`, `readr`, `ggplot2`, `MASS`, `XICOR`, `scales`

### Input data for the R script (Figure S5B)

The R script requires the following two files:
- `allCounts.rds`
- `FACS_bc.csv`

**Important:** Place `allCounts.rds` and `FACS_bc.csv` in the **same folder** as `figS5B_tension_responsive.R` (unless you modify the file paths inside the script).

These data files can be obtained from:
Cougnoux, Antony, et al. "Diffusion Smart-seq3 of breast cancer spheroids to explore spatial tumor biology and test evolutionary principles of tumor heterogeneity." *Scientific Reports* 15.1 (2025): 3811.

### How to run

#### Figure 5B (Mathematica)
1. Download or clone this repository.
2. Open `MEO_spheroid_model_Fig5B.nb` in Mathematica.
3. Evaluate the notebook from top to bottom.
4. The notebook displays the generated plots in the notebook interface.

#### Figure S5B (R)

** in RStudio **
1. Make sure `figS5B_tension_responsive.R`, `allCounts.rds`, and `FACS_bc.csv` are in the same folder.
2. Open `figS5B_tension_responsive.R` in RStudio.
3. Run the script top to bottom.
4. The figure and summary files will be saved to the output folder specified in the script (commonly `results_geneset/`).
