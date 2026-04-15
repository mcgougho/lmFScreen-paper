# lmFScreen-paper

This repository contains R scripts to reproduce all figures for the paper ["Valid F-screening in linear regression"](https://arxiv.org/abs/2505.23113), using a frozen `renv` environment and version 0.2.0 of the **lmFScreen** package.

---

##  Quick start

1. **Clone** this repo with `git clone https://github.com/mcgougho/lmFScreen-paper.git`
in your terminal. 
2. In RStudio, **open** the `lmFScreen-paper.Rproj` file.
3. In the R console, install and restore the project environment:
   - `install.packages("renv")`
   - `renv::activate()`
   - `renv::restore()`

This will recreate the project library using the package versions recorded in `renv.lock`, ensuring that the code runs with the same software environment used for the analysis.
 
It is important to run `renv::restore()` in the R console to install the required packages. This will create a virtual environment with the same package versions used in the original analysis. This ensures that the code runs correctly and produces the same results as in the paper.

4. **Run** the figure scripts in the `figures` folder.

---

## Figures

The `figures/` directory contains R scripts to reproduce all plots from the paper.

### 1. Power and Type I Error (Global and Local Null)
- **Files:** `power_t1error_global_null.R`, `power_t1error_local_null.R`
- **Description:**  
  Demonstrates Type I error control of the selective p-value under the null and compares the power of the selective procedure with a sample-splitting approach. `power_t1error_global_null.R` corresponds to Figure 4 in the paper ["Valid F-screening in linear regression"](https://arxiv.org/abs/2505.23113), `power_t1error_local_null.R` corresponds to Figure C.1 in the appendix.

---

### 2. Confidence Interval Coverage and Width
- **Files:** `CI_globalnull.R`, `CI_localnull.R`
- **Description:**  
  Compares the selective confidence intervals to naive intervals in terms of:
  - Coverage (selective vs. non-selective)
  - Width (one plot with the value of the true coefficient on the horizontal axis and one plot with the number of observations as the horizontal axis)
`CI_globalnull.R` corresponds to Figure 5 in ["Valid F-screening in linear regression"](https://arxiv.org/abs/2505.23113), `CI_localnull.R` corresponds to Figure C.2 in the appendix.
---

### 3. Multiple Testing Adjustments
- **File:** `multiple_corrections.R`
- **Description:**  
  Shows that traditional multiple testing corrections (Bonferroni, Scheffé) do **not** control Type I error after F-screening. `multiple_corrections.R` corresponds to Figure A.1 in the appendix of ["Valid F-screening in linear regression"](https://arxiv.org/abs/2505.23113).

---

### 4. Comparison to [Heller et al.](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12318)
- **File:** `heller_comparison.R`
- **Description:**  
  Compares the performance of the method proposed in the paper ["Valid F-screening in linear regression"](https://arxiv.org/abs/2505.23113) to the method of ["Post-selection estimation and testing following aggregate association tests"](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12318) in the setting of F-screening. `heller_comparison.R` corresponds to Figure B.1 in the appendix of ["Valid F-screening in linear regression"](https://arxiv.org/abs/2505.23113).
  
---

### 5. Selective vs. Naive Point Estimates
- **File:** `point_estimates.R`
- **Description:**  
  Plots the selective vs. non-selective point estimates for the first coefficient, under the global null. `point_estimates.R` corresponds to Figure C.1 in the the appendix of ["Valid F-screening in linear regression"](https://arxiv.org/abs/2505.23113).

---

### 6. Assumption violation sensitivity
- **File:** `sensitivity.R`
- **Description:**  
  Shows the sensitivity of the selective p-value to violations of the assumptions (homoskedasticity, normality). `sensitivity.R` corresponds to Figures D.3, D.4, and D.5 in the appendix of ["Valid F-screening in linear regression"](https://arxiv.org/abs/2505.23113).
  
---

