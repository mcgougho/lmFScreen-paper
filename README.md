# lmFScreen-paper

This repository contains R scripts to reproduce all figures for the “Valid F‑screening in Linear Regression” paper, using a frozen snapshot of the **lmFScreen** package.

---

##  Quick start

1. **Clone** this repo with `git clone https://github.com/mcgougho/lmFScreen-paper.git`
in your terminal. 
2. In RStudio, **open** the `lmFScreen-paper.Rproj` file.
3. **Install** the required packages by running the following command in the R console:
 - `install.packages("renv")`
 - `renv::activate()`
 - `renv::restore()`
 
It is important to run `renv::restore()` in the R console to install the required packages. This will create a virtual environment with the same package versions used in the original analysis. This ensures that the code runs correctly and produces the same results as in the paper.
4. **Run** the figure scripts in the `figures` folder. For example, to run the script for type 1 error and power, run the following command in the R console: `source("power_t1error_global_null.R")`

