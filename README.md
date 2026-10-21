# Replication package

**Digital Transformation and Poverty Reduction in Low-Human-Capital Economies: Panel Evidence with Pre-Existing Network Instruments**

Muhammad Younas Khan Niazi — Department of Business and Commerce, GIFT University, Gujranwala, Pakistan

This repository reproduces every table and figure in the manuscript submitted to *Information Economics and Policy*.

## Quick start

1. Install R (≥ 4.2). The packages `fixest`, `dplyr`, `ggplot2`, `jsonlite` and `scales` are installed automatically if missing; `splines` ships with R.
2. Open R with the working directory set to this folder.
3. Run:

```r
source("run_all.R")
```

Tables are written to `output/tables/` and figures to `output/figures/`. A full run takes roughly 45–60 minutes on a laptop, mostly in the 500-replication bootstraps.

## Repository structure

| Path | Content |
|---|---|
| `run_all.R` | Master script: runs `code/01`–`05` in order |
| `code/functions.R` | Estimators: control-function spline, country-cluster bootstrap, first-stage F, Hansen J, income horse race |
| `code/01_build_data.R` | Builds the instruments and `data/analysis_sample.csv` |
| `code/02_main_results.R` | Tables 1–3, Figures 1–2 (main design) |
| `code/03_instrument_designs.R` | Table 4, Figure 4, Hansen J tests, thesis (year-FE) specification |
| `code/04_robustness.R` | Table 5, Figure 3, leave-one-country-out jackknife, DTI × income horse race |
| `code/05_placebo.R` | Appendix Figure A1 (noise-instrument placebo) |
| `data/raw/` | Input data (see below) |

## Output map

| Manuscript item | File |
|---|---|
| Table 1 | `output/tables/table1_samples.csv` |
| Table 2 | `output/tables/table2_first_stage.csv` |
| Table 3 | `output/tables/table3_poverty_me_by_enrolment.csv` |
| Table 4 | `output/tables/table4_instrument_designs_all_outcomes.csv` |
| Table 5 | `output/tables/table5_robustness_poverty.csv` |
| Income horse race (Section 5.2) | `output/tables/horse_race_income.csv` |
| Figure 1 | `output/figures/fig1_first_stage.png` |
| Figure 2 | `output/figures/fig2_poverty_me_curve.png` |
| Figure 3 | `output/figures/fig3_robustness_poverty.png` |
| Figure 4 | `output/figures/fig4_instrument_designs.png` |
| Figure A1 | `output/figures/figA1_placebo_instruments.png` |

## Method in brief

- **Stage 1:** the Digital Transformation Index (DTI) is regressed on the instruments, institutional quality, trade, FDI and urbanisation, with country and year fixed effects. Instruments: pre-1995 fixed-telephone density × leave-one-out global mean DTI (`z_tel`), and lagged ICT goods imports (`ict_imports_l1`).
- **Stage 2:** each outcome is regressed on the DTI, a natural cubic spline in tertiary enrolment (knots at the 5th/50th/95th percentiles: 10.1/57.8/95.4%), their interaction, controls, and the Stage-1 residual (control function), with country and year fixed effects.
- **Marginal effect:** ME(e) = β + γ′B(e), reported at 11% enrolment (Pakistan) and over 5–100%. Effects are given per DTI unit and per standard deviation of the DTI (σ = 1.47 in the poverty sample).
- **Inference:** 95% percentile intervals from a country-cluster bootstrap (500 replications, seed `20260526`) that resamples whole countries and re-estimates both stages.

## Data

All inputs are public. `data/raw/` contains:

| File | Content | Original source |
|---|---|---|
| `stage1_sample.csv` | Estimation sample: 739 country-years, 37 economies, 2003–2023 | Compiled by the author from the sources below |
| `panel_full.csv` | Full panel: 38 economies, 2000–2023 (DTI sub-components, enrolment) | As above |
| `wdi_IT.MLT.MAIN.P2_fixed_telephone.csv` | Fixed-telephone subscriptions per 100 people, 1975–2000 | World Bank WDI (ITU) |
| `wdi_IT.CEL.SETS.P2_mobile_subscriptions.csv` | Mobile subscriptions per 100 people, 1975–2000 | World Bank WDI (ITU) |

Variable sources:

| Variable | Definition | Source |
|---|---|---|
| `poverty_inv` | 100 − poverty headcount ratio at USD 3.00/day (2021 PPP) | WDI `SI.POV.DDAY` (Poverty and Inequality Platform) |
| `gini_inv` | 100 − Gini index | WDI |
| `gdp_growth` | GDP per capita growth (annual %) | WDI |
| `emp_ratio` | Employment-to-population ratio, 15+ (%) | WDI |
| `fixed_bb`, `mobile_bb` | Fixed / mobile broadband subscriptions per 100 | ITU via WDI |
| `egdi` | E-Government Development Index (biennial 2005–2022, interpolated) | UN DESA E-Government Survey |
| `dfi` | Account ownership at a financial institution or mobile-money provider, % age 15+ (Findex waves 2011/2014/2017/2021, interpolated; held at 2011 value before 2011) | World Bank Global Findex |
| `dti_score` | First principal component of the four standardised sub-components above | Author |
| `tert_enrol` | Gross tertiary enrolment ratio (%) | UNESCO Institute for Statistics |
| `tert_edu_exp` | Expenditure on tertiary education, % of government expenditure on education | WDI `SE.XPD.TERT.ZS` |
| `findev` | IMF Financial Development Index (0–1); 2021–2023 held at the 2020 value | IMF Financial Development Database |
| `inst_composite` | First principal component of the six Worldwide Governance Indicators | WGI |
| `log_gdp_pc_l1` | Lagged log GDP per capita (PPP) | WDI |
| `trade`, `fdi`, `urban`, `life_exp` | Trade openness, FDI net inflows (% GDP), urban share, life expectancy | WDI |
| `ict_imports_l1` | Lagged ICT goods imports (% of total goods imports) | WDI |
| `rq_est_l1` | Lagged regulatory quality | WGI |
| `govcapex_l1` | Lagged **general government total expenditure** (% of GDP); the variable name is historical — it is not capital expenditure | IMF World Economic Outlook `GGX_NGDP` |

Data from the World Bank (WDI, Findex, WGI) are licensed CC BY 4.0. IMF and UN DESA data are used under their terms of use with attribution. Please cite the original sources when reusing the data.

**Note on the DTI scaling.** `panel_full.csv` stores an earlier scaling of the DTI. It is an almost exact linear transformation of the final index used in `stage1_sample.csv` (r = 0.99995). `01_build_data.R` maps it onto the final scale before computing the leave-one-out global DTI path.

## Reproducibility notes

- All bootstraps use fixed seeds (`20260526`; horse race `31`; placebo `4242`), so the results reproduce exactly with the same package versions.
- `output/sessionInfo.txt` records the R and package versions of the run.

## Citation

Niazi, M.Y.K. (2026). Replication package for "Digital Transformation and Poverty Reduction in Low-Human-Capital Economies: Panel Evidence with Pre-Existing Network Instruments". https://github.com/MohsinAli-TCS/Thesis
