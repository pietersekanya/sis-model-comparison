# Comparing a demographic SIS model between Ituri and North Kivu: the 2026 Bundibugyo virus disease outbreak

A MATLAB project for Biomathematics 374. A demographic susceptible–infectious–susceptible (SIS) model is fitted separately to confirmed **cases and deaths** from the 2026 Bundibugyo virus disease outbreak in **Ituri** and **North Kivu** provinces (Democratic Republic of the Congo). The project compares the estimated parameters between the two provinces and compares how well the model predicts future cases and deaths in each, using leave-future-out cross-validation.

**Author:** Kanya Pieterse | **Course:** Biomathematics 374, Stellenbosch University | **Date:** October 2026

---

## Contents

- [Research questions and hypotheses](#research-questions-and-hypotheses)
- [Data](#data)
- [Model](#model)
- [Methods in brief](#methods-in-brief)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [How to run](#how-to-run)
- [Settings you can change](#settings-you-can-change)
- [Outputs](#outputs)
- [Key results](#key-results)
- [Limitations](#limitations)
- [Reproducibility](#reproducibility)
- [Troubleshooting](#troubleshooting)
- [References and data citation](#references-and-data-citation)

---

## Research questions and hypotheses

1. Do the estimated transmission and disease-related mortality parameters, and the quantities derived from them, differ between Ituri and North Kivu?
2. How well does the model predict future cases and deaths in each province, and does predictive performance differ between the provinces?
3. How sensitive are these comparisons to the assumed infectious period?

The hypotheses are:

- **H1:** the disease-related mortality rate, and the model-implied case fatality proportion, differ between the provinces.
- **H2:** the transmission parameter and the initial reproduction measure differ between the provinces, and the initial reproduction measure is greater than one in both.
- **H3:** predictive performance (forecast error relative to the observed counts) differs between the provinces.

## Data

Cumulative confirmed cases and deaths **by health zone and report date** from the INRB-UMIE BDBV2026-Data repository (see [citation](#references-and-data-citation)):

| File (in `data/`) | Columns | Content |
|---|---|---|
| `insp_sitrep__cumulative_confirmed_cases__daily.csv` | `nom`, `date`, `cumulative_confirmed_cases` | Cumulative confirmed cases by health zone |
| `insp_sitrep__cumulative_confirmed_deaths__daily.csv` | `nom`, `date`, `cumulative_confirmed_deaths` | Cumulative confirmed deaths by health zone |

Health zones are assigned to provinces using the lists in `load_province_cases.m` (the same lists are repeated in `load_province_deaths.m`), based on [source of the health-zone lists, e.g. Humanitarian Data Exchange, 2024].

What the loaders do:

- **Cutoff date.** The cutoff is set once in `load_province_cases.m` (currently **12 October 2026**; the latest reports are dated 9 October 2026, giving 118 report dates for cases from 14 May 2026 and 117 for deaths). Do not change it after you have produced your final results.
- **Cleaning.** Variant spellings of health-zone names are matched; non-numeric entries ("ND", "NA") are treated as missing and carried forward; malformed dates and counts are repaired where possible; cases or deaths not assigned to a health zone are excluded (they are returned separately as `Unassigned`).
- **Aggregation.** A health zone missing on a given date keeps its last reported value, and a zone that has not yet reported counts as zero. Province totals that decrease between reports are corrected by keeping the highest value reached so far (the loaders print how many corrections were made).
- **Deaths on the case dates.** Deaths are returned on the same report dates as the cases; where no death report exists for a case date, the latest earlier death count is carried forward.
- **No interpolation.** Only actual reporting dates are used. Days without a report are not assigned zero, and the intervals between reports vary (1 to 7 days).

The data are by **report date**, not symptom-onset date. Some reporting intervals contain unusually large counts (see [Methods in brief](#methods-in-brief)).

## Model

For each province, with $S$ susceptible and $I$ infectious individuals (people), time $t$ in days, and cumulative model cases $C$ and deaths $D$:

$$
\begin{aligned}
\frac{dS}{dt} &= bN - \beta S I + \gamma I - \mu S \\
\frac{dI}{dt} &= \beta S I - \gamma I - \tilde{\mu} I \\
\frac{dC}{dt} &= \beta S I \\
\frac{dD}{dt} &= \tilde{\mu} I
\end{aligned}
$$

| Symbol | Meaning | Value or bounds | Units | Status |
|---|---|---|---|---|
| $\beta N$ | Scaled transmission rate ($\beta = \beta N / N$) | 0.01 to 5.00 | day⁻¹ | Fitted |
| $\tilde{\mu}$ | Disease-related mortality rate | 0 to 0.5 | day⁻¹ | Fitted |
| $I_0$ | Initial infectious individuals | 0.1 to 2,000 | people | Fitted |
| $\gamma$ | Recovery rate | 1/10 (sensitivity: 1/7 and 1/14) | day⁻¹ | Fixed |
| $b$ | Birth rate | $1.120\times10^{-4}$ | day⁻¹ | Fixed |
| $\mu$ | Natural mortality rate | $2.307\times10^{-5}$ | day⁻¹ | Fixed |
| $N$ | Population | Ituri 4,600,499; North Kivu 9,314,854 | people | Fixed |

The same birth and natural mortality rates are used for both provinces. **γ is fixed** because the case and death data cannot separate it from the transmission and mortality parameters: in fits to the case data with β, γ and μ̃ all estimated, γ and μ̃ settled on their bounds and were almost perfectly correlated.

Two derived quantities are reported:

- initial reproduction measure $R_{\mathrm{start}} = \beta N / (\gamma + \tilde{\mu})$ (assumes $S_0 \approx N$ and does not allow for changes in transmission over time);
- model-implied case fatality proportion $\mathrm{CFR}_{\mathrm{model}} = \tilde{\mu} / (\gamma + \tilde{\mu})$ (not the observed clinical case fatality ratio).

## Methods in brief

- **Observation model.** Cumulative counts are converted to counts in each reporting interval, $y_j = Y(t_j) - Y(t_{j-1})$, and compared with the model's change in $C$ and $D$ over the same interval.
- **Fitting.** Parameters are estimated separately for each province by nonlinear least squares (`lsqnonlin`) with the ODEs solved by `ode45` (relative tolerance $10^{-6}$, absolute tolerance $10^{-8}$). The residual for each interval is the square root of the predicted count minus the square root of the observed count; case and death residuals are stacked into one vector. Ten random starting points are used and the best fit is kept. Each province's fitting window starts at the last report before its first confirmed case.
- **Sensitivity analysis.** Both provinces are refitted with infectious periods of 7, 10 and 14 days.
- **Leave-future-out cross-validation.** Forecast origins start at the 40th reporting date and advance by 5 reports (15 origins per province). At each origin the model is refitted to the data up to that date and used to predict the next 5 reporting intervals; accuracy is measured separately for cases and deaths with the mean absolute error (MAE) and root mean squared error (RMSE) per reporting interval.
- **Large reporting intervals.** The script lists the five largest case counts and five largest death counts in a single interval for each province, with the number of days each interval covers. Most follow multi-day intervals, which suggests delayed and batched reporting. No intervals were excluded from the analysis.

## Repository structure

```
.
├── README.md
├── load_province_cases.m      Loads and cleans the case data, sums health zones by province; sets the cutoff date
├── load_province_deaths.m     Loads and cleans the death data on the same report dates
├── province_settings.m        Population, birth and death rates for each province (converted to per-day rates)
├── simulate_sis_cum.m         Demographic SIS model: returns cumulative cases and deaths
├── fit_sis_model.m            Nonlinear least squares fit (lsqnonlin, multiple starting points)
├── run_province_comparison.m  MAIN ANALYSIS: fits both provinces, sensitivity analysis, cross-validation
├── check_province_data.m      Optional: checks and plots the province case series
├── check_province_deaths.m    Optional: checks and plots the death series, saves cases and deaths together
├── cumulativetocase_provinces.m  Optional: counts per reporting interval and exploratory plots
├── data/                      Input data (and files written by the check scripts)
└── results/                   Tables and figures written by the scripts
```

## Requirements

- MATLAB R2026b
- Optimization Toolbox (for `lsqnonlin` and `optimoptions`)

No other toolboxes are needed.

## How to run

1. **Get the code and data.**
   ```
   git clone https://github.com/pietersekanya/sis-model-comparison.git
   ```
   Check that both CSV files are in `data/` with exactly the names in the [Data](#data) table. If they are missing, download them from the BDBV2026-Data repository (the "raw" version of each file) and save them there.
2. **Set the Current Folder in MATLAB** to the top-level repository folder (the one containing `load_province_cases.m` and `data/`). Everything is run from there.
3. **Fill in the province inputs.** `province_settings.m` takes the population, birth rate and death rate for each province (births and deaths per 1000 people per year) and converts them to per-day rates.
4. **Check that the data load** (optional):
   ```matlab
   check_province_data
   check_province_deaths
   ```
5. **Run the main analysis:**
   ```matlab
   run_province_comparison
   ```
   This takes several minutes because of the repeated fitting in the sensitivity analysis and the cross-validation. It prints the estimates and error tables and saves tables and figures to `results/`.

Only the scripts (`run_province_comparison`, `check_province_data`, `check_province_deaths`, `cumulativetocase_provinces`) are meant to be run directly. The other `.m` files are functions called by them.

## Settings you can change

| Setting | Where | Current value |
|---|---|---|
| Cutoff date and data file names | `load_province_cases.m` | 12 October 2026 |
| Health-zone lists and spelling variants | `load_province_cases.m` and `load_province_deaths.m` (keep identical) | — |
| Population, birth rate, death rate | `province_settings.m` | see [Model](#model) |
| Which parameters are fitted, fixed values, bounds | `free`, `p_fixed`, `lb`, `ub` in `run_province_comparison.m` | γ fixed at 1/10 |
| Infectious periods for the sensitivity analysis | `infectious_periods` | 7, 10, 14 days |
| Random starting points | `n_starts` | 10 |
| Cross-validation: first origin, step, horizon, minimum training reports | `first_origin`, `origin_step`, `h`, `min_train` | 40, 5, 5, 15 |

## Outputs

Written to `results/` by `run_province_comparison.m`:

| File | Description |
|---|---|
| `sis_province_parameters.csv` | Fitted parameters, derived quantities and in-sample error metrics for each province |
| `sis_province_parameters_with_se.csv` | Parameter estimates, approximate standard errors and 95% confidence intervals |
| `sis_province_comparison.csv` | Between-province differences, approximate confidence intervals, z-statistics and p-values |
| `sis_province_gamma_sensitivity.csv` | Fitted parameters and derived quantities for assumed infectious periods of 7, 10 and 14 days |
| `sis_province_identifiability.csv` | Condition numbers, maximum absolute parameter correlations and parameters at bounds |
| `sis_province_cv_summary.csv` | Average leave-future-out cross-validation errors for cases and deaths |
| `sis_province_cv_by_origin.csv` | Forecast errors at individual forecast origins |
| `sis_province_largest_intervals.csv` | Five largest reporting intervals for cases and deaths in each province, including interval duration |
| `sis_province_fits.png` | Observed and fitted cases and deaths per reporting interval |
| `sis_province_parameters.png` | Parameter estimates with approximate 95% confidence intervals |
| `sis_province_cv_error.png` | Case and death forecast RMSE across forecast origins |

The optional data-checking scripts also generate intermediate CSV files and exploratory figures in `data/` and `results/`.


## Key results

Latest run, data up to 9 October 2026 (the numbers change if the cutoff date or any setting changes):

| | Ituri | North Kivu |
|---|---|---|
| $\beta N$ (day⁻¹) | 0.163 | 0.249 |
| $\tilde{\mu}$ (day⁻¹) | 0.058 | 0.132 |
| $R_{\mathrm{start}}$ | 1.031 | 1.076 |
| $\mathrm{CFR}_{\mathrm{model}}$ | 36.6% | 56.9% |
| $I_0$ (people) | 158.6 | 8.7 |
| In-sample MAE / RMSE, cases | 24.7 / 45.6 | 6.1 / 8.9 |
| In-sample MAE / RMSE, deaths | 11.6 / 24.3 | 4.7 / 7.1 |
| Cross-validation MAE / RMSE, cases | 56.7 / 69.0 | 11.5 / 13.3 |
| Cross-validation MAE / RMSE, deaths | 22.1 / 29.6 | 7.9 / 8.8 |

**Comparison between provinces**

Approximate two-sided z-tests, comparing Ituri minus North Kivu:

| Parameter | Difference | Approx. 95% CI | p-value |
|---|---:|---:|---:|
| Scaled transmission rate, βN | −0.0869 | −0.1399 to −0.0339 | 0.0013 |
| Disease-related mortality rate, μ̃ | −0.0741 | −0.1272 to −0.0210 | 0.0062 |
| Initial reproduction measure, R_start | −0.0452 | −0.0701 to −0.0204 | 0.00036 |
| Model-implied case fatality proportion | −20.3 percentage points | −32.4 to −8.2 percentage points | 0.0010 |

These are approximate tests conditional on the fixed recovery rate; they should not be interpreted as definitive evidence of differences in the underlying biological processes.

North Kivu had the higher fitted transmission parameter, mortality rate, initial reproduction measure and model-implied case fatality proportion, and this ordering held for infectious periods of 7, 10 and 14 days. Errors are in counts per reporting interval and are larger in Ituri, which reported many more cases and deaths. See the report for the full results and interpretation.

Cross-validation WAPE was 0.84 for cases and 0.83 for deaths in Ituri, compared with 0.57 for cases and 0.78 for deaths in North Kivu. Thus, relative case-forecast errors were larger in Ituri, while relative death-forecast errors were similar. These comparisons are descriptive because the forecast windows overlap and no baseline forecast was evaluated.

## Limitations

- Parameters are assumed constant over time, so the model cannot represent changes in transmission (for example the decline after the August peak, which the fitted curves do not capture).
- The SIS structure ignores immunity among survivors (fewer than 0.2% of each province's population was reported as cases, so the effect within the study window is small).
- The data are by report date: reporting delays, batching and under-reporting are not modelled. Some large counts follow multi-day reporting gaps and were kept in the analysis.
- Cases and deaths not assigned to a health zone are excluded, and the same birth and death rates are used for both provinces.
- γ is fixed, so the estimates are conditional on the assumed infectious period. The scripts calculate approximate standard errors, confidence intervals and z-tests from the Jacobian. These estimates are conditional on the fixed recovery rate and may understate uncertainty because residuals from cumulative data are autocorrelated and the model does not capture all observed features. The fitted transmission and disease-related mortality parameters remain strongly correlated, so they should be interpreted cautiously.
- The cross-validation uses only 15 overlapping forecast windows and no baseline forecast.
- The outbreak is ongoing and later data may differ from those used here.

## Reproducibility

- The cutoff date is fixed in `load_province_cases.m` and the random number generator is seeded (`rng(1)`), so the same data and settings give the same random starting points.
- Small numerical differences are possible across MATLAB versions and operating systems.
- The data are from the repository's `main` branch, which can change; a copy of the files used is stored in `data/` [if the data licence allows redistribution].

## Troubleshooting

| Problem | Likely cause and fix |
|---|---|
| "Unable to find or open ... .csv" | The data files are not in `data/`, the names differ, or the Current Folder is not the top-level folder |
| "Fill in pop, birth_rate and death_rate" | `province_settings.m` has not been completed |
| "Case and death dates do not match." | The loaders were changed separately; the deaths loader takes its dates and cutoff from the cases loader, so keep the zone lists identical |
| Warnings that a parameter is near a bound | Check whether the estimate is genuinely constrained by the bound and whether it is stable across starting points. Cross-validation fits may produce small initial infectious populations. Do not widen bounds automatically; investigate identifiability and assess whether the change is epidemiologically justified. |
| Undefined function or variable | A required `.m` file is missing from the top-level folder |

## References and data citation

INRB-UMIE (2026) *BDBV2026-Data: Data and scripts for epidemiological analysis of the 2026 Bundibugyo Ebola outbreak* [dataset]. GitHub. Available at: https://github.com/INRB-UMIE/BDBV2026-Data (Accessed: 12 October 2026). Please see that repository for the terms of use of the data. 


