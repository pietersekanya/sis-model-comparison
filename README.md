# Basic vs demographic SEIR: does added complexity improve accuracy?

A MATLAB project comparing a basic SEIR model with a demographic SEIR model (births and deaths included) on cumulative confirmed case data from the 2026 Ebola (Bundibugyo) outbreak in the Democratic Republic of the Congo. Models are fitted by nonlinear least squares and compared using leave-future-out cross-validation.

**Author:** [your name] | **Course:** [course name, institution] | **Date:** October 2026

---

## Research question
Does adding demographic terms (births and deaths) to an SEIR model improve its accuracy when fitted to outbreak case data?

The comparison is made on **out-of-sample forecast error**, not on in-sample fit, because extra model structure can lower in-sample error without improving predictions.

## Contents
- [Data](#data)
- [Models](#models)
- [Methods in brief](#methods-in-brief)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [How to run](#how-to-run)
- [Settings you can change](#settings-you-can-change)
- [Outputs](#outputs)
- [Results](#results)
- [Limitations](#limitations)
- [Reproducibility](#reproducibility)
- [References](#references)

## Data
Cumulative confirmed cases for the DRC, by report date, from the INRB-UMIE BDBV2026-Data repository:

- **File:** `data/insp_sitrep/processed/insp_sitrep_national_cumulative_confirmed_cases_daily.csv`
- **Version (tag):** `build-2026-10-01-b077a6c`
- **Accessed:** 3 October 2026

## Models
* Provide a brief overview of the two mathematical models compared here:
  * **Basic SEIR Model:** Standard Susceptible-Exposed-Infectious-Recovered equations without population turnover.
  * **Demographic SEIR Model:** SEIR framework integrated with vital dynamics (birth rates and non-Ebola background death rates).

## Methods in brief
1. **Model Fitting:** Parameters are estimated via nonlinear least squares optimization minimizing the distance to the cumulative confirmed case data.
2. **Validation:** Models are evaluated using leave-future-out cross-validation to simulate real-time forecasting performance and track out-of-sample error accumulation.

## Repository structure
```text
├── data/
│   └── insp_sitrep/
│       └── processed/
│           └── insp_sitrep_national_cumulative_confirmed_cases_daily.csv
├── src/
│   ├── fit_models.m
│   ├── cross_validate.m
│   └── plot_results.m
├── outputs/
├── README.md
└── main.m
```

## Requirements
- MATLAB (Recommended version R2023a or newer)
- Optimization Toolbox

## How to run
1. Open MATLAB and navigate to the project root directory.
2. Run the primary analysis script:
   ```matlab
   run('main.m')
   ```

## Settings you can change
Configuration parameters can be modified inside `main.m` (or your dedicated configuration file):
- `forecast_window`: Number of days to project forward during cross-validation.
- `initial_params`: Starting parameter guesses for the nonlinear least squares solver.

## Outputs
Execution of the pipeline generates:
- Fitted parameter tables saved to the `outputs/` directory.
- Comparative forecasting performance plots showing cumulative trajectories versus actual out-of-sample data points.

## Results
* Add your findings here after running the cross-validation pipeline (e.g., whether the demographic terms reduced out-of-sample forecast error significantly or led to overfitting).*

## Limitations
- Assumes homogenous mixing within the population.
- Relies on cumulative reported data, which may be subject to reporting lags or under-ascertainment during the early stages of the outbreak.

## Reproducibility
To ensure exact replication of the results, use the specific data version tag (`build-2026-10-01-b077a6c`) and keep random number generator seeds fixed if stochastic optimization algorithms are added.

## References
* List any papers, textbook chapters, or data repository links utilized during this research project here.*
