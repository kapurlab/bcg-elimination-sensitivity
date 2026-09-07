# Sensitivity of time to elimination in the Fromsa et al. 2024 BCG model

Formal sensitivity analysis of the within-herd transmission model published in
Fromsa A, Willgert K, Srinivasan S, et al. BCG vaccination reduces bovine
tuberculosis transmission, improving prospects for elimination. Science
2024;383(6690):eadl3962. doi:10.1126/science.adl3962.

Four vaccine properties are varied independently against time to
elimination: direct efficacy (reduction in susceptibility, e_s), indirect
efficacy (reduction in infectiousness of vaccinated animals that become
infected, e_i), duration of protection (D, years) and calfhood coverage (p).

## Source model and data

The authors' code and data are public (Zenodo 10.5281/zenodo.10417489,
GitHub MonkeyMyshkin/BCGCrossover, CC BY-SA 4.0). Three inputs are copied
into `data/`:

| File | Origin | Content |
| --- | --- | --- |
| `R0Estimates.csv` | `R0Estimates/R0Estimates.csv` | 4000 posterior draws of within-herd R0 for 57 Ethiopian dairy herds |
| `herdsize.csv` | derived from `Data/TestAnon.csv` | number of tested animals per herd, used as herd size in the authors' model |
| `dst1_post.csv` | `Posterior/dst1_post.csv` | posterior draws of e_s and e_i from the DST1 chain binomial fit |

Annual mortality is 0.273 (authors' `Demography/MortalityExp.csv`).

## Model

`R/model.R` implements the authors' four-compartment model (S, I, V, IV) as a
mean-field ordinary differential equation system for a set of independent
herds. Herd size is constant, births balance deaths, and infection is
lifelong. Two extensions are needed for the question:

1. Waning of protection, V to S at rate 1/D. D = Inf reproduces the paper's
   assumption of lifelong protection with revaccination as needed. The
   reduced infectiousness of already infected vaccinated animals (IV) is
   assumed to persist; `wane_IV = TRUE` moves IV to I at the same rate.
2. Partial coverage p of newborn calves. The paper compared p = 0 and p = 1.

The reproduction number under vaccination at the disease-free state is

    R_v = R0 [ (1 - f) + f (1 - e_s)(1 - e_i) ],   f = p u / (u + 1/D)

where f is the vaccinated fraction of the herd and u the mortality rate.
Elimination is possible only when R_v < 1. This expression shows the
structure of the sensitivity before any simulation: e_s and e_i enter only
through their product (1 - e_s)(1 - e_i), and D and p enter only through f.

Between-herd cattle movements, which the paper models with ERGM-simulated
networks, are omitted. Movements re-seed herds and make elimination slower,
so the times reported here are optimistic for low-R0 herds.

## Outcome

Time to elimination is the first year in which herd prevalence of infection
falls below 0.1%, starting from the endemic equilibrium at the moment
vaccination begins. Runs are censored at 200 years. The stochastic check uses
the first day with zero infected animals instead.

Herds differ widely in R0 (posterior medians 1.3 to 17.3, median 2.8), so
results are reported for herds at the 25th, 50th, 75th and 90th percentile
of R0 and for the whole 57-herd population.

## Analyses

| Script | Method | Outputs |
| --- | --- | --- |
| `R/run_oat.R` | One parameter at a time, others at the paper's base case (e_s 0.58, e_i 0.74, lifelong, p 1). Local elasticities d ln T / d ln x by central differences. The value of each parameter at which R_v reaches 1. Duration x coverage and direct x indirect heatmaps. | `output/tables/oat_*.csv`, `local_elasticities.csv`, `elimination_thresholds.csv`, figs 1 to 4 |
| `R/run_global.R` | 3000-point Latin hypercube over e_s 0.2 to 0.9, e_i 0 to 0.95, D 2 to 40 years (log-uniform), p 0.3 to 1. Partial rank correlation coefficients with bootstrap intervals. Sobol first-order and total indices (Jansen estimator) for time to elimination and for prevalence at year 50, and a second Sobol run with herd R0 as a fifth factor. | `lhs_samples.csv`, `prcc.csv`, `sobol.csv`, figs 5 to 7 |
| `R/run_stochastic.R` | The authors' SimInf transition list with a waning event added. One-at-a-time sweeps for the median herd with 300 replicates, and the full 57-herd population at the base case with 100 replicates. | `stochastic_*.csv`, fig 8 |

Run each script from the project root with `Rscript`. Required packages:
deSolve, SimInf, lhs, sensitivity, ggplot2, dplyr, tidyr.

## Results

Base case (paper's estimates, lifelong protection, full coverage). Time to
herd prevalence below 0.1%, deterministic model:

| Herd | R0 | R_v | Years to < 0.1% |
| --- | --- | --- | --- |
| 25th percentile | 2.0 | 0.22 | 36 |
| median | 2.8 | 0.31 | 41 |
| 75th percentile | 4.1 | 0.45 | 49 |
| 90th percentile | 9.1 | 1.00 | never |

Across the 57 herds, 89% have R_v below 1, 75% fall below 0.1% within 50
years, and population animal prevalence at year 50 is 4.4%. This matches the
paper's description of a slow path with prevalence still above the
officially-free level at 50 years. The three herds with R0 above 12 never
eliminate under any vaccine property examined, because a vaccine with 89%
total efficacy at full coverage cannot bring R0 = 17 below 1.

Local elasticities, d ln T / d ln x, for the median herd at the base case.
Duration is evaluated at 20 and 30 years because the lifelong base has zero
slope by construction.

| Parameter | Evaluated at | Elasticity |
| --- | --- | --- |
| Coverage | 1.0 | -2.4 |
| Indirect efficacy | 0.74 | -0.75 |
| Duration | 20 y | -0.74 |
| Duration | 30 y | -0.38 |
| Direct efficacy | 0.58 | -0.52 |

Elasticities grow steeply with herd R0 (coverage -1.4 at R0 2.0, -5.0 at
R0 4.1) because every curve has a vertical asymptote where R_v reaches 1.

Value of each parameter, others at base, at which the median herd loses the
possibility of elimination (R_v = 1), and at which time to elimination
doubles from 41 to 81 years:

| Parameter | Doubling | R_v = 1 |
| --- | --- | --- |
| Coverage | 0.82 | 0.72 |
| Duration | 16 y | 9.5 y |
| Indirect efficacy | 0.35 | 0.15 |
| Direct efficacy | never (84 y at e_s = 0) | never |

Global analysis, median herd. Partial rank correlation coefficients with
R_v and with log10 prevalence at year 50 rank coverage first (-0.91 and
-0.91), then duration (-0.87, -0.79), indirect efficacy (-0.76, -0.77) and
direct efficacy (-0.70, -0.81). Over the sampled ranges only 5% of parameter
combinations have R_v below 1, so time to elimination is censored in 96% of
samples and its PRCCs are weak (-0.28 to -0.14) but in the same order. Sobol
total-effect indices for time to elimination give the same ranking: coverage
0.87, duration 0.71, indirect efficacy 0.38, direct efficacy 0.31. First-order
indices are close to zero because the four inputs act through two products,
(1 - e_s)(1 - e_i) and p u / (u + 1/D), and through the R_v = 1 threshold,
so almost all of the variance is interaction. For the uncensored outcome,
log10 prevalence at year 50, the total effects are coverage 0.64, duration
0.44, direct 0.22 and indirect 0.21. When herd R0 is added as a fifth factor
over the observed range it has the largest total effect of all (0.79),
ahead of coverage (0.50).

Stochastic check. In the SimInf model the median herd of 44 animals reaches
zero infected animals in a median of 23 years at the base case (IQR 18 to
27), against 41 years for the deterministic threshold, because 0.1% of 44
animals is far below one animal and chance fade-out removes the last
infected animal earlier. The shape of every one-at-a-time curve, including
the steep rise near R_v = 1, is reproduced. Stochastic fade-out also allows
elimination in small herds with R_v slightly above 1, which softens the
cliff. In 100 replicates of the full 57-herd population no replicate reached
zero infected animals within 200 years, because of the three high-R0 herds.

## Interpretation

Time to elimination is not a linear function of any of the four properties.
It is a hyperbola in R_v with a vertical asymptote at R_v = 1, so the slope
depends entirely on how far the base case sits from that asymptote. For the
median Ethiopian herd the paper's base case sits comfortably below it
(R_v = 0.31), which is why the paper's result is robust to the efficacy
estimates: direct efficacy can fall to zero and the median herd still
eliminates in 84 years, and indirect efficacy can fall to 0.35 before time
doubles. The two properties the paper held fixed, coverage and duration, are
the ones with steep slope. Coverage below 82%, or a mean duration of
protection shorter than 16 years without revaccination, doubles the time to
elimination in the median herd, and coverage below 72% or duration below
9.5 years makes it impossible. The paper's assumption of revaccination as
needed is therefore not a detail. It is the assumption the elimination
result rests on.

Direct and indirect efficacy are interchangeable in this model through
their product, and indirect efficacy carries more weight only because the
base estimate is higher. Duration and coverage are interchangeable through
the vaccinated fraction of the herd, p u / (u + 1/D). With a 3.7-year mean
lifespan, protection lasting 10 years already yields 73% of the lifelong
vaccinated fraction, which is why the duration curve flattens beyond about
15 years.

## Limitations

Herds are independent; the paper's between-herd movement network is
omitted. Herd R0 is fixed at each herd's posterior median. The deterministic
threshold of 0.1% is stricter than one animal for herds under 1000 head.
Waning is exponential, which is the simplest choice, and the reduced
infectiousness of infected vaccinated animals is assumed permanent.
