# Sensitivity of time to elimination in the Fromsa et al. 2024 BCG model

Formal sensitivity analysis of the within-herd transmission model published in
Fromsa A, Willgert K, Srinivasan S, et al. BCG vaccination reduces bovine
tuberculosis transmission, improving prospects for elimination. Science
2024;383(6690):eadl3962. doi:10.1126/science.adl3962.

Four vaccine properties are varied independently against time to
elimination: direct efficacy (reduction in susceptibility, e_s), indirect
efficacy (reduction in infectiousness of vaccinated animals that become
infected, e_i), duration of protection (D, years) and calfhood coverage (p).

## Start here

For a reader who wants the model and the main result without the rest:

1. `R/model.R` is the deterministic version of the S, I, V, IV herd model
   with waning and coverage added; `R/boost.R` adds a birth dose and
   whole-herd booster campaigns with separate durations; `R/onset.R` adds a
   pre-immunity window after birth. Each is under 120 lines.
2. `R/run_memo_figure.R` and `R/run_tiles.R` produce the two figures used in
   the memo (`output/figs/fig14_memo_scenarios.png` and
   `output/figs/fig16_interval_by_duration.png`) in a few minutes each.
3. `translations/` holds English translations of Calmette and Guérin 1920
   and 1924 and Guérin, Richart and Boissière 1927, the papers the
   revaccination schedule goes back to.
4. `sessionInfo.txt` records the R and package versions used.

The authors' original model is reproduced, not modified: the transition
list in `R/run_stochastic.R` is theirs, and `data/` holds three files copied
from their repository (see `LICENSE.md`).

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
| `R/run_onset.R` | Interval from birth to effective immunity (0 to 180 days, covering delay to vaccination and delay to onset) with calf exposure at 1, 3 or 10 times the adult rate, under the paper's regime and under 18-month protection with annual boosting. Model in `R/onset.R`. | `onset*.csv`, figs 12 and 13 |
| `R/run_profile.R` | Protection that rises over 30 days, plateaus, and declines with time since dose, tracked by a 40-stage clock; five profile shapes against four revaccination schedules. Model in `R/profile.R`, which also takes a hazard-ratio curve from a trial (`profile_from_hazard`). | `profile_schedule.csv`, figs 17 and 18 |
| `R/run_nonresponders.R` | Persistent non-responders (never protected by any dose) against random per-dose non-response; analytic ceiling on the tolerable persistent fraction by herd R0. | `nonresponders.csv`, figs 19a and 19b |
| `R/run_import.R` | Infected animals bought into the herd each year, under the paper's regime and under 18-month immunity with annual revaccination. | `imports.csv`, fig 20 |
| `R/run_import_uk.R` | The same import sweep in a low-R0 setting with annual test-and-removal (R0 1.1, infected animals removed at 0.7 per year), against the Ethiopian setting, with and without vaccination. | `imports_uk.csv`, fig 21 |
| `R/run_stochastic_grid.R` | Stochastic (SimInf) elimination on a conceptual grid of within-herd R0 (1.2 to 8) and herd size (5 to 200) under four programmes, 200 replicate herds per cell; probability of a herd with no infected animal by years 10, 20, 50 and 100, and the median time; no-vaccination reference. | `stochastic_grid*.csv`, figs 22 and 23 |
| `R/run_grid_precise.R` | The same R0 by herd-size grid at 2000 replicates per cell, with mean within-herd prevalence over 20 and 50 years as a continuous outcome alongside probability of freedom. Reports the share of infected animal-years averted against no vaccination, with bootstrap intervals, and the share of the lifelong benefit an 18-month vaccine retains. Caches per-replicate results in `output/raw/` so the summaries can be recomputed without rerunning. | `grid_precise_cells.csv`, `grid_precise_comparisons.csv`, figs 35 and 36 |
| `R/run_fig38.R` | Plotting layer only for figure 38: reads the tables written by `run_benefit_timing.R` and draws the 20 strata separately by herd size and R0, with milestone timings in panel C. Writes PNG and vector PDF, and exports the plotted data for each panel. | `fig38_panelA_data.csv`, `fig38_panelB_data.csv`, `fig38_panelC_data.csv`, fig 38 (png and pdf) |
| `R/run_benefit_timing.R` | Follows 20 demographic strata through time under no vaccination and under an 18-month calves-only vaccine, and reports when the burden reduction is achieved against when herds become free of infection. Figure 37 pools the strata. Caches the trajectories so reruns are instant; figure 38 is drawn separately by `R/run_fig38.R`. | `benefit_timing.csv`, `benefit_timing_summary.csv`, figs 37 and 38 |
| `R/run_herd.R` | Age-structured stochastic herd (calves, young stock, adults) with a fixed adult herd and a replacement policy: leavers replaced by home-bred young stock or by purchases from a source population of given prevalence. Adult herd size 5, 20, 100 by R0 1.5, 3, 5 by seven policies, with and without vaccination. Model in `R/herd.R`. | `herd_replacement.csv`, figs 24 and 25 |
| `R/run_alt_views.R` | Alternative views of the stochastic results: difference tiles (programme minus status quo), feasibility frontier contours, archetype slices, time courses of freedom for six representative herds, and prevalence bands under replacement policies. | `survival_curves.csv`, `replacement_trajectories.csv`, figs 26 to 30 |
| `R/run_frontier_fine.R` | Feasibility frontier on an 11 by 10 grid with the 57 survey herds overlaid, by years 20 and 50. | `frontier_fine.csv`, fig 31 |
| `R/run_memo_stochastic.R` | Stochastic versions of the two memo figures for a 100-head herd at the survey-median R0: scenario ladder with status quo, and interval-by-duration grid; median years to no infected animal and share free by year 50. General runner in `R/stochastic.R`, which adds a pre-immunity window and coverage to the SimInf model. | `memo_scenarios_stochastic.csv`, `interval_duration_stochastic.csv`, figs 32 and 33 |
| `R/run_memo_extras.R` | Simplified replacement-policy figure and the archetype table for the memo, from existing tables. | `memo_archetype_table.csv`, fig 34 |
| `R/run_tiles.R` | Booster schedule (none, every 24 months, annual, every 6 months) against birth-dose duration (6 to 24 months), with the booster's own duration equal to, twice, or independent of the birth dose, and a fourth panel for campaign coverage. Model in `R/boost.R`. | `tiles_booster_schedule.csv`, fig 15 |
| `R/run_memo_figure.R` | Scenario ladder: one lever changed per row from the paper's base case. | `memo_scenarios.csv`, fig 14 |
| `R/run_revaccination.R` | Duration of protection from 3 months to lifelong under three strategies: calves at birth only; calves plus an annual campaign vaccinating every unprotected uninfected animal; calves plus an annual campaign vaccinating every uninfected animal, restarting protection in those still protected. Each under gradual (exponential) and near-fixed (Erlang, k = 20) waning. Deterministic times for representative herds and the 57-herd population, time-averaged protected fraction, and a SimInf check. Model in `R/revaccination.R`. | `revacc_*.csv`, figs 9 to 11 |

## Running it in RStudio

Open `Model_Sensitivity_Tradeoffs.Rproj`. That sets the working directory to
the project root, which every script expects. Then either source
`run_all.R`, which installs any missing packages and runs the three analyses
in order, or open the scripts under `R/` and run them one at a time. Each
starts with `R/root.R`, which checks the working directory, and `R/model.R`,
which holds the model functions. Approximate run times on a laptop: one
minute for `run_oat.R`, 15 minutes for `run_global.R`, five minutes for
`run_stochastic.R`. From a terminal, `Rscript run_all.R` does the same.

Packages: deSolve, SimInf, lhs, sensitivity, ggplot2, dplyr, tidyr.

Layout:

    R/model.R           model, R_v, time to elimination, palette
    R/run_oat.R         one-at-a-time sweeps, elasticities, thresholds, heatmaps
    R/run_global.R      Latin hypercube, PRCC, Sobol
    R/run_stochastic.R  SimInf check
    R/revaccination.R   model with protection stages and campaign pulses
    R/run_revaccination.R  duration of protection under annual revaccination
    R/onset.R           model with a pre-protection window after birth
    R/run_onset.R       speed of onset and calf exposure
    R/boost.R           birth dose and booster with separate durations
    R/profile.R         time-since-dose protection profile, non-responders, imports
    R/run_profile.R     profile shape against schedule
    R/run_nonresponders.R  persistent against random non-response
    R/run_import.R      infected purchases
    R/run_import_uk.R   infected purchases with test-and-removal, UK-like R0
    R/run_stochastic_grid.R  stochastic elimination on a conceptual R0 by herd-size grid
    R/herd.R            age-structured stochastic herd with replacement policy
    R/run_herd.R        replacement policy against herd size and R0
    R/run_alt_views.R   difference tiles, frontier, archetype slices, time courses, prevalence bands
    R/run_frontier_fine.R  finer frontier with the survey herds overlaid
    R/stochastic.R      general SimInf runner: coverage, waning, campaigns, onset window
    R/run_memo_stochastic.R  stochastic memo ladder and grid
    R/run_memo_extras.R  memo replacement figure and archetype table
    R/run_tiles.R       booster schedule tile figure
    R/run_memo_figure.R scenario ladder
    data/               three inputs copied from the authors' repository, plus archetypes and the cost template
    output/figs/        fig1 to fig15 (png)
    output/tables/      csv results
    translations/       English translations of the 1920, 1924 and 1927 papers (CC BY 4.0)

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

Revaccination and duration of protection (`run_revaccination.R`). Years to
herd prevalence below 0.1% in the median herd, deterministic model, with all
calves vaccinated at birth and an annual campaign added:

| Duration | Calves only | Annual, unprotected animals | Annual, all uninfected animals |
| --- | --- | --- | --- |
| 6 months, near-fixed | never | never | never |
| 12 months, near-fixed | never | never | 49 |
| 18 months, near-fixed | never | 107 | 40 |
| 24 months, near-fixed | never | 67 | 40 |
| 36 months, near-fixed | never | 52 | 40 |
| 12 months, gradual | never | never | never |
| 24 months, gradual | never | 99 | 99 |
| 36 months, gradual | never | 65 | 65 |
| Lifelong | 41 | 40 | 40 |

Three things decide the outcome. First, an annual campaign that only doses
animals whose protection has lapsed leaves gaps: a 12-month vaccine given
once a year still lets half the campaign cohort and every mid-year calf
lapse before the next round, so the herd averages 71% protected and R_v is
1.03. Second, dosing every uninfected animal at each campaign, protected or
not, closes those gaps once the duration exceeds the 12-month interval,
and from 18 months on it is indistinguishable from lifelong protection.
Third, with gradual (exponential) waning, re-dosing a still-protected
animal changes nothing, because exponential waning has no memory, so the
two campaign strategies coincide and a mean duration of about 24 months is
the shortest that eliminates. Six-month protection fails under every
strategy with an annual interval. The stochastic check for the 44-animal
median herd reproduces the ordering with shorter absolute times; a
12-month near-fixed vaccine with re-dosing of all uninfected animals
reaches zero infected animals in a median of 24 years, against 22 for
lifelong protection. Across the 57 herds, the share eliminated by year 50
under re-dosing of all uninfected animals reaches the lifelong value of 75%
at 18 months, and under lapsed-only dosing at about 20 years.

Speed of onset (`run_onset.R`). Years to elimination in the median herd
when calves are fully susceptible for a window after birth, with exposure
in that window at 1, 3 or 10 times the adult rate, paper regime:

| Window | Exposure x1 | Exposure x3 | Exposure x10 |
| --- | --- | --- | --- |
| 0 days | 41 | 41 | 41 |
| 15 days | 42 | 45 | 65 |
| 30 days | 43 | 51 | never |
| 60 days | 46 | 71 | never |
| 90 days | 50 | 121 | never |
| 180 days | 63 | never | never |

At adult-equivalent exposure a 30-day window costs three years; at
three-fold exposure it costs ten and a 90-day window triples the time; at
ten-fold exposure any window of 30 days or more removes the possibility of
elimination in the median herd. The share of vaccinated calves infected in
the window over the first five years of a programme is 3%, 8% and 27% for a
30-day window at the three exposure levels. Across the 57 herds, the share
eliminated by year 50 falls from 75% with no window to 70%, 49% and 5% for
a 30-day window at the three exposure levels. The 18-month boosted regime
behaves the same, because the window is the same.

Efficacy saturation, lifespan and which effect wanes. With lifelong
protection and full coverage, raising total efficacy from the paper's 89%
to 95%, 98% and 100% shortens the median herd's time from 41 years to 34,
31 and 28. Efficacy is therefore not fully saturated, but its whole range
buys 13 years while duration and coverage decide feasibility. In herds
with longer residence the duration requirement rises: with a six-year mean
residence the duration needed for R_v < 1 in the median herd under
calfhood-only vaccination is 16 years, against 9.6 years at the Ethiopian
3.7-year residence. The model separates two durations. Waning of
susceptibility protection (V to S) is what the sweeps vary. The reduction
in infectiousness of an animal infected while protected is assumed to last
for its remaining life. If that reduction also wanes, elimination in the
median herd under calfhood-only vaccination needs about 20 years of
protection instead of 15.

Booster schedule (`run_tiles.R`). Years to elimination in the median herd
with every calf dosed at birth and whole-herd booster campaigns, booster
lasting as long as the birth dose:

| Schedule | 6 months | 12 months | 18 months | 24 months |
| --- | --- | --- | --- | --- |
| No booster | never | never | never | never |
| Every 24 months | never | never | 84 | 47 |
| Annually | never | 49 | 40 | 40 |
| Every 6 months | 50 | 40 | 39 | 39 |

The rule is that the booster interval must be shorter than the duration
of protection. Where it is, the result equals the paper's lifelong case
within a year or two; where it is longer, the gap between lapse and the
next campaign decides, and 12-month protection with annual boosting
(49 years) is the only combination that eliminates with interval equal to
duration. A booster lasting 5 years makes every schedule work and makes
birth-dose duration nearly irrelevant. Campaign coverage matters as much
as the interval: with 18-month protection and annual boosters, 90% coverage
costs 4 years, 70% costs 24, and 50% removes elimination.

Shape of protection (`run_profile.R`). Replacing all-or-nothing protection
with a curve that rises, holds at 58% and then declines changes the reading
of the schedule grid (figs 17 and 18). Under annual or six-monthly
revaccination the shape barely matters for any profile whose plateau
outlasts the interval: 42 to 46 years in the median herd. Under
revaccination every 24 months the tail is everything: a step to 18 months
gives 97 years, a plateau to 12 months gone by 24 gives 85, the same
plateau with decline stretched to 36 months gives 51, and a plateau to 18
months gone by 36 gives 44. A short profile, rising over 2 months, holding
to 9 and gone by 18, needs annual or six-monthly boosting and then costs
about 20 years (62 to 63), and never eliminates with a 24-month interval.
What the schedule has to match is the protection-years per dose delivered
inside the interval, the area under the trial's hazard-ratio curve, not a
single duration; `profile_from_hazard` in `R/profile.R` converts such a
curve into the model's input.

The credible intervals come from 40 joint posterior draws of direct and
indirect efficacy from the DST1 fit, with the timing of each profile held
fixed. Their lower ends sit 5 to 7 years below the medians. Their upper
ends are "never" in every cell, because the lower tail of the efficacy
posterior (direct 34%, indirect 46%) puts the reproduction number under
vaccination at about 1 even with lifelong protection; the share of draws
that never eliminate is 2 to 8% for the longer profiles under annual or
six-monthly boosting and 12 to 20% for the short profile and the 24-month
interval. The uncertainty in the efficacy estimates therefore matters at
the feasibility boundary, not in the middle of the grid, and timing
uncertainty is not represented at all until the trial supplies it.

Non-responders (`run_nonresponders.R`). If a share of animals never
responds to any dose, elimination stays possible only while the rest of the
herd can carry the threshold alone: the ceiling is 1 minus the required
protected fraction, which is 44% at the 25th-percentile herd, 28% at the
median, 15% at the 75th percentile and zero at the 90th (fig 19a). Under
annual revaccination with 24-month immunity in the median herd, 10% and
20% persistent non-responders give 62 and 137 years, and 30% never
eliminates; the same shares as random per-dose failure, re-drawn at every
campaign, give 47, 55 and 72 years, because a later dose catches most
animals a previous one missed (fig 19b).

Infected purchases (`run_import.R`). Herds in the main analyses are closed.
With infected animals bought in at a constant rate, and herd size held
constant, prevalence stops falling at a floor set by the import rate: one
infected purchase per 100 head every four years holds the median herd at
about 2.5%, one a year at about 9%, and five a year near the unvaccinated
level (fig 20). Both the paper's regime and 18-month immunity with annual
revaccination behave the same, because each bought animal is infected for
life and starts a short chain that vaccination shortens but cannot remove.
Time to elimination is therefore a within-herd quantity; across herds
that trade, the outcome is the prevalence a programme can hold, and it is
set by the sourcing of replacements.

The same sweep in a UK-like herd (`run_import_uk.R`), with infected
animals removed by annual testing at 0.7 per year and within-herd R0 of
1.1 once that removal is counted, gives a different scale but the same
shape (fig 21). Vaccination alone takes such a herd below 0.1% in 7 to 9
years, against 41 to 45 in Ethiopia, because the reproduction number
under vaccination is 0.12 and infected animals last about a year rather
than a lifetime. Imports still set the floor: one infected purchase per
100 head every four years holds prevalence near 0.4%, one a year near
1.6%, and five a year near 8%, which is close to the unvaccinated level of
9%. The floor scales with the import rate divided by the removal rate,
so test-and-removal lowers it about fourfold relative to Ethiopia, and
pre-movement testing, which lowers the import rate itself, does the rest.

Stochastic grid (`run_stochastic_grid.R`). Replacing the survey herds with
a conceptual grid of within-herd R0 (1.2 to 8) and herd size (5 to 200),
and the deterministic threshold with the first day a herd has no infected
animal, gives probabilities rather than times (figs 22 and 23; 200
replicate herds per cell, with no vaccination as the first panel so each
programme reads against its control). Three things stand out. Small herds
fade out by chance: with no vaccination at all, a 5-head herd at R0 2 is
free by year 20 in 74% of replicates and a 10-head herd in 34%, so in smallholder
systems the vaccine's work is less to drive infection out of a herd than
to stop it coming back, which makes purchases and neighbours the binding
constraint there. Large herds need the reproduction number under
vaccination below one: at 100 to 200 head and R0 2 to 3, the paper's
regime and 18-month immunity with an annual booster both reach 96 to 100%
free by year 50, while 18-month immunity without a booster reaches 0 to
7%. And the annual booster reproduces the lifelong assumption cell for
cell, with the six-monthly booster adding little, which is the
deterministic diagonal rule restated with chance included. Probabilities
by year 20 are much lower in large herds (1 to 26% at 100 to 200 head and
R0 2 to 3 under the paper's regime), consistent with the 41-year
deterministic timescale.

Archetypes and cost parameters. `data/archetypes.csv` defines four herd
archetypes, Ethiopia and India by smallholder and commercial dairy, with
the fields the stochastic and economic layers need; the Ethiopian rows
carry survey values marked "to confirm", the Indian rows are placeholders.
`data/cost_parameters_template.csv` is a generic parameter table for the
economic layer with blank values, units and notes on what each entry
should contain.

Replacement policy in an age-structured herd (`run_herd.R`). The herd
now has calves, young stock and adults, a fixed adult number, adults
leaving at 0.2 per year, and each leaver replaced at once, either by a
home-bred heifer or by a purchase from a source population of given
prevalence (figs 24 and 25; 200 replicate herds per cell, 20-year burn-in).
With every replacement home-bred and 18-month immunity with annual
boosters, herds of 5 and 20 adults are free of infection at year 50 in
86 to 100% of replicates across R0 1.5 to 5, and a 100-adult herd in 98%
at R0 1.5, 78% at 3 and 30% at 5. Purchasing changes this more than any
vaccine property. With half the replacements bought from a population at
5% prevalence, the 20-adult herd falls to 31 to 46% and the 100-adult herd
to 0 to 3%; at 20% source prevalence no herd of 20 or more adults is free,
and vaccination holds prevalence at 6 to 17% instead of the 60 to 78% seen
without it. The 5-adult herd is the exception, because chance fade-out and
re-introduction alternate, so it is free at year 50 in 24 to 87% of
replicates even when buying from infected sources. Without vaccination,
only 5-adult herds at R0 1.5 with home-bred replacements are usually free
(93%); every other cell is at or near its endemic level.

The archetype reading is that a smallholder who breeds replacements can
expect a vaccinated herd to clear infection within a generation of cows,
and one who buys from the local market cannot, whatever the vaccine does.
A commercial dairy that buys replacements needs the source to be clean or
tested; vaccination then lowers prevalence severalfold but does not free
the herd. These are the two levers that the economic layer should price
against each other.

Alternative views (`run_alt_views.R`, figs 26 to 30). The difference
tiles subtract the status quo from each programme: the vaccine adds 50 to
75 percentage points to the chance of freedom by year 20 in herds of 20 to
100 head at R0 1.2 to 3, and 5 to 20 points in 5-head herds at low R0,
which clear infection by chance anyway. The feasibility frontier draws,
for each programme, the herd size and R0 at which the chance of freedom by
year 20 is 50% or 80%: the status quo frontier sits at 5 to 25 head, the
paper's regime and the boosted 18-month programmes lift it to 100 to 200
head at R0 1.2 and to about 10 head at R0 8, with the annual and
six-monthly boosters tracking the lifelong assumption and the calves-only
18-month programme falling between them and the status quo. The archetype
slices show five representative herds as dots per programme with the
control alongside. The time courses of freedom show that "free by year
20" is one slice through curves that separate at year 10 to 15 in large
herds and reach 100% by year 40 to 45 under any boosted programme. The
prevalence bands under replacement policies show the floor set by
purchases: with half the replacements bought from a 20% source, vaccination
takes a 100-adult herd from 63% to 10% and holds it there.

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
Waning is exponential in the main analyses, with an Erlang alternative in
the revaccination script, and the reduced infectiousness of infected
vaccinated animals is assumed permanent. Campaigns dose every eligible
animal at once; partial campaign coverage is a parameter (`cv`) that was not
varied.
