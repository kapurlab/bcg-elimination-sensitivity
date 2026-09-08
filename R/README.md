# R scripts

Run every script from the project root (open the `.Rproj` in RStudio or use
`Rscript` from that directory). `run_all.R` at the root runs them in order.

Model files, sourced by the run scripts:

| File | Model |
| --- | --- |
| `model.R` | Deterministic S, I, V, IV herd model of Fromsa et al. (2024) with exponential waning and calfhood coverage added; R_v, time to elimination, palette |
| `revaccination.R` | Erlang-staged protection (near-fixed duration) and annual whole-herd campaign pulses; SimInf stochastic version |
| `boost.R` | Birth dose and booster with separate durations; campaigns at any interval and coverage |
| `onset.R` | Pre-immunity window after birth with a calf exposure multiplier |
| `root.R` | Working-directory guard, creates `output/` |
| `setup.R` | Installs missing packages |

Analysis scripts, each writing to `output/figs` and `output/tables`:

| Script | What it produces | Time |
| --- | --- | --- |
| `run_oat.R` | One-at-a-time sweeps, elasticities, elimination thresholds, heatmaps (figs 1 to 4) | 1 min |
| `run_global.R` | Latin hypercube, PRCC, Sobol indices (figs 5 to 7) | 15 min |
| `run_stochastic.R` | SimInf check of the sweeps using the authors' transition list (fig 8) | 5 min |
| `run_revaccination.R` | Duration of protection under annual campaigns (figs 9 to 11) | 20 min |
| `run_onset.R` | Speed of onset and calf exposure (figs 12 and 13) | 3 min |
| `run_memo_figure.R` | Scenario ladder, one lever changed per row (fig 14) | 2 min |
| `run_tiles.R` | Booster schedule by duration tiles, and campaign coverage (figs 15 and 16) | 3 min |
