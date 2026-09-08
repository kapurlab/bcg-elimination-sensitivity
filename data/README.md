# Input data

Three files copied or derived from the authors' BCGCrossover repository
(https://github.com/MonkeyMyshkin/BCGCrossover, CC BY-SA 4.0). Nothing here
was collected or estimated for this analysis.

| File | Origin | Content | Used for |
| --- | --- | --- | --- |
| `R0Estimates.csv` | `R0Estimates/R0Estimates.csv` | 4000 posterior draws of within-herd R0 for 57 Ethiopian dairy herds, columns `R0.1` to `R0.57` | posterior median per herd defines the herd population; the median across herds (2.8) is the "median herd" |
| `herdsize.csv` | derived from `Data/TestAnon.csv` (animals tested per herd) | herd id and number of animals for the same 57 herds | herd sizes in the population runs |
| `dst1_post.csv` | `Posterior/dst1_post.csv` | posterior draws of transmission rate, direct efficacy `eff_S` and indirect efficacy `eff_I` from the DST1 chain binomial fit; space-separated | the base-case efficacies 0.58 and 0.74 are their medians |

Annual mortality (0.273) is taken from the authors' `Demography/MortalityExp.csv`
and set in `R/model.R`.
