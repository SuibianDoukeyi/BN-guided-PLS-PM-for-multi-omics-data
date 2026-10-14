# BN-guided PLS-PM for multi-omics data

R code for Bayesian network-guided partial least squares path modeling (BN-guided PLS-PM) used in:

> Lin, K., Liu, Y., et al. An endogenous Fenton-like route for enrofloxacin degradation shares an oxidative node with hormesis in a macroalgal holobiont. (Submitted to *Environment International*)

## Repository contents

| File | Description |
|---|---|
| `BN_guided PLS_PM modelling.R` | Full analysis workflow: Level 2 PLS-PM indicator selection, Level 1 BN structure learning, BN-guided PLS-PM estimation, network metrics, dose–latent regressions, and figure generation |
| `data/L2 plspm 82.csv` | Level 2 indicators (transcriptomic, metabolomic, and microbial features). Columns: `type2` (Level 2 module), `var_id` (indicator), followed by 12 samples (N, W, T, E; n = 3 each) |
| `data/L1 latent to L2.csv` | Mapping between Level 1 modules (`type1`) and Level 2 modules (`type2`) |
| `data/l1 sub l101.csv` | Level 1 indicators measured directly (growth, toxicity, antibiotic residues, biodegradation, oxidative response, photosynthetic and other physiological variables). Columns: `type1`, `type2`, followed by the 12 samples |
| `data/DoseM.csv` | ENR dose for each of the 12 samples |

Treatments: N (0 μM), W (0.008 μM), T (8 μM), and E (30 μM) enrofloxacin, 72 h exposure.

Raw sequencing and metabolomics data are available at the National Genomics Data Center (NGDC) under BioProject PRJCA052533.

## Requirements

- R 4.4.2
- R packages: plspm, plspm2, bnlearn, Rgraphviz, igraph, qgraph, sna, tidygraph, ggraph, ggforce, tidyverse (dplyr, tibble, ggplot2, purrr), patchwork, scales, grid

## Usage

1. Download or clone this repository.
2. Set the working directory to the repository root, e.g. `setwd("path/to/BN-guided-PLS-PM-for-multi-omics-data")`.
3. Run `BN_guided PLS_PM modelling.R` from top to bottom.

The script reads all input files from the `data/` folder.

### Scaling of Level 1 data

Lines 110–138 contain two alternative blocks for combining and scaling Level 1 data:
- Block `###1` (active): scales the directly measured Level 1 indicators, then appends the Level 2-derived module scores.
- Block `###2` (commented out): appends first, then scales.

Results reported in the manuscript were obtained with Block `###1`. Do not run both blocks in the same session, as this scales the data twice.

## Troubleshooting

**Error: object `name` not found (Step 10 or Step 12)**

The column name of the PLS-PM outer model output can differ depending on whether `plspm::plspm()` or `plspm2::plspm()` is active in the current R session (`plspm2` is loaded later in the script and masks `plspm`). If this error occurs, replace `name` with `observed` in Step 10 and Step 12:

```r
l2_selected_ids <- as.character(l2_selected_vars$observed)
...
summarise(vars = list(as.character(observed))) %>%
```


## Contact

Jianjun Cui (macroalgae29cui@126.com)
