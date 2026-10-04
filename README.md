# R code for *Wang* & *Xing* "A sampling-theory framework for richness–uniqueness relationships"

## License and citation

Copyright © 2026 Xiaoning Wang and Dingliang Xing.

Software code is released under the MIT License. If you use this code in academic work, please cite the associated publication:

Wang, X. & Xing, D. (2026). A sampling-theory framework for richness–uniqueness relationships. *bioRxiv* preprint.

Input data are stored in `data`, intermediate results in `process data`, and figures in `figure`. LCBD and neutral-model basic functions are defined in `code/01-lcbd-models.R` and loaded by the analysis scripts.

All R scripts and text data use UTF-8 encoding. On Windows, use `source(..., encoding = "UTF-8")` for scripts and `fileEncoding = "UTF-8"` when reading CSV files.

## Model functions

`01-lcbd-models.R` defines the Coverage-deficit and Sørensen LCBD calculations, the parameterized, quenched and annealed neutral expectations, and the model-fitting functions used by the figure scripts.

## Figure 1

Figure 1 compares the parameterized, quenched and annealed expectations across the parameter space. 

```r
source("code/11-simulate-sad-ensembles.R", encoding = "UTF-8")
source("code/12-calculate-parameter-space.R", encoding = "UTF-8")
source("code/13-figure-01-three-baselines.R", encoding = "UTF-8")
```

## Figure 2

Figure 2 compares quenched expectations with spatially explicit and spatially implicit neutral simulations.

```r
source("code/20-simulate-spatially-explicit-community.R", encoding = "UTF-8")
source("code/21-prepare-data.R", encoding = "UTF-8")
source("code/22-figure-02-simulation-quenched.R", encoding = "UTF-8")
```

## Figure 3

Figure 3 validates the annealed expectations using Monte Carlo results.

```r
source("code/31-simulate-annealed.R", encoding = "UTF-8")
source("code/32-figure-03-annealed-pooled-mean.R", encoding = "UTF-8")
```

## Figure 4 and Figure S1

Figure 4 applies the annealed Coverage-deficit expectation to BCI. The same plotting script also exports the corresponding Sørensen analysis as Figure S1. 
Download `bci.tree8.rdata` from the BCI eighth-census Dryad record (https://doi.org/10.15146/5xcp-0d46) and place it directly in `data`.

```r
source("code/40-prepare-bci-data.R", encoding = "UTF-8")
source("code/41-calculate-bci-expectations.R", encoding = "UTF-8")
source("code/42-figure-04-bci-priority-map.R", encoding = "UTF-8")
```

## Figure 5 and Figure S4

Figure 5 applies the annealed Coverage-deficit expectation to BBS. The same plotting script also exports the corresponding Sørensen analysis as Figure S4. 
The WDPCA polygon folders and the original NABCI BCR shapefile must be downloaded directly into `data` before running the preparation and figure scripts.

```r
source("code/51-prepare-bbs-data.R", encoding = "UTF-8")
source("code/52-figure-05-bbs-priority.R", encoding = "UTF-8")
```

## Supplementary figures and tables

Figure S1 is produced with Figure 4, and Figure S4 is produced with Figure 5. Figures S2 and S3 show fitted expectations for all BCRs.

Script T11 produces the BCR species-representation tables.

```r
source("code/S21-prepare-bbs-full-expectations.R", encoding = "UTF-8")
source("code/S22-figures-S2-S3-bbs-annealed-all-bcr.R", encoding = "UTF-8")
source("code/T11-table-bbs-bcr-species-representation.R", encoding = "UTF-8")
```
