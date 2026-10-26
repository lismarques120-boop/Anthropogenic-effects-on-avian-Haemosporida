# Avian Haemosporidian Prevalence and Diversity Analysis

This repository contains the complete, reproducible workflow used to model
prevalence and phylogenetic diversity of avian haemosporidian parasites
(Plasmodium, Parahaemoproteus, and all genera - Plasmodium, Haemoproteus, Leucocytozoon -combined)
as a function of host order (Passeriformes vs. non-Passeriformes), host diversity, climate,
vegetation (NDVI), and human influence (HII), across Neotropical bird communities.

The workflow has two stages:

1. Google Earth Engine ('GEE_scripts' file) - extracts environmental covariates
   (HII, NDVI, WorldClim BIO) within 10 km buffers around each sampling point.
2. R analyses ('Models.R' file) - joins those covariates to the host/parasite dataset
   and fits binomial GLMMs (prevalence) and Gaussian log-link GLMs (diversity).

## Important note on terminology

> The GEE scripts in this repository are written in JavaScript and are
> intended to be pasted into the [Earth Engine Code Editor](https://code.earthengine.google.com/).

## Input asset for GEE: `pontos`

All GEE scripts assume a pre-existing `FeatureCollection` of sampling points:

projects/subtle-seat-379912/assets/Lis/tentativa2

all values contained in `pontos` are already present in `dataset.xlsx` (`community_cluster`, `Latitude`,`Longitude`).

- Geometry: `Point`, CRS `EPSG:4326` (WGS84)
- Required attributes: an identifier column matching
  `community_cluster` in `dataset.xlsx` (plus any metadata you want
  carried through the buffers).
- Created by: uploading a CSV (columns: `id`, `community_cluster`,
  `latitude`, `longitude`) via*Assets ▸ NEW ▸ CSV file in the Code
  Editor, then setting the asset path above.

## Data dictionary - `dataset.xlsx`

One row = one screened bird individual. Columns fall into four groups:

### 1. Host taxonomy and sampling metadata

| Column | Type | Description |
|--------|------|-------------|
| `BirdTree.Order` | text | Avian order (e.g. `Passeriformes`, `Columbiformes`, `Apodiformes`). Used to build `grupo`. |
| `BirdTree.Family` | text | Avian family (e.g. `Tyrannidae`, `Thamnophilidae`). |
| `Host.Latin.Name` | text | Host species, `Genus_species` (underscore-separated). |
| `community_cluster` | integer | Sampling locality / community ID. Primary grouping key used in R aggregation and as the random-effect grouping factor. Matches the GEE `pontos` identifier. |
| `Year` | integer | Year of sampling. Determines which year-specific GEE CSV to join. |
| `Latitude`, `Longitude` | decimal degrees | Individual sampling coordinates (WGS84). |
| `Country` | text | Sampling country (`Brazil`, `Peru`, `Ecuador`). |

### 2. Infection diagnosis and molecular data

| Column | Type | Description |
|--------|------|-------------|
| `Diagnosis` | text | Screening result: `P` = positive (infected by ≥1 haemosporidian genus), otherwise negative. Used to build `inf_bin` in the prevalence models. |
| `Genus.for.all` | text | Parasite genus/genus-combination detected. Codes used by `str_detect` in the R scripts: `PL` = Plasmodium, `PA` = Haemoproteus (Parahaemoproteus), `HA` = Haemoproteus (Haemoproteus), `LE` = Leucocytozoon. Used to build `PL_bin`, `PA_bin`, and the all-genera binary. |
| `GenBank1` … `GenBank4` | text | GenBank accession number for each detected lineage. |
| `Lineage.1.name` … `Lineage.5.name` | text | MalAvi lineage name (e.g. `TARMAJ01`, `LIPVOC01`). Used to compute lineage diversity (q1). |
| `Lineage.1.Cyt.B` … `Lineage.5.Cyt.B` | integer | cytb sequence associated with each lineage. |
| `nome_sequencias1` … `nome_sequencias1` (raw nucleotide strings) | text | Raw cytb sequences for each detected lineage. |

> Note: `Lineage.N.*` columns are stacked per individual - one row can carry up to 5 lineages. Diversity metrics (`q1_*`) summarize these at the community × host-order level.

### 3. Environmental covariates (from GEE)

| Column | Type | Description |
|--------|------|-------------|
| `Latitude_media`, `Longitude_media` | decimal degrees | Community-mean coordinates (averaged across all individuals in `community_cluster`). |
| `NDVI_por_ano` | numeric | Median NDVI within the 10 km buffer for that sampling year. Ranges ~0–1. |
| `HII_por_ano` | numeric | Median Human Influence Index within the 10 km buffer for that year, rescaled to 0–50. |
| `bio01` | numeric (°C) | WorldClim BIO1 - Annual Mean Temperature. |
| `bio04` | numeric | WorldClim BIO4 - Temperature Seasonality (SD × 100). |
| `bio12` | numeric (mm) | WorldClim BIO12 - Annual Precipitation. |

### 4. Diversity indices (`q1_*`)

Diversity is measured with Hill number q = 1 (the exponential of Shannon entropy — "effective number of lineages").
The `_pass` / `_nopass` suffix splits the index by host order, and the unsuffixed columns are the pooled community-level values.

| Column | Type | Description |
|--------|------|-------------|
| `q1_hosts` | numeric | Host (bird) diversity at the community level (q1). |
| `q1_hosts_pass` | numeric | Host diversity restricted to Passeriformes. |
| `q1_hosts_nonpass` | numeric | Host diversity restricted to non-Passeriformes. |
| `q1_all_genera` | numeric | Parasite diversity - all haemosporidian genera combined. |
| `q1_all_genera_pass` | numeric | Same, restricted to Passeriformes. |
| `q1_all_genera_nopass` | numeric | Same, restricted to non-Passeriformes. |
| `q1_plasmodium` | numeric | Parasite diversity - _Plasmodium_ only. |
| `q1_plasmodium_pass` | numeric | Same, restricted to Passeriformes. |
| `q1_plasmodium_nopass` | numeric | Same, restricted to non-Passeriformes. |
| `q1_parahaemoproteus` | numeric | Parasite diversity - _Parahaemoproteus_ only. |
| `q1_parahaemoproteus_pass` | numeric | Same, restricted to Passeriformes. |
| `q1_parahaemoproteus_nopass` | numeric | Same, restricted to non-Passeriformes. |

How the R scripts use these: each script selects the order-appropriate column via `case_when(grupo == "Passeriformes" ~ ..._pass, TRUE ~ ..._nopass)`
to build a single `q1_aves` and `q1_parasite` column for modelling.


## Model overview (R)

| Script | Response | Model | Key covariates |
|--------|----------|-------|----------------|
| `01_prevalence_all_genera.R` | Any parasite (`inf_bin`) | Binomial GLMM, `glmmTMB`, random intercept `(1 | community_cluster)` | `grupo`, NDVI, HII (+ HII²), bio01, bio04, bio12 (+ quadratics), q1_aves, all `grupo:` interactions |
| `02_prevalence_plasmodium.R` | *Plasmodium* (`PL_bin`) | same as above | same |
| `03_prevalence_parahaemoproteus.R` | *Parahaemoproteus* (`PA_bin`) | same as above | same |
| `04_diversity_all_genera.R` | `q1_all_genera` | Gaussian GLM, log link, no random effect | `grupo`, NDVI, HII (+ HII²), q1_aves, `grupo:` interactions |
| `05_diversity_plasmodium.R` | `q1_plasmodium` | same | same |
| `06_diversity_parahaemoproteus.R` | `q1_parahaemoproteus` | same | same |

## Output files (generated at runtime)

Each R script writes one `.xlsx` with the model data and the tidy coefficient table. From `gee/*.js`, each export produces one `.csv` per year (see year-by-year workflow).

| Generated file | Produced by |
|----------------|-------------|
| `HII_principal_escalado_YYYY.csv` | `gee/01_hii.js` |
| `HII_<driver>_escalado_YYYY.csv` (×7 drivers) | `gee/01_hii.js` |
| `NDVI_por_buffer_quadrado_10km_YYYY.csv` | `gee/02_ndvi.js` |
| `WorldClim_BIO_Public_Buffer10km_Celsius.csv` | `gee/03_worldclim_bio.js` |
| `agregated_all_genera.xlsx` | `r/01_prevalence_all_genera.R` |
| `prevalence_all_genera_params.xlsx` | `r/01_prevalence_all_genera.R` |
| `COEFS_model_all_genera.xlsx` | `r/01_prevalence_all_genera.R` |
| `agregated_PL_comunidade_grupo.xlsx` | `r/02_prevalence_plasmodium.R` |
| `3.2-prevalence_PL_params.xlsx` | `r/02_prevalence_plasmodium.R` |
| `COEFS_model_PL.xlsx` | `r/02_prevalence_plasmodium.R` |
| `agregated_PA_comunidade_grupo.xlsx` | `r/03_prevalence_parahaemoproteus.R` |
| `3.2-prevalence_PA_params.xlsx` | `r/03_prevalence_parahaemoproteus.R` |
| `COEFS_model_PA.xlsx` | `r/03_prevalence_parahaemoproteus.R` |
| `model_diversity_all_genera.xlsx` | `r/04_diversity_all_genera.R` |
| `model_diversity_PL.xlsx` | `r/05_diversity_plasmodium.R` |
| `model_diversity_PA.xlsx` | `r/06_diversity_parahaemoproteus.R` |


## Joining GEE outputs to the R dataset

The `community_cluster` column in `dataset.xlsx` is the join key between the GEE CSV exports and the R analysis. If your `pontos` asset uses a
different identifier column (e.g. `id`), rename it in R before joining:

Predictors in `dataset.xlsx` (`NDVI_por_ano`, `HII_por_ano`, `bio01`, `bio04`, `bio12`) are already the joined, year-matched values.

## Reproducibility notes

- R version: ≥ 4.2;
- GEE: JavaScript (V8 runtime), Earth Engine Code Editor;
- File formats: all code is provided as plain text (`.R` /`.txt`) so it can be copied directly;
- Excel locale: `dataset.xlsx` uses comma decimals in some numeric columns (e.g. `Latitude_media` = `-30,20556399`).
  If `readxl` imports these as text, convert with `mutate(across(where(is.character), ~ as.numeric(gsub(",", ".", .x))))` before running the models.
- No separate asset is required: all values contained in `pontos` are already present in `dataset.xlsx` (`community_cluster`, `Latitude`,`Longitude`).
  Users can rebuild an equivalent `FeatureCollection` from those columns and upload it as their own asset.
