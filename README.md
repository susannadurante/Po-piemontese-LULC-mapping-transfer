# Po-piemontese-LULC-mapping-transfer

Hierarchical Land Use / Land Cover (LULC) mapping of the *Parco naturale del Po piemontese* (NW Italy) from multi-temporal Sentinel-2 imagery and Random Forest, with transfer of a pre-existing classification from a reference sector to an unmapped sector in the absence of ground truth.

## Overview

The workflow trains a hierarchical Random Forest classifier on the Vercelli–Alessandria sector, where a 2019 land use/land cover cartography is available, and applies the trained models unchanged to the Turin sector, which has no reference cartography. Classification is organised on two nested levels: **L1** (7 macrotypes) and **L2** (18 classes). Level 2 is split into independent sub-models (`L2_1`, `L2_2`, `L2_4`, `L2_5`), each trained on the members of one L1 macrotype.

The guiding principle is *transfer the model, not the accuracy*: the five fitted models are never retrained on the extended area. Each extended-area script includes a verification block confirming that reference-sector pixels remain byte-identical after mosaicking. Transfer limits are then quantified through independent photo-interpreted validation and an Area of Applicability (AoA) analysis, which flags the Turin floodplain as spectrally/topographically outside the training domain (altimetric extrapolation artefact).

## Repository structure

```
Po-piemontese-LULC-mapping-transfer/
├── R/          # pipeline scripts (numbered 00–26, execution order)
├── Style/      # QGIS style files (.qml) and legend join table
├── outputs/    # example outputs (tables, figures)
├── .gitignore
├── README.md
└── LICENSE
```

All scripts read a single root variable at the top, either edited in place or
supplied through the `POPARK_DATA` environment variable:

```r
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "\n     Set the path at the top of the script, or the POPARK_DATA variable.")
```

`DIR_BASE` is the folder holding the input data and the `output_spettrale*`
directories. No other absolute paths are used. Scripts are designed to be run from a terminal
(`Rscript.exe` on Windows), not from an interactive IDE, to keep memory usage bounded.

## Pipeline

### Block A — Reference sector (Vercelli–Alessandria): training

| Script | Purpose |
|---|---|
| `00_vector_cleaning.R` | Clean and dissolve reference vectors and park boundaries |
| `01_composite_MAM.R` | Spring seasonal composite (Sentinel-2 L2A, 10 bands incl. Red Edge) |
| `02_composite_JJA.R` | Summer seasonal composite |
| `03_composite_SON.R` | Autumn seasonal composite |
| `04_composite_DJF.R` | Winter seasonal composite |
| `05_legend_rasterization.R` | Rasterize the final L1/L2 legend to the composite grid |
| `06_JM_separability.R` | Jeffries–Matusita separability on the 40-band stack (L1 + L2) |
| `07_predictors_extraction.R` | Extract predictors and zonal statistics over training polygons |
| `08_predictors_selection.R` | Per-model predictor selection |
| `09_train_val_split.R` | Polygon-based train/validation split (avoids spatial leakage) |
| `10_RF_train_validate.R` | Parametric RF training + k-fold validation (see note below) |
| `11_gapfill_predictors.R` | Iterative focal-mean gap-filling of NA predictor pixels |
| `12_classification.R` | Apply models to the reference sector |

### Block B — Whole park: model transfer to the Turin sector

| Script | Purpose |
|---|---|
| `13_composite_MAM_whole_park.R` | Spring composite, extended extent |
| `14_composite_JJA_whole_park.R` | Summer composite, extended extent |
| `15_composite_SON_whole_park.R` | Autumn composite, extended extent |
| `16_composite_DJF_whole_park.R` | Winter composite, extended extent |
| `17_indices_whole_park.R` | Spectral indices |
| `18_GLCM_whole_park.R` | GLCM texture predictors |
| `19_topography_riverdist_whole_park.R` | DEM-derived predictors and distance-to-river |
| `20_riverdist_mosaic_whole_park.R` | Mosaic distance-to-river (reference tif + ARPA Po axis) |
| `21_gapfill_predictors_whole_park.R` | Gap-filling of NA predictor pixels |
| `22_fill_GLCM_whole_park.R` | Fill structural GLCM NAs on uniform surfaces |
| `23_classification_whole_park.R` | Apply the five fitted models unchanged to the whole park (Vercelli + Turin) |

### Block C — Validation and applicability

| Script | Purpose |
|---|---|
| `24_validation_metrics_torino.R` | Confusion matrices and accuracy metrics (independent points) |
| `25_export_stack_AoA.R` | Export the predictor stack for the AoA analysis, parametric |
| `26_area_of_applicability.R` | Area of Applicability (Meyer & Pebesma 2021), parametric |

### Note on the parametric scripts

`10_RF_train_validate.R` trains and validates one model per run, selected through
the `MODELLO` variable at the top. Run it five times, setting `MODELLO` to `L1`,
`L2_1`, `L2_2`, `L2_4`, `L2_5` in turn, to reproduce all five fitted models.

`25_export_stack_AoA.R` and `26_area_of_applicability.R` follow the same
convention, through their own `MODELLO` variable. Run them as a pair, with the same
value: `25` writes `stack_predittori_<MODELLO>_intero_parco.tif`, which `26` then
reads. The AoA reported in the thesis is the one for `L2_4`, the forest model
affected by the altimetric extrapolation.

## Requirements

- R ≥ 4.5.3
- R packages: `terra`, `sf`, `gdalcubes`, `raster`, `glcm`, `ranger`, `dplyr`, `data.table`, `tidyr`, `caret`, `Boruta`, `pROC`
- Additional package for the Area of Applicability scripts (`25`, `26`): `CAST`

## Data sources

- Sentinel-2 L2A tiles 32TLQ / 32TLR / 32TMQ / 32TMR — Copernicus Data Space Ecosystem
- DTM Regione Piemonte 5 m
- Reference land use/land cover cartography (Vercelli–Alessandria sector, 2019)
- Park boundaries — Regione Piemonte (`parchi_wgs84.shp`)
- Po river axis — ARPA Piemonte hydrographic network (REST/WFS)
- AGEA 2018 orthophotos — Geoportale Piemonte (WMS), used for photo-interpreted validation

## Reference

This pipeline follows and extends the methodology of:

Richiardi, C., Siniscalco, C., Garbarino, M., Adamo, M.P. (2025).
*Unravelling decades of habitat dynamics in protected areas: A hierarchical approach applied to the Gran Paradiso National Park (NW Italy).*
Environmental Monitoring and Assessment, 197, 1216.
DOI: [10.1007/s10661-025-14669-0](https://doi.org/10.1007/s10661-025-14669-0)

Reference repository: `chiararik/GPNP-habitat-mapping`.

Additional methodological references: Meyer & Pebesma (2021, Area of Applicability);
Olofsson et al. (2014, area-adjusted accuracy); Richards & Jia (2006, JM/Bhattacharyya).

## License

Released under **CC0 1.0 Universal** (public domain dedication): see `LICENSE`.

You are free to use, adapt and redistribute this code without restriction.
If it contributes to your work, a citation of Richiardi et al. (2025) and a link
back to this repository are kindly appreciated but not required.
