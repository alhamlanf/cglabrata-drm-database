# C. glabrata DRM Database — Shiny app

A searchable, publication-linked lookup for **FKS1 / FKS2** mutations in
*Nakaseomyces glabratus* (*Candida glabrata*) and their reported echinocandin
resistance phenotype (anidulafungin, caspofungin, micafungin).

> Reference-lookup tool for research and clinical decision-support — **not a
> diagnostic**. See `PROJECT_BRIEF.md` for scientific context and the locked-in
> design decisions.

## Quick start

```r
# 1. Install dependencies (one time)
Rscript setup.R

# 2. Launch the app
Rscript -e 'shiny::runApp(".")'
# or, from an R session:
shiny::runApp(".")
```

The app opens in your browser (default `http://127.0.0.1:xxxx`).

## What it does

| Feature | Where |
|---|---|
| Free-text search (e.g. `S663P`) + gene / hot-spot / type / variable-phenotype filters | sidebar |
| Results table with per-drug **worst-case** R/I/S and reference count | **Table & detail** tab |
| Detail panel: phenotype pills, variable-phenotype warning banner, legacy aliases, full per-reference phenotype table with PubMed/DOI links | **Table & detail** tab |
| Lollipop / needle plot (FKS1 above, FKS2 below) with shaded hot-spot bands; hover for tooltip, click to select | **Lollipop plot** tab |

Filters in the sidebar apply to **both** the table and the plot. Selecting a
mutation in either the table or the plot fills the detail panel.

### Phenotype aggregation

The table and detail pills show the **worst-case** grade across every reference
for a mutation, using the precedence **R > I > S** (a class string like `R/S` is
split and the most-resistant grade wins). Blank means no phenotype data. The
per-reference table in the detail panel shows each study's raw grade so
discordance is visible.

## Project layout

```
shiny_app/
├── app.R                       # single-file app: data layer + UI + server + lollipop
├── setup.R                     # installs required R packages
├── build_data_from_excel.py    # regenerates data/*.csv from the master Excel
├── README.md                   # this file
├── PROJECT_BRIEF.md            # spec, data model, roadmap, locked decisions
├── Database for C.glabrata DRM_clean5.xlsx   # master source workbook
├── data/
│   ├── mutations.csv                       # 332 rows, one per mutation
│   ├── mutation_reference_phenotype.csv    # 857 rows, one per mutation × reference
│   ├── references.csv                      # 105 references
│   └── mutation_aliases.csv                # legacy-literature aliases
└── data_backups/               # timestamped CSV backups written on each rebuild
```

## Updating the data from the Excel

The four CSVs are **generated** from the master workbook
(`Database for C.glabrata DRM_clean5.xlsx`). After editing the Excel, rebuild with:

```bash
python3 build_data_from_excel.py     # needs openpyxl:  pip install openpyxl
```

It backs up the current CSVs to `data_backups/`, then rewrites all four. Rules
(locked with the project owner):

- **Mutation IDs follow the Excel row**: `mutation_id = MUT_R{row}`. Inserting or
  reordering rows therefore reassigns IDs — IDs are stable only as long as row
  order is.
- **Per-reference discordance phenotypes are carried forward**, not regenerated.
  The Excel's cell comments describe *which drugs a study tested*, not the R/I/S
  grades, so the hand-curated grades can't be rebuilt from the workbook. For a
  mutation whose name still matches the previous CSV, its curated per-reference
  grades are reused; new or renamed mutations get the row-level **aggregate**
  phenotype (the "ANF/CAS/MCF" columns) applied to every cited reference.
- **Reference metadata** for IDs not already in `references.csv` is pulled from
  `References_list_copy_v2.xlsx`; existing reference rows are kept verbatim.
- **Legacy aliases** are re-keyed to the current IDs by the target molecule named
  in each alias note.

## Requirements

- R ≥ 4.5
- Packages (installed by `setup.R`): shiny, bslib, dplyr, readr, stringr, tidyr,
  DT, ggplot2, plotly

## Notes on the data

- **Hot-spot bands** in the plot are drawn by residue position (not the
  `hotspot` column) so placeholder values in the source data don't affect
  shading. Boundaries are locked in `PROJECT_BRIEF.md`.
- **`hotspot` column quirk:** in clean5, ~121 rows carry the gene name
  (`FKS1`/`FKS2`) in the `hotspot` column instead of `HS1`/`HS2`/`Outside`.
  These rows display fine and plot correctly (shading is by position) but are not
  reachable via the sidebar hot-spot filter. Fixing them is an Excel edit; re-run
  `build_data_from_excel.py` afterwards.
- **Aliases** are re-keyed on each rebuild to the current `MUT_R###` id of the
  molecule named in the alias note. Unmatched aliases are simply not shown.
- Combination mutations have no single residue position and therefore appear in
  the table but not on the lollipop plot.

## Roadmap

Stage 2 (lollipop plot) is implemented here. Stage 3 — modularising `app.R` into
`mod_*` files, a 3D structure viewer, export buttons, a combination-network view,
a references-explorer tab, and deployment to shinyapps.io — is described in
`PROJECT_BRIEF.md`.
