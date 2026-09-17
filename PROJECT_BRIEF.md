# C. glabrata DRM Database — Project Brief

## Project Overview

**Name:** Candida glabrata Drug-Resistance Mutation (DRM) Interactive Database
**Purpose:** A searchable, publication-linked web application for looking up FKS1 / FKS2 mutations in C. glabrata (Nakaseomyces glabratus) and their reported echinocandin resistance phenotype, references, and per-study discordance. Intended for clinical decision-support at the reference-lookup level (not diagnostic), and as a citable research resource.

**Scientific context:** C. glabrata is a major cause of invasive candidiasis. Echinocandins (anidulafungin, caspofungin, micafungin) are first-line therapy because the species is intrinsically less susceptible to azoles. Resistance arises from mutations in the glucan synthase catalytic subunits Fks1 and Fks2, concentrated in two "hot spot" regions (HS1 and HS2) of each protein. Clinicians who sequence isolates need an interpretation resource — that is what this app is.

**Reference genome:** Nakaseomyces glabratus strain ATCC 2001 (= CBS138), NCBI assembly GCF_010111755.1 / ASM1011175v1.
- FKS1 → chromosome G (CP048124.1)
- FKS2 → chromosome K (CP048128.1)

**Confirmed hot-spot boundaries (validated against 4+ primary sources):**
- FKS1 HS1 = 625–633 (F-L-I-L-S-L-R-D-P)
- FKS1 HS2 = 1340–1347 (D-W-V-R-R-Y-T-L)
- FKS2 HS1 = 659–667 (F-L-I-L-S-L-R-D-P)
- FKS2 HS2 = 1374–1381 (D-W-L-R-R-Y-T-L)


## Data Model — 4 CSV Schema (in `data/`)

### `mutations.csv` — 303 rows, one per mutation
| column | notes |
|---|---|
| mutation_id | MUT_R### matching source Excel row number |
| excel_row | numeric row for traceability |
| gene | FKS1 / FKS2 |
| hotspot | HS1 / HS2 / Outside / Unknown |
| mutation_name | canonical name (S663P, F659S/S663A/D666E) |
| is_combination | Y / N |
| wt_residue | single letter (blank for combos) |
| position | integer |
| mutant_residue | letter, "del", "*", "fs", or "insX" |
| mutation_type | Missense / Deletion / Premature stop / Frameshift / Insertion / Combination / Position only |
| variable_phenotype | None / ANF / CAS / MCF / "ANF and MCF" / All / etc. |
| component_mutation_ids | for combos: semicolon-separated child IDs |

### `mutation_reference_phenotype.csv` — 824 rows, one per mutation × reference
| column | notes |
|---|---|
| mutation_id | FK → mutations |
| reference_id | FK → references |
| anf_class / cas_class / mcf_class | R / I / S / R/I / R/S / I/S / R/I/S / blank |
| anf_detail / cas_detail / mcf_detail | parsed body of the discordance comment |
| anf_comment / cas_comment / mcf_comment | raw comment text |

Per-reference values override the row-level aggregate for the ~115 rows where the source Excel had a yellow-highlighted cell with a threaded comment.

### `references.csv` — 99 rows
| column | notes |
|---|---|
| reference_id | primary key |
| pmid, doi, pmc | external identifiers |
| first_author, year, journal, title | |
| study_type | Clinical study / WGS / Surveillance / Review / Experimental evolution / etc. |
| full_citation, authors | full APA-style citation + author list |

### `mutation_aliases.csv` — 4 rows (extend as needed)
| column | notes |
|---|---|
| mutation_id | FK to mutations |
| alias_label | the name used in older or discordant literature |
| notes | why the alias exists (typo, C. albicans numbering, etc.) |


## Current Status

**Stage 1 delivered:**
- Data model finalized, all 4 CSVs populated.
- Shiny app skeleton in one heavily-commented `app.R`.
- Search + filter panel + results table + detail panel with:
  - Aggregated phenotype pills (worst-case R/I/S per drug)
  - Variable-phenotype warning banner
  - Legacy alias display
  - Full per-reference phenotype table

**File tree:**
```
shiny_app/
├── app.R
├── setup.R
├── README.md
├── PROJECT_BRIEF.md   ← this file
└── data/
    ├── mutations.csv
    ├── mutation_reference_phenotype.csv
    ├── references.csv
    └── mutation_aliases.csv
```


## Roadmap

### Stage 2 — Lollipop / needle plot visualization
Two aligned panels (FKS1 above, FKS2 below), horizontal protein bars with HS1/HS2 shaded, each mutation drawn as a lollipop.
- Color = mutation type (Missense=blue, Deletion=red, Insertion=purple, Premature stop=black, Frameshift=orange, Combination=green)
- Border thickness = evidence level (1 ref thin, 2–4 medium, 5+ thick)
- Height = number of references
- Hover tooltip shows mutation name, gene, hot spot, type, n refs
- Click selects mutation and populates detail panel
- Sidebar filters apply to plot as well as table
- Tech: ggplot2 + plotly::ggplotly() is the simplest route

### Stage 3 — Polish, modularize, deploy
- Split `app.R` into `mod_search.R`, `mod_lollipop.R`, `mod_detail.R`, `mod_data.R` using `moduleServer()`
- Add 3D structure viewer (NGLVieweR) with AlphaFold models (Fks1 = UniProt Q6FRK7, Fks2 = Q6FQK0)
- Add export buttons (CSV, PNG, BibTeX)
- Add combination network view (visNetwork)
- Add references-explorer tab
- Add version tag and feedback link
- Deploy to shinyapps.io (free) or Posit Connect (institutional)


## Explicit exclusions (per user decision)

- **No sequence-type (ST3, ST7) in search** — ST is a strain lineage, not a mutation.
- **No clinical prediction panel** — this app is a lookup tool, not a diagnostic.


## Design decisions locked in

- IDs are MUT_R### matching Excel row numbers
- Reference genome is ATCC 2001 / GCF_010111755.1 (same as CBS138 under a different catalog code)
- Hot-spot boundaries are as listed above — do not reopen this debate
- Legacy literature errors documented as aliases, not corrections
- Phenotype colors: R=red, I=orange, S=green (universal convention, matches CLSI/EUCAST reports)
- Framework is R + Shiny — chosen for mycology-community familiarity


## Coding conventions

- Use dplyr and the `%>%` or `|>` pipe
- Use `bslib::bs_theme()` for theming
- Use `tags$*` for HTML from R, not raw strings
- Add `# ----` section headers before major blocks
- Snake_case function names, prefixed by domain (get_refs_for_mutation, phenotype_color)
- Inline comment on any line that isn't self-explanatory to a beginner


## Success criteria for Stage 2 + 3

1. User can type "S663P" in the search → the mutation appears in the table.
2. User clicks the row → detail panel fills with phenotype pills, aliases, per-ref table.
3. User sees a lollipop plot with hot-spot bands, hovering shows tooltip, clicking selects.
4. All 4 CSVs load cleanly at startup with no warnings.
5. The app deploys to shinyapps.io and loads in under 3 seconds cold.
