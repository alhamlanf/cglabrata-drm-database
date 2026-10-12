# =============================================================================
# C. glabrata Drug-Resistance Mutation (DRM) Database — Shiny app
# -----------------------------------------------------------------------------
# A searchable, publication-linked lookup for FKS1 / FKS2 mutations in
# Nakaseomyces glabratus (Candida glabrata) and their reported echinocandin
# resistance phenotype.  Reference-lookup tool, NOT a diagnostic.
#
# Single-file app.R:
#   Search + filter + results table + detail panel, plus a Terms of Use tab.
#
# Data lives in ./data as four CSVs (see PROJECT_BRIEF.md for the schema).
# Every non-obvious line is commented for readers new to Shiny/R.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(dplyr)
  library(readr)
  library(stringr)
  library(DT)
})

# ---- Colour conventions ------------------------------------------------------
# Phenotype colours are the universal CLSI/EUCAST convention (locked in brief).
PHENO_COLORS <- c(R = "#c0392b", I = "#e67e22", S = "#27ae60")

# =============================================================================
# DATA LAYER
# =============================================================================

# ---- load_data(): read the four CSVs once at startup ------------------------
# Column types are pinned so the same file always parses identically and we get
# no "guessed type" warnings (Stage 2 success criterion #4).
load_data <- function(data_dir = "data") {

  mutations <- readr::read_csv(
    file.path(data_dir, "mutations.csv"),
    col_types = cols(
      mutation_id           = col_character(),
      excel_row             = col_integer(),
      gene                  = col_character(),
      hotspot               = col_character(),
      mutation_name         = col_character(),
      is_combination        = col_character(),
      wt_residue            = col_character(),
      position              = col_integer(),   # blank for combos -> NA
      mutant_residue        = col_character(),
      mutation_type         = col_character(),
      variable_phenotype    = col_character(),
      component_mutation_ids = col_character()
    )
  )

  phenotypes <- readr::read_csv(
    file.path(data_dir, "mutation_reference_phenotype.csv"),
    col_types = cols(
      mutation_id = col_character(),
      reference_id = col_character(),  # keep as character so it joins cleanly
      anf_class = col_character(), cas_class = col_character(), mcf_class = col_character(),
      anf_detail = col_character(), cas_detail = col_character(), mcf_detail = col_character(),
      anf_comment = col_character(), cas_comment = col_character(), mcf_comment = col_character()
    )
  )

  references <- readr::read_csv(
    file.path(data_dir, "references.csv"),
    col_types = cols(
      reference_id = col_character(),
      pmid = col_character(), doi = col_character(), pmc = col_character(),
      first_author = col_character(), year = col_integer(),
      journal = col_character(), title = col_character(),
      study_type = col_character(), full_citation = col_character(),
      authors = col_character()
    )
  )

  aliases <- readr::read_csv(
    file.path(data_dir, "mutation_aliases.csv"),
    col_types = cols(
      mutation_id = col_character(),
      alias_label = col_character(),
      notes = col_character()
    )
  )

  # -- Normalise alias mutation_id to the canonical MUT_R### form -------------
  # The alias file uses an older "MUT0049" style; mutations use "MUT_R049".
  # Map by the numeric part so the two tables can join.  Best-effort: any alias
  # that still doesn't match simply won't display (it will never crash the app).
  aliases <- aliases %>%
    mutate(
      mutation_id = normalise_mut_id(mutation_id)
    )

  list(
    mutations  = mutations,
    phenotypes = phenotypes,
    references = references,
    aliases    = aliases
  )
}

# ---- normalise_mut_id(): "MUT0049" / "MUT_R49" -> "MUT_R049" -----------------
normalise_mut_id <- function(x) {
  num <- suppressWarnings(as.integer(str_extract(x, "[0-9]+")))  # pull the digits
  ifelse(is.na(num), x, sprintf("MUT_R%03d", num))               # zero-pad to 3
}

# ---- worst_class(): worst-case R/I/S across a set of class strings ----------
# Class strings can be "R", "S", "R/S", "R/I/S", or blank.  We split on "/" and
# take the most-resistant grade seen anywhere.  Precedence: R > I > S.
# Returns "" when there is no data at all (all blank/NA).
worst_class <- function(class_vec) {
  grades <- class_vec %>%
    na.omit() %>%                 # drop NA
    str_split("/") %>%            # "R/S" -> c("R","S")
    unlist() %>%
    str_trim() %>%
    { .[. != ""] }                # drop empties
  if (length(grades) == 0) return("")   # no data
  if ("R" %in% grades) return("R")
  if ("I" %in% grades) return("I")
  if ("S" %in% grades) return("S")
  ""
}

# ---- build_mutation_summary(): one row per mutation for the table + plot -----
# Adds the aggregated worst-case phenotype per drug and the reference count.
build_mutation_summary <- function(data) {

  # Reference count and aggregated phenotype, computed per mutation_id.
  agg <- data$phenotypes %>%
    group_by(mutation_id) %>%
    summarise(
      n_refs   = n_distinct(reference_id),
      anf_worst = worst_class(anf_class),
      cas_worst = worst_class(cas_class),
      mcf_worst = worst_class(mcf_class),
      .groups = "drop"
    )

  data$mutations %>%
    left_join(agg, by = "mutation_id") %>%
    mutate(
      n_refs    = ifelse(is.na(n_refs), 0L, n_refs),
      anf_worst = ifelse(is.na(anf_worst), "", anf_worst),
      cas_worst = ifelse(is.na(cas_worst), "", cas_worst),
      mcf_worst = ifelse(is.na(mcf_worst), "", mcf_worst),
      # Evidence level from the number of supporting references.
      evidence  = case_when(
        n_refs >= 5 ~ "strong (5+ refs)",
        n_refs >= 2 ~ "moderate (2-4 refs)",
        n_refs == 1 ~ "single ref",
        TRUE        ~ "no phenotype refs"
      )
    )
}

# ---- get_refs_for_mutation(): per-reference phenotype rows for detail panel --
get_refs_for_mutation <- function(data, mutation_id) {
  data$phenotypes %>%
    filter(mutation_id == !!mutation_id) %>%
    left_join(data$references, by = "reference_id") %>%
    arrange(as.integer(reference_id))
}

# ---- get_aliases_for_mutation() ---------------------------------------------
get_aliases_for_mutation <- function(data, mutation_id) {
  data$aliases %>% filter(mutation_id == !!mutation_id)
}

# ---- phenotype_pill(): small coloured HTML badge for a drug's worst grade ----
phenotype_pill <- function(drug_label, grade) {
  color <- if (grade %in% names(PHENO_COLORS)) PHENO_COLORS[[grade]] else "#95a5a6"
  shown <- if (grade == "" || is.na(grade)) "no data" else grade
  tags$span(
    style = paste0(
      "display:inline-block;margin:2px 6px 2px 0;padding:4px 10px;border-radius:12px;",
      "font-weight:600;color:white;background:", color, ";"
    ),
    paste0(drug_label, ": ", shown)
  )
}

# ---- variable-phenotype helpers ---------------------------------------------
# Values of `variable_phenotype` that mean "no between-study variability".
NONVARIABLE <- c("None", "No Data", "No phenotypic test", "")

# is_variable(): does this mutation have between-study grade variability at all?
is_variable <- function(vp) !is.na(vp) & !(vp %in% NONVARIABLE)

# drug_varies(): does the named drug ("ANF"/"CAS"/"MCF") vary between studies?
# Driven by the curated `variable_phenotype` field ("All", "CAS and MCF", …),
# the same source the detail-panel warning banner uses.
drug_varies <- function(vp, drug) {
  vp <- ifelse(is.na(vp), "", vp)
  ifelse(vp %in% NONVARIABLE, FALSE, vp == "All" | grepl(drug, vp, fixed = TRUE))
}

# ---- ncbi_link(): turn a PMID / DOI into a clickable link -------------------
ncbi_link <- function(pmid, doi) {
  if (!is.na(pmid) && pmid != "") {
    return(tags$a(href = paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid),
                  target = "_blank", paste0("PMID ", pmid)))
  }
  if (!is.na(doi) && doi != "") {
    return(tags$a(href = paste0("https://doi.org/", doi),
                  target = "_blank", paste0("doi:", doi)))
  }
  tags$span("—")
}

# =============================================================================
# LOAD DATA AT STARTUP  (runs once, shared by all sessions)
# =============================================================================
APP_DATA    <- load_data("data")
MUT_SUMMARY <- build_mutation_summary(APP_DATA)

# =============================================================================
# USER INTERFACE
# =============================================================================
# Flatly's link/accent colour is teal-green (#18bc9c); recolour it blue so links,
# the active tab, and the table's page numbers read blue instead of green.
BLUE      <- "#0072B2"   # primary accent / link colour (recolours flatly's green)
BLUE_DARK <- "#005b8f"

ui <- page_sidebar(
  title = "C. glabrata Echinocandin Resistance — FKS1/FKS2 Mutation Lookup",
  theme = bs_theme(version = 5, bootswatch = "flatly",
                   primary = BLUE, "link-color" = BLUE,
                   "link-hover-color" = BLUE_DARK),

  # Belt-and-suspenders overrides for the few spots that don't follow link-color
  # (DataTables pagination + selected rows).
  tags$head(tags$style(HTML(sprintf("
    a, .nav-link.active { color: %1$s !important; }
    a:hover { color: %2$s !important; }
    /* DataTables pagination renders as Bootstrap pills; flatly colours these
       teal-green — recolour the pill backgrounds blue. */
    .page-link {
      color: #fff !important; background-color: %1$s !important; border-color: %1$s !important;
    }
    .page-item.active .page-link, .page-link:hover {
      color: #fff !important; background-color: %2$s !important; border-color: %2$s !important;
    }
    .page-item.disabled .page-link {
      color: #f1f1f1 !important; background-color: #9dc3d8 !important; border-color: #9dc3d8 !important;
    }
    table.dataTable tbody tr.selected td {
      background-color: %1$s !important; color: #fff !important;
    }
  ", BLUE, BLUE_DARK)))),

  # ---- Sidebar: search + filters (apply to BOTH table and plot) ----
  sidebar = sidebar(
    width = 300,
    textInput("search", "Search mutation (e.g. S663P)", placeholder = "type a name…"),
    # Choices are data-driven so new values (e.g. "multiple" for combinations)
    # appear automatically.
    selectInput("gene", "Gene",
                choices = c("All", sort(unique(MUT_SUMMARY$gene))), selected = "All"),
    selectInput("hotspot", "Hot spot",
                choices = c("All", sort(unique(MUT_SUMMARY$hotspot))), selected = "All"),
    selectInput("mut_type", "Mutation type",
                choices = c("All", sort(unique(MUT_SUMMARY$mutation_type))),
                selected = "All"),
    checkboxInput("variable_only", "Only variable-phenotype mutations", FALSE),
    hr(),
    helpText(
      tags$small(
        "Reference-lookup tool for research and clinical decision-support; ",
        "not a diagnostic. Phenotypes shown are worst-case (most resistant) ",
        "across all references."
      )
    )
  ),

  # ---- Main area: tabs for the table+detail view and Terms of Use ----
  navset_card_tab(
    id = "main_tabs",

    # Tab 1: search results + detail panel
    nav_panel(
      "Table & detail",
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header(textOutput("result_count", inline = TRUE)),
          DTOutput("results"),
          tags$small(
            class = "text-muted",
            HTML("<span style='color:#e67e22;font-weight:700'>*</span> worst-case grade shown; the resistance grade differs between studies for that drug — open the mutation for the per-study breakdown.")
          )
        ),
        card(
          card_header("Mutation detail"),
          uiOutput("detail")
        )
      )
    ),

    # Tab 2: Terms of Use (purpose, disclaimer, citation, contact)
    nav_panel(
      "Terms of Use",
      card(
        card_body(
          div(
            style = "max-width: 820px; margin: 0 auto; line-height: 1.6;",
            tags$h4("Terms of Use"),
            tags$p(class = "text-muted",
                   "C. glabrata Echinocandin Resistance — FKS1/FKS2 Mutation Database"),

            tags$h5("Purpose"),
            tags$p(
              "The ",
              tags$em("C. glabrata Echinocandin Resistance — FKS1/FKS2 Mutation"),
              " database is a curated, openly accessible reference resource that catalogues ",
              tags$em("FKS1"), " and ", tags$em("FKS2"),
              " mutations reported in ", tags$em("Candida glabrata"), " (",
              tags$em("Nakaseomyces glabratus"),
              ") and their association with echinocandin susceptibility (anidulafungin, ",
              "caspofungin, and micafungin). It consolidates mutation–phenotype data from ",
              "the peer-reviewed literature into a single, nomenclature-harmonized resource to ",
              "support research, surveillance, and reference-level interpretation of ",
              tags$em("FKS"), " genotypes."
            ),

            tags$h5("Intended use"),
            tags$p(
              tags$strong("This database is not a diagnostic device."),
              " It is provided as a research and educational reference and a decision-support ",
              "aid for interpreting ", tags$em("FKS"),
              " sequence data, and is not a substitute for laboratory antifungal ",
              "susceptibility testing, current CLSI or EUCAST clinical breakpoints, ",
              "institutional or society treatment guidelines, or the judgment of a qualified ",
              "clinician or clinical microbiologist. Clinical decisions must not be based ",
              "solely on this resource."
            ),

            tags$h5("Data sources and limitations"),
            tags$p(
              "Entries are derived from published studies and reflect the mutations, ",
              "susceptibility results, and interpretations ",
              tags$em("as reported by the original authors"),
              ". Resistance categories (susceptible / intermediate / resistant) depend on the ",
              "testing method and breakpoint version used in each source study and may differ ",
              "from current breakpoints. Where studies disagree, the discordance is displayed ",
              "rather than resolved. The database may contain errors, omissions, or entries ",
              "superseded by newer data, and is provided ",
              tags$strong("“as is,” without warranty"),
              " of accuracy or completeness. Users are responsible for verifying information ",
              "against the primary literature before relying on it."
            ),

            tags$h5("Nomenclature"),
            tags$p(
              "Mutations are standardized to the ", tags$em("C. glabrata"),
              " reference strain ", tags$strong("ATCC 2001 / CBS138"),
              " (assembly GCF_010111755.1). Historical or alternative designations ",
              "(e.g., mutations originally reported using ", tags$em("C. albicans"),
              " numbering) are retained as searchable aliases linked to the corrected entry."
            ),

            tags$h5("Citation"),
            tags$p("If you use this database in a publication, presentation, or report, ",
                   "please cite:"),
            tags$blockquote(
              class = "blockquote",
              tags$p(
                style = "font-size: 0.95rem;",
                "Alhamlan F, et al. ",
                tags$em("C. glabrata Echinocandin Resistance — FKS1/FKS2 Mutation Database"),
                ". Manuscript in preparation; citation to be updated upon publication."
              )
            ),

            tags$h5("Comments and suggestions"),
            tags$p(
              "This resource is maintained as a research project and is updated periodically. ",
              "We welcome corrections, suggestions, and reports of newly published mutations. ",
              "Please direct any comments to the Principal Investigator:"
            ),
            tags$p(
              tags$strong("Fatimah Alhamlan, PhD"), tags$br(),
              "Dr. AlHamlan’s Lab – Research and Innovation", tags$br(),
              "King Faisal Specialist Hospital & Research Centre", tags$br(),
              tags$a(href = "mailto:falhamlan@kfshrc.edu.sa", "falhamlan@kfshrc.edu.sa")
            ),

            tags$hr(),
            tags$p(class = "text-muted", style = "font-size: 0.85rem;",
                   "Database version 1.0 · Last updated October 2026 · ",
                   "Content is subject to revision.")
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # ---- selected_id: the currently-selected mutation, shared across views ----
  selected_id <- reactiveVal(NULL)

  # ---- filtered(): mutation summary after applying sidebar filters ----------
  filtered <- reactive({
    df <- MUT_SUMMARY

    # Free-text search on the mutation name (case-insensitive substring).
    if (!is.null(input$search) && nzchar(input$search)) {
      df <- df %>% filter(str_detect(str_to_lower(mutation_name),
                                     fixed(str_to_lower(input$search))))
    }
    if (input$gene != "All")     df <- df %>% filter(gene == input$gene)
    if (input$hotspot != "All")  df <- df %>% filter(hotspot == input$hotspot)
    if (input$mut_type != "All") df <- df %>% filter(mutation_type == input$mut_type)
    if (isTRUE(input$variable_only)) {
      df <- df %>% filter(is_variable(variable_phenotype))
    }
    df
  })

  # ---- Result count header --------------------------------------------------
  output$result_count <- renderText({
    paste0(nrow(filtered()), " mutation(s)")
  })

  # ---- Results table --------------------------------------------------------
  # Cells show the worst-case grade, with a "*" appended when that drug's grade
  # differs between studies (so "R" and "varies-but-worst-is-R" are not conflated).
  mark_grade <- function(grade, varies) {
    g <- ifelse(is.na(grade) | grade == "", "", grade)
    ifelse(varies & g != "",
           paste0(g, "<sup style='color:#e67e22;font-weight:700'>*</sup>"), g)
  }
  output$results <- renderDT({
    filtered() %>%
      transmute(
        Mutation = mutation_name,
        Gene = gene,
        `Hot spot` = hotspot,
        Type = mutation_type,
        ANF = mark_grade(anf_worst, drug_varies(variable_phenotype, "ANF")),
        CAS = mark_grade(cas_worst, drug_varies(variable_phenotype, "CAS")),
        MCF = mark_grade(mcf_worst, drug_varies(variable_phenotype, "MCF")),
        `n refs` = n_refs
      ) %>%
      datatable(
        selection = "single",
        rownames = FALSE,
        escape = FALSE,   # render the "*" variability marker
        options = list(pageLength = 15, order = list(list(7, "desc")))
      )
  }, server = FALSE)

  # ---- Clicking a table row selects that mutation ---------------------------
  observeEvent(input$results_rows_selected, {
    row <- input$results_rows_selected
    if (length(row)) selected_id(filtered()$mutation_id[row])
  })

  # ---- Detail panel ---------------------------------------------------------
  output$detail <- renderUI({
    id <- selected_id()
    if (is.null(id)) {
      return(tags$em("Select a mutation from the table to see its details."))
    }

    mut <- MUT_SUMMARY %>% filter(mutation_id == id)
    if (nrow(mut) == 0) return(tags$em("Mutation not found."))
    mut <- mut[1, ]

    refs    <- get_refs_for_mutation(APP_DATA, id)
    aliases <- get_aliases_for_mutation(APP_DATA, id)

    # -- Variable-phenotype warning banner (Stage 1 feature) ------------------
    banner <- NULL
    if (is_variable(mut$variable_phenotype)) {
      drugs <- if (mut$variable_phenotype == "All")
        "all three echinocandins" else mut$variable_phenotype
      banner <- div(
        class = "alert alert-warning", role = "alert",
        tags$strong("Variable phenotype: "),
        paste0("resistance grade differs between studies for ", drugs, ".")
      )
    }

    # -- Legacy aliases -------------------------------------------------------
    alias_block <- NULL
    if (nrow(aliases) > 0) {
      alias_block <- div(
        tags$h6("Legacy aliases"),
        tags$ul(lapply(seq_len(nrow(aliases)), function(i) {
          tags$li(tags$strong(aliases$alias_label[i]), " — ",
                  tags$small(aliases$notes[i]))
        }))
      )
    }

    # -- Per-reference phenotype table ----------------------------------------
    ref_rows <- lapply(seq_len(nrow(refs)), function(i) {
      r <- refs[i, ]
      cite <- if (!is.na(r$first_author))
        paste0(r$first_author, " ", r$year, ", ", r$journal) else r$reference_id
      tags$tr(
        tags$td(cite, tags$br(), tags$small(ncbi_link(r$pmid, r$doi))),
        tags$td(style = "text-align:center", r$anf_class %||% ""),
        tags$td(style = "text-align:center", r$cas_class %||% ""),
        tags$td(style = "text-align:center", r$mcf_class %||% "")
      )
    })

    ref_table <- if (nrow(refs) > 0) {
      tags$table(
        class = "table table-sm table-striped",
        tags$thead(tags$tr(
          tags$th("Reference"), tags$th("ANF"), tags$th("CAS"), tags$th("MCF")
        )),
        tags$tbody(ref_rows)
      )
    } else {
      tags$em("No per-reference phenotype records.")
    }

    # -- Assemble the panel ---------------------------------------------------
    tagList(
      tags$h4(mut$mutation_name),
      tags$p(
        tags$span(class = "badge bg-secondary", mut$gene), " ",
        tags$span(class = "badge bg-info text-dark", mut$hotspot), " ",
        tags$span(class = "badge bg-light text-dark", mut$mutation_type)
      ),
      banner,
      tags$h6("Aggregated phenotype (worst case across all references)"),
      div(
        phenotype_pill("Anidulafungin", mut$anf_worst),
        phenotype_pill("Caspofungin",   mut$cas_worst),
        phenotype_pill("Micafungin",    mut$mcf_worst)
      ),
      alias_block,
      tags$h6(paste0("Per-reference phenotype (", nrow(refs), " record(s))")),
      ref_table
    )
  })
}

# `%||%` — small null/NA-coalescing helper used in the reference table.
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# =============================================================================
shinyApp(ui, server)
