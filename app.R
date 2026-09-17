# =============================================================================
# C. glabrata Drug-Resistance Mutation (DRM) Database — Shiny app
# -----------------------------------------------------------------------------
# A searchable, publication-linked lookup for FKS1 / FKS2 mutations in
# Nakaseomyces glabratus (Candida glabrata) and their reported echinocandin
# resistance phenotype.  Reference-lookup tool, NOT a diagnostic.
#
# Single-file app.R holding both:
#   Stage 1 — search + filter + results table + detail panel
#   Stage 2 — lollipop / needle plot (FKS1 above, FKS2 below)
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
  library(ggplot2)
  library(plotly)
})

# ---- Constants: hot-spot boundaries (locked in PROJECT_BRIEF) ---------------
# Position ranges are used to shade the lollipop plot.  We shade purely by
# position (not by the `hotspot` column) because a handful of source rows have
# placeholder values in that column — position is the reliable signal.
HOTSPOTS <- list(
  FKS1 = list(HS1 = c(625, 633), HS2 = c(1340, 1347)),
  FKS2 = list(HS1 = c(659, 667), HS2 = c(1374, 1381))
)

# ---- Colour conventions ------------------------------------------------------
# Phenotype colours are the universal CLSI/EUCAST convention (locked in brief).
PHENO_COLORS <- c(R = "#c0392b", I = "#e67e22", S = "#27ae60")

# Mutation-type colours for the lollipop — the Okabe-Ito palette, the accepted
# colour-blind-safe scheme for scientific figures. Truncating events (premature
# stop) take black by genomics convention; the light Okabe-Ito yellow is omitted
# because it is illegible on a white surface. Assignment is fixed (colour follows
# the category, never its rank), so filtering never repaints the survivors.
TYPE_COLORS <- c(
  "Missense"       = "#0072B2",  # blue
  "Combination"    = "#56B4E9",  # sky blue
  "Deletion"       = "#E69F00",  # orange
  "Insertion"      = "#009E73",  # bluish green
  "Premature stop" = "#000000",  # black (truncating)
  "Frameshift"     = "#D55E00",  # vermillion
  "Position only"  = "#CC79A7"   # reddish purple
)

# Chart chrome & ink — recessive surfaces, hairline grid, muted axis text, so the
# data is the only thing allowed to be loud.
CHROME <- list(
  surface  = "#ffffff",  # chart surface (clean white for a print-like look)
  ink      = "#1a1a1a",  # primary text
  ink2     = "#4a4a4a",  # secondary text (axis titles)
  muted    = "#8a8a8a",  # axis ticks / labels
  grid     = "#ececec",  # hairline gridline
  backbone = "#d9d9d9",  # protein backbone track
  domain   = "#595959",  # hot-spot domain box on the backbone
  stem     = "#bcbcbc",  # lollipop stems
  band     = "#000000"   # hot-spot vertical guide (used at very low alpha)
)

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
      # Evidence level drives lollipop border thickness (brief Stage 2).
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

# Backbone length per gene = furthest observed residue (+ small pad), computed
# once so the protein track stays a fixed size when the plot is filtered.
GENE_LEN <- MUT_SUMMARY %>%
  filter(gene %in% c("FKS1", "FKS2"), !is.na(position)) %>%
  group_by(gene) %>%
  summarise(len = max(position) + 25, .groups = "drop")
GENE_LEN <- setNames(GENE_LEN$len, GENE_LEN$gene)

# =============================================================================
# LOLLIPOP PLOT (Stage 2)
# =============================================================================
# A publication-style needle plot: FKS1 on top, FKS2 below. Each gene is drawn
# as a light protein backbone track with its HS1/HS2 hot-spot domains marked;
# mutations rise as lollipops whose colour = mutation type and height = number
# of reporting studies. The reference axis is square-root-compressed so a single
# heavily-reported residue (e.g. S629P, ~40 studies) doesn't flatten everything
# else. Returns a plotly object carrying mutation_id as customdata for clicks.
make_lollipop <- function(df) {

  # -- Only positioned single mutations in FKS1/FKS2 sit on a residue axis. -----
  plot_df <- df %>%
    filter(gene %in% c("FKS1", "FKS2"), !is.na(position)) %>%
    mutate(
      gene   = factor(gene, levels = c("FKS1", "FKS2")),
      n_refs = pmax(n_refs, 1L),          # 1 so single-study hits are visible
      y      = sqrt(n_refs),              # square-root-compressed stem height
      # Head diameter grows gently with evidence (min ~9px so every head is a target).
      head   = 3.2 + 5.0 * sqrt(n_refs / max(n_refs, na.rm = TRUE)),
      tip = paste0(
        "<b>", mutation_name, "</b>",
        "<br>", gene, "  |  ", hotspot,   # plain ASCII separator (no encoding quirks)
        "<br>Type: ", mutation_type,
        "<br>Studies: ", n_refs
      )
    ) %>%
    droplevels()                          # drop an empty gene panel when filtered

  genes_shown <- levels(plot_df$gene)
  max_y <- sqrt(max(plot_df$n_refs))

  # -- Needle-plot geometry: the backbone track sits below the y=0 baseline; ----
  #    lollipops rise above it; hot-spot labels sit just under the backbone.
  bb_top  <- -0.10 * max_y
  bb_bot  <- -0.30 * max_y
  label_y <- -0.44 * max_y

  # Backbone bar, one per shown gene, spanning residue 1 .. protein length.
  backbone <- data.frame(
    gene = factor(names(GENE_LEN), levels = c("FKS1", "FKS2")),
    xmin = 1, xmax = as.numeric(GENE_LEN)
  ) %>% filter(gene %in% genes_shown) %>% droplevels()

  # Hot-spot domains drawn on the backbone.
  bands <- do.call(rbind, lapply(names(HOTSPOTS), function(g) {
    do.call(rbind, lapply(names(HOTSPOTS[[g]]), function(hs) {
      rng <- HOTSPOTS[[g]][[hs]]
      data.frame(gene = g, hs = hs, xmin = rng[1], xmax = rng[2],
                 xmid = mean(rng), stringsAsFactors = FALSE)
    }))
  }))
  bands <- bands %>% filter(gene %in% genes_shown)
  bands$gene <- factor(bands$gene, levels = genes_shown)

  # Direct-label the most-reported residue in each gene (label the extreme only).
  peaks <- plot_df %>% group_by(gene) %>% slice_max(n_refs, n = 1, with_ties = FALSE)

  # Square-root axis: place ticks at clean study counts.
  brks <- c(1, 2, 5, 10, 20, 40, 60)
  brks <- brks[brks <= max(plot_df$n_refs)]

  p <- ggplot(plot_df) +
    # faint vertical guide over each hot spot (neutral, so it fights no hue)
    geom_rect(data = bands,
              aes(xmin = xmin, xmax = xmax, ymin = bb_bot, ymax = Inf),
              fill = CHROME$band, alpha = 0.035, inherit.aes = FALSE) +
    # protein backbone track
    geom_rect(data = backbone,
              aes(xmin = xmin, xmax = xmax, ymin = bb_bot, ymax = bb_top),
              fill = CHROME$backbone, inherit.aes = FALSE) +
    # hot-spot domain boxes on the backbone, labelled just beneath
    geom_rect(data = bands,
              aes(xmin = xmin, xmax = xmax, ymin = bb_bot, ymax = bb_top),
              fill = CHROME$domain, inherit.aes = FALSE) +
    geom_text(data = bands, aes(x = xmid, y = label_y, label = hs),
              inherit.aes = FALSE, colour = CHROME$ink2, size = 3.1,
              fontface = "bold") +
    # lollipop stems
    geom_segment(aes(x = position, xend = position, y = 0, yend = y),
                 colour = CHROME$stem, linewidth = 0.45) +
    # lollipop heads — white surface ring keeps overlaps legible
    geom_point(aes(x = position, y = y, fill = mutation_type, size = head,
                   text = tip, customdata = mutation_id),
               shape = 21, colour = CHROME$surface, stroke = 1.0) +
    # name the single tallest lollipop per gene
    geom_text(data = peaks, aes(x = position, y = y, label = mutation_name),
              inherit.aes = FALSE, vjust = -1.1, size = 3.2, colour = CHROME$ink) +
    scale_fill_manual(values = TYPE_COLORS, name = NULL, drop = TRUE) +
    scale_size_identity() +
    scale_y_continuous(breaks = sqrt(brks), labels = brks,
                       expand = expansion(mult = c(0.04, 0.14))) +
    facet_wrap(~ gene, ncol = 1, scales = "free_x") +
    labs(x = "Residue position", y = "Reporting studies") +
    theme_minimal(base_size = 14) +
    theme(
      text              = element_text(colour = CHROME$ink2),
      axis.text         = element_text(colour = CHROME$muted, size = 11),
      axis.title        = element_text(colour = CHROME$ink2),
      axis.title.x      = element_text(margin = margin(t = 8)),
      axis.ticks        = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = CHROME$grid, linewidth = 0.4),
      panel.spacing     = unit(26, "pt"),
      strip.text        = element_text(face = "bold", size = 14, hjust = 0,
                                       colour = CHROME$ink),
      plot.background   = element_rect(fill = CHROME$surface, colour = NA),
      panel.background  = element_rect(fill = CHROME$surface, colour = NA),
      legend.position   = "bottom"
    )

  ggplotly(p, tooltip = "text", source = "lollipop") %>%
    layout(
      # Centre the horizontal legend under the plot. (Anchoring it at x = 0 with
      # xanchor "center" put half of it off the left edge, and plotly's autoexpand
      # then inflated the left margin — which is what wasted the plot's width.)
      legend = list(orientation = "h", y = -0.10, x = 0.5, xanchor = "center",
                    title = list(text = ""), font = list(size = 12)),
      paper_bgcolor = CHROME$surface, plot_bgcolor = CHROME$surface,
      margin = list(t = 34, r = 20, l = 60, b = 40)
    ) %>%
    config(displayModeBar = FALSE, responsive = TRUE) %>%
    event_register("plotly_click")
}

# =============================================================================
# USER INTERFACE
# =============================================================================
# Flatly's link/accent colour is teal-green (#18bc9c); recolour it blue so links,
# the active tab, and the table's page numbers read blue instead of green.
BLUE      <- "#0072B2"   # matches the Missense blue in the lollipop plot
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

  # ---- Main area: tabs for table+detail and the lollipop plot ----
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

    # Tab 2: lollipop plot
    nav_panel(
      "Lollipop plot",
      card(
        full_screen = TRUE,
        card_header("FKS1 (top) and FKS2 (bottom) — click a lollipop to select"),
        # fillable = FALSE stops bslib from shrinking the plot to the viewport,
        # so the fixed height is honoured and each gene panel gets real room.
        card_body(
          fillable = FALSE,
          plotlyOutput("lollipop", height = "800px")
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

  # ---- Clicking a lollipop selects that mutation ----------------------------
  observeEvent(event_data("plotly_click", source = "lollipop"), {
    ev <- event_data("plotly_click", source = "lollipop")
    if (!is.null(ev$customdata)) selected_id(ev$customdata)
  })

  # ---- Detail panel ---------------------------------------------------------
  output$detail <- renderUI({
    id <- selected_id()
    if (is.null(id)) {
      return(tags$em("Select a mutation from the table or the lollipop plot."))
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

  # ---- Lollipop plot (respects the same filters) ----------------------------
  output$lollipop <- renderPlotly({
    df <- filtered()
    validate(need(nrow(df) > 0, "No positioned mutations match the filters."))
    make_lollipop(df)
  })
  # Keep the plot live even while the Table tab is in front, so a filter change
  # made there is already applied when the user switches to the plot (otherwise
  # the hidden output is suspended and can show a stale figure).
  outputOptions(output, "lollipop", suspendWhenHidden = FALSE)
}

# `%||%` — small null/NA-coalescing helper used in the reference table.
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# =============================================================================
shinyApp(ui, server)
