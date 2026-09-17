# setup.R — one-time dependency installer for the C. glabrata DRM app
# Run once with:  Rscript setup.R
# ---------------------------------------------------------------------------

# Every package the app.R needs, in one place.
required_packages <- c(
  "shiny",    # web app framework
  "bslib",    # Bootstrap 5 theming (bs_theme)
  "dplyr",    # data wrangling + pipe
  "readr",    # fast, type-safe CSV reading
  "stringr",  # string helpers used when parsing class strings
  "tidyr",    # reshaping (unused today but handy for Stage 3)
  "DT",       # interactive results table
  "ggplot2",  # lollipop / needle plot
  "plotly"    # makes the ggplot interactive (hover + click)
)

# Install only what is actually missing so re-runs are cheap.
missing <- required_packages[!(required_packages %in% rownames(installed.packages()))]

if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed.")
}

# Quick confirmation that everything loads.
invisible(lapply(required_packages, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  message("  loaded: ", p)
}))

message("Setup complete. Launch the app with:  shiny::runApp('.')")
