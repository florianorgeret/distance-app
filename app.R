# ==============================================================================
# File: app.R
# Purpose: Interactive Shiny application for Line Transect Distance Sampling
#          (data cleaning, exploratory analysis, detection function modelling
#          via Distance::ds(), AIC model selection, and CSV export).
# ==============================================================================

library(shiny)
library(tidyverse)
library(sf)
library(Distance)

# Always reload the data-preparation module at startup so edits take effect
source("R/data_prep.R")

# Set a clean, consistent ggplot2 theme across all visualisations
theme_set(theme_bw(base_size = 14))

# ------------------------------------------------------------------------------
# Helper 1: Automatically list all example CSV files available in 'data/'
# ------------------------------------------------------------------------------
get_example_choices <- function() {
  files <- if (dir.exists("data")) {
    list.files("data", pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
  } else {
    character(0)
  }

  if (length(files) == 0) {
    return(c("Aucouturier & Menez" = "AUCOUTURIER_MENEZ_CLEAN.csv"))
  }

  labels <- files |>
    str_remove(regex("_CLEAN\\.csv$|\\.csv$", ignore_case = TRUE)) |>
    str_replace_all("_", " & ") |>
    str_to_title()

  setNames(sort(files), labels[order(files)])
}

# ------------------------------------------------------------------------------
# Helper 2: Resolve file path safely (checks 'data/' first, then project root)
# ------------------------------------------------------------------------------
resolve_example_path <- function(filename) {
  path_in_data <- file.path("data", filename)
  if (file.exists(path_in_data)) return(path_in_data)
  if (file.exists(filename)) return(filename)
  return(NULL)
}

# ------------------------------------------------------------------------------
# Helper 3: Read CSVs robustly (handles ',' vs ';' and UTF-8 vs Windows-1252)
# ------------------------------------------------------------------------------
read_student_csv <- function(filepath) {
  raw_lines <- readLines(filepath, warn = FALSE)
  file_encoding <- if (all(validUTF8(raw_lines))) "UTF-8" else "WINDOWS-1252"
  header_line <- if (length(raw_lines) > 0) raw_lines[1] else ""
  delimiter <- if (str_count(header_line, ";") > str_count(header_line, ",")) ";" else ","

  read_delim(
    filepath,
    delim = delimiter,
    locale = locale(encoding = file_encoding),
    show_col_types = FALSE
  )
}

# ------------------------------------------------------------------------------
# Helper 4: Safely extract SE of average detection probability across mrds versions
# ------------------------------------------------------------------------------
extract_pa_se <- function(mod, sum_mod) {
  candidates <- list(
    sum_mod$dht$individuals$average.p.se,
    mod$dht$individuals$average.p.se,
    tryCatch(summary(mod$ddf)$average.p.se, error = function(e) NULL)
  )
  for (val in candidates) {
    if (!is.null(val) && is.numeric(val) && length(val) > 0 && !is.na(val[1])) {
      return(as.numeric(val[1]))
    }
  }
  return(NA_real_)
}


# ==============================================================================
# 1. USER INTERFACE (UI)
# ==============================================================================
ui <- fluidPage(

  titlePanel("Line Transect Distance Sampling: Exploration, Modelling & Export"),

  sidebarLayout(

    # --------------------------------------------------------------------------
    # SIDEBAR: Step-by-step survey parameters and export buttons
    # --------------------------------------------------------------------------
    sidebarPanel(
      width = 3,

      h4("1. Field Data"),
      radioButtons(
        inputId = "data_source",
        label = "Data source:",
        choices = c(
          "Built-in example datasets" = "example",
          "Upload a CSV file" = "upload"
        ),
        selected = "example"
      ),

      conditionalPanel(
        condition = "input.data_source == 'example'",
        selectInput(
          inputId = "example_file",
          label = "Select an example dataset:",
          choices = get_example_choices()
        )
      ),

      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput(
          inputId = "file_csv",
          label = "Choose raw CSV file:",
          accept = c("text/csv", ".csv"),
          buttonLabel = "Browse...",
          placeholder = "No file selected"
        )
      ),

      checkboxInput(
        inputId = "exclude_test",
        label = "Exclude trial transects ('Test')",
        value = TRUE
      ),

      hr(),

      h4("2. Survey Parameters"),

      numericInput(
        inputId = "h_obs",
        label = "Observer eye height h (m):",
        value = 1.65,
        min = 0.2,
        max = 50,
        step = 0.05
      ),
      helpText("Angle-to-distance formula: d = h * tan((90 - theta) * pi / 180)"),

      numericInput(
        inputId = "area_km2",
        label = "Total study area A (km²):",
        value = 0.01,
        min = 0.0001,
        step = 0.005
      ),
      helpText("Reference: 0.01 km² = 1 hectare = 10,000 m²"),

      numericInput(
        inputId = "default_effort",
        label = "Fallback effort if GPS is missing (km):",
        value = 0.05,
        min = 0.001,
        step = 0.01
      ),

      hr(),

      h4("3. Species & Truncation"),

      selectInput(
        inputId = "selected_species",
        label = "Species / object category:",
        choices = c("All species" = "ALL")
      ),
      helpText("In practice, fit a separate detection function for each species."),

      checkboxInput(
        inputId = "use_truncation",
        label = "Apply right-truncation (w)",
        value = FALSE
      ),

      conditionalPanel(
        condition = "input.use_truncation == true",
        numericInput(
          inputId = "trunc_dist",
          label = "Truncation distance w (m):",
          value = 5.0,
          min = 0.1,
          step = 0.2
        )
      ),

      hr(),

      h4("4. Export Results"),
      p(
        style = "font-size: 0.9em; color: #555;",
        "Download processed tables for your report:"
      ),
      downloadButton(
        outputId = "download_clean_data",
        label = "Cleaned data (.csv)",
        style = "width: 100%; margin-bottom: 8px;"
      ),
      downloadButton(
        outputId = "download_estimates",
        label = "Model estimates (.csv)",
        style = "width: 100%;"
      )
    ),

    # --------------------------------------------------------------------------
    # MAIN PANEL: Automated quality banner + 3 analytical tabs
    # --------------------------------------------------------------------------
    mainPanel(
      width = 9,

      uiOutput("data_quality_banner"),

      tabsetPanel(
        id = "main_tabs",

        # ======================================================================
        # TAB 1: DATA EXPLORATION
        # ======================================================================
        tabPanel(
          title = "1. Exploration",
          br(),

          wellPanel(
            h4("Detection & Effort Summary"),
            uiOutput("summary_metrics")
          ),

          wellPanel(
            fluidRow(
              column(
                width = 8,
                h4("Perpendicular Distance Distribution")
              ),
              column(
                width = 4,
                sliderInput(
                  inputId = "hist_bins",
                  label = "Number of bins:",
                  min = 5,
                  max = 30,
                  value = 12,
                  step = 1
                )
              )
            ),
            plotOutput("plot_hist", height = "320px")
          ),

          wellPanel(
            fluidRow(
              column(
                width = 7,
                h4("Effect of Transect or Covariates on Distances")
              ),
              column(
                width = 5,
                selectInput(
                  inputId = "covar_x",
                  label = "Explanatory variable (X-axis):",
                  choices = c(
                    "Transect (Sample.Label)" = "Sample.Label",
                    "Observer (obs_name)" = "obs_name",
                    "Sea state / Beaufort (beaufort_id)" = "beaufort_id",
                    "Visibility (visibility)" = "visibility"
                  )
                )
              )
            ),
            plotOutput("plot_covar", height = "320px")
          ),

          fluidRow(
            column(
              width = 7,
              wellPanel(
                h4("Effort & Detections by Transect"),
                tableOutput("table_transects")
              )
            ),
            column(
              width = 5,
              wellPanel(
                h4("Cleaned Label Audit"),
                p("Verify that no residual spelling mistakes remain:"),
                tableOutput("table_clean_check")
              )
            )
          )
        ),

        # ======================================================================
        # TAB 2: DETECTION FUNCTION MODELLING (Distance::ds) & AIC COMPARISON
        # ======================================================================
        tabPanel(
          title = "2. Detection Modelling ds()",
          br(),

          uiOutput("species_warning_banner"),

          wellPanel(
            h4("Detection Function Parameters g(x)"),
            fluidRow(
              column(
                width = 5,
                radioButtons(
                  inputId = "ds_key",
                  label = "Key function (key):",
                  choices = c(
                    "Half-normal (hn): smooth shoulder & decay" = "hn",
                    "Hazard-rate (hr): flat shoulder & sharp drop" = "hr"
                  ),
                  selected = "hn"
                )
              ),
              column(
                width = 5,
                selectInput(
                  inputId = "ds_covar",
                  label = "Observation covariate (formula):",
                  choices = c(
                    "None (~ 1)" = "none",
                    "Observer (~ obs_name)" = "obs_name",
                    "Sea state (~ beaufort_id)" = "beaufort_id",
                    "Visibility (~ visibility)" = "visibility"
                  ),
                  selected = "none"
                )
              ),
              column(
                width = 2,
                br(),
                p(strong("Adjustment terms:")),
                code("adjustment = NULL")
              )
            )
          ),

          fluidRow(
            column(
              width = 6,
              wellPanel(
                h4("Fitted Detection Function g(x)"),
                plotOutput("plot_ds_model", height = "350px")
              )
            ),
            column(
              width = 6,
              wellPanel(
                h4("Key Summary Estimates"),
                uiOutput("ds_imred_metrics"),
                hr(),
                h5("Abundance (N) and Density (D) Estimates"),
                tableOutput("table_imred_estimates"),
                helpText("Units: Distance in m, Effort in km, Area A in km², Density D in ind./km².")
              )
            )
          ),

          wellPanel(
            h4("Candidate Model Comparison (AIC)"),
            p(
              "All candidate models below are fitted to the ",
              strong("exact same subset of detections"),
              " (same species, truncation width w, and non-missing covariate rows). ",
              "The top-ranked model has Delta AIC = 0; models with Delta AIC < 2 have comparable support (favour the simpler model for parsimony)."
            ),
            tableOutput("table_aic_comparison"),
            helpText("Reminder: never compare AIC values across different truncation distances w or different sample sizes (n).")
          ),

          wellPanel(
            h4("Full R Output: summary(model)"),
            uiOutput("ds_summary_header"),
            verbatimTextOutput("ds_raw_summary")
          )
        ),

        # ======================================================================
        # TAB 3: PREPARED DATASET PREVIEW (Distance::ds FORMAT)
        # ======================================================================
        tabPanel(
          title = "3. Prepared Data",
          br(),
          p(
            "Dataset generated by ", code("prepare_ds_data()"),
            " formatted for ", code("Distance::ds()"), "."
          ),
          verbatimTextOutput("glimpse_ds"),
          hr(),
          tableOutput("table_preview")
        )
      )
    )
  )
)


# ==============================================================================
# 2. SERVER LOGIC
# ==============================================================================
server <- function(input, output, session) {

  # ----------------------------------------------------------------------------
  # A. Load and validate the raw CSV file
  # ----------------------------------------------------------------------------
  raw_data <- reactive({
    if (input$data_source == "example") {
      req(input$example_file)
      csv_path <- resolve_example_path(input$example_file)
      validate(
        need(
          !is.null(csv_path),
          paste0(
            "Example file not found: '", input$example_file,
            "'. Please place CSV files inside the 'data/' folder of your project."
          )
        )
      )
      df <- tryCatch(read_student_csv(csv_path), error = function(e) e)
    } else {
      req(input$file_csv)
      df <- tryCatch(read_student_csv(input$file_csv$datapath), error = function(e) e)
    }

    validate(
      need(
        !inherits(df, "error") && is.data.frame(df) && nrow(df) > 0,
        "Unable to read the CSV file (file is empty, corrupted, or has an unrecognised format)."
      )
    )

    required_cols <- c(
      "transect_name", "observation_type", "degree",
      "gps_latitude", "gps_longitude"
    )
    missing_cols <- setdiff(required_cols, names(df))

    validate(
      need(
        length(missing_cols) == 0,
        paste0(
          "Missing mandatory columns in CSV: ",
          paste(missing_cols, collapse = ", "),
          ". Please upload the raw field export."
        )
      ),
      need(
        any(c("species_common_name", "species_code") %in% names(df)),
        "No species column found ('species_common_name' or 'species_code')."
      ),
      need(
        any(df$observation_type == "USER", na.rm = TRUE),
        "No rows with observation_type == 'USER' found in this dataset."
      )
    )

    return(df)
  })

  # ----------------------------------------------------------------------------
  # B. Clean data, compute perpendicular distances & GPS effort
  # ----------------------------------------------------------------------------
  prepared_data <- reactive({
    req(raw_data())

    validate(
      need(
        !is.na(input$h_obs) && is.numeric(input$h_obs) && input$h_obs > 0 && input$h_obs <= 50,
        "Observer height h must be strictly positive (between 0.2 m and 50 m)."
      ),
      need(
        !is.na(input$area_km2) && is.numeric(input$area_km2) && input$area_km2 > 0,
        "Total study area A must be strictly greater than 0 km²."
      ),
      need(
        !is.na(input$default_effort) && is.numeric(input$default_effort) && input$default_effort > 0,
        "Fallback effort must be strictly greater than 0 km."
      )
    )

    res <- tryCatch(
      prepare_ds_data(
        df = raw_data(),
        h_obs_m = input$h_obs,
        default_effort_km = input$default_effort,
        exclude_test_transects = isTRUE(input$exclude_test),
        region_label = "Study_Area",
        area_km2 = input$area_km2
      ),
      error = function(e) e
    )

    validate(
      need(
        !inherits(res, "error"),
        paste("Error during data preparation:", if (inherits(res, "error")) res$message else "")
      ),
      need(
        nrow(res) > 0 && sum(res$Effort, na.rm = TRUE) > 0,
        "No transects with valid effort (Effort > 0 km) could be constructed."
      )
    )

    return(res)
  })

  # Update species dropdown and suggest an initial 90th-percentile truncation w
  observeEvent(prepared_data(), {
    df <- prepared_data()

    species_list <- df |>
      filter(!is.na(species)) |>
      count(species, sort = TRUE) |>
      pull(species)

    choices_vec <- c(
      "All species (global exploration)" = "ALL",
      setNames(species_list, species_list)
    )

    current_sel <- if (input$selected_species %in% choices_vec) input$selected_species else "ALL"
    updateSelectInput(session, "selected_species", choices = choices_vec, selected = current_sel)

    valid_dist <- df$distance[!is.na(df$distance)]
    if (length(valid_dist) > 0) {
      suggested_w <- round(as.numeric(quantile(valid_dist, 0.90)), 1)
      updateNumericInput(session, "trunc_dist", value = max(0.5, suggested_w))
    }
  })

  # ----------------------------------------------------------------------------
  # C. Automated data-quality diagnostic banner
  # ----------------------------------------------------------------------------
  output$data_quality_banner <- renderUI({
    df <- prepared_data()
    diag <- attr(df, "diagnostics")
    req(diag)

    alerts <- list()

    if (length(diag$excluded_test_names) > 0) {
      alerts <- append(alerts, list(
        tags$li(
          strong("Excluded trial transect(s): "),
          paste0(
            "transect(s) '", paste(diag$excluded_test_names, collapse = ", "),
            "' were removed from total survey effort (untick the sidebar box to include them)."
          )
        )
      ))
    }

    n_rejected <- diag$n_na_degree + diag$n_le_0 + diag$n_gt_90
    if (n_rejected > 0) {
      details <- c()
      if (diag$n_na_degree > 0) details <- c(details, paste0(diag$n_na_degree, " missing angle(s) NA"))
      if (diag$n_le_0 > 0) details <- c(details, paste0(diag$n_le_0, " angle(s) <= 0°"))
      if (diag$n_gt_90 > 0) details <- c(details, paste0(diag$n_gt_90, " angle(s) > 90°"))

      alerts <- append(alerts, list(
        tags$li(
          strong("Invalid angles removed: "),
          paste0(n_rejected, " detection(s) excluded (", paste(details, collapse = ", "), ").")
        )
      ))
    }

    if (diag$n_eq_90 > 0) {
      alerts <- append(alerts, list(
        tags$li(
          strong("90° angle detected: "),
          paste0(
            diag$n_eq_90,
            " observation(s) with angle = 90° converted to perpendicular distance d = 0 m (on the transect line)."
          )
        )
      ))
    }

    if (diag$n_bad_gps > 0) {
      alerts <- append(alerts, list(
        tags$li(
          strong("Malformed GPS points: "),
          paste0(diag$n_bad_gps, " out-of-bounds or unreadable GPS point(s) ignored when computing effort.")
        )
      ))
    }

    if (length(diag$fallback_transects) > 0) {
      alerts <- append(alerts, list(
        tags$li(
          strong("Fallback effort applied: "),
          paste0(
            "insufficient or stationary GPS track on transect(s): ",
            paste(diag$fallback_transects, collapse = ", "),
            ". Fallback value used: ", input$default_effort, " km."
          )
        )
      ))
    }

    if (diag$n_transects == 1) {
      alerts <- append(alerts, list(
        tags$li(
          strong("Single transect detected (k = 1): "),
          "without spatial replicates, encounter-rate variance across transects is zero and confidence intervals for N̂ are underestimated."
        )
      ))
    }

    if (length(alerts) == 0) return(NULL)

    div(
      class = "alert alert-info",
      style = "margin-bottom: 15px;",
      strong("Field Data Quality Audit:"),
      tags$ul(style = "margin-bottom: 0; margin-top: 5px;", alerts)
    )
  })

  # ----------------------------------------------------------------------------
  # D. Filter dataset by selected species while preserving ALL transects
  # ----------------------------------------------------------------------------
  species_data <- reactive({
    df <- prepared_data()

    all_transects <- df |>
      distinct(Region.Label, Area, Sample.Label, Effort, effort_source)

    detections <- df |>
      filter(!is.na(distance))

    if (input$selected_species != "ALL") {
      detections <- detections |>
        filter(species == input$selected_species)
    }

    all_transects |>
      left_join(
        detections |>
          select(Sample.Label, degree, distance, size, species, obs_name, beaufort_id, visibility),
        by = "Sample.Label"
      ) |>
      mutate(
        object = if_else(!is.na(distance), row_number(), NA_integer_),
        obs_name = as.factor(obs_name),
        beaufort_id = as.factor(beaufort_id),
        visibility = as.factor(visibility)
      ) |>
      as.data.frame()
  })

  valid_truncation <- reactive({
    if (!isTRUE(input$use_truncation)) return(NULL)

    validate(
      need(
        !is.na(input$trunc_dist) && is.numeric(input$trunc_dist) && input$trunc_dist > 0,
        "Truncation distance w must be strictly greater than 0 m."
      )
    )
    return(input$trunc_dist)
  })

  # ----------------------------------------------------------------------------
  # E. Tab 1 Outputs: Exploration metrics, summary tables, and ggplot2 charts
  # ----------------------------------------------------------------------------
  output$summary_metrics <- renderUI({
    df <- species_data()
    trunc_w <- valid_truncation()

    n_transects <- n_distinct(df$Sample.Label)
    total_effort <- df |> distinct(Sample.Label, Effort) |> pull(Effort) |> sum(na.rm = TRUE)

    dets <- df |> filter(!is.na(distance))
    n_total <- nrow(dets)

    n_kept <- if (!is.null(trunc_w)) sum(dets$distance <= trunc_w) else n_total
    n_dropped <- n_total - n_kept

    max_dist <- if (n_total > 0) round(max(dets$distance), 2) else 0
    er <- if (total_effort > 0) round(n_kept / total_effort, 1) else 0

    tagList(
      fluidRow(
        column(3, strong("Transects (k): "), n_transects),
        column(3, strong("Total effort (L): "), paste0(round(total_effort, 3), " km")),
        column(3, strong("Retained detections (n): "), paste0(n_kept, " / ", n_total)),
        column(3, strong("Encounter rate (n/L): "), paste0(er, " / km"))
      ),
      br(),
      fluidRow(
        column(6, strong("Max observed distance: "), paste0(max_dist, " m")),
        column(
          6,
          if (!is.null(trunc_w)) {
            span(
              style = "color: #b22222; font-weight: bold;",
              paste0("Right-truncation active at w = ", trunc_w, " m (", n_dropped, " obs. excluded)")
            )
          } else {
            span(style = "color: #555555;", "No truncation applied (w = max)")
          }
        )
      )
    )
  })

  output$table_transects <- renderTable({
    df <- species_data()
    trunc_w <- valid_truncation()

    df |>
      group_by(Transect = Sample.Label, `Effort (km)` = round(Effort, 4), Source = effort_source) |>
      summarise(
        `Detections (n)` = sum(!is.na(distance) & (is.null(trunc_w) | distance <= trunc_w)),
        `Mean dist. (m)` = round(mean(distance[is.null(trunc_w) | distance <= trunc_w], na.rm = TRUE), 2),
        .groups = "drop"
      ) |>
      mutate(`Mean dist. (m)` = if_else(is.nan(`Mean dist. (m)`), NA_real_, `Mean dist. (m)`))
  }, striped = TRUE, hover = TRUE, na = "-")

  output$table_clean_check <- renderTable({
    prepared_data() |>
      filter(!is.na(distance)) |>
      count(Species = species, Observer = obs_name, name = "Detections") |>
      arrange(desc(Detections))
  }, striped = TRUE, hover = TRUE)

  output$plot_hist <- renderPlot({
    dets <- species_data() |> filter(!is.na(distance))
    trunc_w <- valid_truncation()

    validate(
      need(
        nrow(dets) > 0,
        "No valid detections for this species. Select another species in the sidebar."
      )
    )

    p <- ggplot(dets, aes(x = distance))

    if (input$selected_species == "ALL") {
      p <- p +
        geom_histogram(
          aes(fill = species),
          bins = input$hist_bins,
          boundary = 0,
          colour = "grey20"
        ) +
        scale_fill_brewer(palette = "Pastel1", name = "Species")
    } else {
      p <- p +
        geom_histogram(
          bins = input$hist_bins,
          boundary = 0,
          colour = "grey20",
          fill = "steelblue3"
        )
    }

    if (!is.null(trunc_w)) {
      p <- p +
        geom_vline(
          xintercept = trunc_w,
          linetype = "dashed",
          colour = "red3",
          linewidth = 1.1
        ) +
        annotate(
          "text",
          x = trunc_w,
          y = Inf,
          label = paste0(" Truncation w = ", trunc_w, " m"),
          hjust = -0.05,
          vjust = 1.5,
          colour = "red3",
          fontface = "bold"
        )
    }

    sp_title <- if (input$selected_species == "ALL") "All species" else input$selected_species

    p +
      labs(
        title = paste("Perpendicular Distance Histogram:", sp_title),
        subtitle = "Adjust the bin count to inspect the right tail for extreme values",
        x = "Perpendicular distance (m)",
        y = "Number of detections"
      )
  })

  output$plot_covar <- renderPlot({
    dets <- species_data() |> filter(!is.na(distance))
    trunc_w <- valid_truncation()

    validate(
      need(nrow(dets) > 0, "No valid detections available for this selection.")
    )

    var_x <- input$covar_x
    dets <- dets |>
      mutate(x_group = fct_na_value_to_level(as.factor(.data[[var_x]]), level = "Not recorded"))

    p <- ggplot(dets, aes(x = x_group, y = distance, fill = x_group)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.15, height = 0, size = 2.5, alpha = 0.75, colour = "grey15") +
      scale_fill_brewer(palette = "Pastel1") +
      guides(fill = "none")

    if (!is.null(trunc_w)) {
      p <- p +
        geom_hline(
          yintercept = trunc_w,
          linetype = "dashed",
          colour = "red3",
          linewidth = 0.9
        )
    }

    p +
      labs(
        title = paste("Distance Distribution by:", var_x),
        x = var_x,
        y = "Perpendicular distance (m)"
      )
  })

  # ----------------------------------------------------------------------------
  # F. Tab 2 Outputs: Distance::ds() modelling, summary metrics, and AIC table
  # ----------------------------------------------------------------------------
  output$species_warning_banner <- renderUI({
    if (input$selected_species == "ALL") {
      div(
        class = "alert alert-warning",
        strong("Methodological note: "),
        "you are currently fitting a detection function across all species pooled together. Select each species individually in the sidebar to estimate species-specific detectability and abundance."
      )
    }
  })

  fitted_model <- reactive({
    df <- species_data()
    trunc_val <- valid_truncation()

    validate(
      need(
        all(!is.na(df$Effort) & df$Effort > 0) && sum(df$Effort, na.rm = TRUE) > 0,
        "Zero effort detected (Effort == 0). Check the fallback effort in the sidebar."
      )
    )

    dets <- df |> filter(!is.na(distance))
    if (!is.null(trunc_val)) {
      dets <- dets |> filter(distance <= trunc_val)
    }

    validate(
      need(
        nrow(dets) >= 5,
        paste0(
          "Insufficient detections (n = ", nrow(dets),
          " retained; minimum required = 5) to fit a reliable detection function. ",
          "Increase truncation width w or select a more abundant species."
        )
      ),
      need(
        n_distinct(round(dets$distance, 4)) >= 3,
        "Insufficient variation in distances (< 3 distinct values) to fit a detection curve."
      )
    )

    covar <- input$ds_covar
    if (covar != "none") {
      n_levels <- dets |>
        filter(!is.na(.data[[covar]])) |>
        pull(.data[[covar]]) |>
        droplevels() |>
        nlevels()

      validate(
        need(
          n_levels >= 2,
          paste0(
            "Covariate '", covar,
            "' has only one observed level in this subset (cannot estimate a contrast). ",
            "Select 'None (~ 1)' or another covariate."
          )
        )
      )
    }

    form_str <- if (covar == "none") "~ 1" else paste("~", covar)
    conv_units <- convert_units("meter", "kilometer", "square kilometer")

    # Distance::ds() requires an explicit numeric truncation width rather than NULL
    ds_trunc <- if (!is.null(trunc_val)) trunc_val else max(dets$distance, na.rm = TRUE)

    res <- tryCatch(
      suppressMessages(suppressWarnings(
        ds(
          data = df,
          key = input$ds_key,
          adjustment = NULL,
          formula = as.formula(form_str),
          truncation = ds_trunc,
          convert_units = conv_units
        )
      )),
      error = function(e) e
    )

    validate(
      need(
        !inherits(res, "error"),
        paste0(
          "Model fitting failed in ds() (key = '", input$ds_key, "', formula = ", form_str, "). ",
          "Try the Half-normal ('hn') key function or remove the covariate. Technical detail: ",
          if (inherits(res, "error")) res$message else ""
        )
      )
    )

    return(res)
  })

  output$plot_ds_model <- renderPlot({
    mod <- fitted_model()
    sp_label <- if (input$selected_species == "ALL") "All species" else input$selected_species
    form_label <- if (input$ds_covar == "none") "~ 1" else paste("~", input$ds_covar)

    plot(
      mod,
      main = paste0(sp_label, " | key = '", input$ds_key, "' | formula = ", form_label),
      xlab = "Perpendicular distance (m)",
      ylab = "Detection probability g(x)"
    )
  })

  output$ds_imred_metrics <- renderUI({
    mod <- fitted_model()
    sum_mod <- summary(mod)

    aic_val <- round(as.numeric(AIC(mod)$AIC)[1], 2)
    w_m <- round(as.numeric(mod$ddf$meta.data$width)[1], 2)
    pa <- as.numeric(mod$dht$individuals$average.p)[1]
    pa_se <- extract_pa_se(mod, sum_mod)
    esw_m <- w_m * pa

    se_text <- if (!is.na(pa_se)) paste0(" (SE = ", round(pa_se, 3), ")") else ""

    tagList(
      fluidRow(
        column(
          6,
          p(strong("Model AIC: "), span(style = "font-size: 1.15em; color: #1b4f72;", aic_val)),
          p(strong("Half-width (w): "), paste0(w_m, " m"))
        ),
        column(
          6,
          p(
            strong("Mean detection probability (P̂a): "),
            paste0(round(pa, 3), se_text)
          ),
          p(
            strong("ESW (μ̂ = w × P̂a): "),
            span(
              style = "font-size: 1.15em; color: #145a32; font-weight: bold;",
              paste0(round(esw_m, 2), " m")
            )
          )
        )
      )
    )
  })

  imred_estimates_df <- reactive({
    mod <- fitted_model()
    sum_mod <- summary(mod)

    n_df <- mod$dht$individuals$N
    n_row <- if ("Total" %in% n_df$Label) filter(n_df, Label == "Total") else slice_tail(n_df, n = 1)

    d_df <- mod$dht$individuals$D
    d_row <- if ("Total" %in% d_df$Label) filter(d_df, Label == "Total") else slice_tail(d_df, n = 1)

    w_m <- as.numeric(mod$ddf$meta.data$width)[1]
    pa <- as.numeric(mod$dht$individuals$average.p)[1]
    pa_se <- extract_pa_se(mod, sum_mod)

    tibble(
      Species = if (input$selected_species == "ALL") "All_species" else input$selected_species,
      Key_function = input$ds_key,
      Formula = if (input$ds_covar == "none") "~ 1" else paste("~", input$ds_covar),
      Truncation_w_m = round(w_m, 2),
      AIC = round(as.numeric(AIC(mod)$AIC)[1], 2),
      Pa = round(pa, 3),
      Pa_SE = round(pa_se, 3),
      ESW_m = round(w_m * pa, 2),
      Parameter = c(
        "Total abundance (N)",
        "Density (D, ind./km2)",
        "Density (D, ind./ha)"
      ),
      Estimate = c(
        round(as.numeric(n_row$Estimate)[1], 1),
        round(as.numeric(d_row$Estimate)[1], 1),
        round(as.numeric(d_row$Estimate)[1] / 100, 2)
      ),
      SE = c(
        round(as.numeric(n_row$se)[1], 1),
        round(as.numeric(d_row$se)[1], 1),
        round(as.numeric(d_row$se)[1] / 100, 2)
      ),
      CV_pct = c(
        round(as.numeric(n_row$cv)[1] * 100, 1),
        round(as.numeric(d_row$cv)[1] * 100, 1),
        round(as.numeric(d_row$cv)[1] * 100, 1)
      ),
      CI95_lower = c(
        round(as.numeric(n_row$lcl)[1], 1),
        round(as.numeric(d_row$lcl)[1], 1),
        round(as.numeric(d_row$lcl)[1] / 100, 2)
      ),
      CI95_upper = c(
        round(as.numeric(n_row$ucl)[1], 1),
        round(as.numeric(d_row$ucl)[1], 1),
        round(as.numeric(d_row$ucl)[1] / 100, 2)
      )
    )
  })

  output$table_imred_estimates <- renderTable({
    imred_estimates_df() |>
      select(
        Parameter,
        Estimate,
        SE,
        `CV (%)` = CV_pct,
        `95% CI lower` = CI95_lower,
        `95% CI upper` = CI95_upper
      )
  }, striped = TRUE, hover = TRUE)

  output$table_aic_comparison <- renderTable({
    df <- species_data()
    trunc_val <- valid_truncation()
    covar <- input$ds_covar

    # If a covariate is selected, drop detections where that covariate is NA
    # so that ~ 1 and ~ covar are fitted to the exact same sample size (n)
    if (covar != "none") {
      df <- df |>
        filter(is.na(distance) | !is.na(.data[[covar]]))
    }

    dets <- df |> filter(!is.na(distance))
    if (!is.null(trunc_val)) {
      dets <- dets |> filter(distance <= trunc_val)
    }
    req(nrow(dets) >= 5)

    conv_units <- convert_units("meter", "kilometer", "square kilometer")
    ds_trunc <- if (!is.null(trunc_val)) trunc_val else max(dets$distance, na.rm = TRUE)

    candidates <- list(
      list(label = "Half-normal (no covariate)", key = "hn", formula = "~ 1"),
      list(label = "Hazard-rate (no covariate)", key = "hr", formula = "~ 1")
    )

    if (covar != "none") {
      n_levels <- dets |>
        pull(.data[[covar]]) |>
        droplevels() |>
        nlevels()

      if (n_levels >= 2) {
        candidates <- append(
          candidates,
          list(
            list(label = paste0("Half-normal + ", covar), key = "hn", formula = paste("~", covar)),
            list(label = paste0("Hazard-rate + ", covar), key = "hr", formula = paste("~", covar))
          )
        )
      }
    }

    results <- map_dfr(candidates, function(cand) {
      mod_fit <- tryCatch(
        suppressMessages(suppressWarnings(
          ds(
            data = df,
            key = cand$key,
            adjustment = NULL,
            formula = as.formula(cand$formula),
            truncation = ds_trunc,
            convert_units = conv_units
          )
        )),
        error = function(e) NULL
      )

      if (is.null(mod_fit)) {
        tibble(
          Model = cand$label,
          `Key function` = cand$key,
          Formula = cand$formula,
          n = NA_integer_,
          AIC = NA_real_,
          `P̂a` = NA_real_,
          `ESW (m)` = NA_real_,
          `Abundance N̂` = NA_real_
        )
      } else {
        pa <- as.numeric(mod_fit$dht$individuals$average.p)[1]
        w_m <- as.numeric(mod_fit$ddf$meta.data$width)[1]
        n_df <- mod_fit$dht$individuals$N
        n_est <- if ("Total" %in% n_df$Label) filter(n_df, Label == "Total")$Estimate else tail(n_df$Estimate, 1)

        tibble(
          Model = cand$label,
          `Key function` = cand$key,
          Formula = cand$formula,
          n = nrow(mod_fit$ddf$data),
          AIC = as.numeric(AIC(mod_fit)$AIC)[1],
          `P̂a` = round(pa, 3),
          `ESW (m)` = round(w_m * pa, 2),
          `Abundance N̂` = round(as.numeric(n_est)[1], 1)
        )
      }
    })

    min_aic <- suppressWarnings(min(results$AIC, na.rm = TRUE))

    results |>
      mutate(
        `Delta AIC` = if_else(!is.na(AIC) & is.finite(min_aic), round(AIC - min_aic, 2), NA_real_),
        AIC = round(AIC, 2)
      ) |>
      arrange(is.na(AIC), AIC) |>
      select(Model, `Key function`, Formula, n, AIC, `Delta AIC`, `P̂a`, `ESW (m)`, `Abundance N̂`)
  }, striped = TRUE, hover = TRUE, na = "Did not converge")

  # Header explicitly stating which model is displayed in summary(model)
  output$ds_summary_header <- renderUI({
    mod <- fitted_model()
    sp_label <- if (input$selected_species == "ALL") "All species (pooled)" else input$selected_species
    key_label <- if (input$ds_key == "hn") "Half-normal ('hn')" else "Hazard-rate ('hr')"
    form_str <- if (input$ds_covar == "none") "~ 1" else paste("~", input$ds_covar)
    w_m <- round(as.numeric(mod$ddf$meta.data$width)[1], 2)
    trunc_mode <- if (isTRUE(input$use_truncation)) "user-defined" else "max observed distance"
    n_obs <- nrow(mod$ddf$data)

    r_call <- paste0(
      "ds(data, key = \"", input$ds_key,
      "\", adjustment = NULL, formula = ", form_str,
      ", truncation = ", w_m, ")"
    )

    div(
      style = "background-color: #eef5f9; border-left: 4px solid #2b6cb0; padding: 10px 14px; margin-bottom: 10px;",
      p(
        style = "margin-bottom: 4px;",
        strong("Displayed model: "),
        paste0(
          "Species = ", sp_label,
          " | Key = ", key_label,
          " | Formula = ", form_str,
          " | Truncation w = ", w_m, " m (", trunc_mode, ")",
          " | Detections n = ", n_obs
        )
      ),
      p(
        style = "margin-bottom: 0; font-size: 0.9em;",
        strong("Equivalent R call: "),
        code(r_call)
      )
    )
  })

  output$ds_raw_summary <- renderPrint({
    summary(fitted_model())
  })

  # ----------------------------------------------------------------------------
  # G. Tab 3 Outputs: Prepared dataset preview
  # ----------------------------------------------------------------------------
  output$glimpse_ds <- renderPrint({
    glimpse(species_data())
  })

  output$table_preview <- renderTable({
    head(species_data(), 25)
  }, striped = TRUE, hover = TRUE)

  # ----------------------------------------------------------------------------
  # H. CSV Export Handlers
  # ----------------------------------------------------------------------------
  output$download_clean_data <- downloadHandler(
    filename = function() {
      sp_tag <- if (input$selected_species == "ALL") {
        "all_species"
      } else {
        str_to_lower(str_replace_all(input$selected_species, "\\s+", "_"))
      }
      paste0("cleaned_distance_data_", sp_tag, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      export_df <- species_data() |>
        mutate(
          h_obs_m = input$h_obs,
          distance_m = round(distance, 3),
          Effort_km = round(Effort, 5)
        ) |>
        select(
          Region.Label, Area, Sample.Label, Effort_km, effort_source,
          object, species, obs_name, beaufort_id, visibility,
          size, degree, h_obs_m, distance_m
        )
      write_csv(export_df, file)
    }
  )

  output$download_estimates <- downloadHandler(
    filename = function() {
      sp_tag <- if (input$selected_species == "ALL") {
        "all_species"
      } else {
        str_to_lower(str_replace_all(input$selected_species, "\\s+", "_"))
      }
      paste0("distance_sampling_estimates_", sp_tag, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(imred_estimates_df(), file)
    }
  )
}


# ==============================================================================
# 3. LAUNCH APPLICATION
# ==============================================================================
shinyApp(ui = ui, server = server)