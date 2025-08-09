## Perform PCA on the selected sample IDs
##
## 
## ----------------------------------------------------------------------------

# Function for plotting PCA results with plotly
plotBiplot <- function(obj, x, y, col) {
  # Extract components from the PCA object and bind into a single data.frame
  data <- obj$rotated
  metadata <- obj$metadata
  var_x <- round(obj$variance[x], 1)
  var_y <- round(obj$variance[y], 1)
  d <- cbind(metadata, data)

  vline <- function(x = 0, color = "black") {
    list(
      type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
      line = list(color = color, dash = "dot")
    )
  }

  hline <- function(y = 0, color = "black") {
    list(
      type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
      line = list(color = color, dash = "dot")
    )
  }

  plotly::plot_ly(
    data = d,
    customdata = rownames(d),
    x = ~ get(x),
    y = ~ get(y),
    color = ~ get(col),
    size = 12,
    text = ~ paste(
      "Contrast:", contrast, "\n",
      "BioProject ID:", experiment, "\n",
      "Epigenetic Class:", epigenetic_class, "\n",
      "Tissue:", tissue, "\n",
      "Disease:", disease
    ),
    type = "scatter",
    mode = "markers"
  ) |>
    plotly::layout(
      shapes = list(hline(), vline()),
      dragmode = "lasso",
      xaxis = list(title = paste0(x, " [", var_x, "%]")),
      yaxis = list(title = paste0(y, " [", var_y, "%]")),
      showlegend = FALSE
    ) |>
    plotly::event_register("plotly_selected")
}

# Set up the input parameters for PCA
pcaUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      shinyWidgets::awesomeRadio(
        NS(id, "features"),
        label = "Select features",
        choices = c(
          "Genes" = "gene",
          "Transposable Elements" = "TE",
          "Both" = "both"
        ),
        selected = "gene",
        inline = TRUE
      ),
      shinyWidgets::awesomeRadio(
        NS(id, "dataset"),
        label = "Select data",
        choices = c(
          "z-statistic" = "z",
          "logFC" = "lfc",
          "P-value" = "p"
        ),
        selected = "z",
        inline = TRUE
      ),
      shinyWidgets::prettyCheckbox(
        NS(id, "center"),
        label = "Center data",
        value = TRUE,
        icon = icon("check"),
        status = "success",
        animation = "rotate"
      ),
      shinyWidgets::prettyCheckbox(
        NS(id, "scale"),
        label = "Scale data",
        value = TRUE,
        icon = icon("check"),
        status = "success",
        animation = "rotate"
      ),
      shinyWidgets::prettyCheckbox(
        NS(id, "complete"),
        label = "Use complete cases?",
        value = FALSE,
        icon = icon("check"),
        status = "success",
        animation = "rotate"
      ),
      numericInput(
        NS(id, "rank"),
        label = "Number of components",
        value = 10,
        min = 3,
        max = Inf,
        step = 1
      ),
      numericInput(
        NS(id, "removeVar"),
        label = "Remove proportion low variance features",
        value = 0.9,
        min = 0,
        max = 1,
        step = 0.1
      ),
      shinyWidgets::pickerInput(
        NS(id, "algorithm"),
        label = "PCA algorithm",
        choices = c(
          "FastAuto" = "fast",
          "Irlba" = "irlba",
          "Random" = "random",
          "Exact" = "exact"
        ),
        selected = "fast",
        multiple = FALSE
      ),
      shinyWidgets::actionBttn(
        NS(id, "run"),
        label = "Run PCA",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      shinyWidgets::dropdownButton(
        tags$h3("Biplot parameters:"),
        selectInput(
          NS(id, "x"),
          label = "x",
          choices = paste0("PC", 1:30),
          selected = "PC1"
        ),
        selectInput(
          NS(id, "y"),
          label = "y",
          choices = paste0("PC", 1:30),
          selected = "PC2"
        ),
        selectInput(
          NS(id, "col"),
          label = "Color By",
          choices = c(
            "Epigenetic Class" = "epigenetic_class",
            "Tissue" = "tissue",
            "Drug" = "drug",
            "Disease" = "disease",
            "BioProject" = "experiment"
          )
        ),
        checkboxInput(NS(id, "legend"), label = "Show Legend", value = FALSE),
        status = "danger",
        icon = icon("gear")
      ),
      plotly::plotlyOutput(NS(id, "biplot")),
      gt::gt_output(NS(id, "table"))
    )
  )
}

# PCA server
pcaServer <- function(id, se, keep) {
  moduleServer(id, function(input, output, session) {
    # Perform PCA on 'Run'
    pcaobj <- reactive({
      shinyWidgets::show_alert(
        title = "Performing PCA",
        text = "Please Wait...",
        closeOnClickOutside = FALSE,
        btn_labels = NA,
      )

      # Select the relevant data based on selected IDs and features
      keep_rows <- switch(input$features,
        gene = SummarizedExperiment::rowData(se)$feature_type == "Gene",
        TE = SummarizedExperiment::rowData(se)$feature_type == "TE",
        both = rep(TRUE, nrow(se))
      )

      filtered <- se[keep_rows, keep()]
      df <- data.frame(SummarizedExperiment::colData(filtered))
      m <- switch(input$dataset,
        lfc = SummarizedExperiment::assay(filtered, "logFC"),
        p = SummarizedExperiment::assay(filtered, "P.Value"),
        z = SummarizedExperiment::assay(filtered, "z")
      )

      # Use only complete cases if selected
      if (isTRUE(input$complete)) {
        has_missing <- DelayedMatrixStats::rowAnyNAs(SummarizedExperiment::assay(filtered, "P.Value"))
        m <- m[!has_missing, ]
      }

      # Remove zero-variance features to avoid PCA failing
      if (is.na(input$removeVar) || input$removeVar == 0) {
        m <- m[DelayedMatrixStats::rowVars(m, useNames = FALSE, na.rm = TRUE) != 0, ]
      } else {
        v <- DelayedMatrixStats::rowVars(m, useNames = FALSE, na.rm = TRUE)
        o <- order(v, decreasing=TRUE)
        m <- head(m[o, ], n = max(1, ncol(m) * (1 - input$removeVar)))
      }

      # Select the PCA algorithm
      algo <- switch(input$algorithm,
        fast = BiocSingular::FastAutoParam(),
        irlba = BiocSingular::IrlbaParam(),
        random = BiocSingular::RandomParam(),
        exact = BiocSingular::ExactParam()
      )

      # Attempt PCA computation
      result <- tryCatch(
        {
          PCAtools::pca(
            m,
            metadata = df,
            center = input$center,
            scale = input$scale,
            rank = min(c(input$rank, ncol(m), nrow(m))),
            BSPARAM = algo
          )
        },
        error = function(e) {
          print(e)
          return(NULL)
        },
        warning = function(e) {
          print(e)
          return(NULL)
        }
      )

      if (is.null(result)) {
        shinyWidgets::closeSweetAlert()
        validate("Computation failed! Adjust inputs and try again.")
      }

      shinyWidgets::closeSweetAlert()
      return(result)
    }) |> bindEvent(input$run)

    # Display a biplot of the PCA results
    output$biplot <- plotly::renderPlotly({
      b <- plotBiplot(pcaobj(), x = input$x, y = input$y, col = input$col)
      if (isTRUE(input$legend)) {
        b |> plotly::layout(showlegend = TRUE)
      } else {
        b
      }
    })

    output$table <- gt::render_gt({
      d <- plotly::event_data("plotly_selected")
      df <- pcaobj()$metadata
      if (!is.null(d)) {
        df <- df[d$customdata, ]
      }

      df |>
        gt::gt() |>
        gt::cols_hide(columns = c(id, batch, mutation, comment, desc)) |>
        gt::cols_move(c(tissue, disease), c(cell_line)) |>
        gt::cols_label(
          .list = c(
            "id" = "ID",
            "experiment" = "BioProject ID",
            "contrast" = "Contrast",
            "cell_line" = "Cell Line",
            "drug" = "Drug",
            "dose" = "Dose",
            "time_hr" = "Time (hr)",
            "desc" = "Description",
            "epigenetic_class" = "Epigenetic Class",
            "tissue" = "Tissue",
            "disease" = "Disease"
          )
        ) |>
        gt::cols_width(
          desc ~ px(450),
          contrast ~ px(300),
          experiment ~ px(150)
        ) |>
        gt::tab_header(
          title = gt::md("**Selected PC Data**")
        ) |>
        gt::opt_interactive(use_compact_mode = TRUE)
    })
  })
}
