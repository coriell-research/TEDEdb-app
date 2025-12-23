## Perform UMAP on the PCA results
##
##
## ---------------------------------------------------------------------------

# Function for plotting UMAP results with plotly
plotUmap <- function(df, col) {
  vline <- function(x = 0, color = "black") {
    list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper",
      x0 = x,
      x1 = x,
      line = list(color = color, dash = "dot")
    )
  }

  hline <- function(y = 0, color = "black") {
    list(
      type = "line",
      x0 = 0,
      x1 = 1,
      xref = "paper",
      y0 = y,
      y1 = y,
      line = list(color = color, dash = "dot")
    )
  }

  plotly::plot_ly(
    data = df,
    customdata = rownames(df),
    x = ~V1,
    y = ~V2,
    color = ~ get(col),
    size = 12,
    text = ~ paste(
      "Contrast:",
      contrast,
      "\n",
      "BioProject:",
      experiment,
      "\n",
      "Treatment:",
      treatment,
      "\n",
      "Cell Line:",
      cell_line
    ),
    type = "scatter",
    mode = "markers"
  ) |>
    plotly::layout(
      shapes = list(hline(), vline()),
      dragmode = "lasso",
      xaxis = list(title = "UMAP 1"),
      yaxis = list(title = "UMAP 2"),
      showlegend = FALSE
    ) |>
    plotly::event_register("plotly_selected")
}

umapUI <- function(id) {
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
        label = "Use complete cases",
        value = FALSE,
        icon = icon("check"),
        status = "success",
        animation = "rotate"
      ),
      numericInput(
        NS(id, "removeVar"),
        label = "Remove proportion low variance features",
        value = 0.9,
        min = 0,
        max = 1,
        step = 0.1
      ),
      numericInput(
        NS(id, "neighbors"),
        "N neighbors",
        value = 15,
        min = 2,
        max = Inf,
        step = 1
      ),
      numericInput(
        NS(id, "mindist"),
        "Minimum distance",
        value = 0.01,
        min = 0,
        max = 1,
        step = 0.1
      ),
      shinyWidgets::pickerInput(
        NS(id, "metric"),
        "metric",
        choices = c("Euclidean" = "euclidean", "Manhattan" = "manhattan"),
        selected = "euclidean"
      ),
      shinyWidgets::actionBttn(
        NS(id, "run"),
        label = "Run UMAP",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      shinyWidgets::dropdownButton(
        tags$h3("UMAP parameters:"),
        selectInput(
          NS(id, "col"),
          label = "Color By",
          choices = c(
            "Epigenetic Class (fine)" = "epigenetic_class_fine",
            "Treatment" = "treatment",
            "Mechanism of Action" = "mechanism_of_action",
            "Cell Line" = "cell_line",
            "Cell Lineage" = "oncotree_lineage",
            "Primary Disease" = "oncotree_primary_disease",
            "BioProject" = "experiment"
          )
        ),
        checkboxInput(NS(id, "legend"), label = "Show Legend", value = FALSE),
        status = "danger",
        icon = icon("gear")
      ),
      verbatimTextOutput(NS(id, "test")),
      plotly::plotlyOutput(NS(id, "umap")),
      DT::dataTableOutput(NS(id, "table"))
    )
  )
}

umapServer <- function(id, se, keep) {
  moduleServer(id, function(input, output, session) {
    # Perform UMAP with selected parameters and data
    udata <- reactive({
      shinyWidgets::show_alert(
        title = "Performing UMAP",
        text = "Please Wait...",
        closeOnClickOutside = FALSE,
        btn_labels = NA,
      )

      # Select relevant data and features
      keep_rows <- switch(
        input$features,
        gene = SummarizedExperiment::rowData(se)$feature_type == "Gene",
        TE = SummarizedExperiment::rowData(se)$feature_type == "TE",
        both = rep(TRUE, nrow(se))
      )

      filtered <- se[keep_rows, keep()]
      df <- data.frame(SummarizedExperiment::colData(filtered))
      m <- switch(
        input$dataset,
        lfc = SummarizedExperiment::assay(filtered, "logFC"),
        p = SummarizedExperiment::assay(filtered, "P.Value"),
        z = SummarizedExperiment::assay(filtered, "z")
      )

      # Use only complete cases if selected
      if (isTRUE(input$complete)) {
        has_missing <- DelayedMatrixStats::rowAnyNAs(SummarizedExperiment::assay(
          filtered,
          "P.Value"
        ))
        m <- m[!has_missing, ]
      }

      # Remove zero-variance features to avoid PCA failing
      if (is.na(input$removeVar) || input$removeVar == 0) {
        m <- m[
          DelayedMatrixStats::rowVars(m, useNames = FALSE, na.rm = TRUE) != 0,
        ]
      } else {
        v <- DelayedMatrixStats::rowVars(m, useNames = FALSE, na.rm = TRUE)
        o <- order(v, decreasing = TRUE)
        m <- head(m[o, ], n = max(1, ncol(m) * (1 - input$removeVar)))
      }

      # Attempt UMAP calculation
      result <- tryCatch(
        {
          pca_res <- PCAtools::pca(
            m,
            metadata = df,
            center = input$center,
            scale = input$scale,
            rank = min(c(input$rank, ncol(m), nrow(m))),
            BSPARAM = BiocSingular::FastAutoParam()
          )

          uwot::umap(
            pca_res$rotated,
            n_neighbors = input$neighbors,
            metric = input$metric,
            min_dist = input$mindist
          )
        },
        error = function(e) {
          return(NULL)
        },
        warning = function(e) {
          return(NULL)
        }
      )

      if (is.null(result)) {
        shinyWidgets::closeSweetAlert()
        validate("Computation failed! Adjust inputs and try again.")
      }

      shinyWidgets::closeSweetAlert()

      # Add metadata to the UMAP result
      result <- cbind(as.data.frame(result), df)

      return(result)
    }) |>
      bindEvent(input$run)

    # UMAP plot
    output$umap <- plotly::renderPlotly({
      u <- plotUmap(udata(), col = input$col)
      if (isTRUE(input$legend)) {
        u |> plotly::layout(showlegend = TRUE)
      } else {
        u
      }
    })

    output$table <- DT::renderDataTable({
      d <- plotly::event_data("plotly_selected")

      keep_cols <- c(
        "experiment",
        "contrast",
        "treatment",
        "clinical_phase",
        "mechanism_of_action",
        "targets",
        "epigenetic_class",
        "epigenetic_class_fine",
        "stripped_cell_line",
        "oncotree_lineage",
        "oncotree_primary_disease"
      )

      df <- udata()[, keep_cols]

      if (!is.null(d)) {
        df <- df[d$customdata, ]
      }

      DT::datatable(
        df,
        rownames = FALSE,
        colnames = c(
          "BioProject" = "experiment",
          "Contrast" = "contrast",
          "Treatment" = "treatment",
          "Clinical Phase" = "clinical_phase",
          "Mechanism of Action" = "mechanism_of_action",
          "Target(s)" = "targets",
          "Epigenetic Class" = "epigenetic_class",
          "Epigenetic Class (fine)" = "epigenetic_class_fine",
          "Cell Line" = "stripped_cell_line",
          "Cell Lineage" = "oncotree_lineage",
          "Primary Disease" = "oncotree_primary_disease"
        ),
        lazyRender = FALSE,
        style = "bootstrap4"
      )
    })
  })
}
