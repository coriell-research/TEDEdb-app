## Perform PCA on the selected sample IDs
##
## PCA result object is returned and passed to the UMAP module
## ----------------------------------------------------------------------------

plotBiplot <- function(d, x, y, col) {
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

  plot_ly(
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
    layout(
      shapes = list(hline(), vline()),
      dragmode = "lasso",
      xaxis = list(title = paste(x)),
      yaxis = list(title = paste(y)),
      showlegend = FALSE
    ) |>
    event_register("plotly_selected")
}

# Set up the input parameters for PCA
pcaUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      awesomeRadio(
        NS(id, "features"),
        label = "Select features",
        choices = c(
          "Genes" = "gene", "Transposable Elements" = "TE",
          "Both" = "both"
        ),
        selected = "gene",
        inline = TRUE
      ),
      awesomeRadio(
        NS(id, "dataset"),
        label = "Select data",
        choices = c(
          "logFC" = "lfc", "FDR" = "fdr", "Rank Stat" = "stat",
          "Avg. logCPM" = "lcpm"
        ),
        selected = "lfc",
        inline = TRUE
      ),
      prettyCheckbox(
        NS(id, "center"),
        label = "Center data",
        value = TRUE,
        icon = icon("check"),
        status = "success",
        animation = "rotate"
      ),
      prettyCheckbox(
        NS(id, "scale"),
        label = "Scale data",
        value = TRUE,
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
        value = 0.2,
        min = 0,
        max = 1,
        step = 0.1
      ),
      pickerInput(
        NS(id, "algorithm"),
        label = "PCA algorithm",
        choices = c(
          "FastAuto" = "fast", "Irlba" = "irlba",
          "Random" = "random", "Exact" = "exact"
        ),
        selected = "fast",
        multiple = FALSE
      ),
      actionBttn(
        NS(id, "run"),
        label = "Run PCA",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      dropdownButton(
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
      plotlyOutput(NS(id, "biplot")),
      dataTableOutput(NS(id, "table"))
    )
  )
}

# PCA server
pcaServer <- function(id, se, keep) {
  moduleServer(id, function(input, output, session) {
    pcaobj <- bindEvent(reactive(input$run), {
        keep_rows <- switch(input$features,
          gene = rowData(se)$feature_type == "Gene",
          TE = rowData(se)$feature_type == "TE",
          both = rep(TRUE, nrow(se))
        )

        filtered <- se[keep_rows, keep()]
        df <- data.frame(colData(filtered))

        # Extract the selected data
        M <- switch(input$dataset,
          lfc = assay(filtered, "lfc"),
          fdr = assay(filtered, "fdr"),
          stat = t(impute(t(assay(filtered, "stat")))),
          lcpm = assay(filtered, "lcpm")
        )

        algo <- switch(input$algorithm,
          fast = FastAutoParam(),
          irlba = IrlbaParam(),
          random = RandomParam(),
          exact = ExactParam()
        )

        # perform PCA
        msg <- showNotification("Performing PCA. Please wait...",
          type = "message",
          duration = NULL, closeButton = FALSE
        )
        result <- PCAtools::pca(
          M,
          metadata = df,
          center = input$center,
          scale = input$scale,
          rank = min(input$rank, ncol(M)),
          removeVar = input$removeVar,
          BSPARAM = algo
        )
        removeNotification(msg)
        result
    })

    # PCs must be added to metadata for brushed points to be selected
    df <- cbind(pcaobj()$rotated, pcaobj()$metadata)

    # Display a biplot of the PCA results
    output$biplot <- renderPlotly({
      b <- plotBiplot(df, x = input$x, y = input$y, col = input$col)
      if (isTRUE(input$legend)) {
        b |> layout(showlegend = TRUE)
      } else {
        b
      }
    })

    # Do not plot all PC data - only selected columns and rename in final table
    keep_cols <- c(
      "experiment", "contrast", "tissue", "cell_line", "disease",
      "epigenetic_class"
    )
    new_names <- c(
      "BioProject", "Contrast", "Tissue", "Cell Line", "Disease",
      "Epigenetic Class"
    )
    output$table <- renderDataTable({
      d <- event_data("plotly_selected")
      if (!is.null(d)) {
        df2 <- df[d$customdata, keep_cols]
        rownames(df2) <- NULL
        df2
      }
    })
    pcaobj
  })
}
