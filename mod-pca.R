## Perform PCA on the selected sample IDs
##
##
## ----------------------------------------------------------------------------

# Set up the input parameters for PCA
pcaUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      awesomeRadio(
        NS(id, "features"),
        label = "Select features",
        choices = c("Genes" = "gene", "Transposable Elements" = "TE", 
                    "Both" = "both"),
        selected = "gene",
        inline = TRUE
      ),
      awesomeRadio(
        NS(id, "dataset"),
        label = "Select Data",
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
        choices = c("FastAuto" = "fast", "Irlba" = "irlba", 
                    "Random" = "random", "Exact" = "exact"),
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
        checkboxInput(NS(id, "loadings"), label = "Show Loadings", value = FALSE),
        checkboxInput(NS(id, "legend"), label = "Show Legend", value = FALSE),
        status = "danger",
        icon = icon("gear")
      ),
      plotOutput(NS(id, "biplot"), brush = NS(id, "plot_brush")),
      DT::dataTableOutput(NS(id, "table"))
    )
  )
}

# PCA server
pcaServer <- function(id, se, keep) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$run, {
      keep_rows <- switch(
        input$features,
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
      showModal(modalDialog("Performing PCA. Please Wait..."))
      pcadata <- PCAtools::pca(
        M,
        metadata = df,
        center = input$center,
        scale = input$scale,
        rank = min(input$rank, ncol(M)),
        removeVar = input$removeVar,
        BSPARAM = algo
      )
      removeModal()

      output$biplot <- renderPlot({
        biplot(
          pcadata,
          colby = input$col,
          x = input$x,
          y = input$y,
          showLoadings = input$loadings,
          legendPosition = if (isTRUE(input$legend)) "bottom" else "none",
          lab = NULL,
          hline = 0,
          vline = 0,
          hlineType = 2,
          vlineType = 2,
          sizeLoadingsNames = 4,
          pointSize = 4
        )
      })

      # PCs must be added to metadata for brushed points to be selected
      df <- cbind(pcadata$rotated, pcadata$metadata)
      
      # Do not plot all PC data - only selected columns and rename in final table
      keep_cols <- c(
        "experiment", "contrast", "tissue", "cell_line", "disease",
        "epigenetic_class"
      )
      new_names <- c(
        "BioProject", "Contrast", "Tissue", "Cell Line", "Disease",
        "Epigenetic Class"
      )
      output$table <- DT::renderDataTable({
        pcdf <- brushedPoints(
          df,
          input$plot_brush,
          xvar = input$x,
          yvar = input$y
        )
        pcdf <- pcdf[, keep_cols]
        colnames(pcdf) <- new_names
        rownames(pcdf) <- NULL
        pcdf
      })
    })
  })
}
