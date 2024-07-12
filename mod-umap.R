## Perform UMAP on the PCA results
##
##
## ---------------------------------------------------------------------------

# Function for plotting PCA results with plotly
plotUmap <- function(df, col) {
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
    data = df,
    customdata = rownames(df),
    x = ~UMAP1,
    y = ~UMAP2,
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
      xaxis = list(title = "UMAP 1"),
      yaxis = list(title = "UMAP 2"),
      showlegend = FALSE
    ) |>
    event_register("plotly_selected")
}

umapUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      awesomeRadio(
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
      awesomeRadio(
        NS(id, "dataset"),
        label = "Select data",
        choices = c(
          "z-statistic" = "z",
          "t-statistic" = "t",
          "logFC" = "lfc", 
          "P-value" = "p"
        ),
        selected = "z",
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
      prettyCheckbox(
        NS(id, "complete"),
        label = "Use complete cases?",
        value = FALSE,
        icon = icon("check"),
        status = "success",
        animation = "rotate"
      ),
      numericInput(
        NS(id, "removeVar"),
        label = "Remove proportion low variance features",
        value = NA,
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
        value = 0.1,
        min = 0,
        max = 1,
        step = 0.1
      ),
      numericInput(
        NS(id, "epochs"),
        "Iterations",
        value = 100,
        min = 1,
        max = Inf,
        step = 100
      ),
      pickerInput(
        NS(id, "metric"),
        "metric",
        choices = c("Euclidean" = "euclidean", "Manhattan" = "manhattan"),
        selected = "euclidean"
      ),
      actionBttn(
        NS(id, "run"),
        label = "Run UMAP",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      dropdownButton(
        tags$h3("UMAP parameters:"),
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
      verbatimTextOutput(NS(id, "test")),
      plotlyOutput(NS(id, "umap")),
      gt_output(NS(id, "table"))
    )
  )
}

umapServer <- function(id, se, keep) {
  moduleServer(id, function(input, output, session) {
    
    # Perform UMAP with selected parameters and data
    # Perform PCA on 'Run'
    udata <- reactive({
      
      show_alert(
        title = "Performing UMAP",
        text = "Please Wait...",
        closeOnClickOutside = FALSE,
        btn_labels = NA,
      )
      
      # Select relevant data and features
      keep_rows <- switch(
        input$features,
        gene = rowData(se)$feature_type == "Gene",
        TE = rowData(se)$feature_type == "TE",
        both = rep(TRUE, nrow(se))
      )
      
      filtered <- se[keep_rows, keep()]
      df <- data.frame(colData(filtered))
      m <- switch(
        input$dataset,
        lfc = assay(filtered, "logFC"),
        p = assay(filtered, "P.Value"),
        z = assay(filtered, "z"),
        t = assay(filtered, "t")
      )
      
      # Impute data if not using complete cases
      if (isTRUE(input$complete)) {
        m <- na.omit(m)
      } else {
        val <- switch(
          input$dataset,
          lfc = 0,
          p = 1,
          z = 0,
          t = 0
        )
        m[is.na(m)] <- val
      }
      
      # Remove low/zero-variance features
      if (is.na(input$removeVar) || input$removeVar == 0) {
        m <- m[matrixStats::rowVars(m, useNames = FALSE) != 0, ]
      } else {
        m <- coriell::remove_var(m, input$removeVar)
      }
      
      # Scale and center the data before computation - result is transposed
      m <- as.matrix(scale(t(m), center = input$center, scale = input$scale))
      
      # Attempt UMAP calculation
      result <- tryCatch({
        coriell::UMAP(
          m,
          metadata = df,
          n_neighbors = input$neighbors,
          metric = input$metric,
          min_dist = input$mindist,
          n_epochs = input$epochs
        )}, 
        error = function(e) {
          return(NULL)
        },
        warning = function(e) {
          return(NULL)
        })
      
      if (is.null(result)) {
        closeSweetAlert()
        validate("Computation failed! Adjust inputs and try again.")
      }
      
      closeSweetAlert()
      return(result)
    }) |> bindEvent(input$run)

    # UMAP plot
    output$umap <- renderPlotly({
      u <- plotUmap(udata(), col = input$col)
      if (isTRUE(input$legend)) {
        u |> plotly::layout(showlegend = TRUE)
      } else {
        u
      }
    })

    output$table <- render_gt({
      d <- event_data("plotly_selected")
      df <- udata()
      if (!is.null(d)) {
        df <- df[d$customdata, ]
      }

      df |>
        gt() |>
        cols_hide(columns = c(id, batch, mutation, comment, desc)) |>
        cols_move(c(tissue, disease), c(cell_line)) |>
        cols_label(
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
        cols_width(
          desc ~ px(450),
          contrast ~ px(300),
          experiment ~ px(150)
        ) |>
        tab_header(
          title = gt::md("**Selected UMAP Data**")
        ) |>
        opt_interactive(use_compact_mode = TRUE)
    })
  })
}
