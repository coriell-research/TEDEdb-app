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
    layout(
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
        value = 200,
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
      dataTableOutput(NS(id, "table"))
    )
  )
}

umapServer <- function(id, pcaobj) {
  moduleServer(id, function(input, output, session) {
    # Perform umap with selected parameters
    udata <- reactive({
      msg <- showNotification("Performing UMAP. Please wait...",
        type = "message", duration = FALSE,
        closeButton = FALSE
      )
      u <- coriell::UMAP(
        pcaobj(),
        n_neighbors = input$neighbors,
        metric = input$metric,
        min_dist = input$mindist,
        n_epochs = input$epochs
      )
      removeNotification(msg)
      return(u)
    }) |> bindEvent(input$run)

    # UMAP plot
    output$umap <- renderPlotly({
      u <- plotUmap(udata(), col = input$col)
      if (isTRUE(input$legend)) {
        u |> layout(showlegend = TRUE)
      } else {
        u
      }
    })

    # Selected data for the table
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
      df <- udata()
      if (!is.null(d)) {
        df2 <- df[d$customdata, keep_cols]
        rownames(df2) <- NULL
        df2
      }
    })
  })
}
