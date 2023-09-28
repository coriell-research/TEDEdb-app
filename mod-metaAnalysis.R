## Perform P-value combination and meta-analysis
##
##
## ----------------------------------------------------------------------------

metaUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      pickerInput(
        NS(id, "method"),
        label = "Combination Method",
        choices = c(
          "Berger", "Fisher", "HolmMin", "Pearson", "Simes",
          "Stouffer", "Wilkinson"
        ),
        selected = "Wilkinson"
      ),
      awesomeCheckbox(
        NS(id, "logp"),
        label = "Log Transform P-Values?",
        value = FALSE,
        status = "danger"
      ),
      numericInput(
        NS(id, "min_n"),
        label = "Min N (Holm, Wilkinson)",
        value = 1,
        min = 1,
        max = Inf,
        step = 1
      ),
      numericInput(
        NS(id, "min_prop"),
        label = "Min Prop. (Holm, Wilkinson)",
        value = 0.5,
        min = 0.01,
        max = 1,
        step = 0.1
      ),
      actionBttn(
        NS(id, "run"),
        label = "Run Meta-Analysis",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      dropdownButton(
        tags$h3("Meta-Volcano parameters:"),
        selectInput(
          NS(id, "x"),
          label = "x",
          choices = c("Rep.logFC", "Median.logFC", "Mean.logFC", "Min.logFC", "Max.logFC"),
          selected = "Rep.logFC"
          ),
        selectInput(
          NS(id, "y"),
          label = "y",
          choices = c("Rep.Pval", "Combined.Pval"),
          selected = "Combined.Pval"
        ),
        status = "danger",
        icon = icon("gear")
        ),
      plotlyOutput(NS(id, "metavolcano")),
      gt_output(NS(id, "table"))
      )
  )
}

metaServer <- function(id, se, keep) {
  moduleServer(id, function(input, output, session) {
    results <- reactive({
      msg <- showNotification("Performing Meta-Analysis. Please wait...",
                              type = "message", duration = NULL,
                              closeButton = FALSE
      )
      filtered <- se[, keep()]
      
      selected_assay <- "fdr"
      if (isTRUE(input$logp)) {
        SummarizedExperiment::assay(filtered, "logp") <- log1p(
          SummarizedExperiment::assay(filtered, "fdr")
          )
        selected_assay <- "logp"
      }

      method <- switch(input$method,
                       Berger = metapod::parallelBerger, 
                       Fisher = metapod::parallelFisher, 
                       HolmMin = metapod::parallelHolmMin, 
                       Pearson = metapod::parallelPearson, 
                       Simes = metapod::parallelSimes,
                       Stouffer = metapod::parallelStouffer,
                       Wilkinson = metapod::parallelWilkinson
      )
      
      # Perform P-value combination 
      if (input$method %in% c("HolmMin", "Wilkinson")) {
        res <- coriell::meta_de(filtered, method, fdr = selected_assay, 
                                lfc = "lfc", min.prop = input$min_prop,
                                min.n = input$min_n)
      } else {
        res <- coriell::meta_de(filtered, method, fdr = selected_assay,
                                lfc = "lfc")
      }
      removeNotification(msg)
      return(res)
    }) |> bindEvent(input$run)
    
    # Show the metavolcano
    output$metavolcano <- renderPlotly({
      plot_ly(
        data = results(),
        customdata = results()[, Feature],
        x = ~ get(input$x),
        y = ~ -log10(get(input$y)),
        color = ~ Direction,
        size = 12,
        text = ~ paste("Gene:", Feature),
        type = "scatter",
        mode = "markers",
        colors = c("up" = "red2", "down" = "blue2", "mixed" = "purple4", "none" = "grey40")
      ) |> 
        plotly::layout(
          title = "Meta-Volcano", 
          xaxis = list(title = "logFC"),
          yaxis = list(title = "-log10(PValue)"),
          dragmode = "lasso"
      ) |>
      event_register("plotly_selected") |> 
      toWebGL()
    })
    
    # Show results in table
    output$table <- render_gt({
      d <- event_data("plotly_selected")
      df <- results()
      if (!is.null(d)) {
        df <- df[Feature %chin% d$customdata]
      }
      
      df |> 
        gt() |> 
        fmt_number(columns = c(Rep.logFC), decimals = 2) |> 
        tab_header(
          title = gt::md("**Meta-Significant Features**")
        ) |> 
        opt_interactive(use_compact_mode = TRUE)
    })
  })
}
