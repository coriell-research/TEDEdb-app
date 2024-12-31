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
      ),
      downloadButton(NS(id, "download"))
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
    data <- reactive({
      show_alert(
        title = "Processing Meta-Analysis",
        text = "Please Wait...\nPlots may take additional time to render",
        closeOnClickOutside = FALSE,
        btn_labels = NA,
      )
      filtered <- se[, keep()]
      
      # Impute missing values for testing
      assay(filtered, "P.Value")[is.na(assay(filtered, "P.Value"))] <- 1
      assay(filtered, "logFC")[is.na(assay(filtered, "logFC"))] <- 0

      selected_assay <- "P.Value"
      if (isTRUE(input$logp)) {
        SummarizedExperiment::assay(filtered, "logp") <- log(
          SummarizedExperiment::assay(filtered, "P.Value")
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
      result <- tryCatch(
        {
          if (input$method %in% c("HolmMin", "Wilkinson")) {
            res <- coriell::meta_de(
              filtered,
              method,
              pval = selected_assay,
              min.prop = input$min_prop,
              min.n = input$min_n,
              log.p = isTRUE(input$logp)
            )
          } else {
            res <- coriell::meta_de(
              filtered,
              method,
              pval = selected_assay,
              log.p = isTRUE(input$logp)
            )
          }
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
        closeSweetAlert()
        validate("Computation failed! Adjust inputs and try again.")
      }

      # Rescale log-transformed results back to original scale
      if (isTRUE(input$logp)) {
        res[, Combined.Pval := exp(Combined.Pval)]
      }

      # Drop values that could not be calculated
      before <- nrow(res)
      res <- res[!is.na(Combined.Pval)][order(Combined.Pval)]
      after <- nrow(res)
      showNotification(
        paste("Removing", before - after, "observations where P-values could not be combined"),
        type = "message",
        duration = 10,
        closeButton = TRUE
      )

      # Show this here after creating the results instead of running in reactive below
      showNotification(
        "Creating Meta-Volcano Plot...",
        type = "message",
        duration = 20,
        closeButton = TRUE
      )

      closeSweetAlert()
      return(res)
    }) |> bindEvent(input$run)

    # Show the metavolcano
    output$metavolcano <- renderPlotly({
      plotly::plot_ly(
        colors = c("up" = "red2", "down" = "blue2", "mixed" = "purple2")
      ) |>
        plotly::add_trace(
          data = data()[Direction == "up"],
          x = ~ get(input$x),
          y = ~ -log10(get(input$y)),
          color = ~Direction,
          text = ~ paste("Gene:", Feature),
          customdata = data()[Direction == "up", Feature],
          type = "scatter",
          mode = "markers",
          showlegend = TRUE,
          visible = TRUE
        ) |>
        plotly::add_trace(
          data = data()[Direction == "down"],
          x = ~ get(input$x),
          y = ~ -log10(get(input$y)),
          color = ~Direction,
          text = ~ paste("Gene:", Feature),
          customdata = data()[Direction == "down", Feature],
          type = "scatter",
          mode = "markers",
          showlegend = TRUE,
          visible = TRUE
        ) |>
        plotly::add_trace(
          data = data()[Direction == "mixed"],
          x = ~ get(input$x),
          y = ~ -log10(get(input$y)),
          color = ~Direction,
          text = ~ paste("Gene:", Feature),
          customdata = data()[Direction == "mixed", Feature],
          type = "scatter",
          mode = "markers",
          showlegend = TRUE,
          visible = "legendonly"
        ) |>
        plotly::layout(
          title = "Meta-Volcano",
          xaxis = list(title = "logFC"),
          yaxis = list(title = "-log10(PValue)"),
          dragmode = "lasso"
        ) |>
        plotly::event_register("plotly_selected") |>
        plotly::toWebGL()
    })

    # Show results in table
    output$table <- render_gt({
      d <- event_data("plotly_selected")
      df <- data()
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

    output$download <- downloadHandler(
      filename = function() {
        paste0("meta-combine_", format(Sys.time(), "%Y-%m-%d"), ".zip")
      },
      content = function(file) {
        tmp <- tempdir()
        setwd(tmp)

        data.table::fwrite(data(), "data.tsv", sep = "\t")
        files <- c("data.tsv")

        zip(zipfile = file, files = files)
      },
      contentType = "application/zip"
    )
  })
}
