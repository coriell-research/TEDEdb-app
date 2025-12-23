## GSEA analysis
##
## Performs GSEA on the selected contrast data
##
## ----------------------------------------------------------------------------

gseaUI <- function(id, pathway_names) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput(
        NS(id, "ID"),
        label = "Select Experimental Contrast",
        choices = NULL,
        multiple = FALSE,
        options = list(placeholder = "e.g. Decitabine vs Control")
      ),
      shinyWidgets::pickerInput(
        NS(id, "pathway"),
        label = "MSigDB Gene Set(s)",
        choices = c(
          "HALLMARK (N=50)" = "h",
          "C6 (N=189): Oncogenic Gene Sets" = "c6",
          "C7 (N=5219): Immunologic Signature Gene Sets" = "c7"
        ),
        selected = "h",
        multiple = FALSE,
        options = list(
          title = "Select Experiment",
          size = 10,
          `live-search` = TRUE,
          `actions-box` = TRUE
        )
      ),
      numericInput(
        NS(id, "perm"),
        label = "Permutations",
        value = 1e3,
        min = 10,
        max = Inf,
        step = 100
      ),
      numericInput(
        NS(id, "size"),
        label = "Size of random set",
        value = 101,
        min = 30,
        max = Inf,
        step = 10
      ),
      numericInput(
        NS(id, "min"),
        label = "Minimum gene set size",
        value = 1,
        min = 1,
        max = Inf,
        step = 1
      ),
      numericInput(
        NS(id, "max"),
        label = "Maximum gene set size",
        value = 15000,
        max = Inf,
        min = 10
      ),
      shinyWidgets::actionBttn(
        NS(id, "run"),
        label = "Run GSEA",
        style = "material-flat",
        color = "danger"
      ),
      downloadButton(NS(id, "download"))
    ),
    mainPanel(
      shinyWidgets::dropdownButton(
        tags$h3("GSEA Plot:"),
        shinyWidgets::pickerInput(
          NS(id, "geneset"),
          label = "Select gene set to plot",
          choices = pathway_names,
          selected = NULL,
          multiple = FALSE,
          options = list(
            title = "Select pathway",
            size = 10,
            `live-search` = TRUE,
            `actions-box` = TRUE,
            container = "body"
          )
        ),
        status = "danger",
        icon = icon("gear")
      ),
      plotOutput(NS(id, "plot")),
      DT::dataTableOutput(NS(id, "table"))
    )
  )
}

gseaServer <- function(id, se, pathways, pathway_dt) {
  moduleServer(id, function(input, output, session) {
    choices <- sort(unique(se[["id"]]))
    updateSelectizeInput(
      session,
      "ID",
      choices = choices,
      selected = "PRJNA413957.YB5.HH1.10uM_vs_DMSO_96hr",
      server = TRUE
    )

    data <- reactive({
      req(input$ID)

      shinyWidgets::show_alert(
        title = "Performing GSEA",
        text = "Please Wait...",
        closeOnClickOutside = FALSE,
        btn_labels = NA,
      )

      # Subset for only genes
      filtered <- se[
        SummarizedExperiment::rowData(se)$feature_type == "Gene",
        input$ID
      ]
      z_stats <- SummarizedExperiment::assay(filtered, "z")[, 1]
      names(z_stats) <- rownames(filtered)
      z_stats <- z_stats[!is.na(z_stats)]

      # Run FGSEA
      res <- fgsea::fgsea(
        pathways = pathways[[input$pathway]],
        stats = z_stats,
        sampleSize = input$size,
        minSize = input$min,
        maxSize = input$max
      )
      res <- res[order(padj), ]

      # Make the GSEA plot from the top result
      shinyWidgets::updatePickerInput(
        session,
        "geneset",
        selected = res[1, pathway]
      )
      shinyWidgets::closeSweetAlert()

      return(list(results = res, stats = z_stats))
    }) |>
      bindEvent(input$run)

    # Enrichment Plot
    eplot <- reactive({
      p <- pathway_dt[input$geneset, Pathway]
      fgsea::plotEnrichment(
        pathway = pathways[[p]][[input$geneset]],
        stats = data()[["stats"]]
      ) +
        ggplot2::ggtitle(input$geneset) +
        coriell::theme_coriell()
    })
    output$plot <- renderPlot(eplot())

    # GSEA results table
    output$table <- DT::renderDataTable({
      df <- data()[["results"]]

      validate(
        need(
          is.data.frame(df) && nrow(df) > 0,
          "Please select an experimental contrast to view results."
        )
      )

      cols <- c("pval", "padj", "log2err", "ES", "NES", "size")
      df[, (cols) := lapply(.SD, round, 3), .SDcols = cols]
      df[, leadingEdge := NULL]

      DT::datatable(
        df,
        rownames = FALSE,
        colnames = c(
          "Pathway" = "pathway",
          "P.Value" = "pval",
          "adj.P.Val" = "padj"
        ),
        lazyRender = TRUE,
        style = "bootstrap4",
        class = 'cell-border stripe'
      )
    })

    output$download <- downloadHandler(
      filename = function() {
        paste0("gsea_", format(Sys.time(), "%Y-%m-%d"), ".zip")
      },
      content = function(file) {
        tmp <- tempdir()

        data_path <- file.path(tmp, "data.tsv")
        eplot_path <- file.path(tmp, "enrichment-plot.pdf")

        data.table::fwrite(
          data()[["results"]],
          data_path,
          sep = "\t",
          sep2 = c("", " ", "")
        )
        ggplot2::ggsave(
          eplot_path,
          plot = eplot(),
          device = "pdf",
          width = 8,
          height = 6
        )

        files <- c(data_path, eplot_path)

        zip(zipfile = file, files = files, flags = "-j")
      },
      contentType = "application/zip"
    )
  })
}
