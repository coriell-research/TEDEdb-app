## Differential expression plots for single contrasts
##
##
## ----------------------------------------------------------------------------

# User interface for ID selection
deUI <- function(id) {
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
      selectizeInput(
        NS(id, "ID"),
        label = "Select Experimental Contrast",
        choices = NULL,
        multiple = FALSE,
        options = list(placeholder = "e.g. Decitabine vs Control")
      ),
      numericInput(
        NS(id, "fdr"),
        label = "FDR Cutoff",
        value = 0.05,
        min = 0,
        max = 1,
        step = 0.01
      ),
      numericInput(
        NS(id, "fc"),
        label = "Fold-change Cutoff",
        value = 1.2,
        min = 1,
        max = Inf,
        step = 1
      ),
      downloadButton(NS(id, "download"))
    ),
    mainPanel(
      fluidRow(
        column(6, plotOutput(NS(id, "volcano"))),
        column(6, plotOutput(NS(id, "ma")))
      ),
      fluidRow(
        DT::dataTableOutput(NS(id, "table"))
      )
    )
  )
}

# Display metadata data of selections and return vector of selected IDs
deServer <- function(id, se) {
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

      keep_rows <- switch(
        input$features,
        gene = SummarizedExperiment::rowData(se)$feature_type == "Gene",
        TE = SummarizedExperiment::rowData(se)$feature_type == "TE",
        both = SummarizedExperiment::rowData(se)$feature_type %in%
          c("Gene", "TE")
      )
      keep_col <- input$ID
      filtered <- se[keep_rows, keep_col]

      # Extract assay data as a data.table for plotting
      assay_data <- lapply(names(SummarizedExperiment::assays(filtered)), \(x) {
        SummarizedExperiment::assay(filtered, x)
      })
      df <- as.data.frame(do.call(SummarizedExperiment::cbind, assay_data))
      colnames(df) <- names(SummarizedExperiment::assays(filtered))
      data.table::setDT(df, keep.rownames = "feature_id")
      df <- df[!is.na(adj.P.Val), ]
      setorder(df, adj.P.Val)

      showNotification(
        "Creating Volcano Plot...",
        type = "message",
        duration = 5,
        closeButton = TRUE
      )
      showNotification(
        "Creating MA Plot...",
        type = "message",
        duration = 5,
        closeButton = TRUE
      )

      return(df)
    })

    vplot <- reactive({
      ptitle <- gsub("_vs_", " vs ", gsub("PRJNA[0-9]+\\.", "", input$ID))
      coriell::plot_volcano(
        df = data(),
        y = "adj.P.Val",
        fdr = input$fdr,
        lfc = log2(input$fc),
        up_shape = 16,
        down_shape = 16,
        nonde_shape = "."
      ) +
        ggplot2::ggtitle(ptitle) +
        coriell::theme_coriell()
    })
    output$volcano <- renderPlot(vplot())

    maplot <- reactive({
      ptitle <- gsub("_vs_", " vs ", gsub("PRJNA[0-9]+\\.", "", input$ID))
      coriell::plot_md(
        df = data(),
        x = "AveExpr",
        sig_col = "adj.P.Val",
        fdr = input$fdr,
        lfc = log2(input$fc),
        up_shape = 16,
        down_shape = 16,
        nonde_shape = "."
      ) +
        ggplot2::ggtitle(ptitle) +
        coriell::theme_coriell()
    })
    output$ma <- renderPlot(maplot())

    output$table <- DT::renderDataTable({
      df <- data()

      validate(
        need(
          is.data.frame(df) && nrow(df) > 0,
          "Please select an experimental contrast to view results."
        )
      )

      DT::datatable(df, rownames = FALSE, lazyRender = TRUE, style = "auto")
    })

    output$download <- downloadHandler(
      filename = function() {
        paste0(
          "differential-expression_",
          format(Sys.time(), "%Y-%m-%d"),
          ".zip"
        )
      },
      content = function(file) {
        tmp <- tempdir()

        data_path <- file.path(tmp, "data.tsv")
        volcano_path <- file.path(tmp, "volcano-plot.pdf")
        ma_path <- file.path(tmp, "ma-plot.pdf")

        data.table::fwrite(data(), data_path, sep = "\t")
        ggplot2::ggsave(
          volcano_path,
          plot = vplot(),
          device = "pdf",
          width = 8,
          height = 6
        )
        ggplot2::ggsave(
          ma_path,
          plot = maplot(),
          device = "pdf",
          width = 8,
          height = 6
        )
        files <- c(data_path, volcano_path, ma_path)

        zip(zipfile = file, files = files, flags = "-j")
      },
      contentType = "application/zip"
    )
  })
}
