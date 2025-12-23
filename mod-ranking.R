## Perform ranking of experiments by significant features
##
##
## ----------------------------------------------------------------------------

rankUI <- function(id) {
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
      numericInput(
        NS(id, "fdr"),
        label = "FDR cutoff",
        value = 0.05,
        min = 0,
        max = 1,
        step = 0.01
      ),
      numericInput(
        NS(id, "lfc"),
        label = "logFC cutoff",
        value = 0,
        min = 0,
        max = Inf,
        step = 0.5
      ),
      downloadButton(NS(id, "download"))
    ),
    mainPanel(
      DT::dataTableOutput(NS(id, "table"))
    )
  )
}

rankServer <- function(id, se, keep) {
  moduleServer(id, function(input, output, session) {
    data <- reactive({
      showNotification(
        "Collecting Rank Data...",
        type = "message",
        duration = 5,
        closeButton = TRUE
      )

      keep_rows <- switch(
        input$features,
        gene = SummarizedExperiment::rowData(se)$feature_type == "Gene",
        TE = SummarizedExperiment::rowData(se)$feature_type == "TE",
        both = rep(TRUE, nrow(se))
      )

      filtered <- se[keep_rows, keep()]

      lfc_m <- SummarizedExperiment::assay(filtered, "logFC")
      fdr_m <- SummarizedExperiment::assay(filtered, "adj.P.Val")

      # Calculate the up/down-regulated features
      up_m <- lfc_m > input$lfc & fdr_m < input$fdr
      down_m <- lfc_m < input$lfc & fdr_m < input$fdr

      # Assays contain NA values at same positions and can be used to get total
      # assayed number of features
      total <- DelayedArray::colSums(
        !is.na(SummarizedExperiment::assay(filtered, "AveExpr"))
      )

      # Collect DE results for all contrasts
      dt <- data.table::data.table(
        Total = total,
        N_up = DelayedArray::colSums(up_m, na.rm = TRUE),
        N_down = DelayedArray::colSums(down_m, na.rm = TRUE)
      )

      # Compute percentages
      dt[, `:=`(
        N_non = Total - (N_up + N_down),
        Pct_up = round(N_up / Total * 100, 1),
        Pct_down = round(N_down / Total * 100, 1)
      )][,
        `:=`(Pct_non = round(N_non / Total * 100, 1), Total = NULL)
      ]

      # Add on metadata columns and reorder by most dysregulated
      dt <- data.table::as.data.table(cbind(
        dt,
        data.frame(SummarizedExperiment::colData(filtered))
      ))
      data.table::setorder(dt, Pct_non)
      data.table::setcolorder(
        dt,
        c("experiment", "contrast", "Pct_up", "Pct_down", "Pct_non")
      )

      return(dt)
    })

    output$table <- DT::renderDataTable({
      keep_cols <- c(
        "experiment",
        "contrast",
        "Pct_up",
        "Pct_down",
        "Pct_non",
        "N_up",
        "N_down",
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

      # This is a data.table (..keep_cols needed)
      df <- data()[, ..keep_cols]

      DT::datatable(
        df,
        rownames = FALSE,
        colnames = c(
          "BioProject" = "experiment",
          "Contrast" = "contrast",
          "Percent Up" = "Pct_up",
          "Percent Down" = "Pct_down",
          "Percent Non-DE" = "Pct_non",
          "N Up" = "N_up",
          "N Down" = "N_down",
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

    output$download <- downloadHandler(
      filename = function() {
        paste0("ranking_", format(Sys.time(), "%Y-%m-%d"), ".zip")
      },
      content = function(file) {
        tmp <- tempdir()
        setwd(tmp)

        data.table::fwrite(data(), "data.tsv", sep = "\t")
        files <- "data.tsv"

        zip(zipfile = file, files = files)
      },
      contentType = "application/zip"
    )
  })
}
