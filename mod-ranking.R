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
      gt::gt_output(NS(id, "table"))
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
      
      # Assays contain NA values at same positions
      total <- DelayedArray::colSums(!is.na(lfc_m), na.rm = TRUE)
      
      # Collect DE results for all contrasts
      dt <- data.table::data.table(
        Total = total,
        N_up = DelayedArray::colSums(up_m, na.rm = TRUE),
        N_down = DelayedArray::colSums(down_m, na.rm = TRUE)
      )
      
      # Compute percentages
      dt[, `:=`(N_non = Total - (N_up + N_down),
                Pct_up = round(N_up / Total * 100, 1),
                Pct_down = round(N_down / Total * 100, 1))][,
                `:=`(Pct_non = round(N_non / Total * 100, 1),
                Total = NULL)]
      
      # Add on metadata columns and reorder by most dysregulated
      dt <- data.table::as.data.table(cbind(dt, data.frame(SummarizedExperiment::colData(filtered))))
      data.table::setorder(dt, Pct_non)
      data.table::setcolorder(dt, c("experiment", "contrast", "Pct_up", "Pct_down", "Pct_non"))
      
      return(dt)
    })
    
    output$table <- gt::render_gt({
      data() |> 
        gt::gt() |> 
        gt::cols_hide(columns = c(id, batch, mutation, comment, desc)) |>
        gt::cols_move(c(tissue, disease), c(cell_line)) |> 
        gt::cols_label(
          .list = c(
            "Pct_up" = "Pct. Up",
            "Pct_down" = "Pct. Down",
            "Pct_non" = "Pct. Non",
            "N_up" = "N Up",
            "N_down" = "N Down",
            "N_non" = "N Non",
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
        gt::fmt_percent(
          columns = c("Pct_up", "Pct_down", "Pct_non"), 
          scale_values = FALSE, 
          decimals = 1
        ) |> 
        gt::cols_width(
          contrast ~ px(400),
          experiment ~ px(150),
          epigenetic_class ~ px(150)
        ) |> 
        gt::tab_header(
          title = gt::md("**Experiments Ranked By Count of Significant Features**")
        ) |> 
        gt::opt_interactive(use_compact_mode = TRUE)
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
