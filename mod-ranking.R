## Perform ranking of experiments by significant features
##
##
## ----------------------------------------------------------------------------

rankUI <- function(id) {
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
      gt_output(NS(id, "table"))
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
        gene = rowData(se)$feature_type == "Gene",
        TE = rowData(se)$feature_type == "TE",
        both = rep(TRUE, nrow(se))
      )
      
      filtered <- se[keep_rows, keep()]
      lfc_m <- assay(filtered, "logFC")
      fdr_m <- assay(filtered, "adj.P.Val")
      
      # Calculate the up/down-regulated features
      up_m <- lfc_m > input$lfc & fdr_m < input$fdr
      down_m <- lfc_m < input$lfc & fdr_m < input$fdr
      
      # Assays contain NA values at same positions
      total <- colSums(!is.na(lfc_m), na.rm = TRUE)
      
      # Collect DE results for all contrasts
      dt <- data.table(
        Total = total,
        N_up = colSums(up_m, na.rm = TRUE),
        N_down = colSums(down_m, na.rm = TRUE)
      )
      
      # Compute percentages
      dt[, `:=`(N_non = Total - (N_up + N_down),
                Pct_up = round(N_up / Total * 100, 1),
                Pct_down = round(N_down / Total * 100, 1))][,
                `:=`(Pct_non = round(N_non / Total * 100, 1),
                Total = NULL)]
      
      # Add on metadata columns and reorder by most dysregulated
      dt <- setDT(cbind(dt, data.frame(colData(filtered))))
      setorder(dt, Pct_non)
      setcolorder(dt, c("experiment", "contrast", "Pct_up", "Pct_down", "Pct_non"))
      
      return(dt)
    })
    
    output$table <- render_gt({
      data() |> 
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
          contrast ~ px(400),
          experiment ~ px(150),
          epigenetic_class ~ px(150)
        ) |> 
        tab_header(
          title = gt::md("**Experiments Ranked By Count of Significant Features**")
        ) |> 
        opt_interactive(use_compact_mode = TRUE)
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
