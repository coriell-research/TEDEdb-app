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
      )
    ),
    mainPanel(
      gt_output(NS(id, "table"))
    )
  )
}

rankServer <- function(id, se, keep) {
  moduleServer(id, function(input, output, session) { 
    output$table <- render_gt({
      showNotification(
        "Collecting Rank Data...",
        type = "message", duration = 10,
        closeButton = TRUE
      )
      keep_rows <- switch(input$features,
                          gene = rowData(se)$feature_type == "Gene",
                          TE = rowData(se)$feature_type == "TE",
                          both = rep(TRUE, nrow(se))
      )
      
      filtered <- se[keep_rows, keep()]
      lfc_m <- assay(filtered, "lfc")
      fdr_m <- assay(filtered, "fdr")
      
      # Calculate the up/down-regulated features
      up_m <- lfc_m > input$lfc & fdr_m < input$fdr
      down_m <- lfc_m < input$lfc & fdr_m < input$fdr
      
      dt <- data.table(
        Significant = colSums(up_m + down_m, na.rm = TRUE),
        N_up = colSums(up_m, na.rm = TRUE),
        N_down = colSums(down_m, na.rm = TRUE)
      )
      dt <- setDT(cbind(dt, data.frame(colData(filtered))))
      dt <- dt[order(Significant, decreasing = TRUE)]
      
      dt |> 
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
  })
}
