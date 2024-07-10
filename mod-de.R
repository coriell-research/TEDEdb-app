## Differential expression plots for single contrasts
##
##
## ----------------------------------------------------------------------------

# User interface for ID selection
deUI <- function(id, choice_list) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      awesomeRadio(
        NS(id, "features"),
        label = "Select features",
        choices = c(
          "Genes" = "gene", "Transposable Elements" = "TE",
          "Both" = "both"
        ),
        selected = "gene",
        inline = TRUE
      ),
      pickerInput(
        NS(id, "ID"), 
        label = "Select Experimental Contrast",
        choices = choice_list[["id"]], 
        selected = "PRJNA413957.YB5.HH1.10uM_vs_DMSO_96hr",
        multiple = FALSE,
        options = list(
          title = "Select Experiment",
          size = 10,
          `live-search` = TRUE,
          `actions-box` = TRUE
          )
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
        )
      ),
    mainPanel(
      fluidRow(
        column(6, plotOutput(NS(id, "volcano"))),
        column(6, plotOutput(NS(id, "ma")))
        ),
      fluidRow(
        gt_output(NS(id, "table"))
      )
    )
  )
}

# Display metadata data of selections and return vector of selected IDs
deServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    data <- reactive({
      keep_rows <- switch(input$features,
                          gene = rowData(se)$feature_type == "Gene",
                          TE = rowData(se)$feature_type == "TE",
                          both = rep(TRUE, nrow(se))
      )
      keep_col <- input$ID
      filtered <- se[keep_rows, keep_col]
      
      dt <- data.table(
        feature_id = rownames(filtered),
        FDR = assay(filtered, "adj.P.Val")[, 1],
        logFC = assay(filtered, "logFC")[, 1],
        logCPM = assay(filtered, "AveExpr")[, 1],
        t =  assay(filtered, "t")[, 1],
        z = assay(filtered, "z")[, 1]
      )
      na.omit(dt)
    })
    
    output$volcano <- renderPlot({
      showNotification(
        "Creating Volcano Plot...",
        type = "message", duration = 5,
        closeButton = TRUE
      )
      ptitle <- gsub("_vs_", " vs ", gsub("PRJNA[0-9]+\\.", "", input$ID))
      coriell::plot_volcano(data(), fdr = input$fdr, lfc = log2(input$fc), 
                            up_shape = 16, down_shape = 16, nonde_shape = '.') +
        ggplot2::ggtitle(ptitle) +
        coriell::theme_coriell()
      })
    output$ma <- renderPlot({
      showNotification(
        "Creating MA Plot...",
        type = "message", duration = 5,
        closeButton = TRUE
      )
      ptitle <- gsub("_vs_", " vs ", gsub("PRJNA[0-9]+\\.", "", input$ID))
      coriell::plot_md(data(), fdr = input$fdr, lfc = log2(input$fc),
                       up_shape = 16, down_shape = 16, nonde_shape = '.') +
        ggplot2::ggtitle(ptitle) +
        coriell::theme_coriell()
      })
    
    output$table <- render_gt({ 
      data()[order(FDR)] |> 
        gt() |> 
        tab_header(
          title = gt::md("**Differential Expression Results**")
        ) |> 
        opt_interactive() 
      })
  })
}
