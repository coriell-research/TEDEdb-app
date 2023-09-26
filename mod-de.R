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
        value = 0.1,
        min = 0,
        max = 1,
        step = 0.01
        ),
      numericInput(
        NS(id, "lfc"),
        label = "logFC Cutoff",
        value = 0,
        min = 0,
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
      
      data.frame(
        feature_id = rownames(filtered),
        FDR = assay(filtered, "fdr")[, 1],
        logFC = assay(filtered, "lfc")[, 1],
        logCPM = assay(filtered, "lcpm")[, 1],
        `t` = assay(filtered, "stat")[, 1],
        SE = assay(filtered, "stderr")[, 1]
      )
    })
    
    output$volcano <- renderPlot({
      ptitle <- gsub("_vs_", " vs ", gsub("PRJNA[0-9]+\\.", "", input$ID))
      coriell::plot_volcano(data(), fdr = input$fdr, lfc = input$lfc) +
        ggplot2::ggtitle(ptitle) +
        coriell::theme_coriell()
      })
    output$ma <- renderPlot({
      ptitle <- gsub("_vs_", " vs ", gsub("PRJNA[0-9]+\\.", "", input$ID))
      coriell::plot_md(data(), fdr = input$fdr, lfc = input$lfc) +
        ggplot2::ggtitle(ptitle) +
        coriell::theme_coriell()
      })
    output$table <- render_gt({ 
      data()[order(data()$FDR), ] |> 
        gt() |> 
        tab_header(
          title = gt::md("**Differential Expression Results**")
        ) |> 
        opt_interactive() 
      })
  })
}
