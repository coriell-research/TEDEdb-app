## GSEA analysis
##
## Performs GSEA on the selected contrast data
##
## ----------------------------------------------------------------------------

gseaUI <- function(id, choice_list, pathway_names) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      pickerInput(
        NS(id, "ID"),
        label = "Experimental Contrast",
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
      pickerInput(
        NS(id, "pathway"),
        label = "MSigDB Gene Set(s)",
        choices = c(
          "HALLMARK (N=50)" = "h",
          "C1 (N=300): Positional Gene Sets" = "c1",
          "C2 (N=6495): Curated Gene Sets" = "c2",
          "C3 (N=3713): Regulatory Target Gene Sets" = "c3",
          "C4 (N=858): Computational Gene Sets" = "c4",
          "C5 (N=15,937): Ontology Gene Sets" = "c5",
          "C6 (N=189): Oncogenic Gene Sets" = "c6",
          "C7 (N=5219): Immunologic Signature Gene Sets" = "c7",
          "C8 (N=830): Cell Type Signature Gene Sets" = "c8"
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
      radioGroupButtons(
        NS(id, "score"),
        label = "Score type",
        choices = c("Standard" = "std", "Negative" = "neg", "Positive" = "pos"),
        selected = "std",
        status = "primary"
      ),
      actionBttn(
        NS(id, "run"),
        label = "Run GSEA",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      dropdownButton(
        tags$h3("GSEA Plot:"),
        pickerInput(
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
      gt::gt_output(NS(id, "table"))
    )
  )
}

gseaServer <- function(id, se, pathways, pathway_dt) {
  moduleServer(id, function(input, output, session) {
    data <- reactive({
      show_alert(
        title = "Performing GSEA",
        text = "Please Wait...",
        closeOnClickOutside = FALSE,
        btn_labels = NA,
      )
      
      # Subset for only genes
      filtered <- se[rowData(se)$feature_type == "Gene", ]
      z_stats <- assay(filtered, "stat")[, input$ID]
      names(z_stats) <- rownames(filtered)
      z_stats <- z_stats[!is.na(z_stats)]
      
      # Run FGSEA
      res <- fgsea::fgsea(
        pathways = pathways[[input$pathway]],
        stats = z_stats,
        sampleSize = input$size,
        minSize = input$min,
        maxSize = input$max,
        scoreType = input$score
      )
      res <- res[order(padj)]
      
      # Make the GSEA plot from the top result
      updatePickerInput(session, "geneset", selected = res[1, pathway])

      closeSweetAlert()
      return(list(results = res, stats = z_stats))
    }) |> bindEvent(input$run)
    
    # Enrichment Plot
    output$plot <- renderPlot({
      p <- pathway_dt[input$geneset, Pathway]
      fgsea::plotEnrichment(
        pathway = pathways[[p]][[input$geneset]],
        stats = data()[["stats"]]
        ) +
        ggplot2::ggtitle(input$geneset) +
        coriell::theme_coriell()
    })
    
    # GSEA results table
    output$table <- render_gt({ 
      data()[["results"]] |> 
        gt() |> 
        cols_width(
          pathway ~ px(200),
          leadingEdge ~ px(500)
        ) |> 
        tab_header(
          title = gt::md("**GSEA Results**")
        ) |> 
        opt_interactive()
      })
  
  })
}

