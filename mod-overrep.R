## Perform over-representation testing with clusterProfiler
##
##
## ----------------------------------------------------------------------------

overrepUI <- function(id, choice_list) {
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
      awesomeRadio(
        NS(id, "ontology"),
        label = "Ontology",
        choices = c(
          "Biological Process" = "BP", "Molecular Function" = "MF",
          "Cellular Component" = "CC", "All" = "ALL"
        ),
        selected = "BP"
      ),
      numericInput(
        NS(id, "fdr"),
        label = "FDR cutoff (experiment)",
        value = 0.05,
        min = 0,
        max = 1,
        step = 0.01
      ),
      numericInput(
        NS(id, "pval"),
        label = "P-value cutoff (GO)",
        value = 0.05,
        min = 0,
        max = 1,
        step = 0.01
        ),
      numericInput(
        NS(id, "qval"),
        label = "Q-value cutoff (GO)",
        value = 0.1,
        min = 0,
        max = 1,
        step = 0.01
      ),
      pickerInput(
        NS(id, "adj"),
        label = "P.adjust method",
        choices = c("Holm" = "holm", "Hochberg" = "hochberg", 
                    "Hommel" = "hommel", "Bonferroni" = "bonferroni", 
                    "BH" = "BH", "BY" = "BY", "FDR" = "fdr", "None" = "none"),
        selected = "BH",
        multiple = FALSE
        ),
      numericInput(
        NS(id, "maxss"),
        label = "Max gene set size",
        value = 500,
        min = 1,
        max = Inf,
        step = 10
        ),
      numericInput(
        NS(id, "minss"),
        label = "Min gene set size",
        value = 10,
        min = 1,
        max = Inf,
        step = 10
      ),
      actionBttn(
        NS(id, "run"),
        label = "Run GO",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      fluidRow(
        column(6,
           dropdownButton(
             tags$h3("Dotplot parameters:"),
             selectInput(
               NS(id, "geneset"),
               label = "Gene Set",
               choices = c("Up-regulated" = "up", "Down-regulated" = "down"),
               selected = "up"
             ),
             numericInput(
               NS(id, "n"),
               label = "Number of terms to display",
               value = 10,
               min = 1,
               max = Inf,
               step = 1
             ),
             status = "danger",
             icon = icon("gear")
          ),
          plotOutput(NS(id, "dotplot"))
          ),
        column(6,
           dropdownButton(
             tags$h3("EM plot parameters:"),
             selectInput(
               NS(id, "geneset2"),
               label = "Gene Set",
               choices = c("Up-regulated" = "up", "Down-regulated" = "down"),
               selected = "up"
             ),
             numericInput(
               NS(id, "nodes"),
               label = "Number of terms to display",
               value = 10,
               min = 1,
               max = Inf,
               step = 1
             ),
             numericInput(
               NS(id, "similarity"),
               label = "Minimum similarity threshold",
               value = 0.5,
               min = 0,
               max = 1,
               step = 0.01
             ),
             status = "danger",
             icon = icon("gear")
           ),
          plotOutput(NS(id, "emmap"))
          )
      ),
      gt_output(NS(id, "table"))
    )
  )
}

# Display metadata data of selections and return vector of selected IDs
overrepServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    results <- reactive({
      show_alert(
        title = "Performing Gene Ontology Analysis",
        text = "Please Wait...",
        closeOnClickOutside = FALSE,
        btn_labels = NA,
      )
      
      # Select sample data and remove NA measurements
      filtered <- se[rowData(se)$feature_type == "Gene", input$ID]
      fdr_m <- assay(filtered, "fdr")
      fdr_m <- na.omit(fdr_m)
      lfc_m <- assay(filtered, "lfc")[rownames(fdr_m), ]
      
      # Create vectors of up/down genes
      up <- as.vector(fdr_m < input$fdr & lfc_m > 0)
      down <- as.vector(fdr_m < input$fdr & lfc_m < 0)
      up_genes <- rownames(filtered[up, ])
      down_genes <- rownames(filtered[down, ])
      
      # Perform over-representation analysis
      ego_up <- enrichGO(
        up_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = input$ontology,
        pvalueCutoff = input$pval,
        pAdjustMethod = input$adj,
        qvalueCutoff = input$qval,
        minGSSize = input$minss,
        maxGSSize = input$maxss,
        pool = TRUE,
        readable = TRUE
      )
      
      ego_down <- enrichGO(
        down_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = input$ontology,
        pvalueCutoff = input$pval,
        pAdjustMethod = input$adj,
        qvalueCutoff = input$qval,
        minGSSize = input$minss,
        maxGSSize = input$maxss,
        pool = TRUE,
        readable = TRUE
      )
      
      closeSweetAlert()
      list("down" = ego_down, "up" = ego_up)
    }) |> bindEvent(input$run)
    
    # Table of results
    output$table <- render_gt({
      ego_up <- results()[["up"]]
      ego_down <- results()[["down"]]
      
      if (is.null(ego_up) & is.null(ego_down)) {
        validate("No significant GO terms were found!")
      }
      
      dt <- rbindlist(
        list("Up-regulated" = data.frame(ego_up), 
             "Down-regulated" = data.frame(ego_down)),
             idcol = "Gene Set"
             )
      
      dt |> 
        gt()|> 
        cols_width(
          geneID ~ px(450)
        ) |> 
        tab_header(
          title = gt::md("**Enriched GO terms**")
        ) |> 
        opt_interactive(use_compact_mode = TRUE)
    })
    
    # Dotplot
    output$dotplot <- renderPlot({
      selected <- switch(
        input$geneset,
        up = results()[["up"]],
        down = results()[["down"]]
        )
      
      if (nrow(data.frame(selected)) == 0) {
        validate("No significant GO terms were found!")
      }
      
      enrichplot::dotplot(selected, showCategory = input$n) +
        ggtitle(paste0(tools::toTitleCase(input$geneset), "-regulated GO terms"))
    })
    
    # Dotplot
    output$emmap <- renderPlot({
      selected <- switch(
        input$geneset2,
        up = results()[["up"]],
        down = results()[["down"]]
      )
      
      if (nrow(data.frame(selected)) == 0) {
        validate("No significant GO terms were found!")
      }
      
      res <- enrichplot::pairwise_termsim(selected)
      enrichplot::emapplot(
        res, 
        showCategory = input$nodes, 
        edge.params = list(min = input$similarity)
        ) +
        ggtitle(paste0(tools::toTitleCase(input$geneset2), "-regulated GO terms"))
    })
    
  })
}
