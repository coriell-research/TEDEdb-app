## Module for comparing the samples vs other samples
##
##
## ----------------------------------------------------------------------------

upsetUI <- function(id, choice_list) {
  sidebarLayout(
    sidebarPanel(
      width = 4,
      multiInput(
        NS(id, "ids"),
        label = "Experimental Contrasts",
        choices = choice_list[["id"]],
        selected = c("PRJNA413957.YB5.HH1.10uM_vs_DMSO_96hr", 
                     "PRJNA413957.YB5.HH1.25uM_vs_DMSO_96hr",
                     "PRJNA413957.YB5.HH1.10uM_vs_DMSO_24hr"),
        width = '100%',
        options = list(
          enable_search = TRUE,
          search_placeholder = "Search Experiments...",
          non_selected_header = 'All options',
          selected_header = 'Selected options'
        )
      ),
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
      awesomeRadio(
        NS(id, "mode"),
        label = "Combination set mode",
        choices = c(
          "Distinct" = "distinct", 
          "Intersect" = "intersect",
          "Union" = "union"
        ),
        selected = "intersect",
        inline = TRUE
      ),
      actionBttn(
        NS(id, "run"),
        label = "Run Upset",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      plotOutput(NS(id, "plot"))
    )
  )
}

upsetServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    result <- reactive({
      msg <- showNotification(
        "Performing Overlaps. Please wait...",
        type = "message", duration = NULL,
        closeButton = FALSE
      )
      
      keep_rows <- switch(
        input$features,
        gene = rowData(se)$feature_type == "Gene",
        TE = rowData(se)$feature_type == "TE",
        both = rep(TRUE, nrow(se))
      )
      
      # Extract the selected features
      filtered <- se[keep_rows, se$id %in% input$ids]
      lfc_m <- assay(filtered, "lfc")
      fdr_m <- assay(filtered, "fdr")
      
      # Calculate the up/down-regulated features across contrasts
      up_m <- lfc_m > input$lfc & fdr_m < input$fdr
      down_m <- lfc_m < input$lfc & fdr_m < input$fdr
      
      # Extract lists of names for each selected contrast
      up_features <- apply(up_m, 2, \(x) names(x[!is.na(x) & x == TRUE]), 
        simplify = FALSE
      )
      down_features <- apply(down_m, 2, \(x) names(x[!is.na(x) & x == TRUE]), 
        simplify = FALSE
      )
      
      # Append direction to list element names
      connames <- gsub("PRJNA[0-9]+\\.", "", names(up_features))
      names(up_features) <- paste0(connames, ".up")
      names(down_features) <- paste0(connames, ".down")
      l <- c(up_features, down_features)
      
      removeNotification(msg)
      return(make_comb_mat(l, mode = input$mode))
    }) |> bindEvent(input$run)
    
    output$plot <- renderPlot({
      m <- result()
      m <- m[comb_degree(m) == 2]
      
      UpSet(
        m, 
        comb_order = order(comb_size(m), decreasing = TRUE),
        top_annotation = upset_top_annotation(m, add_numbers = TRUE),
        right_annotation = upset_right_annotation(m, add_numbers = TRUE)
        )
    })
  }) 
}
