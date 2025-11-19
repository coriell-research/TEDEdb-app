## Perform over-representation testing with clusterProfiler
##
##
## ----------------------------------------------------------------------------

overrepUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput(
        NS(id, "ID"),
        label = "Select Experimental Contrast",
        choices = NULL,
        multiple = FALSE,
        options = list(placeholder = "e.g. Decitabine vs Control")
      ),
      shinyWidgets::awesomeRadio(
        NS(id, "ontology"),
        label = "Ontology",
        choices = c(
          "Biological Process" = "BP",
          "Molecular Function" = "MF",
          "Cellular Component" = "CC",
          "All" = "ALL"
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
      shinyWidgets::pickerInput(
        NS(id, "adj"),
        label = "P.adjust method",
        choices = c(
          "Holm" = "holm",
          "Hochberg" = "hochberg",
          "Hommel" = "hommel",
          "Bonferroni" = "bonferroni",
          "BH" = "BH",
          "BY" = "BY",
          "FDR" = "fdr",
          "None" = "none"
        ),
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
      shinyWidgets::actionBttn(
        NS(id, "run"),
        label = "Run GO",
        style = "material-flat",
        color = "danger"
      ),
      downloadButton(NS(id, "download"))
    ),
    mainPanel(
      fluidRow(
        column(
          6,
          shinyWidgets::dropdownButton(
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
        column(
          6,
          shinyWidgets::dropdownButton(
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
      DT::dataTableOutput(NS(id, "table"))
    )
  )
}

# Display metadata data of selections and return vector of selected IDs
overrepServer <- function(id, se) {
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

      shinyWidgets::show_alert(
        title = "Performing Gene Ontology Analysis",
        text = "Please Wait...",
        closeOnClickOutside = FALSE,
        btn_labels = NA,
      )

      # Select sample data and remove NA measurements
      filtered <- se[
        SummarizedExperiment::rowData(se)$feature_type == "Gene",
        input$ID
      ]
      assay_data <- lapply(c("adj.P.Val", "logFC"), \(x) {
        as.matrix(SummarizedExperiment::assay(filtered, x))
      })
      df <- as.data.frame(do.call(cbind, assay_data))
      colnames(df) <- c("adj.P.Val", "logFC")
      data.table::setDT(df, keep.rownames = "feature_id")
      df <- df[!is.na(adj.P.Val)]
      setorder(df, adj.P.Val)

      # Extract the gene sets to test
      up_genes <- df[adj.P.Val < input$fdr & logFC > 0, unique(feature_id)]
      down_genes <- df[adj.P.Val < input$fdr & logFC < 0, unique(feature_id)]
      all_genes <- df[, unique(feature_id)]

      # Perform over-representation analysis
      ego_up <- clusterProfiler::enrichGO(
        gene = up_genes,
        universe = all_genes,
        OrgDb = "org.Hs.eg.db",
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

      ego_down <- clusterProfiler::enrichGO(
        gene = down_genes,
        universe = all_genes,
        OrgDb = "org.Hs.eg.db",
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

      shinyWidgets::closeSweetAlert()
      list("down" = ego_down, "up" = ego_up)
    }) |>
      bindEvent(input$run)

    # Table of results
    output$table <- DT::renderDataTable({
      ego_up <- data()[["up"]]
      ego_down <- data()[["down"]]

      if (is.null(ego_up)) {
        validate("No significant results for up-regulated genes!")
      } else if (is.null(ego_down)) {
        validate("No significant results for down-regulated genes!")
      }

      dt <- rbindlist(
        list(
          "Up-regulated" = data.frame(ego_up),
          "Down-regulated" = data.frame(ego_down)
        ),
        idcol = "Gene Set"
      )

      DT::datatable(dt, rownames = FALSE, lazyRender = TRUE, style = "auto")
    })

    # Dotplot
    dplot <- reactive({
      selected <- switch(
        input$geneset,
        up = data()[["up"]],
        down = data()[["down"]]
      )

      if (nrow(data.frame(selected)) == 0) {
        validate("No significant results for selection! Check other gene set.")
      }

      enrichplot::dotplot(selected, showCategory = input$n) +
        ggplot2::ggtitle(paste0(
          tools::toTitleCase(input$geneset),
          "-regulated GO terms"
        ))
    })
    output$dotplot <- renderPlot(dplot())

    # Network
    eplot <- reactive({
      selected <- switch(
        input$geneset2,
        up = data()[["up"]],
        down = data()[["down"]]
      )

      if (nrow(data.frame(selected)) == 0) {
        validate("No significant results for selection! Check other gene set.")
      }

      res <- enrichplot::pairwise_termsim(selected)
      enrichplot::emapplot(
        res,
        showCategory = input$nodes,
        min_edge = input$similarity
      ) +
        ggplot2::ggtitle(paste0(
          tools::toTitleCase(input$geneset2),
          "-regulated GO terms"
        ))
    })
    output$emmap <- renderPlot(eplot())

    output$download <- downloadHandler(
      filename = function() {
        paste0("over-representation_", format(Sys.time(), "%Y-%m-%d"), ".zip")
      },
      content = function(file) {
        tmp <- tempdir()

        data_up <- data.frame(data()[["up"]])
        data_down <- data.frame(data()[["down"]])

        up_path <- file.path(tmp, "data_up-regulated.tsv")
        down_path <- file.path(tmp, "data_down-regulated.tsv")

        data.table::fwrite(data_up, up_path, sep = "\t")
        data.table::fwrite(data_down, down_path, sep = "\t")
        files <- c(up_path, down_path)

        zip(zipfile = file, files = files, flags = "-j")
      },
      contentType = "application/zip"
    )
  })
}
