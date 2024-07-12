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
      ),
      downloadButton(NS(id, "download"))
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

      # Extract assay data as a data.table for plotting
      assay_data <- lapply(names(assays(filtered)), \(x) assay(filtered, x))
      df <- as.data.frame(do.call(cbind, assay_data))
      colnames(df) <- names(assays(filtered))
      setDT(df, keep.rownames = "feature_id")
      df <- df[!is.na(df$logFC), ]
      setorder(df, adj.P.Val)

      showNotification(
        "Creating Volcano Plot...",
        type = "message", duration = 5,
        closeButton = TRUE
      )
      showNotification(
        "Creating MA Plot...",
        type = "message", duration = 5,
        closeButton = TRUE
      )

      return(df)
    })

    vplot <- reactive({
      ptitle <- gsub("_vs_", " vs ", gsub("PRJNA[0-9]+\\.", "", input$ID))
      coriell::plot_volcano(
        df = data(),
        y = "adj.P.Val",
        fdr = input$fdr,
        lfc = log2(input$fc),
        up_shape = 16,
        down_shape = 16,
        nonde_shape = "."
      ) +
        ggplot2::ggtitle(ptitle) +
        coriell::theme_coriell()
    })
    output$volcano <- renderPlot(vplot())

    maplot <- reactive({
      ptitle <- gsub("_vs_", " vs ", gsub("PRJNA[0-9]+\\.", "", input$ID))
      coriell::plot_md(
        df = data(),
        x = "AveExpr",
        sig_col = "adj.P.Val",
        fdr = input$fdr,
        lfc = log2(input$fc),
        up_shape = 16,
        down_shape = 16,
        nonde_shape = "."
      ) +
        ggplot2::ggtitle(ptitle) +
        coriell::theme_coriell()
    })
    output$ma <- renderPlot(maplot())

    output$table <- render_gt({
      data() |>
        gt() |>
        tab_header(
          title = gt::md("**Differential Expression Results**")
        ) |>
        opt_interactive()
    })

    output$download <- downloadHandler(
      filename = function() {
        paste0("differential-expression_", format(Sys.time(), "%Y-%m-%d"), ".zip")
      },
      content = function(file) {
        tmp <- tempdir()
        setwd(tmp)

        data_file <- data.table::fwrite(data(), "data.tsv", sep = "\t")
        v_plot <- ggsave("volcano-plot.pdf", plot = vplot(), device = "pdf", width = 11, height = 7)
        m_plot <- ggsave("ma-plot.pdf", plot = maplot(), device = "pdf", width = 11, height = 7)
        files <- c("data.tsv", "volcano-plot.pdf", "ma-plot.pdf")

        zip(zipfile = file, files = files)
      },
      contentType = "application/zip"
    )
  })
}
