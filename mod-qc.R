## Display BioProject level results
##
##
## ----------------------------------------------------------------------------

qcPlotsUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput(
        NS(id, "selected_project"),
        label = "Select Experiment",
        choices = NULL,
        multiple = FALSE,
        options = list(placeholder = "e.g. PRJNA12345")
      ),
    ),
    mainPanel(
      uiOutput(NS(id, "qc_plots_grid")),

      tags$hr(style = "margin-top: 30px; margin-bottom: 30px;"),

      tags$div(
        style = "margin-bottom: 40px;",
        tags$h4("Download locus-level quants and metadata"),
        tags$p(
          "Copy and run the snippet below in your terminal to download ",
          "the locus-level data and annotations for this experiment."
        ),
        uiOutput(NS(id, "code_template")),
        tags$p(
          "Then copy and run the snippet below in your local R session to ",
          "create a SummarizedExperiment for the locus-level counts."
        ),
        uiOutput(NS(id, "code_template2"))
      )
    )
  )
}

qcPlotsServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    choices <- sort(unique(se[["experiment"]]))
    updateSelectizeInput(
      session,
      "selected_project",
      choices = choices,
      selected = "PRJNA244098",
      server = TRUE
    )

    output$qc_plots_grid <- renderUI({
      req(input$selected_project)

      plots_info <- list(
        list(file = "abundance-hist.png", title = "Abundance Histogram"),
        list(file = "logCPM-boxplots.png", title = "RLE Boxplots"),
        list(file = "logCPM-density.png", title = "logCPM Density"),
        list(file = "pca-biplot.png", title = "PCA Plot"),
        list(file = "quantro-plot.png", title = "Quantro Plot"),
        list(file = "sa-plot.png", title = "SA Plot")
      )

      plot_tags <- lapply(plots_info, function(plot) {
        img_url <- paste0(
          "https://data.coriell.org/TEDEdb/bioprojects/",
          input$selected_project,
          "/results/figures/",
          plot$file
        )

        column(
          width = 6,
          tags$div(
            style = "text-align: center; margin-bottom: 30px;",
            tags$h5(
              plot$title,
              style = "font-weight: 600; margin-bottom: 10px;"
            ),
            tags$img(
              src = img_url,
              alt = paste(plot$title, "for", input$selected_project),
              style = "max-width: 100%; height: auto; border: 1px solid #dee2e6; border-radius: 4px; padding: 5px; box-shadow: 0 0.125rem 0.25rem rgba(0,0,0,0.075);"
            )
          )
        )
      })

      fluidRow(
        do.call(tagList, plot_tags)
      )
    })

    # wget code for downloading raw data
    output$code_template <- renderUI({
      req(input$selected_project)

      markdown_script <- paste0(
        "```r\n",

        "# Download all locus-level quants for this BioProject:\n",
        "wget -r -np -nH --cut-dirs=3 'https://data.coriell.org/TEDEdb/bioprojects/",
        input$selected_project,
        "/quants/'",

        "\n\n# Download Run annotations for this dataset:\n",
        "wget 'https://data.coriell.org/TEDEdb/bioprojects/",
        input$selected_project,
        "/annotation.tsv'",

        "\n\n# Download annotation mapping each locus MD5sum to gene/repElems:\n",
        "wget https://data.coriell.org/TEDEdb/resources/REdiscoverTE_hg38/tx2gene_REdiscoverTE.tsv.gz\n",

        "```"
      )

      shiny::markdown(markdown_script)
    })

    # Rscript for locus-level analysis
    output$code_template2 <- renderUI({
      req(input$selected_project)

      markdown_script <- paste0(
        '```r\n',
        'library(SummarizedExperiment)\n',
        'library(data.table)\n\n\n',
        '# List paths to all quant files\n',
        'quant_files <- list.files(\n',
        '  path = "quants", \n',
        '  pattern = "quant.sf.gz", \n',
        '  recursive = TRUE, \n',
        '  full.names = TRUE\n',
        ')\n\n',
        '# Name each quant file by its Run ID\n',
        'names(quant_files) <- regmatches(\n',
        '  quant_files, \n',
        '  regexpr("SRR[0-9]+", quant_files)\n',
        ')\n\n',
        '# Read in counts and convert to a matrix\n',
        'counts <- rbindlist(\n',
        '  lapply(quant_files, fread, select = c("Name", "NumReads")), \n',
        '  idcol = "Run"\n',
        ')\n',
        'counts <- dcast(\n',
        '  data = counts, \n',
        '  Name ~ Run, \n',
        '  value.var = "NumReads", \n',
        '  value.fill = 0\n',
        ')\n',
        'counts <- as.matrix(counts, rownames = "Name")\n\n',
        '# Read in the run-level metadata for the project\n',
        'metadata <- fread("annotation.tsv")\n',
        'setDF(metadata, rownames = metadata$Run)\n\n',
        '# Read in the mapping from MD5 hashsum to feature\n',
        'features <- fread("tx2gene_REdiscoverTE.tsv.gz")\n',
        'setDF(features, rownames = features$md5)\n\n',
        '# Organize in SummarizedExperiment for downstream analysis\n',
        'se <- SummarizedExperiment(\n',
        '  assays = list(\n',
        '    counts = counts[match(rownames(features), rownames(counts)), ]\n',
        '  ),\n',
        '  colData = metadata[colnames(counts), ],\n',
        '  rowData = features\n',
        ')\n',
        '```'
      )

      shiny::markdown(markdown_script)
    })
  })
}
