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
        tags$h4("Analyze Raw Data Locally"),
        tags$p(
          "Copy and run the snippet below in your local R environment to download ",
          "and explore the raw SummarizedExperiment object for the selected project."
        ),
        uiOutput(NS(id, "code_template"))
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
      selected = "PRJNA413957",
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

    output$code_template <- renderUI({
      req(input$selected_project)

      markdown_script <- paste0(
        "```r\n",
        "# Load necessary libraries\n",
        "library(SummarizedExperiment)\n\n",
        "# Define the remote FTP path for the selected BioProject\n",
        "data_url <- 'ftp://data.coriell.org/TEDEdb/BioProjects/",
        input$selected_project,
        "/raw_data.rds'\n",
        "dest_file <- '",
        input$selected_project,
        "_raw_data.rds'\n\n",
        "# Download and import the SummarizedExperiment object\n",
        "download.file(url = data_url, destfile = dest_file, mode = 'wb')\n",
        "se <- readRDS(dest_file)\n\n",
        "# Inspect the object metadata\n",
        "print(se)\n",
        "head(colData(se))\n",
        "```"
      )

      shiny::markdown(markdown_script)
    })
  })
}
