## Select IDs to be used in meta-analysis steps
##
##
## ----------------------------------------------------------------------------

selectIdUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput(
        NS(id, "experiment"),
        "Experiment(s)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. 'PRJNA12345'")
      ),
      selectizeInput(
        NS(id, "treatment"),
        "Drug(s)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. Decitabine")
      ),
      selectizeInput(
        NS(id, "clinical_phase"),
        "Clinical Phase",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. Preclinical")
      ),
      selectizeInput(
        NS(id, "epigenetic_class"),
        "Epigenetic Class(es)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. Reader")
      ),
      selectizeInput(
        NS(id, "epigenetic_class_fine"),
        "Epigenetic Class(es) - Fine",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. HDACi")
      ),
      selectizeInput(
        NS(id, "mechanism_of_action"),
        "Mechanism of Action",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. cdk inhibitor")
      ),
      selectizeInput(
        NS(id, "targets"),
        "Target(s)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. CDK12")
      ),
      selectizeInput(
        NS(id, "stripped_cell_line"),
        "Cell Line(s)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. SW48")
      ),
      selectizeInput(
        NS(id, "oncotree_primary_disease"),
        "Primary Disease(s)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. Acute Myeloid Leukemia")
      ),
      selectizeInput(
        NS(id, "oncotree_lineage"),
        "Lineage(s)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. Breast")
      ),
      selectizeInput(
        NS(id, "sample_collection_site"),
        "Collection Site(s)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. Colon")
      ),
      selectizeInput(
        NS(id, "outlier_flags"),
        "Outlier Flag(s)",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. None")
      ),

      tags$hr(),

      # Dynamic selection of additional variables
      selectizeInput(
        NS(id, "additional_vars"),
        "Filter by additional variables:",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "e.g. dose")
      ),

      uiOutput(NS(id, "dynamic_filters_ui")),

      tags$hr(),
      downloadButton(NS(id, "download"))
    ),
    mainPanel(
      DT::dataTableOutput(NS(id, "table"))
    )
  )
}

selectIdServer <- function(id, se, choices) {
  moduleServer(id, function(input, output, session) {
    filter_cols <- c(
      "experiment",
      "treatment",
      "clinical_phase",
      "epigenetic_class",
      "epigenetic_class_fine",
      "mechanism_of_action",
      "targets",
      "stripped_cell_line",
      "oncotree_primary_disease",
      "oncotree_lineage",
      "sample_collection_site",
      "outlier_flags"
    )

    lapply(filter_cols, function(col) {
      updateSelectizeInput(
        session,
        col,
        choices = choices[[col]],
        server = TRUE
      )
    })

    extra_cols <- setdiff(names(choices), filter_cols)
    updateSelectizeInput(
      session,
      "additional_vars",
      choices = extra_cols,
      server = TRUE
    )

    output$dynamic_filters_ui <- renderUI({
      req(input$additional_vars)

      dynamic_inputs <- lapply(input$additional_vars, function(col_name) {
        clean_label <- tools::toTitleCase(gsub("_", " ", col_name))

        selectizeInput(
          inputId = session$ns(col_name),
          label = clean_label,
          choices = choices[[col_name]],
          multiple = TRUE,
          options = list(placeholder = paste("Select", clean_label))
        )
      })

      do.call(tagList, dynamic_inputs)
    })

    filtered_se <- reactive({
      filtered <- se

      all_active_filters <- c(filter_cols, input$additional_vars)

      for (col_name in all_active_filters) {
        selected_values <- input[[col_name]]

        if (!is.null(selected_values) && length(selected_values) > 0) {
          filtered <- filtered[,
            filtered[[col_name]] %in% selected_values,
            drop = FALSE
          ]
        }
      }

      return(filtered)
    })

    selected_df <- reactive({
      se_obj <- filtered_se()
      if (ncol(se_obj) > 0) {
        as.data.frame(SummarizedExperiment::colData(se_obj))
      } else {
        data.frame()
      }
    })

    output$table <- DT::renderDataTable(
      {
        keep_cols <- c(
          "experiment",
          "contrast",
          "treatment",
          "clinical_phase",
          "mechanism_of_action",
          "targets",
          "epigenetic_class",
          "epigenetic_class_fine",
          "stripped_cell_line",
          "oncotree_lineage",
          "oncotree_primary_disease"
        )

        df <- selected_df()[, keep_cols, drop = FALSE]

        req(nrow(df) > 0)

        DT::datatable(
          df,
          rownames = FALSE,
          colnames = c(
            "BioProject" = "experiment",
            "Contrast" = "contrast",
            "Treatment" = "treatment",
            "Clinical Phase" = "clinical_phase",
            "Mechanism of Action" = "mechanism_of_action",
            "Target(s)" = "targets",
            "Epigenetic Class" = "epigenetic_class",
            "Epigenetic Class (fine)" = "epigenetic_class_fine",
            "Cell Line" = "stripped_cell_line",
            "Cell Lineage" = "oncotree_lineage",
            "Primary Disease" = "oncotree_primary_disease"
          ),
          lazyRender = FALSE,
          style = "bootstrap4",
          filter = "top"
        )
      },
      server = TRUE
    )

    output$download <- downloadHandler(
      filename = function() {
        paste0("selected-data_", format(Sys.time(), "%Y-%m-%d"), ".zip")
      },
      content = function(file) {
        df_to_save <- selected_df()
        tmp <- tempdir()
        data.table::fwrite(df_to_save, file.path(tmp, "data.tsv"), sep = "\t")
        zip(zipfile = file, files = file.path(tmp, "data.tsv"), flags = "-j")
      },
      contentType = "application/zip"
    )

    return(reactive(colnames(filtered_se())))
  })
}
