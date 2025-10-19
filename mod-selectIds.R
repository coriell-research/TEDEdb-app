## Select IDs to be used in meta-analysis steps
##
##
## ----------------------------------------------------------------------------

selectIdUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput(NS(id, "experiment"), "Experiment(s)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. 'PRJNA12345'")),
      selectizeInput(NS(id, "contrast"), "Contrast(s)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. Decitabine_vs_DMSO")),
      selectizeInput(NS(id, "cell_line"), "Cell Line(s)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. SW48")),
      selectizeInput(NS(id, "drug"), "Drug(s)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. Decitabine")),
      selectizeInput(NS(id, "epigenetic_class"), "Epigenetic Class(es)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. HDACi")),
      selectizeInput(NS(id, "drug_class"), "Drug Class(es)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. Kinase inhibitor")),
      selectizeInput(NS(id, "mode_of_action"), "Mode of Action", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. CDK9 inhibitor")),
      selectizeInput(NS(id, "target"), "Target", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. DNMT1")),
      selectizeInput(NS(id, "sample_collection_site"), "Collection Site(s)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. Colon")),
      selectizeInput(NS(id, "oncotree_primary_disease"), "Primary Disease(s)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. Acute Myeloid Leukemia")),
      selectizeInput(NS(id, "outlier_flags"), "Outlier Flag(s)", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. None")),
      downloadButton(NS(id, "download"))
    ),
    mainPanel(
      gt::gt_output(NS(id, "table"))
    )
  )
}

selectIdServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    filter_cols <- c(
      "experiment", "contrast", "cell_line", "drug", "epigenetic_class",
      "drug_class", "mode_of_action", "target", "sample_collection_site",
      "oncotree_primary_disease", "outlier_flags"
    )
    
    # Update the selections on the server-side
    lapply(filter_cols, function(col) {
      choices <- sort(unique(se[[col]]))
      updateSelectizeInput(session, col, choices = choices, server = TRUE)
    })
    
    # Iteratively filter the SE object 
    filtered_se <- reactive({
      filtered <- se
      for (col_name in filter_cols) {
        selected_values <- input[[col_name]]
        if (!is.null(selected_values) && length(selected_values) > 0) {
          filtered <- filtered[, filtered[[col_name]] %in% selected_values, drop = FALSE]
        }
      }
      
      return(filtered)
    })
    
    selected_df <- reactive({
      se_obj <- filtered_se()
      if (ncol(se_obj) > 0) {
        as.data.frame(colData(se_obj))
      } else {
        data.frame()
      }
    })
    
    # Render the output table
    output$table <- gt::render_gt({
      df_display <- selected_df()
      
      # Do not render the table if there is no data
      req(nrow(df_display) > 0)
      
      df_display |>
        gt::gt() |>
        gt::cols_hide(columns = any_of(c("id", "batch", "mutation", "comment"))) |>
        gt::cols_move(
          columns = any_of(c("sample_collection_site", "oncotree_primary_disease")), 
          after = any_of("cell_line")
        ) |>
        gt::cols_label(
          .list = c(
            "id" = "ID",
            "experiment" = "BioProject ID",
            "contrast" = "Contrast",
            "cell_line" = "Cell Line",
            "drug" = "Drug",
            "dose" = "Dose",
            "time_hr" = "Time (hr)",
            "description" = "Description",
            "epigenetic_class" = "Epigenetic Class",
            "sample_collection_site" = "Sample Collection Site",
            "oncotree_primary_disease" = "Primary Disease",
            "outlier_flags" = "Outlier Flag(s)"
          )
        ) |>
        gt::cols_width(
          description ~ px(450),
          contrast ~ px(300),
          experiment ~ px(150)
        ) |>
        gt::tab_header(
          title = gt::md("**Selected Data for Meta-Analysis**")
        ) |>
        gt::opt_interactive(use_compact_mode = TRUE)
    })
    
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
