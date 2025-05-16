## Select IDs to be used in meta-analysis steps
##
##
## ----------------------------------------------------------------------------

# generic pickerInput creator
pick <- function(id, uid, ulabel, choice_list, multi = TRUE, size = 10,
                 live = TRUE, action = TRUE) {
  shinyWidgets::pickerInput(
    NS(id, uid),
    ulabel,
    choices = choice_list[[uid]],
    selected = choice_list[[uid]],
    multiple = multi,
    options = list(
      title = paste("Select", ulabel),
      size = size,
      `live-search` = live,
      `actions-box` = action
    )
  )
}

# User interface for ID selection
selectIdUI <- function(id, choice_list) {
  sidebarLayout(
    sidebarPanel(
      width = 3,
      pick(id, "experiment", "Experiment(s)", choice_list),
      pick(id, "contrast", "Contrast(s)", choice_list),
      pick(id, "cell_line", "Cell Line(s)", choice_list),
      pick(id, "drug", "Drug(s)", choice_list),
      pick(id, "epigenetic_class", "Epigenetic Class(es)", choice_list),
      pick(id, "drug_class", "Drug Class(es)", choice_list),
      pick(id, "mode_of_action", "Mode of Action", choice_list),
      pick(id, "target", "Target", choice_list),
      pick(id, "tissue", "Tissue(s)", choice_list),
      pick(id, "disease", "Disease(s)", choice_list),
      pick(id, "outlier_flags", "Outlier Flag(s)", choice_list),
      downloadButton(NS(id, "download"))
    ),
    mainPanel(
      gt::gt_output(NS(id, "table"))
    )
  )
}

# Display metadata data of selections and return vector of selected IDs
selectIdServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    df <- data.frame(SummarizedExperiment::colData(se))
    selected <- reactive({
      subset(
        df,
        experiment %in% input$experiment &
        contrast %in% input$contrast &
        outlier_flags %in% input$outlier_flags &
        cell_line %in% input$cell_line &
        drug %in% input$drug &
        epigenetic_class %in% input$epigenetic_class &
        drug_class %in% input$drug_class &
        mode_of_action %in% input$mode_of_action &
        target %in% input$target &
        tissue %in% input$tissue &
        disease %in% input$disease
      )
    })
    
    output$table <- gt::render_gt({
      selected() |> 
        gt::gt() |> 
        gt::cols_hide(columns = c(id, batch, mutation, comment)) |>
        gt::cols_move(c(tissue, disease), c(cell_line)) |> 
        gt::cols_label(
          .list = c(
            "id" = "ID",
            "experiment" = "BioProject ID",
            "contrast" = "Contrast",
            "cell_line" = "Cell Line",
            "drug" = "Drug",
            "dose" = "Dose",
            "time_hr" = "Time (hr)",
            "desc" = "Description",
            "epigenetic_class" = "Epigenetic Class",
            "tissue" = "Tissue",
            "disease" = "Disease",
            "outlier_flags" = "Outlier Flag(s)"
          )
        ) |> 
        gt::cols_width(
          desc ~ px(450),
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
        tmp <- tempdir()
        setwd(tmp)
        
        data.table::fwrite(selected(), "data.tsv", sep = "\t")
        files <- c("data.tsv")
        
        zip(zipfile = file, files = files)
      },
      contentType = "application/zip"
    )

    # Return the selected IDs
    reactive(rownames(selected()))
  })
}
