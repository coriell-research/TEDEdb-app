## Select IDs to be used in meta-analysis steps
##
##
## ----------------------------------------------------------------------------

# generic pickerInput creator
pick <- function(id, uid, ulabel, choice_list, multi = TRUE, size = 10,
                 live = TRUE, action = TRUE) {
  pickerInput(
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
      pick(id, "experiment", "Experiment(s)", choice_list),
      pick(id, "cell_line", "Cell Line(s)", choice_list),
      pick(id, "drug", "Drug(s)", choice_list),
      pick(id, "epigenetic_class", "Drug Class(es)", choice_list),
      pick(id, "tissue", "Tissue(s)", choice_list),
      pick(id, "disease", "Disease(s)", choice_list)
    ),
    mainPanel(
      dataTableOutput(NS(id, "table"))
    )
  )
}

# Display metadata data of selections and return vector of selected IDs
selectIdServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    df <- data.frame(colData(se))
    selected <- reactive({
      subset(
        df,
        experiment %in% input$experiment &
          cell_line %in% input$cell_line &
          drug %in% input$drug &
          epigenetic_class %in% input$epigenetic_class &
          tissue %in% input$tissue &
          disease %in% input$disease
      )
    })
    output$table <- renderDataTable({
      datatable(
        selected(),
        rownames = FALSE,
        colnames = c(
          "ID" = "id", "BioProject" = "experiment",
          "Contrast" = "contrast", "Cell Line" = "cell_line",
          "Drug" = "drug", "Dose" = "dose", "Time (hr)" = "time_hr",
          "Batch" = "batch", "Mutation" = "mutation",
          "Comment" = "comment", "Description" = "desc",
          "Epigenetic Class" = "epigenetic_class",
          "Tissue" = "tissue", "Disease" = "disease"
        ),
        options = list(autoWidth = TRUE)
      )
    })

    # Return the selected IDs
    reactive(rownames(selected()))
  })
}
