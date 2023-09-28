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
      
    ),
    mainPanel()
  )
}

# Display metadata data of selections and return vector of selected IDs
overrepServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {})
}
