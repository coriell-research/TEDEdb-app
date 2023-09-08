## Perform UMAP on the PCA results
##
##
## ---------------------------------------------------------------------------


umapUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      numericInput(
        NS(id, "neighbors"),
        "N neighbors",
        value = 15,
        min = 2,
        max = Inf,
        step = 1
      ),
      numericInput(
        NS(id, "mindist"),
        "Minimum distance",
        value = 0.1,
        min = 0,
        max = 1,
        step = 0.1
      ),
      numericInput(
        NS(id, "spread"),
        "Spread",
        value = 1,
        min = 0,
        max = Inf,
        step = 0.1
      ),
      numericInput(
        NS(id, "epochs"),
        "Iterations",
        value = 200,
        min = 1,
        max = Inf,
        step = 100
      ),
      pickerInput(
        NS(id, "metric"),
        "metric",
        choices = c("Euclidean" = "euclidean", "Manhattan" = "manhattan"),
        selected = "euclidean"
      ),
      actionBttn(
        NS(id, "run"),
        label = "Run UMAP",
        style = "material-flat",
        color = "danger"
      )
    ),
    mainPanel(
      verbatimTextOutput(NS(id, "test")),
      plotOutput(NS(id, "umap")),
      dataTableOutput(NS(id, "table"))
    )
  )
}

umapServer <- function(id, pcaobj) {
  moduleServer(id, function(input, output, session) {
    udata <- reactive({
      msg <- showNotification("Performing UMAP. Please wait...",
        type = "message", duration = FALSE,
        closeButton = FALSE
      )
      removeNotification(msg)
      sample(100, 1)
    }) |> bindEvent(input$run)

    output$test <- renderText(udata())
  })
}
