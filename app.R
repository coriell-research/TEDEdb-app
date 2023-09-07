suppressPackageStartupMessages(library(coriell))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinythemes))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(BiocSingular))

# Load modules
source("mod-selectIds.R")
source("mod-pca.R")

# Load global data
choices <- readRDS("data/select-inputs.rds")
se <- readRDS("data/se.rds")


# App ---------------------------------------------------------------------


ui <- navbarPage(
  "GEO Cancer RNAseq",
  theme = shinytheme("yeti"),
  tabPanel(
    "Meta-Analysis",
    tabsetPanel(
      tabPanel(
        "1. Data Selection",
        selectIdUI("ids", choices)
      ),
      tabPanel(
        "2. PCA",
        pcaUI("pca")
      ),
      tabPanel("3. UMAP"),
      tabPanel("4. Meta-Combine"),
      tabPanel("5. Ranked Expression")
    )
  ),
  tabPanel("Differential Expression"),
  tabPanel("GSEA"),
  tabPanel("Over-representation"),
  tabPanel("Sample vs. Sample"),
  tabPanel(
    "About",
    htmltools::includeMarkdown("about.md")
  )
)

server <- function(input, output, session) {
  selected_ids <- selectIdServer("ids", se)
  pcaData <- pcaServer("pca", se, selected_ids)
}

shinyApp(ui, server)
