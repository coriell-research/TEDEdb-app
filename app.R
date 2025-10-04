suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(data.table))

# Load modules
source("mod-selectIds.R")
source("mod-pca.R")
source("mod-umap.R")
source("mod-metaAnalysis.R")
source("mod-ranking.R")
source("mod-de.R")
source("mod-gsea.R")
source("mod-overrep.R")

# Load global data
choices <- readRDS("data/select-inputs.rds")
se <- HDF5Array::loadHDF5SummarizedExperiment("data/se_hdf5")
DelayedArray::setAutoBlockSize(250e6)

pathways <- readRDS("data/pathways.rds")
pathway_dt <- data.table::fread(
  "data/pathway_dt.tsv.gz", 
  sep = "\t", 
  colClasses = c("character", "character"), 
  key = "Name"
  )

# App ---------------------------------------------------------------------


ui <- navbarPage(
  "TEDEdb App",
  theme = shinythemes::shinytheme("yeti"),
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
      tabPanel(
        "3. UMAP",
        umapUI("umap")
      ),
      tabPanel(
        "4. Meta-Combine",
        metaUI("meta")
      ),
      tabPanel(
        "5. Ranked Expression",
        rankUI("rank")
      )
    )
  ),
  tabPanel(
    "Differential Expression",
    deUI("de", choices)
  ),
  tabPanel(
    "GSEA",
    gseaUI("gsea", choices, pathway_dt[, Name])
  ),
  tabPanel(
    "Over-representation",
    overrepUI("overrep", choices)
  ),
  tabPanel(
    "About",
    htmltools::includeMarkdown("about.md")
  )
)

server <- function(input, output, session) {
  selected <- selectIdServer("ids", se)
  pcaServer("pca", se, selected)
  umapServer("umap", se, selected)
  metaServer("meta", se, selected)
  rankServer("rank", se, selected)
  deServer("de", se)
  gseaServer("gsea", se, pathways, pathway_dt)
  overrepServer("overrep", se)
}

shinyApp(ui, server)
