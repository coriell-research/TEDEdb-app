suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinythemes))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(coriell))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(BiocSingular))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(metapod))

# Load modules
source("mod-selectIds.R")
source("mod-pca.R")
source("mod-umap.R")
source("mod-metaAnalysis.R")
source("mod-ranking.R")
source("mod-de.R")
source("mod-gsea.R")
source("mod-overrep.R")
source("mod-upset.R")

# Load global data
choices <- readRDS("data/select-inputs.rds")
se <- readRDS("data/se.rds")
pathways <- readRDS("data/pathways.rds")
pathway_dt <- readRDS("data/pathway_dt.rds")
setkey(pathway_dt, Name)

# App ---------------------------------------------------------------------


ui <- navbarPage(
  "ComiTE App",
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
    "Sample vs. Sample",
    upsetUI("upset", choices)
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
  upsetServer("upset", se)
}

shinyApp(ui, server)
