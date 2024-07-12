## GEO RNA-seq Cancer Meta-analysis

This app integrates differential expression (DE) results generated from a uniform
processing pipeline into an interactive meta-analysis framework. The app allows 
for clustering of DE results, dimensionality reduction, exploration of volcano 
and MA plots, gene set enrichment, and over-representation analysis. 

## Meta-analysis

The meta-analysis tab enables the exploration of DE results across many 
experiments. Each of it's sub-tabs allows for different levels of across 
aggregated differential expression results. 

### 1. Data Selection

The data selection tab provides an interface for selecting contrasts 
(i.e. Treatment vs Control comparisons) that will be included in any page that 
falls under the "meta-analysis" tab. Users can use the dropdowns to dynamically 
select subsets of contrasts by their BioProject ID, drug class, tissue type, etc. 
By default, every contrast in the database is used. 

### 2. PCA

PCA is performed on data from the selected contrasts. The user can opt to 
perform PCA on either the genes alone, TEs alone, or both. The input data values 
for PCA are one of the values computed from the differential expression analysis 
for a particular feature, e.g. logFC, t-statistic, P-value, or z-statistic. 
Using these values for dimensionality reduction allows for the exploration of
patterns of differential expression across disparate treatments, cell lines, 
experiments, etc. 

The available options for PCA are:

- **Center data?**: Whether or not the input data should be mean-centered prior to performing PCA
- **Scale data?**: Whether or not the input data should be scaled prior to performing PCA
- **Use complete cases?**: Should only features present across all contrasts be used in 
the PCA computation? If complete cases is not selected then default values 
are imputed for each data type (logFC=0, P.Value=1, z-statistic=0, t-statistic=0).
- **Number of components**: How many components should be calculated when performing PCA
- **Remove proportion low variance features**: Removes this proportion of features 
in the selected data by unsupervised variance across all contrasts. If left blank 
then only features with constant variance are removed prior to downstream 
PCA computation.
- **PCA algorithm**: A `BiocSingularParam` used to specify the algorithm to use for 
SVD.

After performing PCA, a biplot of the results is displayed to the user. Clicking 
the gear icon at the top left allows the user to select which components to 
use on each axis as well as options for coloring data points by metadata 
information. Users can use the lasso selector on the plot to select specific 
data points. Metadata pertaining to the selected points is displayed below the
plot.

### 3. UMAP

Like PCA, UMAP is performed on data from the selected contrasts. The user can opt 
to perform UMAP on either the genes alone, TEs alone, or both. The input data values 
for UMAP are one of the values computed from the differential expression analysis 
for a particular feature, e.g. logFC, t-statistic, P-value, or z-statistic. 

The available options for UMAP are:

- **Center data?**: Whether or not the input data should be mean-centered prior to performing UMAP
- **Scale data?**: Whether or not the input data should be scaled prior to performing UMAP
- **Use complete cases?**: Should only features present across all contrasts be used in 
the UMAP computation? If complete cases is not selected then default values 
are imputed for each data type (logFC=0, P.Value=1, z-statistic=0, t-statistic=0).
- **Remove proportion low variance features**: Removes this proportion of features 
in the selected data by unsupervised variance across all contrasts. If left blank 
then only features with constant variance are removed prior to downstream 
UMAP computation.
- **N neighbors**: Number of nearest neighbors
- **Minimum distance**: Determines how close points are in the final layout
- **Iterations**: the number of iterations used during layout optimization
- **Metric**: Determines how distances between points are computed

After performing UMAP, a scatter plot of the computed embeddings is displayed 
to the user. Clicking the gear icon at the top left allows the user to select 
options for coloring data points by metadata information. Users can use the 
lasso selector on the plot to select specific data points. Metadata pertaining 
to the selected points is displayed below the plot.

### 4. Meta-Combine

The meta-combine tab allows for the exploration of global patterns of 
dysregulation across all selected contrasts by implementing p-value combination 
methods on the differential expression results via the `metapod` R package 
(see package [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/metapod/inst/doc/overview.html#4_Summarizing_the_direction) for more details). Users can select one of the following p-value combination techniques to be performed on the raw p-values and logFC values from each contrast:

The following yield a significant result where any individual tests are significant

- **Fisher**: Most sensitive to smallest p-value
- **Pearson**: Most sensitive to largest p-value
- **Stouffer**: Compromise between Fisher's and Pearson's methods
- **Simes**: More conservative than the above but functional in the presence of dependencies

The following yield significant results where *some* of the individual tests are
significant (a minimum number of tests)

- **Wilkinson**: Requires independence. User most supply minimum proportion of 
tests or minimum number.
- **Holm**: Does not require independence. User most supply minimum proportion of 
tests or minimum number.
- **Berger**: Will yield a significant result for groups where *all* of the 
individual tests are significant

After calculating the meta-differential expression results a meta-volcano plot is
displayed. The meta-volcano can be used to examine the representative logFC vs 
the combined p-values for all features of the selected contrasts. Users can use the 
lasso selector on the plot to select specific data points to be displayed in the 
results table below.

### 5. Ranked Expression

The ranked expression tab simply allows users to examine a table of all selected
contrasts ranked by the number of differentially expressed features (genes, TEs, 
or both) per contrast. This table is useful for quickly finding and comparing 
those contrasts where expression is most dysregulated.

## Differential Expression

The differential expression tab enables the exploration of individual 
experiment-level DE results as volcano or MA plots. The user can select features
(genes, TEs, or both) and create volcano plots of the differential expression
results for any contrast. Users can adjsut the FDR cutoff and logFC thresholds 
used in the plots. A table of differential expression results for the selected
contrast is displayed below the plots.

## GSEA

The GSEA tab enables interactive gene-set enrichment analysis for any selected 
contrast. All MSigDB gene sets are available for analysis. GSEA is performed 
using the `fgsea::fgsea()` function. After performing GSEA, an enrichment plot
is displayed for the most significant result by default. Users can use the gear 
icon in the top left of the plot to create an enrichment plot for any available 
contrast. A results table is displayed below the plot. 

NOTE: GSEA is performed only on the genes from the differential expression 
analysis. The ranking stat used for GSEA is the z-statistic which is computed 
from the moderated t-statistic of the differential expression results.

## Over-representation

The over-representation tab enables interactive over-representation analysis 
between sets of selected features using the `clusterProfiler::enrichGO()` 
function for the selected contrast. Users can select the ontology set to 
perform GO analysis on and an FDR cutoff for defining significant genes in the
study. After performing GO analysis, a dotplot and network graph of the 
significantly enriched terms is displayed for each gene set 
(up or down-regulated). 

NOTE: GO is performed only on the genes from the differential expression 
analysis. the background gene set used is comprised of all genes observed in 
that particular study.

## Sample vs Sample

The Sample-vs-Sample tab allows users to directly compare the number and 
direction of shared dysregulated features across collections of contrasts. Users 
can select the FDR cutoff used to determine significance. Significant feature for
each contrast are then displayed as an UpSet plot via `ComplexHeatmap::UpSet()`.
The number of shared features depending on the combination set mode for each 
pairwise combination is displayed on the plot. 
