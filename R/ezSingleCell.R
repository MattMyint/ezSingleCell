#' ezSingleCell: Interactive single-cell data analysis using the Seurat pipeline
#' 
#' This package was made to provide a GUI for carrying out single cell data analysis on Seurat, as 
#' opposed to running the entire pipeline from the command line. Analysis using the package starts after alignment
#' and quantification of counts from raw reads, taking a sparse matrix of expression values.
#' 
#' Data loading
#'   Upon loading count data, counts are log normalised, and filtered based on user input of minumum cut-offs,
#'   as well as expression thresholds. Users can also set cell population identities based on the formatting
#'   of cell names in their expression table.
#'  
#' Quality check plots
#'   Immediately after counts are loaded, plots visualising metadata such as nUMI and percentage of mitochondrial
#'   genes are generated to allow the user to determine if any cells are low quality and need to be filtered out
#'   from downstream analysis.
#' 
#' Variable gene identification
#'   Variable genes are identified with user-selected dispersion and mean cut-offs. Clicking "Find Variable
#'   genes" will return the number of variable genes identified for downstream analysis. This is sufficient to
#'   proceed with later steps (i.e. plotting the graph is not necessary). It should be noted that this step
#'   needs to be done so that dimension reduction and clustering can be done later on. 
#' 
#' PCA 
#'   PCA will be run on identified variable genes, and users can visualise 2D plots of selected PCs.
#'   This can provide better visualisation of any outliers. After PCA is done, users can search for clusters.
#'   
#' Diagnostic analysis of PCs with Jackstraw and Elbow plots
#'   Seurat's tSNE clustering outcomes are notably dependent on the PCs used, so users should take time to 
#'   determine which PCs are significant for more accurate results.
#' 
#' tSNE dimension reduction
#'   tSNE will collapse the chosen PCs into lower dimensions and provide additional visualisation of clustering
#'   of cell populations.
#' 
#' Differentially Expressed Gene (DEG) analysis
#'   Users can identify differentially expressed genes across groups, either as a one-vs-all comparison, or a
#'   1-on-1 comparison between 2 selected groups.
#' 
#' The generated figures can be saved to PDF and CSV files, and if the user wishes to export their current analysis,
#' they can use the "Save Data" button to save the Seurat object as an .RObj file.
#' 
#' @author Matthew Myint
#' @references placeholder
#' @docType package
#' @name ezSingleCell-package
#' 
NULL

#' Launch shinyApp
#'
#' Use this function to run the shinyApp
#' @keywords shiny
#' @return Launches the shiny app
#' @import shiny
#' @author Matthew Myint
#' @export
#' @examples
#' ezSingleCell()

ezSingleCell <- function(){
  appDir <- system.file('shiny', package = "ezSingleCell")
  shiny::runApp(appDir, display.mode = "normal")
}