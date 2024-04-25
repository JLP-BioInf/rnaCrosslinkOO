#' rnaCrosslinkOO: A package for analysing a rnaCrosslink dataset
#'
#' The package consists of 1 Object:
#'    rnaCrosslinkDataSet
#'
#' @section \code{\link{rnaCrosslinkDataSet}}
#'
#'
#' @docType package
#' @name rnaCrosslinkOO
#' @keywords internal
"_PACKAGE"


globalVariables(c( "count",   "ID" ,"PCa" ,"PCb" ,"clusteredCds",  "dgs" ,"k", 
                      "sampleTable2" ,"sd" ,"value","stopCluster","variable"))

#' @importFrom seqinr read.fasta
#' @importFrom  stats aggregate median complete.cases prcomp reorder
#' @import ggplot2
#' @import reshape2
#' @rawNamespace import("MASS", except = area)
#' @import mixtools
#' @import utils
#' @importFrom doParallel stopImplicitCluster registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom igraph cluster_walktrap graph_from_adjacency_matrix membership
#' @import R4RNA
#' @import RColorBrewer
#' @import heatmap3 
#' @import TopDom
#' @import tidyverse
#' @import RRNA
#' @import patchwork 
#' @import ggrepel
#' @import foreach
#' @import GenomicRanges
#' @importFrom ClassDiscovery aspectHeatmap
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges IRanges
#' @importFrom methods new slot
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @aliases rnaCrosslinkOOPackage

## usethis namespace: start
## usethis namespace: end
NULL
