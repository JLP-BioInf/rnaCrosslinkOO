
#' comradesOO: A package for analysing a COMRADES dataset
#'
#' The package consists of 1 Object:
#'    comradesDataSet
#'
#' @section \code{\link{comradesDataSet}}
#'
#'
#' @docType package
#' @name comradesOO
#' @keywords internal
"_PACKAGE"


globalVariables(c(    "ID" ,"PCa" ,"PCb" ,"clusteredCds",  "dgs" ,"k", 
                      "sampleTable2" ,"sd" ,"value","stopCluster"))

#' @importFrom seqinr read.fasta
#' @importFrom  stats aggregate median complete.cases prcomp reorder
#' @import ggplot2
#' @import reshape2
#' @import MASS
#' @import mixtools
#' @import utils
#' @importFrom doParallel stopImplicitCluster registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom igraph cluster_walktrap graph_from_adjacency_matrix membership
#' @import R4RNA
#' @import RColorBrewer
#' @import heatmap3 
#' @import TopDom
#' @import tidyverse
#' @import RRNA
#' @import ggrepel
#' @import foreach
#' @import GenomicRanges
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges IRanges
#' @importFrom methods new slot
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @aliases comradesOOPackage

## usethis namespace: start
## usethis namespace: end
NULL
