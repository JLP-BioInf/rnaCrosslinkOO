#' @include rnaCrosslinkOO.R


#' rnaCrosslinkDataSet
#'
#' An S4 class to represent a rnaCrosslink dataset
#' @rdname rnaCrosslinkDataSet
#' @export
#'
setClass(
  "rnaCrosslinkDataSet",
  slots = c(
    rnas = "character",
    rnaSize = "numeric",
    sampleTable = "data.frame",
    # meta data
    InputFiles = "list",
    # data tables
    matrixList = "list",
    clusterGrangesList = "list",
    clusterTableList = "list",
    clusterTableFolded = "data.frame",
    interactionTable = "data.frame",
    viennaStructures = "list",
    dgs = "list"
  ),
  prototype = list()
)

setValidity("rnaCrosslinkDataSet", function(object) {
  
})




#' rnaCrosslinkDataSet
#'
#' \code{rnaCrosslinkDataSet} objects are used to store the input meta-data, data and
#' create a framework for the storage of results. Whilst creating the object,
#' the original Input files are also filtered for the RNA of interest. Check the 
#' package vignette for more information.
#'
#'
#' @param rnas vector - The names of the RNA interest, these must be displayed
#' the same way as in the input Input Files.
#' @param rnaSize named list - The sizes (nt) of the RNAs of interest, the list
#'  elements must have same names as the \code{rnas} vector and each each contain
#'  one numeric value.
#' @param sampleTable string - The address of the sample table, the sample table
#'  must have 4 columns, fileName (the full path and file name of the input
#'  Input file for each sample ), group ("s" - sample or "c" - control),
#'  sample (1,2,3, etc), sampleName (must be unique).
#'
#' @return A rnaCrosslinkDataSet object.
#'
#'
#'
#' @slot clusterTableFolded table - a table similar to the \code{clusterTableList}
#' it contains coordinates of the clusters along with vienna format fold and
#' RNA sequences for each cluster
#' @slot clusterTableList List - Follows the pattern for list slots of rnaCrosslinkDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. contains a table
#' with coordinates and information about the clusters identified
#' @slot clusterGrangesList List - Follows the pattern for list slots of rnaCrosslinkDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. contains GRanges
#' objects of the original duplexes with their cluster membership
#' @slot sampleTable table - Column names; fileName, group (s or c),
#'  sample (1,2,3, etc), sampleName (must be unique)
#' @slot rnas string - a single RNA to analyse - must be present in \code{rnas(cdsObject)}
#' @slot rnaSize if set to 0 this will be calculated 
#' @slot matrixList List - Follows the pattern for list slots of rnaCrosslinkDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. Contains a set
#' of contact matrices, each cell contains the number of duplexes identified
#' for position x,y.
#' @slot InputFiles List - Follows the pattern for list slots of rnaCrosslinkDataSet
#' objects, \code{InputFiles(cds)[[rna]][[type]][[sample]]}. Contains a set of
#' tables, these are the original Input files that were read in.
#' @slot interactionTable Table of interactions discovered in step1 of the folding
#' @slot viennaStructures List of vienna format structures from final prediction
#' @slot dgs List of free energies
#'
#' @name rnaCrosslinkDataSet
#' @docType class
#' @rdname rnaCrosslinkDataSet
#'
#' @examples 
#' # make example input
#' cds = makeExamplernaCrosslinkDataSet()
#'
#'cds
#'
#' @export
rnaCrosslinkDataSet <- function(rnas,
                            rnaSize = 0 ,
                            sampleTable) {
  ###########################################################
  # Read in the sample table
  ###########################################################
  # check the inputs here, stop if wrong
  message(" ******************************************** ")
  message(" *****            rnaCrosslink-OO          ****** ")
  message(" ******************************************** ")
  message(" *****-------*******************-------****** ")
  message(" *****       Reading SampleTable       ****** ")
  
  
  # Read in sample table

  #check for more than two samples
  # if( nrow(sampleTable) < 2 ){
  #        stop( "The sample Table must contain at least 1
  #              sample and 1 control" )
  #    }
  message(paste(" *****       Detected ",
                nrow(sampleTable), " Samples      ****** "))
  
  
  
  ###########################################################
  #check column names of sampleTable
  colnamesST = c("file", "group", "sample", "sampleName")
  if (all(colnames(sampleTable) != colnamesST)) {
    stop("Column names of metaData table should be :
              file, group, sample, sampleNames")
  }
  
  # Get the comparison groups check group has the c and s
  if (!(unique(as.character(sampleTable$group))[1] %in% c("c", "s") &
        unique(as.character(sampleTable$group))[2] %in% c("c", "s"))) {
    stop("Groups should be c and s")
  }
  
  
  # Make group into a list with control and sample
  group = sampleTable[, "group"]
  group2 = list()
  group2[["c"]] = which(group == "c")
  group2[["s"]] = which(group == "s")
  
  group = group2
  message(paste(
    " *****     detected group c::",
    paste(group[["c"]],
          collapse = " ") ,
    paste(rep(" ",
              (
                length(group[["c"]]) * (3 - length(group[["c"]]))) * 2),
          collapse = ""),
    "   ***** "
  ))
  message(paste(
    " *****     detected group s::",
    paste(group[["s"]],
          collapse = " ") ,
    paste(rep(" ",
              (
                length(group[["s"]]) * (3 - (length(group[["s"]]))) * 2)),
          collapse = ""),
    "   ***** "
  ))
  
  
  
  ###########################################################
  # Get the sampleNames
  sampleNames = c()
  if (is.null(sampleTable$sampleName)) {
    stop("The sample Table must have a column named sampleName")
  } else if (length(unique(sampleTable$sampleName)) !=
             length(sampleTable$sampleName)) {
    stop("Sample names must be unique")
  } else{
    sampleNames = as.character(sampleTable$sampleName)
    spaces =  (length(sampleNames) * (3 - length(sampleNames))) * 2 
    if(spaces < 0 ){spaces = 0}
    message(paste(
      " ****  ",
      paste(rep(" ",
                spaces,
            collapse = ""),
      " Sample Names: ",
      paste(sampleNames, collapse = " "),
      paste(rep(" ",
                spaces,
            collapse = ""),
      " **** "
    ))))
  }
  
  
  
  ###########################################################
  # Read in the  Input files
  ###########################################################
  #load the files into a list
  message(" *****         Reading Input Files        ***** ")
  
  InputFiles = list()
  InputFiles[["all"]] = list()
  InputFiles[["all"]][["all"]] = list()
  
  
  
  # read in the tables
  inputs <- lapply(as.character(sampleTable$file),
                   function(file)
                     read.table(file,colClasses = c("character",
                                                    "character",
                                                    "character",
                                                    "character",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "character",
                                                    "character",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "character")))
  
  #check column names
  if (all(sapply(inputs, function(file)
    ! (identical(
      colnames(file),
      c(
        "V1",
        "V2",
        "V3",
        "V4",
        "V5",
        "V6",
        "V7",
        "V8",
        "V9",
        "V10",
        "V11",
        "V12",
        "V13",
        "V14",
        "V15"
      )
    ))))) {
    stop(
      " The input files do not look they are produced with the
                 nextflow pipeline, please check the documentation. "
    )
  }
  
  
  checkRNA = function(inputs, rna) {
    rnas %in% c(inputs$V4, inputs$V10)
    
  }  
  rnaCheck = sapply(inputs, checkRNA, rna = rnas)
  if(!any(rnaCheck)){
    stop(
      "No match for chosen RNA exists in this dataSet"
    )
  }
  
  inputs = lapply(inputs, function(x) x[complete.cases(x),] )
  inputs = lapply(inputs, function(x) x[!is.infinite(x$V8),] )
  inputs = lapply(inputs, function(x) x[!is.infinite(x$V14),] )

  
  InputFiles[["all"]][["all"]] = inputs
  names(InputFiles[["all"]][["all"]]) = sampleNames
  
  
  #get the maxpos
  getMax = function(file, rna) {
    return(max(file[file$V10 == rnas, "V14"],
               max(file[file$V4 == rnas, "V8"])))

  }
 
  
  maxPos = max(sapply(inputs, getMax, rna = rnas))
  

  
  
  
  
  ###########################################################
  message(" *****     Getting RNAs of Interest    ****** ")
  message(" *****    RNA of interest + Host RNA    ***** ")
  InputFiles[[rnas]][["original"]] = swapInputs(InputList = InputFiles[["all"]][["all"]],
                                            rna = rnas)
  InputFiles[[rnas]][["host"]] = swapInputs3(InputList = InputFiles[["all"]][["all"]],
                                         rna = rnas)
  message(" *****      RNA of interest Alone       ***** ")
  InputFiles[[rnas]][["noHost"]] = swapInputs2(InputList = InputFiles[["all"]][["all"]],
                                           rna = rnas)
  
  
  
  ###########################################################
  # Make matrices of the specific RNA without host
  ###########################################################
  message(" *****         Making Matrices         ****** ")
  
  matrixList = list()
  matrixList[[rnas]] = list()
  if (rnaSize > 0) {
    rnaSize2 =   rnaSize
  } else{
    rnaSize2 = maxPos
  }
  
  message(paste(" *****          RNA Size: ",
                rnaSize2, "        ***** "))
  matrixList[[rnas]][["noHost"]] = list()
  matrixList[[rnas]][["noHost"]] = getMatrices(InputFiles[[rnas]][["noHost"]],
                                               rnas, rnaSize2)
  names(matrixList[[rnas]][["noHost"]]) = sampleNames
  matrixList[[rnas]][["original"]] = list()
  matrixList[[rnas]][["original"]] = getMatrices(InputFiles[[rnas]][["original"]],
                                                 rnas, rnaSize2)
  names(matrixList[[rnas]][["original"]]) = sampleNames
  
  
  
  
  ###########################################################
  message(" *****         Creating object          ***** ")
  message(" *****-------*******************-------****** ")
  message(" ******************************************** ")
  message(" ******************************************** ")
  #create rnaCrosslink dataset object
  object  = new(
    "rnaCrosslinkDataSet",
    rnas = rnas,
    rnaSize = rnaSize2,
    sampleTable = sampleTable,
    InputFiles = InputFiles,
    matrixList = matrixList,
    clusterGrangesList = list(),
    clusterTableList = list(),
    clusterTableFolded = data.frame(),
    interactionTable = data.frame(),
    viennaStructures = list(),
    dgs = list()
  )
  
  return(object)
  
  
}
