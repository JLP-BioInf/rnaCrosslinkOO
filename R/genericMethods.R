#' @include  rnaCrosslinkDataSet.R 
NULL


##################################################
###############      METHODS       ###############
##################################################

##################################################
###############      Show          ###############
##################################################

 
setMethod("show", "rnaCrosslinkDataSet", function(object) {
    cat("rnaCrosslinkDataSet Object \n")
    cat("RNAs Analysed                -              - ",rnas(object), "\n")
    cat("Samples Analysed             -              - ",sampleNames(object), "\n")
    types = c()
    for(i in names(InputFiles(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("InputFiles           - Raw data               - ", types, "\n") 
    
    types = c()
    for(i in names(matrixList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("matrixList         - Contact Matricies      - ", types, "\n")
     
    types = c()
    for(i in names(clusterTableList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("clusterTableList   - Cluster positions      - ", types, "\n")
    
    types = c()
    for(i in names(clusterGrangesList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("clusterGrangesList - Granges of clusters    - ", types, "\n")
    
    cat("interactionTable   - Predictied interaction - ", dim(object@interactionTable)[1], "\n")
    
    cat("viennaStructures   - Predicted Structures   - ",
         length(object@viennaStructures), sep = "", 
        "\n")
    

    
    
    
})


##################################################
###############  Accessors         ###############
##################################################



#' getData
#' 
#' Get data is more generic method for retrieving data from the object
#' and returns a list, the number of entries in the list is number of
#' samples in the dataset and the list contain entries of the data type
#' and analysis stage you select.
#' 
#' @param x A rnaCrosslinkDataSet object
#' @param data The data type to return <InputFiles | matrixList | clusterGrangesList | clusterTableList>
#' @param type The analysis stage <original | noHost | originalClusters | trimmedClusters> 
#' @name getData
#' @docType methods
#' @rdname getData
#' @aliases getData,rnaCrosslinkDataSet-method
#' @return A list of the chosen data type - one entry for each sample
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' getData(cds, 'matrixList','original')
#' 
#' @export
setGeneric("getData", function(x,data, type ) standardGeneric("getData"))
setMethod("getData", "rnaCrosslinkDataSet", function( x,data, type )  
    slot(x, data)[[rnas(x)]][[type]] )


###############  metaData        ###############

#' sampleNames
#' 
#' Extract the sample names for the instance
#' 
#' @param x A rnaCrosslinkDataSet object
#' @name sampleNames
#' @docType methods
#' @rdname sampleNames
#' @aliases sampleNames,rnaCrosslinkDataSet-method
#' @return A character vector - the sample names
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' sampleNames(cds)
#' @export
setGeneric("sampleNames", function(x ) standardGeneric("sampleNames"))
setMethod("sampleNames", "rnaCrosslinkDataSet", function( x)  
    as.character( sampleTable(x)$sampleName ) )

#' sampleTable
#' 
#' Extract the sample table for the instance
#' @param x A rnaCrosslinkDataSet object
#' @name sampleTable
#' @docType methods
#' @aliases sampleTable,rnaCrosslinkDataSet-method
#' @rdname sampleTable
#' @return A data frame - The orginal meta-data table 
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' sampleTable(cds)
#' @export
setGeneric("sampleTable", function(x) standardGeneric("sampleTable"))
setMethod("sampleTable", "rnaCrosslinkDataSet", function(x)   
  x@sampleTable)

#' group
#' 
#' Extract the indeces for each group for the instance
#' @param x A rnaCrosslinkDataSet object
#' @aliases group,rnaCrosslinkDataSet-method
#' @return A list - The indices of the sample in the control and sample groups
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' group(cds)
#' @export
setGeneric("group", function(x) standardGeneric("group"))
setMethod("group", "rnaCrosslinkDataSet", function(x)
   list("c" =  which(sampleTable(x)$group == "c"),
        "s" = which(sampleTable(x)$group== "s") )
)


#' rnas
#' 
#' Extract the rna ID for the instance
#' @param x A rnaCrosslinkDataSet object
#' @name rnas
#' @docType methods
#' @rdname rnas
#' @aliases rnas,rnaCrosslinkDataSet-method
#' @return A character - the ID of the RNA
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' rnas(cds)
#' @export
setGeneric("rnas", function(x) standardGeneric("rnas"))
setMethod("rnas", "rnaCrosslinkDataSet", function(x) 
  x@rnas)



#' rnaSize
#' 
#' Extract the size of the RNA for the instance
#' @param x A rnaCrosslinkDataSet object
#' @name rnaSize
#' @docType methods
#' @rdname rnaSize
#' @aliases rnaSize,rnaCrosslinkDataSet-method
#' @return A numeric - the size of the RNA (nucleotides)
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' rnaSize(cds)
#' @export
setGeneric("rnaSize", function(x) standardGeneric("rnaSize"))
setMethod("rnaSize", "rnaCrosslinkDataSet", function(x) 
  x@rnaSize)

#' InputFiles
#' 
#' Extract the data in original format 
#' @param x A rnaCrosslinkDataSet object
#' @name InputFiles
#' @docType methods
#' @rdname InputFiles
#' @aliases InputFiles,rnaCrosslinkDataSet-method
#' @return A list of tables in the original input format, one entry for each sample
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' InputFiles(cds)
#' @export
setGeneric("InputFiles", function(x) standardGeneric("InputFiles"))
setMethod("InputFiles", "rnaCrosslinkDataSet", function(x)   x@InputFiles)


#' matrixList
#' 
#' Extract the contact matrices
#' @param x A rnaCrosslinkDataSet object
#' @name matrixList
#' @docType methods
#' @rdname matrixList
#' @aliases matrixList,rnaCrosslinkDataSet-method
#' @return A list of contract matrices, one entry for each sample
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' matrixList(cds)
#' @export
setGeneric("matrixList", function(x) standardGeneric("matrixList"))
setMethod("matrixList", "rnaCrosslinkDataSet", function(x)   x@matrixList)






#' clusterGrangesList
#' 
#' Extract the cluster coordinates in granges format
#' @param x A rnaCrosslinkDataSet object
#' @name clusterGrangesList
#' @docType methods
#' @rdname clusterGrangesList
#' @aliases clusterGrangesList,rnaCrosslinkDataSet-method
#' @return A list of Granges objects showing the positions of each cluster, one entry for each sample
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusterGrangesList(cds)
#' @export
setGeneric("clusterGrangesList", function(x) standardGeneric("clusterGrangesList"))
setMethod("clusterGrangesList", "rnaCrosslinkDataSet", function(x)  x@clusterGrangesList)



#' clusterTableList
#' 
#' Extract the cluster coordinates in data frame format
#' @param x A rnaCrosslinkDataSet object
#' @name clusterTableList
#' @docType methods
#' @rdname clusterTableList
#' @aliases clusterTableList,rnaCrosslinkDataSet-method
#' @return A list of tables showing the vienna structures of each cluster
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusterTableList(cds)
#' @export
setGeneric("clusterTableList", function(x) standardGeneric("clusterTableList"))
setMethod("clusterTableList", "rnaCrosslinkDataSet", function(x)  x@clusterTableList)


#' clusterTableFolded
#' 
#' Extract the cluster coordinates with fold prediciton in data frame format
#' @param x A rnaCrosslinkDataSet object
#' @name clusterTableFolded
#' @docType methods
#' @rdname clusterTableFolded
#' @aliases clusterTableFolded,rnaCrosslinkDataSet-method
#' @return A table showing the vienna structures of each cluster
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusterTableFolded(cds)
#' @export
setGeneric("clusterTableFolded", function(x) standardGeneric("clusterTableFolded"))
setMethod("clusterTableFolded", "rnaCrosslinkDataSet", function(x)  x@clusterTableFolded)



##################################################
###############  Setters           ###############
##################################################


#' matrixList
#' 
#' Set new matrixList slot
#' @param x A rnaCrosslinkDataSet object
#' @param value A replacement
#' @name matrixList<-
#' @docType methods
#' @rdname matrixList_set
#' @aliases matrixList<-,rnaCrosslinkDataSet-method
#' @return No return - Sets a new matrixList slot
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' newMatrixList <- matrixList(cds)
#' matrixList(cds) <- newMatrixList
#' @export
setGeneric("matrixList<-", function(x, value) standardGeneric("matrixList<-"))
setMethod("matrixList<-", "rnaCrosslinkDataSet", function(x, value) {
    x@matrixList  = value
})

#' clusterGrangesList<-
#' 
#' Set new clusterGrangesList slot
#' @param x A rnaCrosslinkDataSet object
#' @param value A replacement
#' @name clusterGrangesList<-
#' @docType methods
#' @rdname clusterGrangesList_set
#' @aliases clusterGrangesList<-,rnaCrosslinkDataSet-method
#' @return No return - Sets a new clusterGrangesList slot
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' newclusterGrangesList <- clusterGrangesList(cds)
#' clusterGrangesList(cds) <- newclusterGrangesList
#' @export
setGeneric("clusterGrangesList<-", function(x, value) standardGeneric("clusterGrangesList<-"))
setMethod("clusterGrangesList<-", "rnaCrosslinkDataSet", function(x, value) {
    x@clusterGrangesList  = value
})

#' clusterTableList<-
#' 
#' Set new clusterTableList slot
#' 
#' @param x A rnaCrosslinkDataSet object
#' @param value A replacement
#' @name clusterTableList<-
#' @docType methods
#' @rdname clusterTableList_set
#' @aliases clusterTableList<-,rnaCrosslinkDataSet-method
#' @return No return - Sets a new clusterTableList slot
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' newclusterGrangesList <- clusterTableList(cds)
#' clusterTableList(cds) <- newclusterGrangesList
#' @export
setGeneric("clusterTableList<-", function(x, value) standardGeneric("clusterTableList<-"))
setMethod("clusterTableList<-", "rnaCrosslinkDataSet", function(x, value) {
    x@clusterTableList  = value
})









