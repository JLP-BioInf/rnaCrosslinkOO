#' @include  comradesDataSet.R 
NULL


##################################################
###############      METHODS       ###############
##################################################

##################################################
###############      Show          ###############
##################################################

 
setMethod("show", "comradesDataSet", function(object) {
    cat("comradesDataSet Object \n")
    cat("RNAs Analysed       - ",rnas(object), "\n")
    cat("Samples Analysed    - ",sampleNames(object), "\n")
    types = c()
    for(i in names(hybFiles(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Raw data            - ", types, "\n") 
    
    types = c()
    for(i in names(matrixList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Matrix Types        - ", types, "\n")
    
    types = c()
    for(i in names(clusterTableList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Cluster Types       - ", types, "\n")
    
    types = c()
    for(i in names(clusterGrangesList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Granges Types       - ", types, "\n")
    
    cat("Interactions        - ", dim(object@interactionTable), "\n")
    
    cat("Vienna Structures   - ", length(object@viennaStructures), "\n")
    

    
    
    
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
#' @param x A comradesDataSet object
#' @param data The data type to return <hybFiles | matrixList | clusterGrangesList | clusterTableList>
#' @param type The analysis stage <original | noHost | originalClusters | trimmedClusters> 
#' @name getData
#' @docType methods
#' @rdname getData
#' @aliases getData,comradesDataSet-method
#' @return A list of the chosen data type - one entry for each sample
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' getData(cds, 'matrixList','original')
#' 
#' @export
setGeneric("getData", function(x,data, type ) standardGeneric("getData"))
setMethod("getData", "comradesDataSet", function( x,data, type )  
    slot(x, data)[[rnas(x)]][[type]] )


###############  metaData        ###############

#' sampleNames
#' 
#' Extract the sample names for the instance
#' 
#' @param x A comradesDataSet object
#' @name sampleNames
#' @docType methods
#' @rdname sampleNames
#' @aliases sampleNames,comradesDataSet-method
#' @return A character vector - the sample names
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' sampleNames(cds)
#' @export
setGeneric("sampleNames", function(x ) standardGeneric("sampleNames"))
setMethod("sampleNames", "comradesDataSet", function( x)  
    as.character( sampleTable(x)$sampleName ) )

#' sampleTable
#' 
#' Extract the sample table for the instance
#' @param x A comradesDataSet object
#' @name sampleTable
#' @docType methods
#' @aliases sampleTable,comradesDataSet-method
#' @rdname sampleTable
#' @return A data frame - The orginal meta-data table 
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' sampleTable(cds)
#' @export
setGeneric("sampleTable", function(x) standardGeneric("sampleTable"))
setMethod("sampleTable", "comradesDataSet", function(x)   x@sampleTable)

#' group
#' 
#' Extract the indeces for each group for the instance
#' @param x A comradesDataSet object
#' @aliases group,comradesDataSet-method
#' @return A list - The indices of the sample in the control and sample groups
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' group(cds)
#' @export
setGeneric("group", function(x) standardGeneric("group"))
setMethod("group", "comradesDataSet", function(x)
    
   list("c" =  which(sampleTable(x)$group == "c"),
        "s" = which(sampleTable(x)$group== "s") )


)



#' rnas
#' 
#' Extract the rna ID for the instance
#' @param x A comradesDataSet object
#' @name rnas
#' @docType methods
#' @rdname rnas
#' @aliases rnas,comradesDataSet-method
#' @return A character - the ID of the RNA
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' rnas(cds)
#' @export
setGeneric("rnas", function(x) standardGeneric("rnas"))
setMethod("rnas", "comradesDataSet", function(x)  x@rnas)



#' rnaSize
#' 
#' Extract the size of the RNA for the instance
#' @param x A comradesDataSet object
#' @name rnaSize
#' @docType methods
#' @rdname rnaSize
#' @aliases rnaSize,comradesDataSet-method
#' @return A numeric - the size of the RNA (nucleotides)
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' rnaSize(cds)
#' @export
setGeneric("rnaSize", function(x) standardGeneric("rnaSize"))
setMethod("rnaSize", "comradesDataSet", function(x)  x@rnaSize)

#' hybFiles
#' 
#' Extract the data in original format 
#' @param x A comradesDataSet object
#' @name hybFiles
#' @docType methods
#' @rdname hybFiles
#' @aliases hybFiles,comradesDataSet-method
#' @return A list of tables in the original input format, one entry for each sample
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' hybFiles(cds)
#' @export
setGeneric("hybFiles", function(x) standardGeneric("hybFiles"))
setMethod("hybFiles", "comradesDataSet", function(x)   x@hybFiles)


#' matrixList
#' 
#' Extract the contact matrices
#' @param x A comradesDataSet object
#' @name matrixList
#' @docType methods
#' @rdname matrixList
#' @aliases matrixList,comradesDataSet-method
#' @return A list of contract matrices, one entry for each sample
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' matrixList(cds)
#' @export
setGeneric("matrixList", function(x) standardGeneric("matrixList"))
setMethod("matrixList", "comradesDataSet", function(x)   x@matrixList)






#' clusterGrangesList
#' 
#' Extract the cluster coordinates in granges format
#' @param x A comradesDataSet object
#' @name clusterGrangesList
#' @docType methods
#' @rdname clusterGrangesList
#' @aliases clusterGrangesList,comradesDataSet-method
#' @return A list of Granges objects showing the positions of each cluster, one entry for each sample
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' clusterGrangesList(cds)
#' @export
setGeneric("clusterGrangesList", function(x) standardGeneric("clusterGrangesList"))
setMethod("clusterGrangesList", "comradesDataSet", function(x)  x@clusterGrangesList)



#' clusterTableList
#' 
#' Extract the cluster coordinates in data frame format
#' @param x A comradesDataSet object
#' @name clusterTableList
#' @docType methods
#' @rdname clusterTableList
#' @aliases clusterTableList,comradesDataSet-method
#' @return A list of tables showing the vienna structures of each cluster
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' clusterTableList(cds)
#' @export
setGeneric("clusterTableList", function(x) standardGeneric("clusterTableList"))
setMethod("clusterTableList", "comradesDataSet", function(x)  x@clusterTableList)


#' clusterTableFolded
#' 
#' Extract the cluster coordinates with fold prediciton in data frame format
#' @param x A comradesDataSet object
#' @name clusterTableFolded
#' @docType methods
#' @rdname clusterTableFolded
#' @aliases clusterTableFolded,comradesDataSet-method
#' @return A table showing the vienna structures of each cluster
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' clusterTableFolded(cds)
#' @export
setGeneric("clusterTableFolded", function(x) standardGeneric("clusterTableFolded"))
setMethod("clusterTableFolded", "comradesDataSet", function(x)  x@clusterTableFolded)



##################################################
###############  Setters           ###############
##################################################


#' matrixList
#' 
#' Set new matrixList slot
#' @param x A comradesDataSet object
#' @param value A replacement
#' @name matrixList<-
#' @docType methods
#' @rdname matrixList_set
#' @aliases matrixList<-,comradesDataSet-method
#' @return No return - Sets a new matrixList slot
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' newMatrixList <- matrixList(cds)
#' matrixList(cds) <- newMatrixList
#' @export
setGeneric("matrixList<-", function(x, value) standardGeneric("matrixList<-"))
setMethod("matrixList<-", "comradesDataSet", function(x, value) {
    x@matrixList  = value
})

#' clusterGrangesList<-
#' 
#' Set new clusterGrangesList slot
#' @param x A comradesDataSet object
#' @param value A replacement
#' @name clusterGrangesList<-
#' @docType methods
#' @rdname clusterGrangesList_set
#' @aliases clusterGrangesList<-,comradesDataSet-method
#' @return No return - Sets a new clusterGrangesList slot
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' newclusterGrangesList <- clusterGrangesList(cds)
#' clusterGrangesList(cds) <- newclusterGrangesList
#' @export
setGeneric("clusterGrangesList<-", function(x, value) standardGeneric("clusterGrangesList<-"))
setMethod("clusterGrangesList<-", "comradesDataSet", function(x, value) {
    x@clusterGrangesList  = value
})

#' clusterTableList<-
#' 
#' Set new clusterTableList slot
#' 
#' @param x A comradesDataSet object
#' @param value A replacement
#' @name clusterTableList<-
#' @docType methods
#' @rdname clusterTableList_set
#' @aliases clusterTableList<-,comradesDataSet-method
#' @return No return - Sets a new clusterTableList slot
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' newclusterGrangesList <- clusterTableList(cds)
#' clusterTableList(cds) <- newclusterGrangesList
#' @export
setGeneric("clusterTableList<-", function(x, value) standardGeneric("clusterTableList<-"))
setMethod("clusterTableList<-", "comradesDataSet", function(x, value) {
    x@clusterTableList  = value
})









