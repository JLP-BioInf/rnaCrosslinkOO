#' @include  comradesDataSet.R
NULL


#' clusterNumbers
#'
#' This method prints a table showing the number of clusters in each step 
#' of the analysis
#'
#' @param knowClusteredCds A comradesDataSet object after clustering has been performed
#' @param rna The RNA ID of interest - use rna(cdsObject).
#' @name clusterNumbers
#' @docType methods
#' @rdname clusterNumbers
#' @aliases clusterNumbers,comradesDataSet-method
#' @return A data.frame shoing the number of clusters for each sample
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' clusteredCds = clusterComrades(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 1)
#' clusterNumbers(clusteredCds)
#' @export
setGeneric("clusterNumbers",
           function(knowClusteredCds,
                    rna)
             standardGeneric("clusterNumbers"))

setMethod("clusterNumbers", "comradesDataSet", function(knowClusteredCds,
                                                        rna)  {
  sampleTable = sampleTable(knowClusteredCds)
  clusters = clusterTableList(knowClusteredCds)
  
  
  
  
  for (rna in rnas(knowClusteredCds)) {
    sampleTable$rna = rna
    for (type in 1:length(clusters[[rna]])) {
      typeID = names(clusters[[rna]])[type]
      sampleTable[, typeID] = 0
      
      v = c()
      for (sample in 1:length(clusters[[rna]][[type]])) {
        sampleID = sampleNames(knowClusteredCds)[sample]
        
        v = c(v, as.numeric(nrow(clusters[[rna]][[type]][[sample]])))
        
        
      }
      sampleTable[, typeID] = v
      
    }
  }
  return(sampleTable)
  
})





#' readNumbers
#'
#'
#' This method prints a table showing the number of duplexes in 
#' the clusters in each step of the analysis
#'
#' @param knowClusteredCds A comradesDataSet object after clustering has been performed
#' @param rna The RNA ID of interest - use rna(cdsObject).
#'
#' @name readNumbers
#' @docType methods
#' @rdname readNumbers
#' @aliases readNumbers,comradesDataSet-method
#' @return A data.frame shoing the number of reads in clusters for each sample
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' clusteredCds = clusterComrades(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 1)
#' readNumbers(clusteredCds)
#' @export
#'
setGeneric("readNumbers",
           function(knowClusteredCds,
                    rna)
             standardGeneric("readNumbers"))

setMethod("readNumbers", "comradesDataSet", function(knowClusteredCds,
                                                     rna)  {
  sampleTable = sampleTable(knowClusteredCds)
  clusters = hybFiles(knowClusteredCds)
  #get numbers of original hyb files
  for (rna in c(rnas(knowClusteredCds), "all")) {
    sampleTable$rna = rnas(knowClusteredCds)
    for (type in 1:length(clusters[[rna]])) {
      typeID = names(clusters[[rna]])[type]
      sampleTable[, typeID] = 0
      
      v = c()
      for (sample in 1:length(clusters[[rna]][[type]])) {
        sampleID = sampleNames(knowClusteredCds)[sample]
        
        v = c(v, as.numeric(nrow(clusters[[rna]][[type]][[sample]])))
        
        
      }
      sampleTable[, typeID] = v
      
    }
  }
  
  if (class(knowClusteredCds)[1] == "comradesDataSet") {
    clusters = clusterGrangesList(knowClusteredCds)
    
    for (rna in rnas(knowClusteredCds)) {
      sampleTable$rna = rnas(knowClusteredCds)
      for (type in 1:length(clusters[[rna]])) {
        typeID = names(clusters[[rna]])[type]
        sampleTable[, typeID] = 0
        
        v = c()
        for (sample in 1:length(clusters[[rna]][[type]])) {
          sampleID = sampleNames(knowClusteredCds)[sample]
          if (length(clusters[[rna]][[type]][[sample]]) == 2) {
            v = c(v, as.numeric(length(unlist(
              clusters[[rna]][[type]][[sample]]
            ))))
          } else{
            v = c(v, as.numeric(length((
              clusters[[rna]][[type]][[sample]]
            ))))
          }
          
        }
        sampleTable[, typeID] = v
        
      }
    }
  }
  return(sampleTable)
  
})















#  plotMatrices
#'
#'
#' Plots a number of contact maps to file of each sample for a stage in the analysis
#'
#' @param cds A comradesDataSet object 
#' @param type The analysis stage to plot 
#' @param directory An output directory for the heatmap (use 0 for no output)
#' @param a To make a subsetted plot (left value on x)
#' @param b To make a subsetted plot (right value on x)
#' @param c To make a subsetted plot (left value on y)
#' @param d To make a subsetted plot (right value on y)
#' @param h Height of image (inches) ( only useful if plotting)
#' @name plotMatrices
#' @docType methods
#' @rdname plotMatrices
#' @aliases plotMatrices,comradesDataSet-method
#' @return A heatmap of the reads in the analysis stage chosen
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' plotMatrices(cds,
#'             b = rnaSize(cds),
#'             d = rnaSize(cds))
#' @export
setGeneric("plotMatrices", function(cds,
                                    type = 'original', 
                                    directory = 0,
                                    a = 1,
                                    b = 50,
                                    c = 1,
                                    d = 50,
                                    h= 3)
  standardGeneric("plotMatrices"))

setMethod("plotMatrices", "comradesDataSet", function(cds, 
                                                      type = 'original', 
                                                      directory = 0, 
                                                      a = 1,
                                                      b = 50,
                                                      c = 1,
                                                      d = 50,
                                                      h= 3)  {
  hybMatList = matrixList(cds)
  rnaS = rnas(cds)
  sampleNames = names(hybMatList[[1]][[type]])
  if (is.null(sampleNames)) {
    sampleNames = 1:length(hybMatList[[1]][[type]])
  }
  
  
  
  for (sample in  1:length(sampleNames)) {
    for (rna in c(rnaS)) {
      cols = log2(max(hybMatList[[rna]][[type]][[sample]][a:b, c:d] + 1))
      
      myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols -
                                                                        1))
      
      if (cols > 14) {
        cols = log2(max(hybMatList[[rna]][[type]][[sample]][a:b, c:d][hybMatList[[rna]][[type]][[sample]][a:b, c:d] < 30000] +
                          1))
        myCol = c("black",
                  colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1),
                  rep("white", (14 - cols) + 1))
      }
      
      if (directory == 0) {
        heatmap3((log2(t(
          hybMatList[[rna]][[type]][[sample]][a:b, c:d] + 1
        ))),
        col = myCol,
        scale = "none" ,
        Rowv = NA,
        Colv = NA,
        useRaster = T
        )
      } else{
        pdf(
          paste(
            directory,
            "/",
            rna ,
            "_",
            sampleNames[sample],
            "-",
            type ,
            ".pdf",
            sep = ""
          ),
          height = h,
          width = h
        )
        heatmap3((log2(t(
          hybMatList[[rna]][[type]][[sample]][a:b, c:d] + 1
        ))),
        col = myCol,
        scale = "none" ,
        Rowv = NA,
        Colv = NA,
        useRaster = T
        )
        dev.off()
      }
      
    }
  }
  
  
})



#'  plotMatricesAverage
#'
#'
#' Plots a contact map to file of all samples for a stage in the analysis
#'
#' @param cds A comradesDataSet object 
#' @param type The analysis stage to plot 
#' @param directory An output directory for the heatmap (use 0 for no output)
#' @param a To make a subsetted plot (left value on x)
#' @param b To make a subsetted plot (right value on x)
#' @param c To make a subsetted plot (left value on y)
#' @param d To make a subsetted plot (right value on y)
#' @param h Height of image (inches) ( only useful if plotting)
#' 
#' @name plotMatricesAverage
#' @docType methods
#' @rdname plotMatricesAverage
#' @aliases plotMatricesAverage,comradesDataSet-method
#' @return A heatmap of the reads in the analysis stage chosen
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' plotMatricesAverage(cds,
#'             b = rnaSize(cds),
#'             d = rnaSize(cds))
#'
#' @export
setGeneric("plotMatricesAverage", function(cds, type, directory, a, b, c, d, h)
  standardGeneric("plotMatricesAverage"))

setMethod("plotMatricesAverage", "comradesDataSet", function(cds, type, directory, a, b, c, d, h)  {
  for (rna in rnas(cds)) {
    hybMatList = matrixList(cds)
    hybMatList2 = hybMatList
    for (i in c("c", "s")) {
      c = 1
      for (j in group(cds)[[i]]) {
        if (length(group(cds)[[i]]) < 2 | c == 1) {
          sum(hybMatList[[rna]][[type]][[j]])
          hybMatList2[[rna]][[type]][[i]] =   hybMatList[[rna]][[type]][[j]]
          
          
          # hybMatList2[[rna]][[type]][["s"]] =   hybMatList[[rna]][[type]][[ j ]]
        } else{
          hybMatList2[[rna]][[type]][[i]] =  hybMatList2[[rna]][[type]][[i]] +
            hybMatList[[rna]][[type]][[j]]
          
          
          # hybMatList2[[rna]][[type]][["s"]] =  Reduce('+', hybMatList[[rna]][[type]][[ j ]] )
        }
        c = c + 1
      }
      
    }
    
    sampleNames = c("s", "c")
    #c = 1
    for (sample in  sampleNames) {
      cols = log2(max(hybMatList2[[rna]][[type]][[sample]][a:b, c:d] + 1))
      
      myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols -
                                                                        1))
      
      if (cols > 14) {
        cols = log2(max(hybMatList2[[rna]][[type]][[sample]][a:b, c:d][hybMatList2[[rna]][[type]][[sample]][a:b, c:d] < 30000] +
                          1))
        myCol = c("black",
                  colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1),
                  rep("white", (15 - cols) + 2))
      }
      
      # myCol = myCol[1:cols]
      
      
      if (directory == 0) {
        heatmap3((log2(t(
          hybMatList2[[rna]][[type]][[sample]][a:b, c:d] + 1
        ))),
        col = myCol,
        scale = "none" ,
        Rowv = NA,
        Colv = NA,
        useRaster = T
        )
        
      } else{
        pdf(
          paste(directory, "/", rna , "_", sample , "-", type , ".pdf", sep = ""),
          height = h,
          width = h
        )
        heatmap3((log2(t(
          hybMatList2[[rna]][[type]][[sample]][a:b, c:d] + 1
        ))),
        col = myCol,
        scale = "none" ,
        Rowv = NA,
        Colv = NA,
        useRaster = T
        )
        dev.off()
      }
      
      
      
      #  c = c + 1
    }
  }
  
  
})


