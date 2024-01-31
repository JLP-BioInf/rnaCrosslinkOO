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
#' Plots a contact map of all samples for two chosen stages in the analysis, with each chosen stage on separate halves of the contact map
#'
#' @param cds A comradesDataSet object 
#' @param type1 The analysis stage to plot on the upper half of the heatmap (use 'blank' for to leave this half blank)
#' @param type2 The analysis stage to plot on the lower half of the heatmap (use 'blank' for to leave this half blank)
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
#' @return A heatmap of the reads in the two analysis stages chosen, with each chosen stage on a separate half of the heatmap
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' plotMatricesAverage(cds,
#'             b = rnaSize(cds),
#'             d = rnaSize(cds))
#'
#' @export
setGeneric("plotMatricesAverage", function(cds,
                                           type1 = 'original',
                                           type2 = 'blank',
                                           directory = 0,
                                           a = 1,
                                           b = 50,
                                           c = 1,
                                           d = 50,
                                           h= 3)
  standardGeneric("plotMatricesAverage"))

setMethod("plotMatricesAverage", "comradesDataSet", function(cds,
                                                             type1 = 'original',
                                                             type2 = 'blank',
                                                             directory = 0,
                                                             a = 1,
                                                             b = 50,
                                                             c = 1,
                                                             d = 50,
                                                             h= 3)  {
  if (type2 != 'blank') {
    typeCombinationName = paste(type1, type2, sep = "-")
  } else {
    typeCombinationName = type1
  }
  for (rna in rnas(cds)) {
    hybMatList = matrixList(cds)
    hybMatList2 = hybMatList
    for (i in c("c", "s")) {
      for (matrixhalfNumber in 1:2){
        if (matrixhalfNumber == 1){
          type = type1
        } else{
          type = type2
        }
        if (type != "blank"){
          c = 1
          for (j in group(cds)[[i]]) {
            if (length(group(cds)[[i]]) < 2 | c == 1) {
              sum(hybMatList[[rna]][[type]][[j]])
              hybMatList2[[rna]][[type]][[i]] = hybMatList[[rna]][[type]][[j]]
            } else{
              hybMatList2[[rna]][[type]][[i]] = hybMatList2[[rna]][[type]][[i]] + hybMatList[[rna]][[type]][[j]]
            }
            c = c + 1
          }
        }
      }
      if (type2 == 'blank'){
        hybMatList2[[rna]][[typeCombinationName]][[i]] = hybMatList2[[rna]][[type1]][[i]]
      } else{
        if (type1 == 'blank') {
          hybMatList2[[rna]][[typeCombinationName]][[i]] = t(hybMatList2[[rna]][[type2]][[i]])
        } else{
          hybMatList2[[rna]][[typeCombinationName]][[i]] = hybMatList2[[rna]][[type1]][[i]] + t(hybMatList2[[rna]][[type2]][[i]])
        }
      }
    }
    sampleNames = c("s", "c")
    for (sample in  sampleNames) {
      matrixToPlot = hybMatList2[[rna]][[typeCombinationName]][[sample]][a:b, c:d]
      cols = log2(max(matrixToPlot + 1))
      
      myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols -
                                                                        1))
      
      if (cols > 14) {
        cols = log2(max(matrixToPlot[matrixToPlot < 30000] +
                          1))
        myCol = c("black",
                  colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1),
                  rep("white", (15 - cols) + 2))
      }
      
      if (directory == 0) {
        heatmap3((log2(t(
          matrixToPlot + 1
        ))),
        col = myCol,
        scale = "none" ,
        Rowv = NA,
        Colv = NA,
        useRaster = T
        )
        
      } else{
        pdf(
          paste(directory, "/", rna , "_", sample , "-", typeCombinationName , ".pdf", sep = ""),
          height = h,
          width = h
        )
        heatmap3((log2(t(
          matrixToPlot + 1
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

#  plotCombinedMatrix
#'
#'
#' Plots a contact map of two chosen samples for chosen stages in the analysis, with each chosen sample on separate halves of the contact map 
#'
#' @param cds A comradesDataSet object 
#' @param type1 The analysis stage to plot on the upper half of the heatmap
#' @param type2 The analysis stage to plot on the lower half of the heatmap
#' @param sample1 The sample number to plot on the upper half of the heatmap
#' @param sample2 The sample number to plot on the upper half of the heatmap
#' @param directory An output directory for the heatmap (use 0 for no output)
#' @param a To make a subsetted plot (left value on x)
#' @param b To make a subsetted plot (right value on x)
#' @param c To make a subsetted plot (left value on y)
#' @param d To make a subsetted plot (right value on y)
#' @param h Height of image (inches) (only useful if plotting)
#' @name plotCombinedMatrix
#' @docType methods
#' @rdname plotCombinedMatrix
#' @aliases plotCombinedMatrix,comradesDataSet-method
#' @return A heatmap of the reads of the chosen sample numbers, in the analysis stages chosen, with each chosen sample on a separate half of the heatmap 
#' @examples 
#' cds = makeExampleComradesDataSet()
#' 
#' plotCombinedMatrix(cds,
#'             type1 = "original",
#'             type2 = "noHost",
#'             b = rnaSize(cds),
#'             d = rnaSize(cds))
#' @export
setGeneric("plotCombinedMatrix", function(cds,
                                          type1 = 'original',
                                          type2 = 'original',
                                          sample1 = 1,
                                          sample2 = 1,
                                          directory = 0,
                                          a = 1,
                                          b = 50,
                                          c = 1,
                                          d = 50,
                                          h= 3)
  standardGeneric("plotCombinedMatrix"))

setMethod("plotCombinedMatrix", "comradesDataSet", function(cds, 
                                                            type1 = 'original',
                                                            type2 = 'original',
                                                            sample1 = 1,
                                                            sample2 = 1,
                                                            directory = 0, 
                                                            a = 1,
                                                            b = 50,
                                                            c = 1,
                                                            d = 50,
                                                            h= 3)  {
  hybMatList = matrixList(cds)
  rnaS = rnas(cds)
  sampleNames = names(hybMatList[[1]][[type1]])
  if (is.null(sampleNames)) {
    sampleNames = 1:length(hybMatList[[1]][[type1]])
  }
  for (rna in c(rnaS)) {
    sumOfUpper = getData(cds, data = "matrixList", type = type1)[[sample1]]
    sumOfLower = getData(cds, data = "matrixList", type = type2)[[sample2]]
    matrixToPlot = sumOfUpper + t(sumOfLower)
    matrixToPlot = matrixToPlot[a:b, c:d]
    
    # choose colour pallet
    cols = log2(max(matrixToPlot + 1))
    myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1))
    if (cols > 14) {
      cols = log2(max(matrixToPlot[matrixToPlot < 30000] + 1))
      myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1),  rep("white", (14 - cols) + 1))
    }
    
    # plot the heatmap
    if (directory == 0) {
      heatmap3((log2(t(
        matrixToPlot + 1
      ))),
      col = myCol,
      scale = "none" ,
      Rowv = NA,
      Colv = NA,
      useRaster = T
      )
    } else{
      pdf(
        paste(directory, "/", rna, "_", sampleNames[sample1], "-", type1 , "-", sampleNames[sample2], "-", type2 , ".pdf", sep = ""),
        height = h,
        width = h
      )
      heatmap3((log2(t(
        matrixToPlot + 1
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
})