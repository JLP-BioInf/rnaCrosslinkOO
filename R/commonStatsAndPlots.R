#' @include  rnaCrosslinkDataSet.R
NULL


#' clusterNumbers
#'
#' This method prints a table showing the number of clusters in each step 
#' of the analysis
#'
#' @param knowClusteredCds A rnaCrosslinkDataSet object after clustering has been performed
#' @param rna The RNA ID of interest - use rna(cdsObject).
#' @name clusterNumbers
#' @docType methods
#' @rdname clusterNumbers
#' @aliases clusterNumbers,rnaCrosslinkDataSet-method
#' @return A data.frame shoing the number of clusters for each sample
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusteredCds = clusterrnaCrosslink(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 1)
#' clusterNumbers(clusteredCds)
#' @export
setGeneric("clusterNumbers",
           function(knowClusteredCds,
                    rna)
             standardGeneric("clusterNumbers"))

setMethod("clusterNumbers", "rnaCrosslinkDataSet", function(knowClusteredCds,
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
#' @param knowClusteredCds A rnaCrosslinkDataSet object after clustering has been performed
#' @param rna The RNA ID of interest - use rna(cdsObject).
#'
#' @name readNumbers
#' @docType methods
#' @rdname readNumbers
#' @aliases readNumbers,rnaCrosslinkDataSet-method
#' @return A data.frame shoing the number of reads in clusters for each sample
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusteredCds = clusterrnaCrosslink(cds,
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

setMethod("readNumbers", "rnaCrosslinkDataSet", function(knowClusteredCds,
                                                     rna)  {
  sampleTable = sampleTable(knowClusteredCds)
  clusters = InputFiles(knowClusteredCds)
  #get numbers of original Input files
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
  
  if (class(knowClusteredCds)[1] == "rnaCrosslinkDataSet") {
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
#' @param cds A rnaCrosslinkDataSet object 
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
#' @aliases plotMatrices,rnaCrosslinkDataSet-method
#' @return A heatmap of the reads in the analysis stage chosen
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
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

setMethod("plotMatrices", "rnaCrosslinkDataSet", function(cds, 
                                                      type = 'original', 
                                                      directory = 0, 
                                                      a = 1,
                                                      b = 50,
                                                      c = 1,
                                                      d = 50,
                                                      h= 3)  {
  InputMatList = matrixList(cds)
  rnaS = rnas(cds)
  sampleNames = names(InputMatList[[1]][[type]])
  if (is.null(sampleNames)) {
    sampleNames = 1:length(InputMatList[[1]][[type]])
  }
  
  
  
  for (sample in  1:length(sampleNames)) {
    for (rna in c(rnaS)) {
      cols = log2(max(InputMatList[[rna]][[type]][[sample]][a:b, c:d] + 1))
      
      myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols -
                                                                        1))
      
      if (cols > 14) {
        cols = log2(max(InputMatList[[rna]][[type]][[sample]][a:b, c:d][InputMatList[[rna]][[type]][[sample]][a:b, c:d] < 30000] +
                          1))
        myCol = c("black",
                  colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1),
                  rep("white", (14 - cols) + 1))
      }
      
      if (directory == 0) {
        heatmap3((log2(t(
          InputMatList[[rna]][[type]][[sample]][a:b, c:d] + 1
        ))),
        col = myCol,
        scale = "none" ,
        Rowv = NA,
        Colv = NA,
        useRaster = TRUE
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
          InputMatList[[rna]][[type]][[sample]][a:b, c:d] + 1
        ))),
        col = myCol,
        scale = "none" ,
        Rowv = NA,
        Colv = NA,
        useRaster = TRUE
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
#' @param cds A rnaCrosslinkDataSet object 
#' @param type1 The analysis stage to plot on the upper half of the heatmap (use 'blank' to leave this half blank)
#' @param type2 The analysis stage to plot on the lower half of the heatmap (use 'blank' to leave this half blank)
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
#' @aliases plotMatricesAverage,rnaCrosslinkDataSet-method
#' @return A heatmap of the reads in the two analysis stages chosen, with each chosen stage on a separate half of the heatmap
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
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

setMethod("plotMatricesAverage", "rnaCrosslinkDataSet", function(cds,
                                                             type1 = 'noHost',
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
    InputMatList = matrixList(cds)
    InputMatList2 = InputMatList
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
              sum(InputMatList[[rna]][[type]][[j]])
              InputMatList2[[rna]][[type]][[i]] = InputMatList[[rna]][[type]][[j]]
            } else{
              InputMatList2[[rna]][[type]][[i]] = InputMatList2[[rna]][[type]][[i]] + InputMatList[[rna]][[type]][[j]]
            }
            c = c + 1
          }
        }
      }
      if (type2 == 'blank'){
        InputMatList2[[rna]][[typeCombinationName]][[i]] = InputMatList2[[rna]][[type1]][[i]]
      } else{
        if (type1 == 'blank') {
          InputMatList2[[rna]][[typeCombinationName]][[i]] = t(InputMatList2[[rna]][[type2]][[i]])
        } else{
          InputMatList2[[rna]][[typeCombinationName]][[i]] = InputMatList2[[rna]][[type1]][[i]] + t(InputMatList2[[rna]][[type2]][[i]])
        }
      }
    }
    sampleNames = c("s", "c")
    for (sample in  sampleNames) {
      matrixToPlot = InputMatList2[[rna]][[typeCombinationName]][[sample]][a:b, c:d]
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
        useRaster = TRUE
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
        useRaster = TRUE
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
#' @param cds A rnaCrosslinkDataSet object 
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
#' @aliases plotCombinedMatrix,rnaCrosslinkDataSet-method
#' @return A heatmap of the reads of the chosen sample numbers, in the analysis stages chosen, with each chosen sample on a separate half of the heatmap 
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
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

setMethod("plotCombinedMatrix", "rnaCrosslinkDataSet", function(cds, 
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
  InputMatList = matrixList(cds)
  rnaS = rnas(cds)
  sampleNames = names(InputMatList[[1]][[type1]])
  if (is.null(sampleNames)) {
    sampleNames = 1:length(InputMatList[[1]][[type1]])
  }
  #sample1 = which(sampleTable(cds)[sampleTable(cds)$sampleName == sample1])
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
      useRaster = TRUE
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
      useRaster = TRUE
      )
      dev.off()
    }
  }
})

#  plotInteractions
#'
#'
#' Plots a contact map of interactions of each sample of an RNA (interactor) on the RNA of interest
#'
#' @param cds A rnaCrosslinkDataSet object 
#' @param rna The RNA of interest
#' @param interactor The RNA to show interactions with
#' @param directory An output directory for the heatmap (use 0 for no output)
#' @param a To make a subsetted plot (left value on x)
#' @param b To make a subsetted plot (right value on x) (use 'max' to plot the whole RNA strand length)
#' @param c To make a subsetted plot (left value on y)
#' @param d To make a subsetted plot (right value on y) (use 'max' to plot the whole RNA strand length)
#' @param h Height of image (inches) (only useful if plotting)
#' @name plotInteractions
#' @docType methods
#' @rdname plotInteractions
#' @aliases plotInteractions,rnaCrosslinkDataSet-method
#' @return A heatmap of interactions of the RNA (interactor) on the RNA of interest
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' plotInteractions(cds,
#'             rna = "transcript1",
#'             interactor = "transcript2",
#'             b = "max",
#'             d = "max")
#' @export
setGeneric("plotInteractions", function(cds,
                                        rna,
                                        interactor,
                                        directory = 0,
                                        a = 1,
                                        b = 50,
                                        c = 1,
                                        d = 50,
                                        h= 3)
  standardGeneric("plotInteractions"))

setMethod("plotInteractions", "rnaCrosslinkDataSet", function(cds, 
                                                          rna,
                                                          interactor,
                                                          directory = 0, 
                                                          a = 1,
                                                          b = 50,
                                                          c = 1,
                                                          d = 50,
                                                          h= 3)  {
  for (i in names(InputFiles(cds)[[rnas(cds)]][["host"]])) {
    table = InputFiles(cds)[[rnas(cds)]][["host"]][[i]]
    InputOutput =  table[as.character(table$V4) == rna & as.character(table$V10) == interactor,] 
    InputOutput = unique(InputOutput)
    startsends = InputOutput[,c(7,8,13,14)]
    maxX = max(InputOutput[,8])
    maxY = max(InputOutput[,14])
    matrixToPlot = matrix(0, nrow = maxX, ncol = maxY)
    for(rowNumber in 1:nrow(startsends)){
      data = startsends[rowNumber,]
      xData = seq(data$V7, data$V8)
      yData = seq(data$V14, data$V13)
      matrixToPlot[xData,yData] = matrixToPlot[xData,yData] + 1
    }

    bActual = b
    if (b == "max"){
      bActual = maxX
    }
    dActual = d
    if (d == "max"){
      dActual = maxY
    }
    matrixToPlot = matrixToPlot[a:bActual, c:dActual]
    
    # choose colour pallet
    cols = log2(max(matrixToPlot + 1))
    myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1))
    if (cols > 14) {
      cols = log2(max(matrixToPlot[matrixToPlot < 30000] + 1))
      myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1),  rep("white", (14 - cols) + 1))
    }
    
    # plot the heatmap
    if (directory == 0) {
      aspectHeatmap((log2(t(
        matrixToPlot + 1
      ))),
      col = myCol,
      scale = "none" ,
      Rowv = NA,
      Colv = NA,
      hExp = dActual/bActual
      )
    } else{
      pdf(
        paste(directory, "/", rna, "_", i, "-", interactor, ".pdf", sep = ""),
        height = h,
        width = h
      )
      aspectHeatmap((log2(t(
        matrixToPlot + 1
      ))),
      col = myCol,
      scale = "none" ,
      Rowv = NA,
      Colv = NA,
      hExp = dActual/bActual
      )
      dev.off()
    }
  }
})

#  plotInteractionsAverage
#'
#'
#' Plots a contact map of interactions of all samples of an RNA (interactor) on the RNA of interest
#'
#' @param cds A rnaCrosslinkDataSet object 
#' @param rna The RNA of interest
#' @param interactor The RNA to show interactions with
#' @param directory An output directory for the heatmap (use 0 for no output)
#' @param a To make a subsetted plot (left value on x)
#' @param b To make a subsetted plot (right value on x) (use 'max' to plot the whole RNA strand length)
#' @param c To make a subsetted plot (left value on y)
#' @param d To make a subsetted plot (right value on y) (use 'max' to plot the whole RNA strand length)
#' @param h Height of image (inches) (only useful if plotting)
#' @name plotInteractionsAverage
#' @docType methods
#' @rdname plotInteractionsAverage
#' @aliases plotInteractionsAverage,rnaCrosslinkDataSet-method
#' @return A heatmap of interactions of all samples of the RNA (interactor) on the RNA of interest
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' plotInteractionsAverage(cds,
#'             rna = "transcript1",
#'             interactor = "transcript2",
#'             b = "max",
#'             d = "max")
#' @export
setGeneric("plotInteractionsAverage", function(cds,
                                               rna,
                                               interactor,
                                               directory = 0,
                                               a = 1,
                                               b = 50,
                                               c = 1,
                                               d = 50,
                                               h= 3)
  standardGeneric("plotInteractionsAverage"))

setMethod("plotInteractionsAverage", "rnaCrosslinkDataSet", function(cds, 
                                                                     rna,
                                                                     interactor,
                                                                     directory = 0, 
                                                                     a = 1,
                                                                     b = 50,
                                                                     c = 1,
                                                                     d = 50,
                                                                     h= 3)  {
  
  for (cors in c("s", "c")) {
    counter=1
    maxXAll=0
    maxYAll=0
    matrices=list()
    for (sampleName in group(cds)[[cors]]) {
      table = InputFiles(cds)[[rnas(cds)]][["host"]][[sampleName]]
      InputOutput =  table[as.character(table$V4) == rna & as.character(table$V10) == interactor,] 
      InputOutput = unique(InputOutput)
      startsends = InputOutput[,c(7,8,13,14)]
      maxX = max(InputOutput[,8])
      if (maxXAll<maxX){
        maxXAll=maxX
      }
      maxY = max(InputOutput[,14])
      if (maxYAll<maxY){
        maxYAll=maxY
      }
      matrixToAdd = matrix(0, nrow = maxX, ncol = maxY)
      for(rowNumber in 1:nrow(startsends)){
        data = startsends[rowNumber,]
        xData = seq(data$V7, data$V8)
        yData = seq(data$V14, data$V13)
        matrixToAdd[xData,yData] = matrixToAdd[xData,yData] + 1
      }
      matrices[[counter]]=matrixToAdd
      counter=counter+1
    }
    matrixToPlot = matrix(0, nrow = maxXAll, ncol = maxYAll)
    for (mat in matrices){
      matrixToPlot[1:dim(mat)[[1]],1:dim(mat)[[2]]] = matrixToPlot[1:dim(mat)[[1]],1:dim(mat)[[2]]] + mat
    }
    bActual = b
    if (b == "max"){
      bActual = maxXAll
    }
    dActual = d
    if (d == "max"){
      dActual = maxYAll
    }
    matrixToPlot = matrixToPlot[a:bActual, c:dActual]
    
    # choose colour pallet
    cols = log2(max(matrixToPlot + 1))
    myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1))
    if (cols > 14) {
      cols = log2(max(matrixToPlot[matrixToPlot < 30000] + 1))
      myCol = c("black", colorRampPalette(c(brewer.pal(9, "YlOrRd")))(cols - 1),  rep("white", (14 - cols) + 1))
    }
    
    # plot the heatmap
    if (directory == 0) {
      aspectHeatmap((log2(t(
        matrixToPlot + 1
      ))),
      col = myCol,
      scale = "none" ,
      Rowv = NA,
      Colv = NA,
      hExp = dActual/bActual
      )
    } else{
      pdf(
        paste(directory, "/", rna, "_", cors, "-", interactor, ".pdf", sep = ""),
        height = h,
        width = h
      )
      aspectHeatmap((log2(t(
        matrixToPlot + 1
      ))),
      col = myCol,
      scale = "none" ,
      Rowv = NA,
      Colv = NA,
      hExp = dActual/bActual
      )
      dev.off()
    }
  }
})
