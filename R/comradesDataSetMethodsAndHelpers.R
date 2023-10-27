#' @include  comradesDataSet.R
NULL



################################################################################
# Methods functions and methods relating to creation and use of
# comradesDataSets
################################################################################


#' featureInfo
#'
#' Produces a list list of 2 elemnts 'transcript' and 'family'
#' Each element contains a table with the counts for each RNA in each sample.
#'
#' @param cds a \code{comradesDataSet} object
#'
#' 
#' @name featureInfo
#' @docType methods
#' @aliases featureInfo,comradesDataSet-method
#' @rdname featureInfo
#' @return A list - Feature level and transcript level counts for each sample
#' @examples 
#' cds = makeExampleComradesDataSet()
#' featureInfo(cds)
#'
#' @export
setGeneric("featureInfo",
           function(cds)
             standardGeneric("featureInfo"))

setMethod("featureInfo",
          "comradesDataSet",
          function(cds)  {
            alteredHybList = list()
            TE = rnas(cds)
            hybList = getData(x = cds, 
                              data = "hybFiles", 
                              type = "original")
            
            for (hyb in 1:length(hybList)) {
              controlHyb = hybList[[hyb]]
              #get the right columns
              controlHyb  = as.data.frame(cbind(
                as.character(controlHyb$V4),
                as.character(controlHyb$V10)
              ))
              
              #get only those lines with the RNA
              controlHyb = controlHyb
              
              #remove the factor levels
              controlHyb$V1 = as.character(controlHyb$V1)
              controlHyb$V2 = as.character(controlHyb$V2)
              
              #Swap the columns around 
              controlHybtmp1 = controlHyb[controlHyb$V1 == TE, ]
              controlHybtmp2 = controlHyb[!(controlHyb$V1 == TE), ]
              controlHybtmp2 = rev(controlHybtmp2)
              colnames(controlHybtmp2) = c("V1", "V2")
              controlHyb = as.data.frame(rbind(controlHybtmp1, controlHybtmp2))
              
              #add to list
              alteredHybList[[TE]][[hyb]] = controlHyb
              
            }
            
            ##############################
            # Now Combine the two
            ##############################
            
            df = list()
            df[[TE]] = data.frame()
            aggList = list()
            totalNames = c()
            for (hyb in 1:length(alteredHybList[[TE]])) {
              
              
              sampleHyb = alteredHybList[[TE]][[hyb]]
              sData = sampleHyb[sampleHyb$V1 == TE, ]
              freqSample2 = aggregate(sData$V1, by = list(sData$V2), FUN = length)
              freqSample = freqSample2$x
              names(freqSample) = freqSample2$Group.1
              aggList[[TE]][[hyb]] = freqSample 
              # Get the total features that exist in the dataset
              totalNames = unique(sort(c(totalNames, names(freqSample))))
            }
            
            
            
            # geta. matrix for plotting
            tmpMat = matrix(0, nrow = length(totalNames), ncol = length(hybList))
            row.names(tmpMat) = totalNames
            colnames(tmpMat) = sampleNames(cds)
            for (i in 1:nrow(tmpMat)) {
              for (j in 1:ncol(tmpMat)) {
                tmpMat[i, j] = aggList[[TE]][[j]][row.names(tmpMat)[i]]
              }
            }
            tmpMat[is.na(tmpMat)] <- 0.00001
            df[[TE]] = as.data.frame(tmpMat)
            
            
            
            
            # now make a stats list
            statsList = list()
            statsList[[TE]] = df[[TE]]
            statsList[[TE]]$ID = sapply(row.names(statsList[[TE]]), function(x)
              tail(strsplit(x, "_")[[1]],n=1), USE.NAMES = FALSE)
            statsList[[TE]]
            
            
            # Now get the stats from aggregate
            aggList = list()
            aggList2 = list()
            
            
            aggList[[TE]] =   aggregate(. ~ ID , statsList[[TE]], sum)
            
            
            #get sample Table
            st = sampleTable(cds)
            
            #find out which samples have controls
            samples = st[which(duplicated(st$sample)), "sample"]
            
            
            
            for (i in samples) {
              # now get the control and sample index
              control =  as.numeric(which(st$sample == i &
                                            st$group == "c")) + 1
              sample = as.numeric(which(st$sample == i &
                                          st$group == "s")) + 1
              
              
              factor1 = sum(aggList[[TE]][, sample]) / sum(aggList[[TE]][, control])
              
              aggList[[TE]][, paste("norm", i, sep = "_")] = aggList[[TE]][, sample] / (aggList[[TE]][, control] * factor1)
              
              
            }
            
            data = melt(aggList[[TE]], id.vars = c("ID"))
            
            samples = paste("norm", samples, sep = "_")
            
            
            
            plot(ggplot() +
                   geom_boxplot(data = data[data$variable %in% samples, ], aes(
                     x = reorder(data[data$variable %in% samples, ]$ID, 
                                 data[data$variable %in% samples, ]$value, 
                                 FUN = median), y = log2(value)
                   )) +
                   geom_point(data = data[data$variable %in% samples, ], aes(
                     x = reorder(data[data$variable %in% samples, ]$ID, 
                                 data[data$variable %in% samples, ]$value, 
                                 FUN = median),
                     y = log2(value),
                     fill = ID
                   )) +
                   geom_hline(yintercept = 0, colour = "darkred") +
                   theme_classic() +
                   theme(axis.text.x = element_text(angle = 90, hjust = 1)) )
            
            
            featureStats = list()
            featureStats[["transcript"]] = df[[TE]]
            featureStats[["family"]] = aggList[[TE]]
            return(featureStats)
          })





#' topTranscripts
#'
#' This method prints the top transcripts that have the most duplexes
#' assigned
#'
#' @param cds a \code{comradesDataSet} object
#' @param ntop the number of entries to display
#'
#' @name topTranscripts
#' @docType methods
#' @rdname topTranscripts
#' @aliases topTranscripts,comradesDataSet-method
#' @return A table, the number of counts per sample per transcript
#' @examples 
#' cds = makeExampleComradesDataSet()
#' topTranscripts(cds)
#' @export
#'
setGeneric("topTranscripts",
           function(cds,
                    ntop = 10)
             standardGeneric("topTranscripts"))

setMethod("topTranscripts",
          "comradesDataSet",
          function(cds,
                   ntop = 10)  {
            c = group(cds)[["s"]]
            vect = c()
            for (i in c) {
              vect = c(vect,
                       hybFiles(cds)[["all"]][["all"]][[i]]$V4,
                       hybFiles(cds)[["all"]][["all"]][[i]]$V10)
            }
            x = table(vect)[order(table(vect), decreasing = T)]
            
            
            c = group(cds)[["c"]]
            vect = c()
            for (i in c) {
              vect = c(vect,
                       hybFiles(cds)[["all"]][["all"]][[i]]$V4,
                       
                       hybFiles(cds)[["all"]][["all"]][[i]]$V10)
            }
            y = table(vect)[order(table(vect), decreasing = T)]
            
            
            y = y[names(x[1:ntop])]
            x = x[1:ntop]
            
            
            if(length(which(complete.cases(c(x)))) == 1){
              x2 = data.frame(
                t = names(x), 
                s = x,
                c = y)
              colnames(x2) = c("RNA", "Samples", "Control")
              x2$enrichment = x2$Samples /  x2$Control
              
            }else{
            
            x2 = as.data.frame(x)
            x2$control = y
            colnames(x2) = c("RNA", "Samples", "Control")
            x2$enrichment = x2$Samples /  x2$Control
            
            }
            
            for (i in names(hybFiles(cds)[["all"]][["all"]])) {
              t = c(hybFiles(cds)[["all"]][["all"]][[i]]$V4,
                    hybFiles(cds)[["all"]][["all"]][[i]]$V10)
              t = table(t)
              
              t = t[names(x[1:ntop])]
              x2[, i] = t
            }
            x2 = x2[, c(1, 5:ncol(x2), 2, 3, 4)]
            x2
          })





#' topInteractions
#'
#' This method prints the top transcript interactions that have the most duplexes
#' assigned
#'
#' @param cds a \code{comradesDataSet} object
#' @param ntop the number of entries to display
#'
#' @name topInteractions
#' @docType methods
#' @rdname topInteractions
#' @aliases topInteractions,comradesDataSet-method
#' @return A table, the number of counts per sample per interaction
#' @examples 
#' cds = makeExampleComradesDataSet()
#' topInteractions(cds)
#'
#' @export
setGeneric("topInteractions",
           function(cds,
                    ntop = 10)
             standardGeneric("topInteractions"))

setMethod("topInteractions",
          "comradesDataSet",
          function(cds,
                   ntop = 10)  {
            c = group(cds)[["s"]]
            vect = c()
            for (i in c) {
              vect = c(vect,
                       paste(hybFiles(cds)[["all"]][["all"]][[i]]$V4,
                             hybFiles(cds)[["all"]][["all"]][[i]]$V10, sep = "::"))
            }
            x = table(vect)[order(table(vect), decreasing = T)]
            
            
            c = group(cds)[["c"]]
            vect = c()
            for (i in c) {
              vect = c(vect,
                       paste(hybFiles(cds)[["all"]][["all"]][[i]]$V4,
                             hybFiles(cds)[["all"]][["all"]][[i]]$V10, sep = "::"))
            }
            y = table(vect)[order(table(vect), decreasing = T)]
            
            
            y = y[names(x[1:ntop])]
            x = x[1:ntop]
            
            if(length(which(complete.cases(c(x)))) == 1){
              x2 = data.frame(
                t = names(x), 
                s = x,
                c = y)
              colnames(x2) = c("RNA", "Samples", "Control")
              x2$enrichment = x2$Samples /  x2$Control
              
            }else{
            
            x2 = as.data.frame(x)
            x2$control = y
            colnames(x2) = c("RNA", "Samples", "Control")
            x2$enrichment = x2$Samples /  x2$Control
            
            }
            
            for (i in names(hybFiles(cds)[["all"]][["all"]])) {
              t =   paste(hybFiles(cds)[["all"]][["all"]][[i]]$V4,
                          hybFiles(cds)[["all"]][["all"]][[i]]$V10, sep = "::")
              t = table(t)
              
              t = t[names(x[1:ntop])]
              x2[, i] = t
            }
            x2 = x2[, c(1, 5:ncol(x2), 2, 3, 4)]
            x2
            
          })





#' topInteracters
#'
#' This method prints the top transcripts that have the most duplexes
#' assigned that interact with the transcript of interest 
#'
#' @param cds a \code{comradesDataSet} object
#' @param ntop the number of entries to display
#'
#' @name topInteracters
#' @docType methods
#' @rdname topInteracters
#' @aliases topInteracters,comradesDataSet-method
#' @return A table, the number of counts per sample per interacting transcript
#' @examples 
#' cds = makeExampleComradesDataSet()
#' topInteracters(cds)
#'
#' @export
setGeneric("topInteracters",
           function(cds,
                    ntop = 10)
             standardGeneric("topInteracters"))

setMethod("topInteracters",
          "comradesDataSet",
          function(cds,
                   ntop = 10)  {
            c = group(cds)[["s"]]
            vect = c()
            for (i in c) {
              vecto = hybFiles(cds)[[2]][["host"]][[i]]
              vecto = vecto[vecto$V4 == rnas(cds), ]
              vect = c(vect, vecto$V10)
            }
            x = table(vect)[order(table(vect), decreasing = T)]
            
            
            c = group(cds)[["c"]]
            vect = c()
            for (i in c) {
              vecto = hybFiles(cds)[[2]][["host"]][[i]]
              vecto = vecto[vecto$V4 == rnas(cds), ]
              vect = c(vect, vecto$V10)
            }
            y = table(vect)[order(table(vect), decreasing = T)]
            
            
            y = y[names(x[1:ntop])]
            x = x[1:ntop]
            
            if(length(which(complete.cases(c(x)))) == 1){
              x2 = data.frame(
                t = names(x), 
                s = x,
                c = y)
              colnames(x2) = c("RNA", "Samples", "Control")
              x2$enrichment = x2$Samples /  x2$Control
              
            }else{
              
              x2 = as.data.frame(x)
              x2$control = y
              colnames(x2) = c("RNA", "Samples", "Control")
              x2$enrichment = x2$Samples /  x2$Control
            }
            
            for (i in names(hybFiles(cds)[[rnas(cds)]][["host"]])) {
              t =   hybFiles(cds)[[rnas(cds)]][["host"]][[i]]$V10
              t = table(t)
              
              t = t[names(x[1:ntop])]
              x2[, i] = t
            }
            x2 = x2[, c(1, 5:ncol(x2), 2, 3, 4)]
            x2
            
          })



#' getInteractions
#'
#' This method returns a table interactions of an RNA (interactor) on the RNA of
#' interest use topInteracters. 
#'
#' @param cds a \code{comradesDataSet} object
#' @param interactor The rna to show interactions with
#'
#'
#' @name getInteractions
#' @docType methods
#' @rdname getInteractions
#' @return A table showign the read coverage of the interacting RNA 
#' @aliases getInteractions,comradesDataSet-method
#' @examples
#' cds = makeExampleComradesDataSet()
#' getInteractions(cds, 'transcript2')
#'
#' @export
setGeneric("getInteractions",
           function(cds,
                    interactor)
             standardGeneric("getInteractions"))

setMethod("getInteractions",
          "comradesDataSet",
          function(cds,
                   interactor)  {
            seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
            
            # subset the hyb files based on the interacter of choice
            table = data.frame()
            for (i in names(hybFiles(cds)[[rnas(cds)]][["host"]])) {
              x = hybFiles(cds)[[rnas(cds)]][["host"]][[i]]
              x = x[x$V10 == interactor, ]
              if (nrow(x) == 0) {
                v = as.data.frame(0)
                v$rna = interactor
                v$sample = i
                colnames(v) = c("Position", "rna", "sample")
                table = rbind.data.frame(table, v)
              } else {
                starts = x$V7
                ends = x$V8
                v = as.data.frame(unlist(seq2(from = starts, to = ends)))
                v$rna = interactor
                v$sample = i
                colnames(v) = c("Position", "rna", "sample")
                table = rbind.data.frame(table, v)
              }
              
            }
            
            colnames(table) = c("Position", "rna", "sample")
            table = aggregate(
              table$Position,
              by = list(table$rna,
                        table$sample,
                        table$Position),
              FUN = length
            )
            colnames(table) = c("rna", "sample", "Position", "depth")
            table
          })




#' getReverseInteractions
#'
#' This method prints interactions of the RNA of interest on another RNA
#' transcript.
#'
#' @param cds a \code{comradesDataSet} object
#' @param interactor The rna to show interactions with
#'
#' @name getReverseInteractions
#' @docType methods
#' @rdname getReverseInteractions
#' @aliases getReverseInteractions,comradesDataSet-method
#' @return A long format table shoing the read coverage of chosen RNA
#' @examples
#' cds = makeExampleComradesDataSet()
#' getReverseInteractions(cds, 'transcript2')
#'
#' @export
setGeneric("getReverseInteractions",
           function(cds,
                    interactor)
             standardGeneric("getReverseInteractions"))

setMethod("getReverseInteractions",
          "comradesDataSet",
          function(cds,
                   interactor)  {
            seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
            
            # subset the hyb files based on the interacter of choice
            table = data.frame()
            for (i in names(hybFiles(cds)[[rnas(cds)]][["host"]])) {
              x = hybFiles(cds)[[rnas(cds)]][["host"]][[i]]
              x = x[x$V10 == interactor, ]
              if (nrow(x) == 0) {
                v = as.data.frame(0)
                v$rna = interactor
                v$sample = i
                colnames(v) = c("Position", "rna", "sample")
                table = rbind.data.frame(table, v)
              } else {
                starts = x$V13
                ends = x$V14
                v = as.data.frame(unlist(seq2(from = starts, to = ends)))
                v$rna = interactor
                v$sample = i
                colnames(v) = c("Position", "rna", "sample")
                table = rbind.data.frame(table, v)
              }
              
            }
            
            colnames(table) = c("Position", "rna", "sample")
            table = aggregate(
              table$Position,
              by = list(table$rna,
                        table$sample,
                        table$Position),
              FUN = length
            )
            colnames(table) = c("rna", "sample", "Position", "depth")
            table
          })





#' swapHybs
#'
#' Swap the table to ensure that 3 prime most duplex side is on he left of the
#' table used to make one sides heatmaps and other reasons where having the
#' left of the table coming after the right side is a problem. Different from
#' swapHybs as it ensure that BOTH duplex sides originate from the RNA of
#' interest.
#'
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#'
#' @name swapHybs
#' @docType methods
#' @rdname swapHybs
#' @return A list of "swapped" hyb datas
#'
swapHybs = function(hybList, 
                    rna) {
  hybListSwap = hybList
  for (hyb in 1:length(hybList)) {
    hybListS  =  hybList[[hyb]][hybList[[hyb]]$V4 == rna |
                                  hybList[[hyb]]$V10 == rna, ]
    
    
    
    comb = hybListS
    
    
    
    
    hybListSwap[[hyb]] = comb
  }
  
  return(hybListSwap)
}



#' swapHybs2
#'
#' Swap the table to ensure that 3 prime most duplex side is on the left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Different from swapHybs as it
#' ensure that BOTH duplex sides originate from the RNA of interest.
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#' @name swapHybs2
#' @docType methods
#' @rdname swapHybs2
#' @return A list of "swapped" hyb data
#'
swapHybs2 = function(hybList, 
                     rna) {
  for (hyb in 1:length(hybList)) {
    hybList18S = hybList[[hyb]][as.character(hybList[[hyb]]$V4) == rna &
                                  as.character(hybList[[hyb]]$V10) == rna, ]
    tmp1 = hybList18S[hybList18S$V7 < hybList18S$V13, ]
    tmp2 = hybList18S[hybList18S$V7 > hybList18S$V13, ]
    tmp2 = tmp2[, c(
      "V1",
      "V2",
      "V3",
      "V4",
      "V11",
      "V12",
      "V13",
      "V14",
      "V15",
      "V10",
      "V5",
      "V6",
      "V7",
      "V8",
      "V9"
    )]
    colnames(tmp2) = colnames(tmp1)
    comb = rbind.data.frame(tmp1, tmp2)
    hybList[[hyb]] = comb
  }
  
  return(hybList)
}




#' swapHybs3
#'
#' Swap the table to ensure that 3 prime most duplex side is ont he left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Different from swapHybs as it
#' ensure that BOTH duplex sides originate from the RNA of interest.
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#'
#' @return A list of "swapped" hyb datas
#' @name swapHybs3
#' @docType methods
#' @rdname swapHybs3
swapHybs3 = function(hybList, rna) {
  hybListSwap = hybList
  for (hyb in 1:length(hybList)) {
    hybListS  =  hybList[[hyb]][hybList[[hyb]]$V4 == rna |
                                  hybList[[hyb]]$V10 == rna, ]
    tmp1 = hybListS[hybListS$V4 == rna, ]
    tmp2 = hybListS[hybListS$V4 != rna, ]
    tmp2 = tmp2[, c(
      "V1",
      "V2",
      "V3",
      "V10",
      "V11",
      "V12",
      "V13",
      "V14",
      "V15",
      "V4",
      "V5",
      "V6",
      "V7",
      "V8",
      "V9"
    )]
    colnames(tmp2) = colnames(tmp1)
    comb = rbind.data.frame(tmp1, tmp2, stringsAsFactors = F)
    hybListSwap[[hyb]] = comb
  }
  return(hybListSwap)
}









#' makeExampleComradesDataSet
#'
#' Creat a minimal example comradesdataSetObject
#'
#' @return An example comradesDataSet objext
#' @name makeExampleComradesDataSet
#' @rdname makeExampleComradesDataSet
#' @examples
#' cds = makeExampleComradesDataSet()
#' @export
makeExampleComradesDataSet = function() {
  
  c4 = c(rep("transcript1",100),rep("transcript2",100) )
  c10 = c(rep("transcript1",200) )
  c1 = 1:200
  c2 = rep(paste(rep("A", 40), collapse = ""),200)
  c3 = rep(".",200)
  c9 = rep(".",200)
  c15 = rep(".",200)
  c5 = rep(1,200)
  c11 = rep(21,200)
  c6 = rep(20,200)
  c12= rep(40,200)
  # short distance 50
  c7 = sample(1:5, 50, replace = T)
  c8 = sample(20:25, 50, replace = T)
  c13 = sample(20:25, 50, replace = T)
  c14 = sample(40:45, 50, replace = T)
  # long distance 50
  c7 = c(c7,sample(1:5, 50, replace = T))
  c8 = c(c8,sample(20:25, 50, replace = T))
  c13 = c(c13,sample(60:70, 50, replace = T))
  c14 = c(c14,sample(80:83, 50, replace = T))
  # inter RNA 100
  c7 = c(c7,sample(1:5, 100, replace = T))
  c8 = c(c8,sample(20:25, 100, replace = T))
  c13 = c(c13,sample(1:5, 100, replace = T))
  c14 = c(c14,sample(20:25, 100, replace = T))
  
  exampleInput = data.frame(V1 = c1,
                            V2 = c2,
                            V3 = c3,
                            V4 = c4,
                            V5 = as.numeric(c5),
                            V6 = as.numeric(c6),
                            V7 = as.numeric(c7),
                            V8 = as.numeric(c8),
                            V9 = c9,
                            V10 = c10,
                            V11 = as.numeric(c11),
                            V12 = as.numeric(c12),
                            V13 = as.numeric(c13),
                            V14 = as.numeric(c14),
                            V15 = c15)
  
  
  file = tempfile()
  write.table(exampleInput,file = file, quote = F, row.names = F, sep = "\t", col.names = F)
  
  
  
  
  c4 = c(rep("transcript1",55),rep("transcript2",90) )
  c10 = c(rep("transcript1",145) )
  c1 = 1:145
  c2 = rep(paste(rep("A", 40), collapse = ""),145)
  c3 = rep(".",145)
  c9 = rep(".",145)
  c15 = rep(".",145)
  c5 = rep(1,145)
  c11 = rep(21,145)
  c6 = rep(20,145)
  c12= rep(40,145)
  # short distance 55
  c7 = sample(1:5, 55, replace = T)
  c8 = sample(20:25, 55, replace = T)
  c13 = sample(20:25, 55, replace = T)
  c14 = sample(40:45, 55, replace = T)

  # inter RNA 100
  c7 = c(c7,sample(1:40, 90, replace = T))
  c8 = c(c8,sample(20:75, 90, replace = T))
  c13 = c(c13,sample(1:40, 90, replace = T))
  c14 = c(c14,sample(20:75, 90, replace = T))
  
  exampleInput = data.frame(V1 = c1,
                            V2 = c2,
                            V3 = c3,
                            V4 = c4,
                            V5 = as.numeric(c5),
                            V6 = as.numeric(c6),
                            V7 = as.numeric(c7),
                            V8 = as.numeric(c8),
                            V9 = c9,
                            V10 = c10,
                            V11 = as.numeric(c11),
                            V12 = as.numeric(c12),
                            V13 = as.numeric(c13),
                            V14 = as.numeric(c14),
                            V15 = c15)
  
  
  file2 = tempfile()
  write.table(exampleInput,file = file2, quote = F, row.names = F, sep = "\t", col.names = F)
  
  
  
  
  
  
  
  # Set up the sample table. ----
  sampleTabler1 = c(file, "s", "1", "s1")
  sampleTabler2 = c(file2, "c", "1", "c1")
  
  # make the sample table 
  sampleTable2 = rbind.data.frame(sampleTabler1, sampleTabler2)
  
  # add the column names 
  colnames(sampleTable2) = c("file", "group", "sample", "sampleName")
  
  # Choose RNA and set up the object ----
  rna = c("transcript1")
  # load the object
  cds = comradesDataSet(rnas = rna,
                        rnaSize = 0,
                        sampleTable = sampleTable2)
  
  return(cds)
}








