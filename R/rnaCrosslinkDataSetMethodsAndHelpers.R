#' @include  rnaCrosslinkDataSet.R
NULL



################################################################################
# Methods functions and methods relating to creation and use of
# rnaCrosslinkDataSets
################################################################################


#' featureInfo
#'
#' Produces a list list of 2 elemnts 'transcript' and 'family'
#' Each element contains a table with the counts for each RNA in each sample that 
#' interact with the target RNA
#'
#' @param cds a \code{rnaCrosslinkDataSet} object
#'
#' 
#' @name featureInfo
#' @docType methods
#' @aliases featureInfo,rnaCrosslinkDataSet-method
#' @rdname featureInfo
#' @return A list - Feature level and transcript level counts for each sample
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' featureInfo(cds)
#'
#' @export
setGeneric("featureInfo",
           function(cds)
             standardGeneric("featureInfo"))

setMethod("featureInfo",
          "rnaCrosslinkDataSet",
          function(cds)  {
            alteredInputList = list()
            TE = rnas(cds)
            InputList = getData(x = cds, 
                                data = "InputFiles", 
                                type = "original")
            
            for (Input in 1:length(InputList)) {
              controlInput = InputList[[Input]]
              #get the right columns
              controlInput  = as.data.frame(cbind(
                as.character(controlInput$V4),
                as.character(controlInput$V10)
              ))
              
              #get only those lines with the RNA
              controlInput = controlInput
              
              #remove the factor levels
              controlInput$V1 = as.character(controlInput$V1)
              controlInput$V2 = as.character(controlInput$V2)
              
              #Swap the columns around 
              controlInputtmp1 = controlInput[controlInput$V1 == TE, ]
              controlInputtmp2 = controlInput[!(controlInput$V1 == TE), ]
              controlInputtmp2 = rev(controlInputtmp2)
              colnames(controlInputtmp2) = c("V1", "V2")
              controlInput = as.data.frame(rbind(controlInputtmp1, controlInputtmp2))
              
              #add to list
              alteredInputList[[TE]][[Input]] = controlInput
              
            }
            
            ##############################
            # Now Combine the two
            ##############################
            
            df = list()
            df[[TE]] = data.frame()
            aggList = list()
            totalNames = c()
            for (Input in 1:length(alteredInputList[[TE]])) {
              
              
              sampleInput = alteredInputList[[TE]][[Input]]
              sData = sampleInput[sampleInput$V1 == TE, ]
              freqSample2 = aggregate(sData$V1, by = list(sData$V2), FUN = length)
              freqSample = freqSample2$x
              names(freqSample) = freqSample2$Group.1
              aggList[[TE]][[Input]] = freqSample 
              # Get the total features that exist in the dataset
              totalNames = unique(sort(c(totalNames, names(freqSample))))
            }
            
            
            
            # geta. matrix for plotting
            tmpMat = matrix(0, nrow = length(totalNames), ncol = length(InputList))
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




#' topInteractions
#'
#' This method prints the top transcript interactions that have the most duplexes
#' assigned
#'
#' @param cds a \code{rnaCrosslinkDataSet} object
#' @param ntop the number of entries to display
#'
#' @name topInteractions
#' @docType methods
#' @rdname topInteractions
#' @aliases topInteractions,rnaCrosslinkDataSet-method
#' @return A table, the number of counts per sample per interaction
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' topInteractions(cds)
#'
#' @export
setGeneric("topInteractions",
           function(cds,
                    ntop = 10)
             standardGeneric("topInteractions"))

setMethod("topInteractions",
          "rnaCrosslinkDataSet",
          function(cds,
                   ntop = 10)  {
            ci = group(cds)[["c"]]
            si = group(cds)[["s"]]
            
            
            
            
            
            df= data.frame()
            
            for (i in names(InputFiles(cds)[["all"]][["all"]])) {
              
              
              same = which(InputFiles(cds)[["all"]][["all"]][[i]]$V4 ==
                             InputFiles(cds)[["all"]][["all"]][[i]]$V10 )
              notSame = which(InputFiles(cds)[["all"]][["all"]][[i]]$V4 !=
                                InputFiles(cds)[["all"]][["all"]][[i]]$V10)
              
              
              s =   paste(InputFiles(cds)[["all"]][["all"]][[i]]$V4[same],
                          InputFiles(cds)[["all"]][["all"]][[i]]$V10[same], sep = "::")
              #s = unlist(lapply(s, function(x)  
              #  unlist(strsplit(x,split = "::"))[1] ))
              s = as.data.frame(table(s))
              s$type="intra"
              s$sample = i
              colnames(s) = c("RNA","reads","type","sample")
              
              
              
              t =   paste(InputFiles(cds)[["all"]][["all"]][[i]]$V4[notSame],
                          InputFiles(cds)[["all"]][["all"]][[i]]$V10[notSame], sep = "::")
              
              #t1 = unlist(lapply(t, function(x)  
              #  unlist(strsplit(x,split = "::"))[1] ))
              #t2 = unlist(lapply(t, function(x)  
              #  unlist(strsplit(x,split = "::"))[2] ))
              #t = c(t1,t2)
              t = as.data.frame(table(t))
              t$type="inter"
              t$sample = i
              colnames(t) = c("RNA","reads","type","sample")
              
              df =rbind.data.frame(df,t,s)
            }
            
            df = dcast(df,formula = RNA + type ~ sample,   value.var = "reads")
            
            
            
            ci = sampleTable(cds)[ci,"sampleName"]
            si = sampleTable(cds)[si,"sampleName"]
            
            if(length(si) == 1){
              df$samples = df[,si]
              df$control = df[,ci]
              df$enrichment = df$samples / df$control
              
            }else{
              
              df$samples = rowSums(df[,si],na.rm = TRUE)
              df$control = rowSums(df[,ci],na.rm = TRUE)
              df$enrichment = df$samples / df$control
              
            }
            
            df = df[order(df$samples,decreasing = TRUE),]
            df = df[1:ntop,]
            
            
            df
          })





















#' topTranscripts
#'
#' This method prints the top transcripts that have the most duplexes
#' assigned
#'
#' @param cds a \code{rnaCrosslinkDataSet} object
#' @param ntop the number of entries to display
#'
#' @name topTranscripts
#' @docType methods
#' @rdname topTranscripts
#' @aliases topTranscripts,rnaCrosslinkDataSet-method
#' @return A table, the number of counts per sample per transcript
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' topTranscripts(cds)
#' @export
#'
setGeneric("topTranscripts",
           function(cds,
                    ntop = 10)
             standardGeneric("topTranscripts"))

setMethod("topTranscripts",
          "rnaCrosslinkDataSet",
          function(cds,
                   ntop = 10)  {
            ci = group(cds)[["c"]]
            si = group(cds)[["s"]]
            
            
            
            
            
            df= data.frame()
            
            for (i in names(InputFiles(cds)[["all"]][["all"]])) {
              
              
              same = which(InputFiles(cds)[["all"]][["all"]][[i]]$V4 ==
                             InputFiles(cds)[["all"]][["all"]][[i]]$V10 )
              notSame = which(InputFiles(cds)[["all"]][["all"]][[i]]$V4 !=
                                InputFiles(cds)[["all"]][["all"]][[i]]$V10)
              
              
              s =   paste(InputFiles(cds)[["all"]][["all"]][[i]]$V4[same],
                          InputFiles(cds)[["all"]][["all"]][[i]]$V10[same], sep = "::")
              s = unlist(lapply(s, function(x)  
                unlist(strsplit(x,split = "::"))[1] ))
              s = as.data.frame(table(s))
              s$type="intra"
              s$sample = i
              colnames(s) = c("RNA","reads","type","sample")
              
              
              
              t =   paste(InputFiles(cds)[["all"]][["all"]][[i]]$V4[notSame],
                          InputFiles(cds)[["all"]][["all"]][[i]]$V10[notSame], sep = "::")
              
              t1 = unlist(lapply(t, function(x)  
                unlist(strsplit(x,split = "::"))[1] ))
              t2 = unlist(lapply(t, function(x)  
                unlist(strsplit(x,split = "::"))[2] ))
              t = c(t1,t2)
              t = as.data.frame(table(t))
              t$type="inter"
              t$sample = i
              colnames(t) = c("RNA","reads","type","sample")
              
              df =rbind.data.frame(df,t,s)
            }
            
            df = dcast(df,formula = RNA + type ~ sample,   value.var = "reads")
            
            
            
            ci = sampleTable(cds)[ci,"sampleName"]
            si = sampleTable(cds)[si,"sampleName"]
            
            if(length(si) == 1){
              df$samples = df[,si]
              df$control = df[,ci]
              df$enrichment = df$samples / df$control
              
            }else{
              
              df$samples = rowSums(df[,si],na.rm = TRUE)
              df$control = rowSums(df[,ci],na.rm = TRUE)
              df$enrichment = df$samples / df$control
              
            }
            
            df = df[order(df$samples,decreasing = TRUE),]
            df = df[1:ntop,]
            
            
            df
            
            
            
          })





#' topInteracters
#'
#' This method prints the top transcripts that have the most duplexes
#' assigned that interact with the transcript of interest 
#'
#' @param cds a \code{rnaCrosslinkDataSet} object
#' @param ntop the number of entries to display
#' @param sds known bug, doesn't work for small data sets fix incoming
#'
#' @name topInteracters
#' @docType methods
#' @rdname topInteracters
#' @aliases topInteracters,rnaCrosslinkDataSet-method
#' @return A table, the number of counts per sample per interacting transcript
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' topInteracters(cds, sds = TRUE)
#'
#' @export
setGeneric("topInteracters",
           function(cds,
                    ntop = 10,
                    sds = TRUE)
             standardGeneric("topInteracters"))

setMethod("topInteracters",
          "rnaCrosslinkDataSet",
          function(cds,
                   ntop = 10,
                   sds = TRUE)  {
            if(sds != TRUE){
            ci = group(cds)[["c"]]
            si = group(cds)[["s"]]
            
            
            
            
            
            df= data.frame()
            
            for (i in names(InputFiles(cds)[[rnas(cds)]][["noHost"]])) {
              iFile = InputFiles(cds)[[rnas(cds)]][["noHost"]][[i]]
              
              same = which(iFile$V4 == iFile$V10 )
              notSame = which(iFile$V4 != iFile$V10 )
              
              
              s =   iFile$V4[same]

              s = as.data.frame(table(s))
              s$type="intra"
              s$sample = i
              colnames(s) = c("RNA","reads","type","sample")
              
              
              v4rna = which(iFile$V4[notSame] == rnas(cds))
              t1 = iFile$V4[notSame][-v4rna]
              t2 = iFile$V10[notSame][v4rna]
              t = c(t1,t2)
              
              
              if(length(notSame) == 0){           
                df =rbind.data.frame(df,t,c(NA,NA, "inter", i),stringsAsFactors = F)
              }else{
                
                t = as.data.frame(table(t),stringsAsFactors = F)
                t$type="inter"
                t$sample = i
                colnames(t) = c("RNA","reads","type","sample")
              }
              if(length(same) == 0){           df =rbind.data.frame(df,s,c(NA,NA, "inter", i),stringsAsFactors = F)
              }else{
                df =rbind.data.frame(df,t,s,stringsAsFactors = F)
              }
            }
            
            df$reads = as.numeric(df$reads)
            
            
            df = dcast(df,formula = RNA + type ~ sample,   value.var = "reads")
            
            
            
            ci = sampleTable(cds)[ci,"sampleName"]
            si = sampleTable(cds)[si,"sampleName"]
            
            if(length(si) == 1){
              df$samples = df[,si]
              df$control = df[,ci]
              df$enrichment = df$samples / df$control
              
            }else{
              
              df$samples = rowSums(df[,si],na.rm = TRUE)
              df$control = rowSums(df[,ci],na.rm = TRUE)
              df$enrichment = df$samples / df$control
              
            }
            
            df = df[order(df$samples,decreasing = TRUE),]
            df = df[1:ntop,]
            
            
            df
            }
          })



#' getInteractions
#'
#' This method returns a table of interactions of an RNA (interactor) on the RNA of interest. 
#'
#' @param cds a \code{rnaCrosslinkDataSet} object
#' @param interactors A vector containing the names of RNAs to show interactions with
#'
#'
#' @name getInteractions
#' @docType methods
#' @rdname getInteractions
#' @return A table showing the read coverage of the specified interacting RNAs 
#' @aliases getInteractions,rnaCrosslinkDataSet-method
#' @examples
#' cds = makeExamplernaCrosslinkDataSet()
#' getInteractions(cds, c("transcript1","transcript2"))
#'
#' @export
setGeneric("getInteractions",
           function(cds,
                    interactors)
             standardGeneric("getInteractions"))

setMethod("getInteractions",
          "rnaCrosslinkDataSet",
          function(cds,
                   interactors)  {
            seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
            # subset the Input files based on the interactor of choice
            table = data.frame()
            for (interactor in interactors){
              for (i in names(InputFiles(cds)[[rnas(cds)]][["host"]])) {
                x = InputFiles(cds)[[rnas(cds)]][["host"]][[i]]
                x = x[x$V10 == interactor, ]
                if (nrow(x) == 0) {
                  v = as.data.frame(0)
                } 
                else {
                  starts = x$V7
                  ends = x$V8
                  v = as.data.frame(unlist(seq2(from = starts, to = ends)))
                  if(ncol(v) > 1){
                    next
                  }
                }
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
#' @param cds a \code{rnaCrosslinkDataSet} object
#' @param interactor The rna to show interactions with
#'
#' @name getReverseInteractions
#' @docType methods
#' @rdname getReverseInteractions
#' @aliases getReverseInteractions,rnaCrosslinkDataSet-method
#' @return A long format table shoing the read coverage of chosen RNA
#' @examples
#' cds = makeExamplernaCrosslinkDataSet()
#' getReverseInteractions(cds, 'transcript2')
#'
#' @export
setGeneric("getReverseInteractions",
           function(cds,
                    interactor)
             standardGeneric("getReverseInteractions"))

setMethod("getReverseInteractions",
          "rnaCrosslinkDataSet",
          function(cds,
                   interactor)  {
            seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
            
            # subset the Input files based on the interacter of choice
            table = data.frame()
            for (i in names(InputFiles(cds)[[rnas(cds)]][["host"]])) {
              x = InputFiles(cds)[[rnas(cds)]][["host"]][[i]]
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





#' swapInputs
#'
#' Swap the table to ensure that 3 prime most duplex side is on he left of the
#' table used to make one sides heatmaps and other reasons where having the
#' left of the table coming after the right side is a problem. Different from
#' swapInputs as it ensure that BOTH duplex sides originate from the RNA of
#' interest.
#'
#' @param InputList the original InputList created with readInputFiles or subsetInputList
#' @param rna The rna of interest
#'
#' @name swapInputs
#' @docType methods
#' @rdname swapInputs
#' @return A list of "swapped" Input datas
#'
swapInputs = function(InputList, 
                      rna) {
  InputListSwap = InputList
  for (Input in 1:length(InputList)) {
    InputListS  =  InputList[[Input]][InputList[[Input]]$V4 == rna |
                                        InputList[[Input]]$V10 == rna, ]
    
    
    
    comb = InputListS
    
    
    
    
    InputListSwap[[Input]] = comb
  }
  
  return(InputListSwap)
}



#' swapInputs2
#'
#' Swap the table to ensure that 3 prime most duplex side is on the left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Different from swapInputs as it
#' ensure that BOTH duplex sides originate from the RNA of interest.
#' @param InputList the original InputList created with readInputFiles or subsetInputList
#' @param rna The rna of interest
#' @name swapInputs2
#' @docType methods
#' @rdname swapInputs2
#' @return A list of "swapped" Input data
#'
swapInputs2 = function(InputList, 
                       rna) {
  for (Input in 1:length(InputList)) {
    InputList18S = InputList[[Input]][as.character(InputList[[Input]]$V4) == rna &
                                        as.character(InputList[[Input]]$V10) == rna, ]
    tmp1 = InputList18S[InputList18S$V7 < InputList18S$V13, ]
    tmp2 = InputList18S[InputList18S$V7 > InputList18S$V13, ]
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
    InputList[[Input]] = comb
  }
  
  return(InputList)
}




#' swapInputs3
#'
#' Swap the table to ensure that 3 prime most duplex side is ont he left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Different from swapInputs as it
#' ensure that BOTH duplex sides originate from the RNA of interest.
#' @param InputList the original InputList created with readInputFiles or subsetInputList
#' @param rna The rna of interest
#'
#' @return A list of "swapped" Input datas
#' @name swapInputs3
#' @docType methods
#' @rdname swapInputs3
swapInputs3 = function(InputList, rna) {
  InputListSwap = InputList
  for (Input in 1:length(InputList)) {
    InputListS  =  InputList[[Input]][InputList[[Input]]$V4 == rna |
                                        InputList[[Input]]$V10 == rna, ]
    tmp1 = InputListS[InputListS$V4 == rna, ]
    tmp2 = InputListS[InputListS$V4 != rna, ]
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
    InputListSwap[[Input]] = comb
  }
  return(InputListSwap)
}









#' makeExamplernaCrosslinkDataSet
#'
#' Creat a minimal example rnaCrosslinkdataSetObject
#'
#' @return An example rnaCrosslinkDataSet objext
#' @name makeExamplernaCrosslinkDataSet
#' @rdname makeExamplernaCrosslinkDataSet
#' @examples
#' cds = makeExamplernaCrosslinkDataSet()
#' @export
makeExamplernaCrosslinkDataSet = function() {
  
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
  c7 = sample(1:5, 50, replace = TRUE)
  c8 = sample(20:25, 50, replace = TRUE)
  c13 = sample(20:25, 50, replace = TRUE)
  c14 = sample(40:45, 50, replace = TRUE)
  # long distance 50
  c7 = c(c7,sample(1:5, 50, replace = TRUE))
  c8 = c(c8,sample(20:25, 50, replace = TRUE))
  c13 = c(c13,sample(60:70, 50, replace = TRUE))
  c14 = c(c14,sample(80:83, 50, replace = TRUE))
  # inter RNA 100
  c7 = c(c7,sample(1:5, 100, replace = TRUE))
  c8 = c(c8,sample(20:25, 100, replace = TRUE))
  c13 = c(c13,sample(1:5, 100, replace = TRUE))
  c14 = c(c14,sample(20:25, 100, replace = TRUE))
  
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
  write.table(exampleInput,
              file = file, 
              quote = FALSE,
              row.names = FALSE, 
              sep = "\t", col.names = F)
  
  
  
  
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
  c7 = sample(1:5, 55, replace = TRUE)
  c8 = sample(20:25, 55, replace = TRUE)
  c13 = sample(20:25, 55, replace = TRUE)
  c14 = sample(40:45, 55, replace = TRUE)
  
  # inter RNA 100
  c7 = c(c7,sample(1:40, 90, replace = TRUE))
  c8 = c(c8,sample(20:75, 90, replace = TRUE))
  c13 = c(c13,sample(1:40, 90, replace = TRUE))
  c14 = c(c14,sample(20:75, 90, replace = TRUE))
  
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
  write.table(exampleInput,
              file = file2, 
              quote = FALSE, 
              row.names = FALSE, 
              sep = "\t",
              col.names = F)
  
  
  
  
  
  
  
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
  cds = rnaCrosslinkDataSet(rnas = rna,
                            rnaSize = 0,
                            sampleTable = sampleTable2)
  
  return(cds)
}








