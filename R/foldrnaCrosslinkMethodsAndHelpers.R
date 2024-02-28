#' @include  rnaCrosslinkDataSet.R
NULL


#' findBasePairsRNAcoFold2
#'
#' Folds the clusters using Vienna RNAfold
#'
#' @param  startPos1  Start of the cluster side x
#' @param  endPos1  End of the cluster side x
#' @param  seq1 Sequence of x
#' @param  startPos2 Start of the cluster side y
#' @param  endPos2 End of the cluster side y
#' @param  seq2 Sequence of y
#' @param  fasta \code{rnaRefs}
#' @param  shape shape reactivities
#'
#' @return A table of clusters and coordinates with folds
#' @name findBasePairsRNAcoFold2
#' @docType methods
#' @rdname findBasePairsRNAcoFold2
#'
findBasePairsRNAcoFold2 = function(startPos1,
                                   endPos1,
                                   seq1,
                                   startPos2,
                                   endPos2,
                                   seq2,
                                   fasta,
                                   shape) {
  if (length(shape) == 1) {
    table = data.frame(stringsAsFactors = FALSE)
    
    As = 20
    #make a new sequence with AAA as the middle
    seq = paste(seq1, paste(rep("A", As), collapse = ""), seq2, sep = "")
    
    #make the constraints, you want the 40 A's not to be paired
    start = nchar(seq1) + 1
    
    #P i 0 k
    row = c("P", start, 0, As)
    
    csfile = tempfile()
    write.table(
      t(as.data.frame(row)),
      file = csfile,
      quote = F,
      row.names = F,
      col.names = F
    )
    
    
    command = paste(
      "echo \">",
      startPos1,
      "-",
      endPos1,
      "\n",
      seq,
      "\" | RNAfold  --noPS --constraint=",csfile,
      sep = ""
    )
    
    
    #create vienna command and run capturing output
     x = system(command, intern = T)
    
    
    #extract the vienna and the start locations
    vienna = sub(" .*", "", x[3])
    
    helix = viennaToHelix(vienna)
    
    # to anything below a certain number we want to add startPos1 -1
    # and anything over a certain number we want to add startpos2 -2
    cutoff = nchar(seq1) + As
    helix2 = helix
    for (i in c(1, 2)) {
      for (j in 1:nrow(helix)) {
        if (helix[j, i] > cutoff) {
          helix2[j, i] = (helix2[j, i] +  startPos2) - (As + 1 + nchar(seq1))
        } else{
          helix2[j, i] = helix2[j, i] + startPos1 - 1
        }
        
      }
    }
    helix = helix2
    helix = helix[, -c(3, 4)]
    helix$rl = helix$i
    helix$rr = helix$j
    helix$pl = fasta[[1]][[1]][helix$rl]
    helix$pr = fasta[[1]][[1]][helix$rr]
    helix$veinna = vienna
    helix$seq1new = seq1
    helix$seq2new = seq2
    colnames(helix) = c("l",
                        "r",
                        "rl",
                        "rr",
                        "pl",
                        "pr",
                        "vienna",
                        "seq1new",
                        "seq2new")
    helix$check = 1
    aggTable = aggregate(
      helix$check,
      by = list(
        helix$rl,
        helix$rr,
        helix$pl,
        helix$pr,
        helix$vienna,
        helix$seq1new,
        helix$seq2new
      ),
      FUN = sum
    )
    return(aggTable)
  }   else{
    As = 15
    shape1 = shape[startPos1:endPos1]
    shape2 = shape[startPos2:endPos2]
    length1 = (endPos1 - startPos1) + 1
    length2 = (endPos2 - startPos2) + 1
    

    shapeTable = data.frame("x" = c(1:(length1), (length1 + As + 1):((length1 +
                                                                        As) + length2)),
                            "y" = c(shape1, shape2))
    shapeTable = shapeTable[complete.cases(shapeTable), ]
    

    sfile = tempfile()
    write.table(
      shapeTable,
      sfile,
      quote = F,
      row.names = F,
      col.names = F
    )
    
    table = data.frame(stringsAsFactors = FALSE)
    
    
    #make a new sequence with AAA as the middle
    seq = paste(seq1, paste(rep("A", As), collapse = ""), seq2, sep = "")
    
    #make the constraints, you want the 40 A's not to be paired
    start = nchar(seq1) + 1
    
    #P i 0 k
    row = c("P", start, 0, As)
    csfile = tempfile()
    write.table(
      t(as.data.frame(row)),
      file = csfile,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    
    command = paste(
      "echo \">",
      startPos1,
      "-",
      endPos1,
      "\n",
      seq,
      "\" | RNAfold  --noPS --constraint=",csfile," --shape=",sfile,
      sep = ""
    )
    
    
    #create vienna command and run capturing output
    #command = paste("echo \">",startPos1,"-", endPos1,"\n",seq1,"\n>",startPos2,"-",endPos2,"\n",seq2,"\" | RNAduplex", sep = "")
    x = system(command, intern = TRUE)
    
    
    #extract the vienna and the start locations
    vienna = sub(" .*", "", x[3])
    
    helix = viennaToHelix(vienna)
    
    # to anything below a certain number we want to add startPos1 -1
    # and anything over a certain number we want to add startpos2 -2
    cutoff = nchar(seq1) + As
    helix2 = helix
    for (i in c(1, 2)) {
      for (j in 1:nrow(helix)) {
        if (helix[j, i] > cutoff) {
          helix2[j, i] = (helix2[j, i] +  startPos2) - (As + 1 + nchar(seq1))
        } else{
          helix2[j, i] = helix2[j, i] + startPos1 - 1
        }
        
      }
    }
    helix = helix2
    helix = helix[, -c(3, 4)]
    helix$rl = helix$i
    helix$rr = helix$j
    helix$pl = fasta[[1]][[1]][helix$rl]
    helix$pr = fasta[[1]][[1]][helix$rr]
    helix$veinna = vienna
    helix$seq1new = seq1
    helix$seq2new = seq2
    colnames(helix) = c("l",
                        "r",
                        "rl",
                        "rr",
                        "pl",
                        "pr",
                        "vienna",
                        "seq1new",
                        "seq2new")
    helix$check = 1
    aggTable = aggregate(
      helix$check,
      by = list(
        helix$rl,
        helix$rr,
        helix$pl,
        helix$pr,
        helix$vienna,
        helix$seq1new,
        helix$seq2new
      ),
      FUN = sum
    )
    return(aggTable)
  }
  
}











#' findBasePairsRNAfold
#'
#' Folds the clusters using Vienna RNA duplex
#'
#' @param  startPos  Start of the cluster side x
#' @param  endPos  End of the cluster side x
#' @param  seqs Sequence of x
#' @param  fasta \code{rnaRefs}
#' @param  shape shape reactivities
#'
#' @return A table of clusters and coordinates with folds
#' @name findBasePairsRNAfold
#' @docType methods
#' @rdname findBasePairsRNAfold
#'
findBasePairsRNAfold = function(startPos, 
                                endPos, 
                                seqs, 
                                fasta, 
                                shape) {
  if (length(shape) == 1) {
    table = data.frame(stringsAsFactors = FALSE)
    
    
    #  for(i in 1:10){
    #command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold --ImFeelingLucky", sep = "")
    command = paste("echo \">",
                    startPos,
                    "-",
                    endPos,
                    "\n",
                    seqs,
                    "\" | RNAfold   --noPS",
                    sep = "")
    x = system(command, intern = TRUE)
    
    #extract the vienna and the start locations
    vienna = sub(" .*", "", x[3])
    
    helix = viennaToHelix(vienna)
    
    helix = helix[, -c(3, 4)]
    helix$rl = helix$i + startPos - 1
    helix$rr = helix$j + startPos  - 1
    helix$pl = fasta[[1]][[1]][helix$rl]
    helix$pr = fasta[[1]][[1]][helix$rr]
    helix$veinna = vienna
    helix$seq1new = seqs
    helix$seq2new = ""
    colnames(helix) = c("l",
                        "r",
                        "rl",
                        "rr",
                        "pl",
                        "pr",
                        "vienna",
                        "seq1new",
                        "seq2new")
    helix$check = 1
    aggTable = aggregate(
      helix$check,
      by = list(
        helix$rl,
        helix$rr,
        helix$pl,
        helix$pr,
        helix$vienna,
        helix$seq1new,
        helix$seq2new
      ),
      FUN = sum
    )
    return(aggTable)
    
  } else{
    shape2 = shape[startPos:endPos]
    length = (endPos - startPos) + 1
    shapeTable = data.frame("x" = 1:length,
                            "y" = shape2)
    shapeTable = shapeTable[complete.cases(shapeTable), ]

    sfile = tempfile()
    write.table(
      shapeTable,
      sfile,
      quote = F,
      row.names = F,
      col.names = F
    )
    table = data.frame(stringsAsFactors = TRUE)
    
    
    #  for(i in 1:10){
    #command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold --ImFeelingLucky", sep = "")
    command = paste(
      "echo \">",
      startPos,
      "-",
      endPos,
      "\n",
      seqs,
      "\" | RNAfold   --noPS --shape=",sfile ,
      sep = ""
    )
    x = system(command, intern = TRUE)
    
    #extract the vienna and the start locations
    vienna = sub(" .*", "", x[3])
    
    helix = viennaToHelix(vienna)
    
    helix = helix[, -c(3, 4)]
    helix$rl = helix$i + startPos - 1
    helix$rr = helix$j + startPos  - 1
    helix$pl = fasta[[1]]$NR_003286.4[helix$rl]
    helix$pr = fasta[[1]]$NR_003286.4[helix$rr]
    helix$veinna = vienna
    helix$seq1new = seqs
    helix$seq2new = ""
    colnames(helix) = c("l",
                        "r",
                        "rl",
                        "rr",
                        "pl",
                        "pr",
                        "vienna",
                        "seq1new",
                        "seq2new")
    helix$check = 1
    aggTable = aggregate(
      helix$check,
      by = list(
        helix$rl,
        helix$rr,
        helix$pl,
        helix$pr,
        helix$vienna,
        helix$seq1new,
        helix$seq2new
      ),
      FUN = sum
    )
    return(aggTable)
    
  }
}




#' getClusterClusterShortRangeWhole
#'
#' Decides if a cluster is long or short range
#' then either grabs the whole sequence or the sequence of the two sides of the
#' interaction separately.
#' 
#' @param  cluster cluster positions
#' @param  seq sequence of transcript
#' @name getClusterClusterShortRangeWhole
#' @docType methods
#' @rdname getClusterClusterShortRangeWhole
#' @return The same table with an extra column
#'
getClusterClusterShortRangeWhole = function(cluster, 
                                            seq) {
  fasta = seq
  if (cluster$rs - cluster$le <= 30) {
    coords = list()
    coords[[1]] = cluster$ls
    coords[[2]] = cluster$re
    seqs1 = paste(fasta[[1]][coords[[1]]:coords[[2]]], collapse = "")
    type = "short"
    seqs2 = ""
    return(list(type, seqs1, seqs2, cluster))
  } else{
    coords = list()
    coords[[1]] = cluster$ls
    coords[[2]] = cluster$le
    coords[[3]] = cluster$rs
    coords[[4]] = cluster$re
    seqs1 = paste(fasta[[1]][coords[[1]]:coords[[2]]], collapse = "")
    seqs2 = paste(fasta[[1]][coords[[3]]:coords[[4]]], collapse = "")
    type = "long"
    return(list(type, seqs1, seqs2, cluster))
  }
}






#' findBasePairsRNAfold2
#'
#' Folds the clusters using Vienna RNA duplex
#'
#' @param  startPos  Start of the cluster side x
#' @param  endPos  End of the cluster side x
#' @param  seqs Sequence of x
#' @param  fasta \code{rnaRefs}
#' @name findBasePairsRNAfold2
#' @docType methods
#' @rdname findBasePairsRNAfold2
#' @return A table of clusters and coordinates with folds
#'
findBasePairsRNAfold2 = function(startPos, 
                                 endPos, 
                                 seqs, 
                                 fasta) {
  table = data.frame(stringsAsFactors = FALSE)
  #  for(i in 1:10){
  #command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold --ImFeelingLucky", sep = "")
  command = paste("echo \">",
                  startPos,
                  "-",
                  endPos,
                  "\n",
                  seqs,
                  "\" | RNAfold  --noPS",
                  sep = "")
  x = system(command, intern = T)
  
  #extract the vienna and the start locations
  vienna = sub(" .*", "", x[3])
  
  helix = viennaToHelix(vienna)
  
  helix = helix[, -c(3, 4)]
  helix$rl = helix$i + startPos - 1
  helix$rr = helix$j + startPos  - 1
  helix$pl = fasta[[1]]$NR_003286.4[helix$rl]
  helix$pr = fasta[[1]]$NR_003286.4[helix$rr]
  helix$veinna = vienna
  helix$seq1new = seqs
  helix$seq2new = ""
  colnames(helix) = c("l", "r", "rl", "rr", "pl", "pr", "vienna", "seq1new", "seq2new")
  helix$check = 1
  aggTable = aggregate(
    helix$check,
    by = list(
      helix$rl,
      helix$rr,
      helix$pl,
      helix$pr,
      helix$vienna,
      helix$seq1new,
      helix$seq2new
    ),
    FUN = sum
  )
  return(aggTable)
}








#' compareKnownStructures
#'
#' This method compares the predicted structures to a set of known interactions
#' 
#'
#' @param foldedCds a \code{rnaCrosslinkDataSet} object
#'
#' @param  foldedCds  \code{rnaCrosslinkDataSet} after running foldrnaCrosslink
#' @param  file a two column file with column header i and j with numeric values showing
#' nucleoide i binds to nucleotide j
#'
#'
#' @return Returns a dataframe
#' @name compareKnownStructures
#' @docType methods
#' @rdname compareKnownStructures
#' @aliases compareKnownStructures,rnaCrosslinkDataSet-method
#' @return a tables showing the number of predicted interactions and their agreement
#' @examples
#' \dontrun{
#' cds = makeExamplernaCrosslinkDataSet()
#' clusteredCds = clusterrnaCrosslink(cds = cds,
#'                                cores = 3,
#'                                stepCount = 2,
#'                                clusterCutoff = 1)
#' 
#' 
#' trimmedClusters = trimClusters(clusteredCds = clusteredCds,trimFactor = 1, clusterCutoff = 1)
#' 
#' 
#' 
#' fasta = paste(c(rep('A',25),
#'                 rep('T',25),
#'                 rep('A',10),
#'                 rep('T',23)),collapse = "")
#' 
#' header = '>transcript1'
#' 
#' 
#' fastaFile = tempfile()
#' writeLines(paste(header,fasta,sep = "\n"),con = fastaFile)
#' 
#' 
#' rnaRefs = list()
#' rnaRefs[[rnas(cds)]] = read.fasta(fastaFile)
#' rnaRefs
#' 
#' 
#' 
#' foldedCds = foldrnaCrosslink(trimmedClusters,
#'                          rnaRefs = rnaRefs,
#'                          start = 1,
#'                          end = 83,
#'                          shape = 0,
#'                          ensembl = 5,
#'                          constraintNumber  = 1,
#'                          evCutoff = 1)
#' 
#' 
#' # make an example table of "know" interactions
#' file = data.frame(V1 = c(6), 
#'                   V2 = c(80))
#' compareKnownStructures(foldedCds,file)
#' }
#' @export
#'
setGeneric("compareKnownStructures",
           function(foldedCds, 
                    file)
             standardGeneric("compareKnownStructures"))

setMethod("compareKnownStructures",
          "rnaCrosslinkDataSet",
          function(foldedCds, 
                   file)  {
            known18S = file
            
            
            start = as.numeric(sub(":.*", "", names(foldedCds@viennaStructures)))
            end = as.numeric(sub(".*:", "", names(foldedCds@viennaStructures)))
            
            
            viennas = unlist(foldedCds@viennaStructures[[1]])
            
            length = end - start + 1
            
            
            ## make a helix variable as well.
            helixes = list()
            for (i in 1:length(viennas)) {
              helixes[[i]] = viennaToHelix(viennas[i])
              helixes[[i]]$i = helixes[[i]]$i + start - 1
              helixes[[i]]$j = helixes[[i]]$j + start - 1
              
            }
            
            
            # for each helix retrieve the sensitvity and specificity
            
            known18S = known18S[known18S$V1 > start &
                                  known18S$V1 < end &
                                  known18S$V2 > start &
                                  known18S$V2 < end, ]

            known18SID = paste(known18S$V1, known18S$V2, sep = "-")
            viennaScores = data.frame(stringsAsFactors = FALSE)
            for (i in 1:length(viennas)) {
              structureIDs = paste(helixes[[i]]$i, helixes[[i]]$j, sep = "-")
              row = c(length(which(structureIDs %in% known18SID)), length(structureIDs),
                      length(known18SID))
              viennaScores = rbind.data.frame(viennaScores, row)
              
            }
            
            colnames(viennaScores) = c("positive", "totalStructure", "totalKnown")
            
            # FP FN - TN TP
            viennaScores$FP = viennaScores$totalStructure - viennaScores$positive
            viennaScores$FN = viennaScores$totalKnown - viennaScores$positive
            
            viennaScores$TN = ((length * length) - ((viennaScores$totalKnown -  viennaScores$positive) + viennaScores$totalStructure
            ))
            viennaScores$TP = viennaScores$positive
            
            
            viennaScores$sensitivity = viennaScores$TP / (viennaScores$TP + viennaScores$FN)
            viennaScores$precision = viennaScores$positive / (viennaScores$TP + viennaScores$FP)
            viennaScores$structureID = names(unlist(foldedCds@viennaStructures[[1]]))
            viennaScores
            
          })





#' plotEnsemblePCA
#'
#' This method plots a PCA of the ensembl
#'
#'
#' @param  foldedCds  \code{rnaCrosslinkDataSet} after running foldrnaCrosslink
#' @param labels plot with labels or not (TRUE/FALSE)
#' @param split split the plot using facets based on the samples  (TRUE/FALSE)
#'
#' @name plotEnsemblePCA
#' @docType methods
#' @rdname plotEnsemblePCA
#' @aliases plotEnsemblePCA,rnaCrosslinkDataSet-method
#' @return a PCA plot of the ensembl
#' @examples
#' \dontrun{
#' cds = makeExamplernaCrosslinkDataSet()
#' clusteredCds = clusterrnaCrosslink(cds = cds,
#'                                cores = 3,
#'                                stepCount = 2,
#'                                clusterCutoff = 1)
#' 
#' 
#' trimmedClusters = trimClusters(clusteredCds = clusteredCds,trimFactor = 1, clusterCutoff = 1)
#' 
#' 
#' 
#' fasta = paste(c(rep('A',25),
#'                 rep('T',25),
#'                 rep('A',10),
#'                 rep('T',23)),collapse = "")
#' 
#' header = '>transcript1'
#' 
#' 
#' fastaFile = tempfile()
#' writeLines(paste(header,fasta,sep = "\n"),con = fastaFile)
#' 
#' 
#' rnaRefs = list()
#' rnaRefs[[rnas(cds)]] = read.fasta(fastaFile)
#' rnaRefs
#' 
#' 
#' 
#' foldedCds = foldrnaCrosslink(trimmedClusters,
#'                          rnaRefs = rnaRefs,
#'                          start = 1,
#'                          end = 83,
#'                          shape = 0,
#'                          ensembl = 5,
#'                          constraintNumber  = 1,
#'                          evCutoff = 1)
#' 
#' 
#' plotEnsemblePCA(foldedCds)
#' }
#' @export


setGeneric("plotEnsemblePCA",
           function(foldedCds,
                    labels = TRUE,
                    split = TRUE)
             standardGeneric("plotEnsemblePCA"))

setMethod("plotEnsemblePCA",
          "rnaCrosslinkDataSet",
          function(foldedCds,
                   labels = TRUE,
                   split = TRUE)  {
            discordMat = matrix(0, ncol = length(unlist(foldedCds@viennaStructures[[1]])), nrow = length(unlist(foldedCds@viennaStructures[[1]])))
            
            for (i in 1:length(unlist(foldedCds@viennaStructures[[1]]))) {
              for (j in 1:length(unlist(foldedCds@viennaStructures[[1]]))) {
                a = viennaToHelix(unlist(foldedCds@viennaStructures[[1]])[i],-31.5)
                b = viennaToHelix(unlist(foldedCds@viennaStructures[[1]])[j],-31.5)
                
                disI = length(which(!(paste(a$i, a$j) %in% paste(b$i, b$j))))
                disJ = length(which(!(paste(b$i, b$j) %in% paste(a$i, a$j))))
                
                dis = disI + disJ
                discordMat[i, j] = dis
                
                
              }
            }
            row.names(discordMat) = names(unlist(foldedCds@viennaStructures[[1]]))
            colnames(discordMat) = names(unlist(foldedCds@viennaStructures[[1]]))
            
            
            # perform a PCA on the data in assay(x) for the selected genes
            pca <- prcomp(t(discordMat))
            
            # the contribution to the total variance for each component
            percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
            
            
            # assembly the data for the plot
            d <- data.frame(PCa = pca$x[, 1], PCb = pca$x[, 2])
            d$sample = substr(names(unlist(foldedCds@viennaStructures[[1]])), 1, 2)
            d$dgs = as.numeric(as.character(unlist(foldedCds@dgs[[1]])))
            
            
            
            p1 = ggplot(data = d, aes(x = PCa, y = PCb, colour = sample)) +
              geom_point(mapping = aes(size = dgs), alpha = 0.2) +
              xlab(paste0("PC", 1, ": ", round(percentVar[1] * 100), "% variance")) +
              ylab(paste0("PC", 2, ": ", round(percentVar[2] * 100), "% variance")) +
              theme_classic() 
              
              if (split == T) {
                p1 = p1 +   facet_grid(sample ~ .)
              }
            
            if (labels == T) {
              p1 = p1 +    geom_text_repel(
                aes(label = row.names(d)),
                color = "gray20",
                data = d,
                force = 10,max.overlaps = 400
              )
              
            }
            plot(p1)
            
})






#' plotComparisonArc
#'
#' This method plots two structures chosen from the
#' plotEnsemblePCA method
#'
#' @param  foldedCds  \code{rnaCrosslinkDataSet} after running foldrnaCrosslink
#' @param  s1 sample of structure 1
#' @param  s2 sample of structure 2
#' @param  n1 number of structure from first sample
#' @param  n2 number of structure from first sample
#'
#' @name plotComparisonArc
#' @docType methods
#' @rdname plotComparisonArc
#' @aliases plotComparisonArc,rnaCrosslinkDataSet-method
#' @return an ark diagram of the two strcutures
#' @examples
#' \dontrun{
#' cds = makeExamplernaCrosslinkDataSet()
#' clusteredCds = clusterrnaCrosslink(cds = cds,
#'                                cores = 3,
#'                                stepCount = 2,
#'                                clusterCutoff = 1)
#' 
#' 
#' trimmedClusters = trimClusters(clusteredCds = clusteredCds,trimFactor = 1, clusterCutoff = 1)
#' 
#' 
#' 
#' fasta = paste(c(rep('A',25),
#'                 rep('T',25),
#'                 rep('A',10),
#'                 rep('T',23)),collapse = "")
#' 
#' header = '>transcript1'
#' 
#' 
#' fastaFile = tempfile()
#' writeLines(paste(header,fasta,sep = "\n"),con = fastaFile)
#' 
#' 
#' rnaRefs = list()
#' rnaRefs[[rnas(cds)]] = read.fasta(fastaFile)
#' rnaRefs
#' 
#' 
#' 
#' foldedCds = foldrnaCrosslink(trimmedClusters,
#'                          rnaRefs = rnaRefs,
#'                          start = 1,
#'                          end = 83,
#'                          shape = 0,
#'                          ensembl = 5,
#'                          constraintNumber  = 1,
#'                          evCutoff = 1)
#' 
#' 
#' plotComparisonArc(foldedCds,"s1","s1",1,3)
#' }
#' @export
setGeneric("plotComparisonArc",
           function(foldedCds,
                    s1 = "s1",
                    s2 = "s2",
                    n1 = 1,
                    n2 = 2)
             standardGeneric("plotComparisonArc"))

setMethod("plotComparisonArc",
          "rnaCrosslinkDataSet",
          function(foldedCds,
                   s1 = "s1",
                   s2 = "s2",
                   n1 = 1,
                   n2 = 2)  {
            #now try and plot two different structures in ark
            
            a = viennaToHelix(foldedCds@viennaStructures[[1]][[s1]][n1])
            b = viennaToHelix(foldedCds@viennaStructures[[1]][[s2]][n2])
            
            a$value = as.numeric(as.character(foldedCds@dgs[[s1]][n1]))
            b$value = as.numeric(as.character(foldedCds@dgs[[s2]][n2]))
            
            a = expandHelix(a)
            b = expandHelix(b)
            
            
            plotDoubleHelix(a, b, line = TRUE, arrow = TRUE)
            
            
          })






#' plotStructure
#'
#' This method plots a structures chosen from the
#' plotEnsemblePCA method
#'
#' @param  foldedCds  \code{rnaCrosslinkDataSet} after running foldrnaCrosslink
#' @param  rnaRefs A fasta of the transcript (made with seqinr::read.fasta)
#' @param  s sample of structure
#' @param  n number of structure
#'
#'
#'
#'
#' @name plotStructure
#' @docType methods
#' @rdname plotStructure
#' @aliases plotStructure,rnaCrosslinkDataSet-method
#' @return a diagram of the predicted structure
#' @examples
#' \dontrun{
#' cds = makeExamplernaCrosslinkDataSet()
#' clusteredCds = clusterrnaCrosslink(cds = cds,
#'                                cores = 3,
#'                                stepCount = 2,
#'                                clusterCutoff = 1)
#' 
#' 
#' trimmedClusters = trimClusters(clusteredCds = clusteredCds,trimFactor = 1, clusterCutoff = 1)
#' 
#' 
#' 
#' fasta = paste(c(rep('A',25),
#'                 rep('T',25),
#'                 rep('A',10),
#'                 rep('T',23)),collapse = "")
#' 
#' header = '>transcript1'
#' 
#' 
#' fastaFile = tempfile()
#' writeLines(paste(header,fasta,sep = "\n"),con = fastaFile)
#' 
#' 
#' rnaRefs = list()
#' rnaRefs[[rnas(cds)]] = read.fasta(fastaFile)
#' rnaRefs
#' 
#' 
#' 
#' foldedCds = foldrnaCrosslink(trimmedClusters,
#'                          rnaRefs = rnaRefs,
#'                          start = 1,
#'                          end = 83,
#'                          shape = 0,
#'                          ensembl = 5,
#'                          constraintNumber  = 1,
#'                          evCutoff = 1)
#' 
#' 
#' plotStructure(foldedCds,rnaRefs,"s1",3)
#' }
#' @export


setGeneric("plotStructure",
           function(foldedCds,
                    rnaRefs,
                    s = "s1",
                    n = 1)
             standardGeneric("plotStructure"))

setMethod("plotStructure",
          "rnaCrosslinkDataSet",
          function(foldedCds,
                   rnaRefs,
                   s = "s1",
                   n = 1)  {
            ct = makeCt(foldedCds@viennaStructures[[1]][[s]][n],
                        paste(rnaRefs[[1]][[1]][as.numeric(as.character(sub(
                          ":.*", "", names(foldedCds@viennaStructures)
                        ))):as.numeric(as.character(sub(
                          ".*:", "", names(foldedCds@viennaStructures)
                        )))], collapse = ""))
            
            coord = ct2coord(ct)
            
            RNAPlot(
              coord,
              nt = TRUE,
              seqTF = TRUE,
              labTF = TRUE,
              tsize = 1
            )
            
            
          })
