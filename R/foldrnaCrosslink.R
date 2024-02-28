#' @include  rnaCrosslinkDataSet.R
NULL



#' foldrnaCrosslink
#'
#'
#' This methods folds an ensebl of structures for the whole RNA or chosen region
#' of the RNA. See \code{rnaCrosslinkDataSet} for slot information.
#'
#' @param cdsObject rnaCrosslinkDataSet object created with rnaCrosslinkDataSet
#' @param rnaRefs named List - a list with named elements that correspond to the
#'     .RNA of interest. The element of the list must be a fasta file that has
#'     been read with \code{seqinr::read.fasta()}
#' @param start Start of segmnent to fold
#' @param end End of segmnent to fold
#' @param evCutoff Mininum number of read support for contraint to be included in folding
#' @param ensembl Number of structures to Nake
#' @param constraintNumber Number of constraints to add to each final fold
#' @param shape shape reactivities (0 for no constraints)
#' @name foldrnaCrosslink
#' @docType methods
#' @rdname foldrnaCrosslink
#' @aliases foldrnaCrosslink,rnaCrosslinkDataSet-method
#' @return a rnaCrosslinkDataSet object
#' @examples 
#' \dontrun{
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusteredCds = clusterrnaCrosslink(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 0)
#'                 
#' trimmedClusters = trimClusters(clusteredCds = clusteredCds,
#'              trimFactor = 1, 
#'              clusterCutoff = 0)
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
#' foldedCds
#' }
#' @export
setGeneric("foldrnaCrosslink",
           function(cdsObject,
                    rnaRefs,
                    start,
                    end,
                    evCutoff = 1,
                    ensembl = 50,
                    constraintNumber = 20,
                    shape = 0)
             standardGeneric("foldrnaCrosslink"))

setMethod("foldrnaCrosslink",
          "rnaCrosslinkDataSet",
          function(cdsObject,
                   rnaRefs,
                   start,
                   end,
                   evCutoff = 1,
                   ensembl = 50,
                   constraintNumber = 20,
                   shape = 0)  {
            ###########################################################
            # Fold each cluster and add vienna to the table
            ###########################################################
            
            rna = rnas(cdsObject)
            ##############################
            # get trimmed cluster tables
            
            clusterPositionsListTrimmed = clusterTableList(cdsObject)[[rna]][["trimmedClusters"]]
            
            
            ##############################
            #make combined tables for the samples
            clusterPositionsListTrimmedSarsCombined = list()
            for (j in c(group(cdsObject)[["s"]],group(cdsObject)[["c"]])) {
              clusterPositionsListTrimmed[[j]]$sample = sampleTable(cdsObject)[j, "sampleName"]
              clusterPositionsListTrimmedSarsCombined = rbind.data.frame(
                clusterPositionsListTrimmedSarsCombined,
                clusterPositionsListTrimmed[[j]],
                stringsAsFactors = FALSE
              )
            }
            
            ##############################
            # add the sequences to the table
            # seq1 = left seq2 = right, type = short or long (range interaction)
            colnames(clusterPositionsListTrimmedSarsCombined)
            clusterPositionsListTrimmedSarsCombined$type = ""
            clusterPositionsListTrimmedSarsCombined$seq1 = ""
            clusterPositionsListTrimmedSarsCombined$seq2 = ""
            # add the sequences tot he table
            for (i in 1:nrow(clusterPositionsListTrimmedSarsCombined)) {
              x = getClusterClusterShortRangeWhole(clusterPositionsListTrimmedSarsCombined[i, ],
                                                   rnaRefs[[rna]])
              clusterPositionsListTrimmedSarsCombined$type[i] = x[[1]]
              clusterPositionsListTrimmedSarsCombined$seq2[i] = x[[3]]
              clusterPositionsListTrimmedSarsCombined$seq1[i] = x[[2]]
              
            }
            
            
            ##############################
            # now Fold
            
            tableAll = data.frame()
            
            
            for (i in 1:nrow(clusterPositionsListTrimmedSarsCombined)) {
              table = c()
              
              if (clusterPositionsListTrimmedSarsCombined[i, "type"] == "short") {
                table = findBasePairsRNAfold(
                  startPos = clusterPositionsListTrimmedSarsCombined[i, "ls"],
                  endPos = clusterPositionsListTrimmedSarsCombined[i, "re"],
                  seqs = clusterPositionsListTrimmedSarsCombined[i, "seq1"],
                  fasta = rnaRefs,
                  shape = shape
                )
              } else{
                table = findBasePairsRNAcoFold2(
                  startPos1 = clusterPositionsListTrimmedSarsCombined[i, "ls"],
                  endPos1 = clusterPositionsListTrimmedSarsCombined[i, "le"],
                  seq1 = clusterPositionsListTrimmedSarsCombined[i, "seq1"],
                  startPos2 = clusterPositionsListTrimmedSarsCombined[i, "rs"],
                  endPos2 = clusterPositionsListTrimmedSarsCombined[i, "re"],
                  seq2 = clusterPositionsListTrimmedSarsCombined[i, "seq2"],
                  fasta = rnaRefs,
                  shape = shape
                )
              }
              
              table$cluster = clusterPositionsListTrimmedSarsCombined[i, "id"]
              table$evidence = clusterPositionsListTrimmedSarsCombined[i, "size.x"]
              table$sample = clusterPositionsListTrimmedSarsCombined[i, "sample"]
              table$type = clusterPositionsListTrimmedSarsCombined[i, "type"]
              clusterPositionsListTrimmedSarsCombined$vienna[i] = table$Group.5[1]
              clusterPositionsListTrimmedSarsCombined$seq1new[i] = table$Group.6[1]
              clusterPositionsListTrimmedSarsCombined$seq2new[i] = table$Group.7[1]
              tableAll = rbind.data.frame(table, tableAll)
              
            }
            
            
            colnames(tableAll) = c(
              "p1",
              "p2",
              "nt1",
              "nt2",
              "vienna",
              "seq",
              "x1",
              "x2",
              "cluster",
              "evidence",
              "sample",
              "type"
            )
            clusterPositionsListTrimmedSarsCombinedWithStructures = clusterPositionsListTrimmedSarsCombined
            
            
            
            
            ############################################################################
            # Fold the whole molecule
            ############################################################################
            
            # get interactions
            interactionTable = tableAll
            #interactionTable = interactionTable[interactionTable$type == "long",]
            
            
            # Just get the columns needed
            interactionTable = interactionTable[, c(1, 2, 3, 4, 9, 10, 11, 12)]
            
            
            # Aggregate the table to combine interactions by sample
            interactionTable2 = aggregate(
              interactionTable$evidence,
              by = list(
                interactionTable$p1,
                interactionTable$p2,
                interactionTable$nt1 ,
                interactionTable$nt2,
                interactionTable$sample
              ),
              FUN = sum
            )
            colnames(interactionTable2) = c("p1", "p2", "nt1", "nt2", "sample", "evidence")
            interactionTable3 = aggregate(
              interactionTable2$sample,
              by = list(
                interactionTable2$p1,
                interactionTable2$p2,
                interactionTable2$nt1 ,
                interactionTable2$nt2
              ),
              FUN = paste,
              collapse = ","
            )
            interactionTable4 = aggregate(
              interactionTable2$evidence,
              by = list(
                interactionTable2$p1,
                interactionTable2$p2,
                interactionTable2$nt1 ,
                interactionTable2$nt2
              ),
              FUN = paste,
              collapse = ","
            )
            
            interactionTable5 = aggregate(
              interactionTable2$evidence,
              by = list(
                interactionTable2$p1,
                interactionTable2$p2,
                interactionTable2$nt1 ,
                interactionTable2$nt2
              ),
              FUN = sum
            )
            
            
            interactionTable3$evidence = interactionTable5$x
            interactionTable3$evidence2 = interactionTable4$x
            
            colnames(interactionTable3) = c("p1", "p2", "nt1", "nt2", "samples" , "evidence", "evidence2")
            
            
            
            
            
            #remove the need for evidence as we use it later
            interactionTable3_sub = interactionTable3[interactionTable3$evidence > evCutoff ,]
            
            
            
            #interactionTable3_sub = interactionTable3
            
            unique(interactionTable3_sub$samples)
            
            
            
            
            #head(interactionTable3_sub)
            #nrow(interactionTable3_sub)
            
            # find the probability of each constraint
            
            
            
            
            # now get the constraints and subset based on the start and end
            interactionTable3_sub = interactionTable3_sub[interactionTable3_sub$p1  > start  &
                                                            interactionTable3_sub$p1  < end &
                                                            interactionTable3_sub$p2  > start  &
                                                            interactionTable3_sub$p2  < end  ,]
            
            
            interactionTable3_sub$p1 = interactionTable3_sub$p1 - start + 1
            interactionTable3_sub$p2 = interactionTable3_sub$p2 - start + 1
            
            
            
            # this is how you sample the "bag"
            #sample(1:nrow(interactionTable3_sub),1,prob = normalized_evidence)
            
            # just gte the at staggered
            # #get the staggered Us
            #bpMatrix = outer(rnaRefs[[1]][[1]], rnaRefs[[1]][[1]], paste, sep=".")
            #diagi = c()
            #diagj = c()
            #for(i in 1:(ncol(bpMatrix)-1)){
            #    for(j in 1:(ncol(bpMatrix)-1)){
            #        if(bpMatrix[i,j] == "a.t"){
            #            print("paired Nucl detected")
            #            if(bpMatrix[i+1,j+1] == "t.a"){
            #                #print("Diaganal of pair detected")
            #                diagi = c(diagi,i)
            #                diagj = c(diagj,j)
            #            }
            #        }else if(bpMatrix[i,j] == "t.a"){
            #            print("paired Nucl detected")
            #            if(bpMatrix[i+1,j+1] == "a.t"){
            #                #print("Diaganal of pair detected")
            #                diagi = c(diagi,i)
            #                diagj = c(diagj,j)
            #            }
            #        }
            #    }
            #}
            
            x = interactionTable3_sub
            
            
            #index = c(paste(diagi,diagj),paste(diagj,diagi) )
            #this subsets ofr staggered U's
            #interactionTable3_sub = interactionTable3_sub[paste(interactionTable3_sub$p1,interactionTable3_sub$p2) %in% index,]
            # end get staggered
            
            
            
            normalized_evidence = interactionTable3_sub$evidence / sum(interactionTable3_sub$evidence)
            
            
            samples = sampleTable(cdsObject)$sampleName
            
            # run the structures
            
            prevConstraints = 0
            viennas = list()
            dgs = list()
            for (sample in samples) {
              interactionTable3_subSample = interactionTable3_sub[grepl(sample, interactionTable3_sub$sample),]
              normalized_evidenceSample = normalized_evidence[grepl(sample, interactionTable3_sub$sample)]
              
              dgSample = c()
              viennaSample = c()
              for (j in 1:ensembl) {
                goodvienna = ""
                prevConstraints = c()
                
                
                for (i in 1:constraintNumber) {
                  # pull constraints and re pick if constraint breaks the structure
                  constraints = c(
                    sample(
                      1:nrow(interactionTable3_subSample),
                      1,
                      prob = normalized_evidenceSample,
                      replace = FALSE
                    ),
                    prevConstraints
                  )
                  
                  
                  prevConstraints = constraints
                  #write contraint to file
                  
                  constraints = unique(constraints)
                  
                  constraintFile = interactionTable3_subSample[constraints, c("p1", "p2")]
                  #F i j k
                  constraintFile$F = "F"
                  constraintFile$k = 1
                  constraintFile = constraintFile[, c(3, 1, 2, 4)]
                  
                  
                  tmpFile = tempfile()
                  write.table(
                    constraintFile,
                    file = tmpFile,
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE
                  )
                  
                  
                  
                  
                  
                  table = data.frame(stringsAsFactors = FALSE)
                  
                  if (length(shape) == 1) {
                    command = paste(
                      "echo \">",
                      rna,
                      "\n",
                      paste(rnaRefs[[1]][[1]][start:end], collapse = ""),
                      "\" | RNAfold  --noPS --constraint=",tmpFile,
                      sep = ""
                    )
                    x = system(command, intern = TRUE)
                  } else{
                    shape2 = shape[start:end]
                    length = (end - start) + 1
                    shapeTable = data.frame("x" = 1:length,
                                            "y" = shape2)
                    shapeTable = shapeTable[complete.cases(shapeTable), ]
                    
                    tmpFile2 = tempfile()
                    write.table(
                      shapeTable,
                      tmpFile2,
                      quote = FALSE,
                      row.names = FALSE,
                      col.names = FALSE
                    )
                    command = paste(
                      "echo \">",
                      rna,
                      "\n",
                      paste(rnaRefs[[1]][[1]][start:end], collapse = ""),
                      "\" | RNAfold  --noPS --constraint=",tmpFile," --shape=",tmpFile2,
                      sep = ""
                    )
                    x = system(command, intern = TRUE)
                  }
                  if (!(grepl("\\(\\(", x[3]))) {
                    prevConstraints = prevConstraints[-1]
                    next
                    
                  } else{
                    goodvienna = x[3]
                  }
                  
                  
                  
                  
                }
                viennaSample = c(viennaSample, sub("\\s.*", "", goodvienna))
                
                dg = sub(".*\\s", "", goodvienna)
                dg2 = sub("\\(", "", dg)
                dgSample = c(dgSample, sub("\\)", "", dg2))
                #print("NEW STRUCTURE")
                
                
              }
              viennas[[sample]] = viennaSample
              dgs[[sample]] = dgSample
            }
            
            
            
            
            
            
            
            viennas2 = list()
            
            viennas2[[paste(start, ":", end, sep = "")]] = viennas
            
            viennas = viennas2
            
            
            
            
            
            ###########################################################
            # Make object
            ###########################################################
            message(" *** Creating object ***")
            #create rnaCrosslink dataset object
            object  = new(
              "rnaCrosslinkDataSet",
              rnas = rnas(cdsObject),
              rnaSize = rnaSize(cdsObject),
              sampleTable = sampleTable(cdsObject),
              InputFiles = InputFiles(cdsObject),
              matrixList = matrixList(cdsObject),
              clusterTableList = clusterTableList(cdsObject),
              clusterGrangesList = clusterGrangesList(cdsObject),
              clusterTableFolded = clusterPositionsListTrimmedSarsCombinedWithStructures,
              interactionTable = tableAll,
              viennaStructures = viennas,
              dgs = dgs
            )
            
            return(object)
            
          })
