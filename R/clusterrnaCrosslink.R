#' @include rnaCrosslinkOO.R rnaCrosslinkDataSet.R
NULL




#' clusterrnaCrosslink
#'
#'  This method clusters the duplexes. 
#'
#'
#' @param cds rnaCrosslinkDataSet object created with rnaCrosslinkDataSet
#' @param cores numeric - The number of cores to use 
#' @param stepCount Stringency for clustering 
#' @param clusterCutoff The minimum number of reads a cluster requires 
#'
#'
#' @name clusterrnaCrosslink
#' @aliases clusterrnaCrosslink,rnaCrosslinkDataSet-method
#' @docType methods
#' @rdname clusterrnaCrosslink
#' @return A \code{rnaCrosslinkDataSet} object
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusterrnaCrosslink(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 0)
#' @export

setGeneric("clusterrnaCrosslink",
           function(cds,
                    cores = 3,
                    stepCount = 2,
                    clusterCutoff = 20) standardGeneric("clusterrnaCrosslink" ) )

setMethod("clusterrnaCrosslink",
          "rnaCrosslinkDataSet",
          function(cds,
                   cores = 3,
                   stepCount = 2,
                   clusterCutoff = 20)  {
            message("********************************************")
            # Set up variables
            clusters = list()
            clusterTables = list()
            for(rna in rnas(cds)){
              
              # Set up variables
              clusters[[rna]] = list() #will contain the clusterGranges
              clusterTables[[rna]] = list() #will contain the cluster tables
              
              rnaSize = rnaSize(cds)
              message(paste("**** ", rna, " *****"))
              message(paste("****              ", rnaSize," nt ", "           ****"))
              
   
              ##############################
              message(paste("****       Assessing Long Range         ****"))
              longDistInput = subsetInputList2(InputFiles(cds)[[rna]][[ "noHost" ]],
                                           11,rnaSize,length = 800)
              chimeraList = InputToGRanges(longDistInput,rna)
              message(paste("****        Sampling Long Range         ****"))
              chimeraListSampled = sampleChimeras(chimeraList)
              
              matrixList = matrixList(cds)
              
              cl <- makeCluster(cores)
              registerDoParallel(cl)
              

              plottingList = list()
              for(i in 1:length(sampleNames(cds))){
                plottingList[[i]] = list()
                
                if(any(is.na(chimeraListSampled[[i]][["left"]]))){
                  plottingList[[i]][[k]] = GRanges()
                }else{
                  
                  foreach (k=1:length( chimeraListSampled[[i]][["gap"]])) %do% {
                    adjacancyMat = getAdjacancyMat(chimeraListSampled[[i]][["gap"]][[k]],
                                                   "nucleotide",
                                                   15)
                    net = graph_from_adjacency_matrix(adjacancyMat,
                                                      mode = "undirected",
                                                      weighted = T)
                    clustering = cluster_walktrap(net,steps = stepCount)
                    #remove low evidence clusters
                    highest_clusters = 
                      names(table(membership(clustering))[table(membership(clustering)) >
                                                            clusterCutoff])
                    plottingList[[i]][[k]]  = 
                      printClustersFast(dir,clustering, 
                                        highest_clusters,
                                        chimeraListSampled[[i]][["left"]][[k]],
                                        chimeraListSampled[[i]][["right"]][[k]])
                  }
                }
                
              }
              
              longRange = plottingList
              
              ##############################
              # short Long interactions
              message(paste("****      Assessing Short Range         ****"))
              longDistInput = subsetInputList2(InputFiles(cds)[[rna]][[ "noHost" ]],
                                           1,10,length = 800)
              chimeraList = InputToGRanges(longDistInput,rna)
              
              message(paste("****      Sampling short Range          ****"))
              chimeraListSampled = sampleChimeras(chimeraList)
              
              plottingList = list()
              for(i in 1:length(sampleNames(cds))){
                plottingList[[i]] = list()
                
                if(any(is.na(chimeraListSampled[[i]][["left"]]))){
                  plottingList[[i]][[k]] = GRanges()
                }else{
                  
                  foreach (k=1:length( chimeraListSampled[[i]][["gap"]])) %do% {
                    adjacancyMat = getAdjacancyMat(chimeraListSampled[[i]][["gap"]][[k]],"nucleotide", 15)
                    net = graph_from_adjacency_matrix(adjacancyMat, mode = "undirected", weighted = T)
                    clustering = cluster_walktrap(net,steps = 2)
                    highest_clusters = names(table(membership(clustering))[table(membership(clustering)) > clusterCutoff])
                    plottingList[[i]][[k]]  = printClustersFast(dir,clustering, highest_clusters,
                                                                chimeraListSampled[[i]][["left"]][[k]],
                                                                chimeraListSampled[[i]][["right"]][[k]])
                    message(paste("*****        done ",sampleNames(cds)[i] , k, " /" , length( chimeraListSampled[[i]][["gap"]]), "          *****"))
                    
                    
                  }
                }
                
              }
              
              
              shortRange = plottingList
              
              
              
              ##############################
              #make the IDs for the table and combine
              combinedPlotting = shortRange
              for(i in 1:length(shortRange)){
                plotting = GRangesList(shortRange[[i]])
                for(j in 1:length(plotting)){
                  if(length(plotting[[j]]) == 0){next}
                  plotting[[j]]$k = paste(plotting[[j]]$cluster,"binShort",j, sep = ":")
                }
                combinedPlotting[[i]] =  plotting
              }
              
              longRange
              for(i in 1:length(longRange)){
                plotting = GRangesList(longRange[[i]])
                for(j in 1:length(plotting)){
                  if(length(plotting[[j]]) == 0){next}
                  plotting[[j]]$k = paste(plotting[[j]]$cluster,"binLong",j, sep = ":")
                
                combinedPlotting[[i]] =  c(combinedPlotting[[i]],plotting)
                }
              }
              
              
              clusters[[rna]][["original"]] = combinedPlotting
              
              
              
              
              ##############################
              # Make a matrix and table for the clusters
              clusterPositionsList = list()
              matList = list()
              for(j in 1:length(sampleNames(cds))){
                
                
                plotting =  unlist(combinedPlotting[[j]])
                
                lengths = aggregate(mcols(plotting)$cluster, by = list(mcols(plotting)$k), FUN = length)
                row.names(lengths) = lengths$Group.1
                #for each cluster get the min start and max end
                plottingSplit = split(plotting, paste(mcols(plotting)$k, mcols(plotting)$type))
                minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
                maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
                clusterPositionsList[[j]] = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                                                       "ls" = minStarts[seq(1,length(minStarts),2)],
                                                       "le" = maxEnd[seq(1,length(maxEnd),2)],
                                                       "rs" = minStarts[seq(2,length(minStarts),2)],
                                                       "re" = maxEnd[seq(2,length(maxEnd),2)],
                                                       "size" = lengths[sub("\\s.*",
                                                                            "",
                                                                            names(maxEnd)[seq(1,length(minStarts),2)]),])
                
                
                
                
                matList[[j]] = matrix(0,nrow = rnaSize, ncol = rnaSize)
                
                clusterPositionsList[[j]] = clusterPositionsList[[j]][!(is.na(clusterPositionsList[[j]]$size.x)),]
                for(i in 1:nrow(clusterPositionsList[[j]])){
                  matList[[j]][clusterPositionsList[[j]][i,"ls"]:clusterPositionsList[[j]][i,"le"],
                               clusterPositionsList[[j]][i,"rs"]:clusterPositionsList[[j]][i,"re"]] =
                    matList[[j]][clusterPositionsList[[j]][i,"ls"]:clusterPositionsList[[j]][i,"le"],
                                 clusterPositionsList[[j]][i,"rs"]:clusterPositionsList[[j]][i,"re"]] +
                    clusterPositionsList[[j]][i, "size.x"]
                  
                }
              }
              names(matList) = sampleNames(cds)
              names(clusterPositionsList) = sampleNames(cds)
              
              clusterTables[[rna]][["original"]] =clusterPositionsList
              matrixList[[rna]][["originalClusters"]] = matList
              clusterGranges = clusters
            }
            
            
            
            ###########################################################
            # Make object
            ###########################################################
            message("*****          Creating object         *****")
            message("******************************************** ")
            #create rnaCrosslink dataset object
            object  = new("rnaCrosslinkDataSet",
                          rnas = rnas(cds),
                          rnaSize = rnaSize(cds),
                          sampleTable = sampleTable(cds),
                          InputFiles = InputFiles(cds),
                          matrixList = matrixList,
                          clusterTableList = clusterTables,
                          clusterGrangesList = clusterGranges,
                          clusterTableFolded = data.frame(),
                          interactionTable = data.frame(),
                          viennaStructures = list(),
                          dgs = list()
            )
            stopCluster(cl)
            return(object)
            
            
          })
