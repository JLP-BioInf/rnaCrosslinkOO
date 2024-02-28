#' @include  rnaCrosslinkDataSet.R 
NULL

#' trimClusters
#'
#' Trimming of the clusters removes redundant information derived from random
#' fragmentation of the reads during library preparation. This method takes
#' a \code{rnaCrosslinkDataSet} object where clustering has been performed with 
#' the clusterrnaCrosslink method and trims the clusters according to the 
#' trimFactor argument.
#' 
#' The 3 attributes; matrixList, clusterTableList and clusterGrangesList 
#' will gain the \code{types} "superClusters" and "trimmedClusters"
#'
#' @param clusteredCds a \code{rnaCrosslinkDataSet} object
#' @param trimFactor a positive value that defines how much the clusters will 
#' @param clusterCutoff Minimum number of reads before discarding cluster
#' be trimmed = mean + ( sd * trimFactor )
#' 
#' @return Returns a \code{rnaCrosslinkDataSet} object
#' @name trimClusters
#' @docType methods
#' @rdname trimClusters
#' @aliases trimClusters,rnaCrosslinkDataSet-method
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusteredCds = clusterrnaCrosslink(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 0)
#'                 
#' trimClusters(clusteredCds = clusteredCds,
#'              trimFactor = 1, 
#'              clusterCutoff = 0)
#' @export
#' 
setGeneric("trimClusters",
           function(clusteredCds,
                    trimFactor = 2.5, 
                    clusterCutoff = 1) 
             standardGeneric("trimClusters" ) )

setMethod("trimClusters",
          "rnaCrosslinkDataSet",
          function(clusteredCds,
                   trimFactor = 2.5, 
                   clusterCutoff = 1)  {
              message("********************************************")
              ##############################
              # set up variables
              allChimerasForSuperClustersPlotting = list()
              goodClusters = list()
              for(rna in rnas(clusteredCds)){  ## for each RNA
                  
                  # size of rna
                  rnaSize = rnaSize(clusteredCds)
                  # original clusters
                  originalClusters =  clusterTableList(clusteredCds)[[rna]][["original"]]
                  
                  
                  ##############################
                  # Now cluster the clusters
                  #get original tables
                  clusterPositionsList = clusterTableList(clusteredCds)[[rna]][["original"]]
                  #get original gRanges
                  combinedPlotting     = clusterGrangesList(clusteredCds)[[rna]][["original"]]
                  
                  ##############################
                  # Set up new tables, matric and granges lists
                  superclustersPoisitonList = list()
                  superclustersPlotting = list()
                  matList = list()
                  
                  for(z in 1:length(sampleNames(clusteredCds))){
         
                      clusterPositions = clusterPositionsList[[z]]
                      ##############################
                      # changes coordinates of clusters where the 2 sides overlap
                      clusterPositions2 = clusterPositions
                      for(i in 1:nrow(clusterPositions )){
                          if(clusterPositions$le[i] > clusterPositions$rs[i]){
                              clusterPositions2[i,"rs"] =     clusterPositions[i,"le"] +1
                          }
                      }
                      ##############################
                      #make Granges left right and gap
                      left = GRanges(seqnames=rna,
                                     IRanges(start=clusterPositions2$ls,
                                             end=clusterPositions2$le))
                      names(left) <- clusterPositions2$id
                      right= GRanges(seqnames=rna,
                                     IRanges(start=clusterPositions2$rs,
                                             end=clusterPositions2$re))
                      names(right) <- clusterPositions2$id
                      distances = GRanges(seqnames=rna,
                                          IRanges(start=clusterPositions2$le,
                                                  end=clusterPositions2$rs))
                      names(distances) <- clusterPositions2$id
                      
                      
                      ##############################
                      # Now make super clusters
                      ##############################
                      
                      # from the gaps, make a adjacancy matrix
                      adjacancyMat = getAdjacancyMat(distances,"nucleotide", 35)
                      # create Graph
                      cluster = c()
                      if(!any(is.na(adjacancyMat))){
                          net = graph_from_adjacency_matrix(adjacancyMat,
                                                            mode = "undirected",
                                                            weighted = T)
                          # clusterGraph
                          clustering = cluster_walktrap(net,steps = 1)
                          
                          ##############################
                          # Store Super-clusters
                          highest_clusters = names(table(membership(clustering)))
                          # printClustersFast function creates a the standard clustering
                          # table from the iGraph output
                          superclustersPlotting[[z]]  = printClustersFast(tempfile(),clustering, highest_clusters, left, right)
                          plottingListFull = superclustersPlotting[[z]]
                          
                          ##############################
                          # Identify orphan clusters that missed with superclustering
                          missing = as.character(clusterPositions$id[which( !(as.character(clusterPositions$id) %in% unique(names(plottingListFull)) ) )])
                          clusterPositionsmissing = clusterPositions[clusterPositions$id %in% missing,]
                          
                          ##############################
                          #  Get super cluster and cluster identity
                          cluster = mcols(plottingListFull)$cluster
                          names(cluster)= names(plottingListFull)
                          #  add this super cluster membership to the clustering table
                          clusterPositions$superCluster = cluster[as.character((clusterPositions$id))]
                          clusterPositions = clusterPositions[!is.na(clusterPositions$superCluster),]
                          
                          
                          # no find the number of chimeras in each supercluster
                          lengths = aggregate(clusterPositions$size.x, by = list(clusterPositions$superCluster), FUN = sum)
                          row.names(lengths) = lengths$Group.1
                          
                          
                          # subset by the cluster membersip
                          #goodClusters = lengths[lengths$x > clusterCutoff, "Group.1"]
                          #clusterPositions = clusterPositions[clusterPositions$superCluster %in% goodClusters,]
                          
                          
                          ##############################
                          # make Table
                          #for each cluster get the min start and max end
                          plottingSplit = split(plottingListFull, paste(mcols(plottingListFull)$cluster, mcols(plottingListFull)$type))
                          
                          #returns the min start and max ends of each cluster
                          minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
                          maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
                          # Make the clustering table
                          clusterPositionsCombined = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                                                                "ls" = minStarts[seq(1,length(minStarts),2)],
                                                                "le" = maxEnd[seq(1,length(maxEnd),2)],
                                                                "rs" = minStarts[seq(2,length(minStarts),2)],
                                                                "re" = maxEnd[seq(2,length(maxEnd),2)],
                                                                "size" = lengths[as.numeric(sub("\\s.*","",names(maxEnd)[seq(1,length(minStarts),2)])),])
                          
                          # Make the clustering table
                          superclustersPoisitonList[[z]] = rbind.data.frame(clusterPositionsmissing,
                                                                            clusterPositionsCombined,
                                                                            stringsAsFactors = FALSE)
                          
                          goodClusters[[z]] = superclustersPoisitonList[[z]][superclustersPoisitonList[[z]]$size.x > clusterCutoff,]                          
                          goodClusters[[z]] = sub("\\s.*","",row.names(superclustersPoisitonList[[z]]))
   
                          ##############################
                          # Make matrices of superclusters
                          clusterPositions = superclustersPoisitonList[[z]]
                          
                      }else {
                          clusterPositions = clusterPositions
                          superclustersPoisitonList[[z]] = clusterPositions
                      }
                      
                      
                      
                      mat = matrix(0,nrow = rnaSize, ncol = rnaSize)
                      for(i in 1:nrow(clusterPositions)){
                          mat[clusterPositions[i,"ls"]:clusterPositions[i,"le"],
                              clusterPositions[i,"rs"]:clusterPositions[i,"re"]] =    mat[clusterPositions[i,"ls"]:clusterPositions[i,"le"],
                                                                                          clusterPositions[i,"rs"]:clusterPositions[i,"re"]] + clusterPositions[i, "size.x"]
                          
                      }
                      matList[[z]] = mat
                  }
                  ##############################
                  # Save new Table and matrix list -  super clusters
                  message(  "******        Trimming Clusters       ******")
                  message(  "******            Saving              ******")
                  message(  "******        Saving mat list         ******")
                  ml = matrixList(clusteredCds)
                  ml[[rna]][["superClusters"]] = matList
                  
                  message("******       Saving table list        ******")
                  ctl = clusterTableList(clusteredCds)
                  ctl[[rna]][["superClusters"]] =   superclustersPoisitonList
              }
              
              
              ##############################
              # Get Granges List for the super clusters (containing original duplexes)
              allChimerasForSuperClustersPlotting = list()
              combinedPlottingSplit = list()
              combinedPlottingUnlist = list()
              
              for(b in 1:length(sampleNames(clusteredCds))){
                  
                  combinedPlottingUnlist = unlist(combinedPlotting[[b]])
                  combinedPlottingSplit = split(combinedPlottingUnlist,
                                                mcols(combinedPlottingUnlist)$cluster)
                  
                  
                  
                  superClusterArray =  sub("\\s.*","",names(superclustersPlotting[[b]][duplicated(names(superclustersPlotting[[b]]))]))
                  names(superClusterArray) = superclustersPlotting[[b]][duplicated(names(superclustersPlotting[[b]]))]$cluster
                  
                  x = superclustersPoisitonList[[b]][ grep("bin", row.names(superclustersPoisitonList[[b]])),]
                  names = c(names(superClusterArray), row.names(x))
                  
                  superClusterArray = c(superClusterArray,row.names(x))
                  
                  names(superClusterArray) = names
                  combinedPlottingUnlist$superCluster = "X"
                  combinedPlottingUnlist = combinedPlottingUnlist[which(!(is.na(combinedPlottingUnlist$k ))),]
                  for( z in 1:length(superClusterArray)){
                      supercluster = names(superClusterArray)[z]
                      cluster = unique(sub("\\s.*","",superClusterArray[z]))
                      combinedPlottingUnlist[combinedPlottingUnlist$k == cluster,]$superCluster = supercluster
                  }
                  allChimerasForSuperClustersPlotting[[b]] = combinedPlottingUnlist
              }
              
              
              
              
              
              
              ##############################
              # Save new Granges super clusters
   
              cgr = clusterGrangesList(clusteredCds)
              cgr[[rna]][["superClusters"]]   =  allChimerasForSuperClustersPlotting
              
              
              
              
              
              
              ##############################
              # Now Trim the clusters
              # the new granges list
              allChimerasForSuperClustersPlottingTrimmed = list()
              # for each sample
              for(i in 1:length(sampleNames(clusteredCds))){
                  allChimerasForSuperClustersPlottingTrimmed[[i]] = GRanges()
                  #for each cluster
                  for(cluster in unique(allChimerasForSuperClustersPlotting[[i]]$superCluster)[unique(allChimerasForSuperClustersPlotting[[i]]$superCluster) %in% 
                                                                                               c(goodClusters[[i]], paste(goodClusters[[i]],"left"))] ){

                    cluster2 = sub("\\s.*","" ,cluster)
                      lefty = list()
                      
                      #for the left and right sides for each cluster
                      # cut the ends based mean and sd of evidence
                      for(l in c("left","right")){
                          clusterrange = allChimerasForSuperClustersPlotting[[i]][allChimerasForSuperClustersPlotting[[i]]$superCluster == cluster & allChimerasForSuperClustersPlotting[[i]]$type == l  ,]
                          min = min(start(clusterrange[clusterrange$superCluster == cluster,]))
                          max = max(end(clusterrange[clusterrange$superCluster == cluster,]))
                          s = (start(clusterrange[clusterrange$superCluster == cluster &clusterrange$type == l  ,]))
                          e = (end(clusterrange[clusterrange$superCluster == cluster &clusterrange$type == l  ,]))
                          
                          # function that vectorises Seq
                          seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
                          # get a vector of each cluster and side from start to end
                          # the more eveidence the more times a number appears
                          x = unlist(c(seq2(from = (s), to = e)))
                          xt = table(x)
                          x1 = x
                          # find rhe mean + and - one sd and
                          # make a GRanges with this value
                          removal =  mean(x) + sd(x)*trimFactor
                          included  = GRanges(seqnames=rna,
                                              IRanges(
                                                  start=rep(min, length(clusterrange)),
                                                  end=rep(removal, length(clusterrange))
                                              ))
                          if( l == "left"){
                              removal =  mean(x) -  sd(x)*trimFactor
                              included = GRanges(seqnames=rna,
                                                 IRanges(
                                                     start=rep(removal, length(clusterrange)),
                                                     end=rep(max, length(clusterrange))
                                                 ))
                          }
                          
                          # now remove any bases that overlap with that
                          t = pintersect( clusterrange,included)
                          
                          allChimerasForSuperClustersPlottingTrimmed[[i]] = c(  allChimerasForSuperClustersPlottingTrimmed[[i]],t)
                          # un comment to print the views of the trimming
                          # s = (start(t))
                          #    e = (end(t))
                          #    x = unlist(c(seq2(from = s, to = e)))
                          #    if( l == "left"){
                          #        lefty[[1]] = x
                          #            lefty[[2]] = x1
                          #        }
                          # }
                          
                          #print(cluster )
                          #tbl1 = data.frame(table(c(x1,lefty[[2]])))
                          #tbl2 = data.frame(table(c(x,lefty[[1]])))
                          #plot(ggplot(mapping =  aes(x = Var1, y = as.numeric(as.character(Freq))))+
                          #       geom_bar(data = tbl1, stat = "identity")+
                          #       geom_bar(data = tbl2, stat = "identity", colour = "firebrick") +
                          #       theme_classic())
                      }
                  }
              }
              
              
              
              
              
              ##############################
              # From the trimmed Granges make a cluster table
              # foir the trimmed super clusters
              matListTrimmed = list()
              clusterPositionsListTrimmed = list()
              for(j in 1:length(sampleNames(clusteredCds))){
                  
                  plotting  = allChimerasForSuperClustersPlottingTrimmed[[j]]
                  lengths = aggregate(mcols(plotting)$superCluster, by = list(mcols(plotting)$superCluster), FUN = length)
                  row.names(lengths) = lengths$Group.1
                  plottingSplit = split(plotting, paste(mcols(plotting)$superCluster, mcols(plotting)$type))
                  
                  minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
                  maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
                  x = sub("\\sleft","",names(maxEnd)[seq(1,length(minStarts),2)])
                  x = sub("\\sright","",x,2)
                  clusterPositionsListTrimmed[[j]] = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                                                                "ls" = minStarts[seq(1,length(minStarts),2)],
                                                                "le" = maxEnd[seq(1,length(maxEnd),2)],
                                                                "rs" = minStarts[seq(2,length(minStarts),2)],
                                                                "re" = maxEnd[seq(2,length(maxEnd),2)],
                                                                "size" = lengths[x,])
                  
                  ###################################
                  # make the matrices
                  matListTrimmed[[j]] = InputFiles(clusteredCds)[[rna]][["noHost"]][[j]][InputFiles(clusteredCds)[[rna]][["noHost"]][[j]]$V1 %in% names( allChimerasForSuperClustersPlottingTrimmed[[j]]),]

              }
              
              
              ###################################
              # And save
              message("******           Saving  End          ******")
              message("******      Saving mat list  End      ******")

              ml[[rna]][["trimmedClusters"]] =     getMatrices(matListTrimmed,
                                                               rna, rnaSize)
              names(ml[[rna]][["trimmedClusters"]]) = sampleNames(clusteredCds)
              
              message("******     Saving granges list        ******")
              cgr[[rna]][["trimmedClusters"]]   =  allChimerasForSuperClustersPlottingTrimmed
              
              message("******     Saving table list  End     ******")
              ctl[[rna]][["trimmedClusters"]] =   clusterPositionsListTrimmed
              message("********************************************")
              
              ###################################
              # Re-make the object
              object  = new("rnaCrosslinkDataSet",
                            rnas = rnas(clusteredCds),
                            rnaSize = rnaSize(clusteredCds),
                            sampleTable = sampleTable(clusteredCds),
                            InputFiles = InputFiles(clusteredCds),
                            matrixList = ml,
                            clusterTableList = ctl,
                            clusterGrangesList = cgr
              )
              return(object)
          })








#' InputToGRanges
#'
#' This function is useful to turn a list of Input data into lists of GRanges
#' It creates a list for each sample one for the left side one for the right
#' side and one for the gap in the middle.
#' 
#' @param InputList the original InputList created with readInputFiles or subsetInputList
#' @param rna The rna of interest
#' 
#' @return A list of GRanges data in Input format
#' @name InputToGRanges
#' @aliases InputToGRanges
#' @docType methods
#' @rdname InputToGRanges
#'  
InputToGRanges = function(InputList, 
                        rna){
    seqName = rna
    InputOutput2 = InputList
    gList = list()
    for(i in 1:length(InputOutput2)){
      if(!is.data.frame(InputOutput2[[i]])){
        gList[[i]] = list()
        gList[[i]][["left"]] = NA
        gList[[i]][["right"]] = NA
        gList[[i]][["gap"]] = NA
      }else{
        InputOutput = InputOutput2[[i]]
        gList[[i]] = GRangesList()
        #make a Granges from the left
        left <- GRanges(seqnames=seqName,
                        IRanges(
                            start=InputOutput$V7,
                            end=InputOutput$V8
                        ))
        names(left) <- InputOutput$V1
        
        
        
        
        #make a GRanges from the right
        right <-GRanges(seqnames=seqName,
                        IRanges(
                            start=InputOutput$V13,
                            end=InputOutput$V14
                        ))
        names(right) <- InputOutput$V1
        
        
        
        distances = GRanges(seqnames=seqName,
                            IRanges(
                                start=end(left),
                                end=start(right)
                            ))
        names(distances) <- InputOutput$V1
        
        gList[[i]][["left"]] = left
        gList[[i]][["right"]] = right
        gList[[i]][["gap"]] = distances
      }
    }
    return(gList)
}




#' sampleChimeras
#'
#' This function samples chimeras into smaller chunks so that clustering is 
#' quicker 
#' 
#' @param chimeraList list of chimeras
#' @name sampleChimeras
#' @aliases sampleChimeras
#' @rdname sampleChimeras
sampleChimeras = function(chimeraList){
    
    chimeraListSampled =list()
    
    for(i in 1:length(chimeraList)){
        max =  length(chimeraList[[i]][["left"]])
        seq = c(1,max)
        
        if(max > 10000){
            seq = seq(1,max,by = 3000)
            seq = c(seq,max)
            
        }
        
        chimeraListSampled[[i]] = list()
        
        for(j in c("left","right","gap")){
            chimeraListSampled[[i]][[j]] = list()
            for(k in 1:(length(seq)-1)){
                sample = seq[k]:seq[k+1]
                chimeraListSampled[[i]][[j]][[k]] = chimeraList[[i]][[j]][sample]
            }
        }
        
    }
    return(chimeraListSampled)
}







#' compareKnown
#'
#' This method compares the current object to a know structure.run 
#' \code{trimClusters()} on the  \code{rnaCrosslinkDataSet} first
#'
#' @param trimmedClusters a \code{rnaCrosslinkDataSet} object, 
#' run \code{trimClusters()} on the  \code{rnaCrosslinkDataSet} first
#' 
#' @param knownMat Matrix - A marix(ncol = lengthRNA,nrow = lengthRNA) where a
#' value in matrix[x,y] would indicate a known interation between nucleotide 
#' x and nucleotide y 
#' @param type string - the Analysis stage of clusters you would like to compare you can find 
#' available types by just running the objects name
#' 
#' 
#' @return Returns a \code{rnaCrosslinkClusteredDataSet} object
#' 
#' The 3 attributes matrixList, clusterTableList and clusterGrangesList 
#' will gain the \code{types} "known" and "novel" and "knownAndNovel"
#' @name compareKnown
#' @docType methods
#' @rdname compareKnown
#' @aliases compareKnown,rnaCrosslinkDataSet-method
#' @examples 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' clusteredCds = clusterrnaCrosslink(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 0)
#' knownMat = matrix(0, ncol = rnaSize(cds), nrow = rnaSize(cds))
#' knownMat[7,27] = 1
#' # use compare known to gett he known and not know clusters
#' knowClusteredCds = compareKnown(clusteredCds,
#'                                 knownMat,
#'                                 "original")
#' clusterNumbers(knowClusteredCds)                 
#'                 
#'
#' @export
setGeneric("compareKnown", 
           function(trimmedClusters, 
                    knownMat,
                    type) standardGeneric("compareKnown"))

setMethod("compareKnown", "rnaCrosslinkDataSet", function(trimmedClusters, 
                                                      knownMat,
                                                      type)  {
    
    ###################################
    # Inputs
    rna = rnas(trimmedClusters)
    k18Smat = knownMat
    type = type
    trimmedClusters = trimmedClusters
    sampleNames = sampleNames(trimmedClusters)
    ml = matrixList(trimmedClusters)
    rnaSize = rnaSize(trimmedClusters)
    
    ###################################
    # set up variables
    ml[[rna]][["KnownAndNovel"]] = list()
    novelClusters = list()
    novelClustersMat = list()
    novelClustersMat2 = list()
    cannonicalClusters = list()
    cannonicalClustersMat = list()
    cannonicalClustersMat2 = list()
    
    
    clusterPositionsListTrimmed = clusterTableList(trimmedClusters)[[rna]][[type]]
    # for each sample
    for(i in 1:length(clusterPositionsListTrimmed)){
        
        # set up the matrix for this sample and the cluster positions
        clusters = clusterPositionsListTrimmed[[i]]
        
        ###################################
        # test each cluster against the known interactions
        #this a matrix of known positions:
        k18Smat2 = k18Smat
        tf = c()
        for(j in 1:nrow(clusters)){
            #For each cluster, make a individual matrix
            clusterMat = matrix(0, nrow = rnaSize, ncol = rnaSize)
            # add 10 to the positions of this cluster
            clusterMat[ clusters$ls[j]:clusters$le[j] ,  clusters$rs[j]:clusters$re[j] ] =
                clusterMat[  clusters$ls[j]:clusters$le[j] ,  clusters$rs[j]:clusters$re[j] ] + 10
            
            # Add that to the matric of cannonical interactions
            k18Smat2 = k18Smat + clusterMat
            # find those when the two annotations overlap.
            tf = c(tf, all(k18Smat2 < 11))
        }
        
        ###################################
        # Get novel and cannonical tables and matrices
        #print(which((tf == F)))
        novelClusters[[i]] = clusterPositionsListTrimmed[[i]][which((tf == TRUE)),]
        novelClustersMat[[i]] = matrix(0, nrow = rnaSize, ncol = rnaSize)
        if(nrow(novelClusters[[i]]) > 0){
            for(j in 1:nrow(novelClusters[[i]])){
                
                
                
                
                novelClustersMat[[i]][novelClusters[[i]]$ls[j]:novelClusters[[i]]$le[j] ,
                                      novelClusters[[i]]$rs[j]:novelClusters[[i]]$re[j] ] =
                    novelClustersMat[[i]][ novelClusters[[i]]$ls[j]:novelClusters[[i]]$le[j] ,
                                           novelClusters[[i]]$rs[j]:novelClusters[[i]]$re[j] ] +  novelClusters[[i]]$size.x[j]
                
            }
        }
        
        #print(which((tf == T)))
        cannonicalClusters[[i]] = clusterPositionsListTrimmed[[i]][which((tf == FALSE)),]
        cannonicalClustersMat[[i]] = matrix(0, nrow = rnaSize, ncol = rnaSize)
        
        
        
        for(j in 1:nrow(cannonicalClusters[[i]])){
            
            cannonicalClustersMat[[i]][cannonicalClusters[[i]]$ls[j]:cannonicalClusters[[i]]$le[j] ,
                                       cannonicalClusters[[i]]$rs[j]:cannonicalClusters[[i]]$re[j] ] =
                cannonicalClustersMat[[i]][ cannonicalClusters[[i]]$ls[j]:cannonicalClusters[[i]]$le[j] ,
                                            cannonicalClusters[[i]]$rs[j]:cannonicalClusters[[i]]$re[j] ] + cannonicalClusters[[i]]$size.x[j]
        }
        
        ###################################
        # Add the known interactions to the known and novel matrices
        cannonicalClustersMat2[[i]] = cannonicalClustersMat[[i]] + knownMat*30000
        novelClustersMat2[[i]] = novelClustersMat[[i]] + knownMat*30000
        ml[[rna]][["KnownAndNovel"]][[i]] = ml[[rna]][[type]][[i]] + knownMat*30000
    }
    
    
    
    ###################################
    # add to the lists for the object
    ml[[rna]][["novel"]] = novelClustersMat2
    ml[[rna]][["known"]] = cannonicalClustersMat2
    ctl = clusterTableList(trimmedClusters)
    
    ctl[[rna]][["novel"]] = novelClusters
    ctl[[rna]][["known"]] = cannonicalClusters
    cgl = clusterGrangesList(trimmedClusters)
    
    
    
    ###################################
    # create object
    object  = new("rnaCrosslinkDataSet",
                  rnas = rnas(clusteredCds),
                  rnaSize = rnaSize(clusteredCds),
                  sampleTable = sampleTable(clusteredCds),
                  InputFiles = InputFiles(clusteredCds),
                  matrixList = ml,
                  clusterTableList = ctl,
                  clusterGrangesList = cgl,   
                  clusterTableFolded = data.frame(),
                  interactionTable = data.frame(),
                  viennaStructures = list(),
                  dgs = list()
    
                  
    )
    
    return(object)
    
})






#' printClustersFast
#' 
#' Makes a table with the coordinates of the clusters
#'
#' Does the same as printClusters but is a lot faster and does not create plots
#' of each cluster
#' @param  dir the directory that contains the *Inputrids.Input files
#' @param  clustering The output from the iGraph function cluster_walktrap for the (made with adjacency matrix input)
#' @param  highest_clusters The cluster you are interested in keeping
#' @param  left list created with InputToGRanges (but just the left section of the list)
#' @param  right list created with InputToGRanges (but just the right section of the list)
#' @return A table of clusters and coordinates
#' @name printClustersFast
#' @docType methods
#' @rdname printClustersFast
printClustersFast = function(dir, 
                             clustering, 
                             highest_clusters, 
                             left, 
                             right ){
    
    plotting = GRanges()
    for(i in highest_clusters ){
        c1c1 = names(membership(clustering)[membership(clustering) == i])
        plotting2 = GRanges()
        plotting2 = addCluster(left, c1c1, plotting2,i  ,"left")
        plotting2 = addCluster(right, c1c1, plotting2,i ,"right" )
        plotting = addCluster(left, c1c1, plotting,i  ,"left")
        plotting = addCluster(right, c1c1, plotting,i ,"right" )
    }
    return(plotting)
}



#helper function for printClusters
addCluster = function(granges, indexes, prev, cluster, type){
    x = granges[as.numeric(indexes)]
    x$cluster = cluster
    x$type = type
    prev= c(prev, x)
}








#' subsetInputList2
#' 
#' Subset a list of Input files
#'
#' Function used to subset a list of Input data created by readInputFiles
#' This function produces the same size list as before but 
#' it returns ONLY the rna of interest and also
#' Choose duplexes where the nt difference in position between the
#' one side and other side of an interaction is between min and max
#' @param InputList the original InputList created with readInputFiles
#' @param min the rna of interest that you want to subset
#' @param max The number of randomly subsetted chimeric reads you need
#' @param length The number of randomly subsetted chimeric reads you need
#' @return A list of subsetted Input files
#' @name subsetInputList2
#' @docType methods
#' @rdname subsetInputList2
subsetInputList2 = function(InputList, 
                          min, 
                          max, 
                          length){
    longDistInput = list()
    for (i in 1:length(InputList)){
        InputList[[i]]$dist = InputList[[i]]$V13 -  InputList[[i]]$V8
        longDistInput[[i]] = InputList[[i]][InputList[[i]]$dist < max & InputList[[i]]$dist >= min,]
        leftLength = longDistInput[[i]]$V6 - longDistInput[[i]]$V5
        rightLength = longDistInput[[i]]$V12 - longDistInput[[i]]$V11
        longDistInput[[i]] = longDistInput[[i]][leftLength < length & rightLength < length,]
    }
    for (i in 1:length(InputList)){
      if(nrow(longDistInput[[i]]) == 0 ){
        longDistInput[[i]] = NA
      }
    }
    
    return(longDistInput)
}












#  plotClusterAgreementHeat
#'
#'
#' Plot a heatmap that plots the agreements between replicates 
#' after clusterrnaCrosslink has been performed
#'
#' @param cds A rnaCrosslinkDataSet object 
#' @param analysisStage The stage of the analysis to plot
#' @name plotClusterAgreementHeat
#' @docType methods
#' @rdname plotClusterAgreementHeat
#' @aliases plotClusterAgreementHeat,rnaCrosslinkDataSet-method
#' @return A heatmap of the agreement between replicates in the analysis stage chosen
#' @examples 
#' 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' 
#' clusteredCds = clusterrnaCrosslink(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 0)
#' 
#' 
#' plotClusterAgreementHeat(cds)
#' 
#' 
#'
#' @export
setGeneric("plotClusterAgreementHeat", 
           function(cds,
                    analysisStage = 'originalClusters')
             standardGeneric("plotClusterAgreementHeat") )

setMethod("plotClusterAgreementHeat", 
          "rnaCrosslinkDataSet", function(cds,
                                      analysisStage = 'originalClusters')  {
            
            # Get the cluster matrices
            samples= group(cds)$s
            if(length(samples) < 2){message("try using more than one replicate")
              }else{
            
            # Make an emtpy
            mat = matrix(0, 
                         nrow = nrow(getData(cds,
                                             "matrixList", 
                                             analysisStage)[[samples[1]]]), 
                         ncol = ncol(getData(cds,
                                             "matrixList", 
                                             analysisStage)[[samples[1]]]))
            # Add overlap
            df = data.frame()
            for(i in 1:length(samples)){
              x = getData(cds,"matrixList", analysisStage)[[samples[i]]]
              x[x>0] = 1
              mat = mat + x
              df = rbind.data.frame(df,table(x))
            }
            # plot
            myCol = c("darkgrey", "#52934E","#3F7B39","#244420")
            heatmap3(t(mat),
                     scale = "none" ,col = myCol,
                     Rowv = NA,
                     Colv = NA,
                     useRaster = TRUE
            )
}
          })





#  plotClusterAgreement
#'
#'
#' Plot a heatmap that plots the agreements between replicates 
#' after clusterrnaCrosslink has been performed
#'
#' @param cds A rnaCrosslinkDataSet object 
#' @param analysisStage The stage of the analysis to plot
#' @name plotClusterAgreement
#' @docType methods
#' @rdname plotClusterAgreement
#' @aliases plotClusterAgreement,rnaCrosslinkDataSet-method
#' @return A heatmap of the agreement between replicates in the analysis stage chosen
#' @examples 
#' 
#' cds = makeExamplernaCrosslinkDataSet()
#' 
#' 
#' clusteredCds = clusterrnaCrosslink(cds,
#'                 cores = 1,
#'                 stepCount = 1,
#'                 clusterCutoff = 0)
#' 
#' 
#' plotClusterAgreement(cds)
#' 
#' 
#' @export
setGeneric("plotClusterAgreement", 
           function(cds,
                    analysisStage = 'originalClusters')
             standardGeneric("plotClusterAgreement") )

setMethod("plotClusterAgreement", 
          "rnaCrosslinkDataSet", function(cds,
                                      analysisStage = 'trimmedClusters')  {
            
            
            samples= group(cds)$s
            samples= group(cds)$s
            if(length(samples) < 2){message("try using more than one replicate")
            }else{
            
            mat = matrix(0, 
                         nrow = nrow(getData(cds,
                                             "matrixList", 
                                             analysisStage)[[samples[1]]]), 
                         ncol = ncol(getData(cds,
                                             "matrixList", 
                                             analysisStage)[[samples[1]]]))
            
            df = data.frame()
            for(i in 1:length(samples)){
              x = getData(cds,"matrixList", analysisStage)[[samples[i]]]
              x[x>0] = 1
              mat = mat + x
              df = rbind.data.frame(df,table(x))
            }
            
            df$sample = sampleTable(cds)$sampleName[samples]
            colnames(df) = c(0,1,"sample")
            df = melt(df, id.vars = list("sample"))
            df2 = as.data.frame(table(mat))
            df2$sample = "combinedMat"
            colnames(df2) = c("variable", "value","sample")
            df = rbind.data.frame(df,df2[,c(3,1,2)])
            
            a = ggplot(df[df$variable!=0,]) +
              geom_bar(aes(x = sample,y = value,fill = variable), stat = "identity") +
              theme_classic()
            
            
            # get the cluster agreement 
            clusterDF = data.frame()
            for(sample in samples){
              
              clusters = cds@clusterTableList[[rnas(cds)]][[analysisStage]][[sample]]
              numberOfClusters = nrow(clusters)
              highestValue = rep(0,nrow(clusters))
              for(i in 1:nrow(clusters)){
                highestValue[i] = (max(mat[clusters$ls[i]:clusters$le[i],clusters$rs[i]:clusters$re[i]]))
              }
              d = data.frame(sample = sampleTable(cds)$sampleName[sample],
                             highestValue = highestValue,
                             numberOfClusters = numberOfClusters) 
              clusterDF = rbind.data.frame(clusterDF,d)
            }
            clusterDF
            
            c = ggplot(clusterDF) +
              geom_bar(aes(x = sample, fill = as.factor(highestValue)),stat = "count",
                       position = "fill")+
              theme_classic()
            
            a+c
            }
          })
















