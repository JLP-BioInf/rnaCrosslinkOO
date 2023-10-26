#' @include  comradesDataSet.R 
NULL


#' getMatrices
#' 
#' Make a matrix of contact interactions
#'
#' Function used to create a list of matrices for plotting with
#' plotMatrixList or plotMatrixListFull, the output list will be same as the
#' input except for an extra list layer for the specific RNA
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna the RNA of interest that you want to subset
#' @param size The size of the RNA
#' @return A list of matrices
#' @name getMatrices
#' @docType methods
#' @rdname getMatrices
getMatrices = function(hybList, 
                       rna, 
                       size){
    hybMatList = list()
    for(hyb in 1:length(hybList)){

        hybOutputO = hybList[[hyb]]
        #print(hybOutputO)
        
        
        hybOutput =  hybOutputO[as.character(hybOutputO$V4) == rna & as.character(hybOutputO$V10) == rna,]
        hybOutput = unique(hybOutput)
        
        startsends = hybOutput[,c(7,8,13,14)]
        

        
        mat = matrix(0, ncol = size, nrow =size)
        for(i in 1:nrow(startsends)){
            data = startsends[i,]
            xData = seq(data$V7, data$V8)
            yData = seq(data$V14, data$V13)
            mat[xData,yData] = mat[xData,yData] +1
            
            
        }
        
        hybMatList[[hyb]] = mat
        
    }
    return(hybMatList)
}






#' getAdjacancyMat
#' 
#' Makes and adjacency matrix list (for clustering)
#'
#' Makes and adjacency matrix list (for clustering)
#' @param hybGranges list created with hybToGRanges (but just the gap section of the list)
#' @param  nucletideOrPerc measure difference by percentage or nucleotides
#' @param  cutoff The maximum difference before giving these two gaps 0
#' @return A list of Adjacancy matrices
#' @name getAdjacancyMat
#' @docType methods
#' @rdname getAdjacancyMat
getAdjacancyMat = function(hybGranges, nucletideOrPerc, cutoff){
    distances = hybGranges
    max = max(width(distances))
    
    #get overlapping
    hits <- findOverlaps(distances, drop.self=T, drop.redundant=F)
    # get the relative overlap for the weights
    x <- distances[queryHits(hits)]
    y <- distances[subjectHits(hits)]
    
    relative_overlap <-  width(pintersect(x, y)) / pmax(width(x), width(y))
    
    hitsWithOverlap = hits
    # parameter for relative overlap
    
    #print(length(hitsWithOverlap))
    
    if(nucletideOrPerc == "none"){
        hitsWithOverlap = hits
    } else if(nucletideOrPerc == "nucleotide"){
        relative_overlap =     pmax(width(x), width(y)) - width(pintersect(x, y))
        relative_overlap = cutoff - relative_overlap
        
        
        hitsWithOverlap = hits[relative_overlap <= cutoff & relative_overlap >= 0 ]
        relative_overlap = relative_overlap[relative_overlap <= cutoff  & relative_overlap >= 0 ]
    } else if(nucletideOrPerc == "perc"){
        
        relative_overlap = (1- (width(pintersect(x, y)) / max))
        
        hitsWithOverlap = hits[relative_overlap >= cutoff]
        if(length(hitsWithOverlap ) ==0 ){return(0)}else{}
        # print(length(hitsWithOverlap))
        relative_overlap = relative_overlap[relative_overlap >= cutoff]
        
    }
    if(length(hitsWithOverlap) == 0){
        return(NA)
    }else{
        
        
        hitsMat = as.data.frame(hitsWithOverlap)
        hitsMat$weight = relative_overlap
        
        
        testLong = dcast(hitsMat, queryHits ~ subjectHits, value.var = "weight")
        row.names(testLong) = testLong$queryHits
        testLong = testLong[,-1]
        testLong = as.matrix(testLong)
        testLong[is.na(testLong)] =0
        
        return(testLong)
    }
    
}




