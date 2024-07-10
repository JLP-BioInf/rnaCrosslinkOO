#' @include rnaCrosslinkOO.R
NULL


#' rnaCrosslinkQC
#' 
#' get a plot fo the read lengths and transcripts in the dataset
#' The fucntion will make 1 pdf and 2 text file in the directory provided
#'
#' @param sampleTable string - The address of the sample table, the sample table
#'  must have 4 columns, fileName (the full path and file name of the input
#'  Input file for each sample ), group ("s" - sample or "c" - control),
#'  sample (1,2,3, etc), sampleName (must be unique).
#' @param directory A directory address to write the files
#' @param topTranscripts If FALSE a table of top trandscirpts will not be written to file
#' 
#' @return ggplot and txt file
#' @name rnaCrosslinkQC
#' @docType methods
#' @rdname rnaCrosslinkQC
#' @examples 
#'   c4 = c(rep("transcript1",100),rep("transcript2",100) )
#'  c10 = c(rep("transcript1",200) )
#'  c1 = 1:200
#'  c2 = rep(paste(rep("A", 40), collapse = ""),200)
#'  c3 = rep(".",200)
#'  c9 = c3
#'  c15 = c3
#'  c5 = rep(1,200)
#'  c11 = rep(21,200)
#'  c6 = rep(20,200)
#'  c12= rep(40,200)
#'  # short distance 50
#'  
#'  # long distance 50
#'  c7 = sample(1:5, 100, replace = TRUE)
#'  c8 = sample(20:25, 100, replace = TRUE)
#'  
#'  # inter RNA 100
#'  c7 = c(c7,sample(1:5, 100, replace = TRUE))
#'  c8 = c(c8,sample(20:25, 100, replace = TRUE))
#'  c13 = c7
#'  c14 = c8
#'  
#'  exampleInput = data.frame(V1 = c1,
#'                            V2 = c2,
#'                            V3 = c3,
#'                            V4 = c4,
#'                            V5 = as.numeric(c5),
#'                            V6 = as.numeric(c6),
#'                            V7 = as.numeric(c7),
#'                            V8 = as.numeric(c8),
#'                            V9 = c9,
#'                            V10 = c10,
#'                            V11 = as.numeric(c11),
#'                            V12 = as.numeric(c12),
#'                            V13 = as.numeric(c13),
#'                            V14 = as.numeric(c14),
#'                            V15 = c15)
#'  
#'  
#'  file = tempfile()
#'  write.table(exampleInput,
#'              file = file, 
#'              quote = FALSE,
#'              row.names = FALSE, 
#'              sep = "\t", col.names = FALSE)
#'
#'  
#'  # Set up the sample table. ----
#'  sampleTabler1 = c(file, "s", "1", "s1")
#'  sampleTabler2 = c(file, "c", "1", "c1")
#'  # make the sample table 
#'  sampleTable2 = rbind.data.frame(sampleTabler1, sampleTabler2)
#'  # add the column names 
#'  colnames(sampleTable2) = c("file", "group", "sample", "sampleName")
#' 
#' rnaCrosslinkQC(sampleTable2,tempdir())
#' @export
rnaCrosslinkQC = function(sampleTable, 
                          directory,
                          topTranscripts = TRUE){
  if(topTranscripts == TRUE){
  message(" ******************************************** ")
  message(" *****            Collect Metrics      ****** ")
  message(" ******************************************** ")
  message(" *****-------*******************-------****** ")
  message(" *****       Reading SampleTable       ****** ")
  
  
  # Read in sample table
  
  #check for more than two samples
  # if( nrow(sampleTable) < 2 ){
  #        stop( "The sample Table must contain at least 1
  #              sample and 1 control" )
  #    }
  message(paste(" *****       Detected ",
                nrow(sampleTable), " Samples      ****** "))
  
  
  
  ###########################################################
  #check column names of sampleTable
  colnamesST = c("file", "group", "sample", "sampleName")
  if (all(colnames(sampleTable) != colnamesST)) {
    stop("Column names of metaData table should be :
              file, group, sample, sampleNames")
  }
  
  # Get the comparison groups check group has the c and s
  if (!(unique(as.character(sampleTable$group))[1] %in% c("c", "s") &
        unique(as.character(sampleTable$group))[2] %in% c("c", "s"))) {
    stop("Groups should be c and s")
  }
  
  
  # Make group into a list with control and sample
  group = sampleTable[, "group"]
  group2 = list()
  group2[["c"]] = which(group == "c")
  group2[["s"]] = which(group == "s")
  
  group = group2
  message(paste(
    " *****     detected group c::",
    paste(group[["c"]],
          collapse = " ") ,
    paste(rep(" ",
              
              2),
          collapse = ""),
    "   ***** "
  ))
  message(paste(
    " *****     detected group s::",
    paste(group[["s"]],
          collapse = " ") ,
    paste(rep(" ",
              (
                2)),
          collapse = ""),
    "   ***** "
  ))
  
  
  
  ###########################################################
  # Get the sampleNames
  sampleNames = c()
  if (is.null(sampleTable$sampleName)) {
    stop("The sample Table must have a column named sampleName")
  } else if (length(unique(sampleTable$sampleName)) !=
             length(sampleTable$sampleName)) {
    stop("Sample names must be unique")
  } else{
    sampleNames = as.character(sampleTable$sampleName)
    spaces =  (length(sampleNames) * (3 - length(sampleNames))) * 2 
    if(spaces < 0 ){spaces = 0}
    message(paste(
      " ****  ",
      paste(rep(" ",
                spaces,
                collapse = ""),
            " Sample Names: ",
            paste(sampleNames, collapse = " "),
            paste(rep(" ",
                      spaces,
                      collapse = ""),
                  " **** "
            ))))
  }
  
  
  
  ###########################################################
  # Read in the  Input files
  ###########################################################
  #load the files into a list
  message(" *****         Reading Input Files        ***** ")
  
  InputFiles = list()
  InputFiles[["all"]] = list()
  InputFiles = list()
  
  
  
  # read in the tables
  inputs <- lapply(as.character(sampleTable$file),
                   function(file)
                     read.table(file,colClasses = c("character",
                                                    "character",
                                                    "character",
                                                    "character",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "character",
                                                    "character",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "character")))
  
  #check column names
  if (all(sapply(inputs, function(file)
    ! (identical(
      colnames(file),
      c(
        "V1",
        "V2",
        "V3",
        "V4",
        "V5",
        "V6",
        "V7",
        "V8",
        "V9",
        "V10",
        "V11",
        "V12",
        "V13",
        "V14",
        "V15"
      )
    ))))) {
    stop(
      " The input files do not look they are produced with the
                 nextflow pipeline, please check the documentation. "
    )
  }
  
  
  # now get the most interacting partners and transcripts as well as the 
  # read size distriubtion for each side.
  
  widthLeft <- lapply(inputs, function(x) x$V6-x$V5)
  widthRight <- lapply(inputs, function(x) x$V12-x$V11)
  
  
  
  widthLeft <- lapply(1:nrow(sampleTable), 
                      function(x) { 
                        x = data.frame(width = widthLeft[[x]],
                                      type = "left",
                                      sample =  sampleTable[x,"sampleName"] )
                      })
  
  
  
  widthRight <- lapply(1:nrow(sampleTable), 
                      function(x) { 
                        x = data.frame(width = widthRight[[x]],
                                       type = "right",
                                       sample =  sampleTable[x,"sampleName"] )
                      })
  
  
  
  widthdf = rbind.data.frame(do.call(rbind,widthLeft ),do.call(rbind,widthRight ))
  
  
  
  widthdf = aggregate(widthdf$width, 
                      by = list(widthdf$width, 
                      widthdf$type,widthdf$sample),
                      FUN = length )
  colnames(widthdf) = c("width","type","sample","count")

  plot(ggplot(widthdf) +
    geom_bar(mapping = aes(x = width, y = count), stat = "identity") +
    facet_grid(sample~type) +
    theme_bw())
  

  
  
  # now get the transcript sin the dataset
  

  
  ci = group[["c"]]
  si = group[["s"]]
  
  
  
  
  
  df= data.frame()
  
  for (i in 1:length(inputs)) {
    
    
    same = which(inputs[[i]]$V4 ==
                   inputs[[i]]$V10 )
    notSame = which(inputs[[i]]$V4 !=
                      inputs[[i]]$V10)
    
    
    s =   paste(inputs[[i]]$V4[same],
                inputs[[i]]$V10[same], sep = "::")
    s = unlist(lapply(s, function(x)  
      unlist(strsplit(x,split = "::"))[1] ))
    s = as.data.frame(table(s))
    s$type="intra"
    s$sample = i
    colnames(s) = c("RNA","reads","type","sample")
    
    
    
    t =   paste(inputs[[i]]$V4[notSame],
                inputs[[i]]$V10[notSame], sep = "::")
    
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
  
  
  
  ci = sampleTable[ci,"sampleName"]
  si = sampleTable[si,"sampleName"]
  
  colnames(df) = c("rna","type",sampleTable$sampleName)
  
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

  
  
  df
  
  
  
  
  
  # now write these to a file 
  
  write.table(df, 
              file = paste(directory,"/topTranscripts_all.txt", sep = ""),
              quote = FALSE, row.names = FALSE)
  }

  
  
}


















