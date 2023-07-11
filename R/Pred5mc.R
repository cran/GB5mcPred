#' @title Pred5mc
#' @description Prediction of sequences with 5mc sites.
#' @param Fastafile Sequence file (.fasta format)
#' @return MethStatus: Sequences with their methylation state (methylated or non-methylated)
#' @import stats devtools tidyverse seqinr Biostrings splitstackshape entropy party stringr tibble doParallel parallel e1071 caret randomForest gbm foreach iterators ftrCOOL
#' @export
#'
#' @examples
#' \donttest{
#' library(GB5mcPred)
#' data<-system.file("exdata/test.fasta", package = "GB5mcPred")
#' pred<-Pred5mc(Fastafile=data)
#' }
#' @references Chen, W., Lv, H., Nie, F., & Lin, H. (2019). i6mA-Pred: identifying DNA N6-methyladenine 	sites in the rice genome. Bioinformatics, 35(16), 2796-2800.
requireNamespace("Biostrings","caret","devtools","doParallel","e1071","entropy","foreach","ftrCOOL","gbm","iterators","parallel","party","randomForest","seqinr","splitstackshape","stats","stringr","tibble","tidyverse")

Pred5mc<- function(Fastafile){
  utils::globalVariables("gb_model")
  gb_model<-NULL
  gb_model<-readRDS("inst/exdata/gradient_boost.rda")
  FastaToTabular <- function (filename){

    #read fasta file

    file1 <- readLines(filename)

    #find the genename location by grepping >

    location <- which((str_sub(file1,1,1))==">")

    #start an empty vector to collect name and sequence

    name=c()
    sequence =c()



    #number of genes= number of loops
    #extract name first
    for ( i in 1:length(location)){
      name_line = location[i]
      name1 = file1[name_line]
      name=c(name,name1)
      #extract sequence between the names
      #the last sequence will be missed using this strategy
      #so, we are using if condition to extract last sequence
      start= location[i]+1
      end = location[i+1]-1
      if ( i < length (location)){

        end=end

      } else {

        end=length(file1)
      }

      lines = start:end
      sequence1= as.character(paste(file1[lines],collapse = ""))
      sequence =c(sequence,sequence1)
    }

    #now create table using name and sequence vector

    data <- tibble(name,sequence)



    #finally export the file
    #before that remove preexisting file
    unlink(c("dna_table.csv"),force=TRUE)
    as.matrix(data,"dna_table.csv")

    #function ends
  }
  #########################alphabetcheck###########################
  alphabetCheck<-function (sequences, alphabet = "aa", label = c())
  {
    if (length(sequences) == 0) {
      stop("ERROR: sequence parameter is empty")
    }
    if (length(label) != 0 && length(label) != length(sequences)) {
      stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
    }
    if (alphabet == "rna") {
      alphabet <- c("A", "C", "G", "U")
    }
    else if (alphabet == "dna") {
      alphabet <- c("A", "C", "G", "T")
    }
    else if (alphabet == "aa") {
      alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I",
                    "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                    "W", "Y")
    }
    else {
      stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
    }
    alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                               split = "")[[1]] %in% alphabet))
    flag = 0
    if (length(label) == length(sequences)) {
      flag = 1
      label = label[alphabetCheck]
    }
    else if (length(label) > 0 && length(label) != length(sequences)) {
      stop("ERROR: The number of labels is not equal to the number of sequences!")
    }
    if (is.null(names(sequences))) {
      names(sequences) <- as.character(1:length(sequences))
    }
    nonstanSeq <- names(sequences)[!alphabetCheck]
    if (length(nonstanSeq) != 0) {
      nonstanSeq <- toString(nonstanSeq)
      warMessage <- paste("The sequences (", nonstanSeq, ") were deleted. They contained non-standard alphabets")
      message(warMessage)
    }
    sequences = sequences[alphabetCheck]
    if (length(sequences) == 0) {
      stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
    }
    if (flag == 1) {
      names(label) = names(sequences)
    }
    seq_lab <- list(sequences = sequences, Lab = label)
    return(seq_lab)
  }
  #################################NCP_DNA############################

  ncp_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                     label = c())
  {
    if (length(seqs) == 1 && file.exists(seqs)) {
      seqs <- fa.read(seqs, alphabet = "dna")
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else if (is.vector(seqs)) {
      seqs <- sapply(seqs, toupper)
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else {
      stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
    }
    lenSeqs <- sapply(seqs, nchar)
    nucs <- list(A = c(1, 1, 1), C = c(0, 0, 1), G = c(1, 0,
                                                       0), T = c(0, 1, 0), U = c(0, 1, 0))
    numSeqs <- length(seqs)
    if (outFormat == "mat") {
      if (length(unique(lenSeqs)) > 1) {
        stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
      }
      if (binaryType == "strBin") {
        nucs <- c(A = "111", C = "001", G = "100", T = "010",
                  U = "010")
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        colnames(featureMatrix) <- paste("pos", 1:lenSeqs[1],
                                         sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "logicBin") {
        nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                    TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                   FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("pos", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "numBin") {
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("pos", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else {
        stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
      }
      return(featureMatrix)
    }
    else if (outFormat == "txt") {
      nucs <- c(A = "111", C = "001", G = "100", T = "010",
                U = "010")
      counter <- 0
      namesSeqs <- names(seqs)
      codes <- lapply(seqs, function(x) {
        counter <- counter + 1
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        namecods <- namesSeqs[counter]
        cods <- unlist(cods)
        cods <- c(namecods, cods)
        temp <- paste(cods, collapse = "\t")
        write(temp, outputFileDist, append = TRUE)
      })
    }
    else {
      stop("ERROR: outFormat should be 'mat' or 'txt' ")
    }
  }



  ###########################GC_Content##########################
  GC.content <- function(fasta_file_train){
    x <- read.fasta(file=fasta_file_train)
    tt<-function(x){
      res<-GC(x)
      val=round(res,4)
      return(val)
    }

    f_res<-lapply(x,tt)
    s=data.frame(f_res)

    rownames(s) <- c("GC-content")

    w=t(s)
    return(w)
  }
  ################################ONF#################################
  oligo.freq <- function(fasta_file_train,f){
    x<- readDNAStringSet(fasta_file_train)
    y <- oligonucleotideFrequency(x,width = f)
    z <- data.frame(y)
    rownames(z) <- names(x)

    return(z)
  }


  #######################mononucleotide_binary_encoding##################################

  FastaToTabular <- function (filename){

    #read fasta file

    file1 <- readLines(filename)

    #find the genename location by grepping >

    location <- which((str_sub(file1,1,1))==">")

    #start an empty vector to collect name and sequence

    name=c()
    sequence =c()



    #number of genes= number of loops
    #extract name first
    for ( i in 1:length(location)){
      name_line = location[i]
      name1 = file1[name_line]
      name=c(name,name1)
      #extract sequence between the names
      #the last sequence will be missed using this strategy
      #so, we are using if condition to extract last sequence
      start= location[i]+1
      end = location[i+1]-1
      if ( i < length (location)){

        end=end

      } else {

        end=length(file1)
      }

      lines = start:end
      sequence1= as.character(paste(file1[lines],collapse = ""))
      sequence =c(sequence,sequence1)
    }

    #now create table using name and sequence vector

    data <- tibble(name,sequence)



    #finally export the file
    #before that remove preexisting file
    unlink(c("dna_table.csv"),force=TRUE)
    as.matrix(data,"dna_table.csv")

    #function ends
  }
  #########################alphabetcheck###########################
  alphabetCheck<-function (sequences, alphabet = "aa", label = c())
  {
    if (length(sequences) == 0) {
      stop("ERROR: sequence parameter is empty")
    }
    if (length(label) != 0 && length(label) != length(sequences)) {
      stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
    }
    if (alphabet == "rna") {
      alphabet <- c("A", "C", "G", "U")
    }
    else if (alphabet == "dna") {
      alphabet <- c("A", "C", "G", "T")
    }
    else if (alphabet == "aa") {
      alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I",
                    "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                    "W", "Y")
    }
    else {
      stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
    }
    alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                               split = "")[[1]] %in% alphabet))
    flag = 0
    if (length(label) == length(sequences)) {
      flag = 1
      label = label[alphabetCheck]
    }
    else if (length(label) > 0 && length(label) != length(sequences)) {
      stop("ERROR: The number of labels is not equal to the number of sequences!")
    }
    if (is.null(names(sequences))) {
      names(sequences) <- as.character(1:length(sequences))
    }
    nonstanSeq <- names(sequences)[!alphabetCheck]
    if (length(nonstanSeq) != 0) {
      nonstanSeq <- toString(nonstanSeq)
      warMessage <- paste("The sequences (", nonstanSeq, ") were deleted. They contained non-standard alphabets")
      message(warMessage)
    }
    sequences = sequences[alphabetCheck]
    if (length(sequences) == 0) {
      stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
    }
    if (flag == 1) {
      names(label) = names(sequences)
    }
    seq_lab <- list(sequences = sequences, Lab = label)
    return(seq_lab)
  }
  #################################MBE_DNA############################

  mbe_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                     label = c())
  {
    if (length(seqs) == 1 && file.exists(seqs)) {
      seqs <- fa.read(seqs, alphabet = "dna")
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else if (is.vector(seqs)) {
      seqs <- sapply(seqs, toupper)
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else {
      stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
    }
    lenSeqs <- sapply(seqs, nchar)
    nucs <- list(A = c(1, 0, 0, 0), C = c(0, 1, 0, 0), G = c(0, 0, 1, 0), T = c(0, 0, 0, 1), U = c(0, 0, 0, 1))
    numSeqs <- length(seqs)
    if (outFormat == "mat") {
      if (length(unique(lenSeqs)) > 1) {
        stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
      }
      if (binaryType == "strBin") {
        nucs <- c(A = "1000", C = "0100", G = "0010", T = "0001",
                  U = "0001")
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        colnames(featureMatrix) <- paste("pos_mbe", 1:lenSeqs[1],
                                         sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "logicBin") {
        nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                    TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                   FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("pos_mbe", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "numBin") {
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("pos_mbe", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else {
        stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
      }
      return(featureMatrix)
    }
    else if (outFormat == "txt") {
      nucs <- c(A = "1000", C = "0100", G = "0010", T = "0001",
                U = "0001")
      counter <- 0
      namesSeqs <- names(seqs)
      codes <- lapply(seqs, function(x) {
        counter <- counter + 1
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        namecods <- namesSeqs[counter]
        cods <- unlist(cods)
        cods <- c(namecods, cods)
        temp <- paste(cods, collapse = "\t")
        write(temp, outputFileDist, append = TRUE)
      })
    }
    else {
      stop("ERROR: outFormat should be 'mat' or 'txt' ")
    }
  }

  ######################Test data######################
  fasta_file_test<-Fastafile
  res<-FastaToTabular(fasta_file_test)
  data<-as.vector(res[,2])
  mat<-as.matrix(ncp_dna(seqs = data,binaryType="strBin",outFormat="mat"))
  sequence<-rownames(mat)
  seq_id<-res[,1]
  ncp<-cbind(seq_id,sequence,mat)
  rownames(ncp)<-seq_id
  ncp_temp<-data.frame(ncp[,-1], stringsAsFactors = FALSE)
  ncp_final<-as.data.frame(apply(ncp_temp[,-1], 2, as.numeric))
  log_gc_temp<-log((GC.content(fasta_file_test))*100, base = exp(1))
  log_gc<-as.data.frame(as.numeric((ifelse(log_gc_temp>0,log_gc_temp,'0'))))
  onf<-oligo.freq(fasta_file_test, 2)
  res_temp_mbe<-FastaToTabular(fasta_file_test)
  data_mbe<-as.vector(res_temp_mbe[,2])
  res_mbe<-as.matrix(mbe_dna(seqs = data_mbe,binaryType="strBin",outFormat="mat"))
  mbe_temp<-data.frame(res_mbe[,-1], stringsAsFactors = FALSE)
  mbe_final<-as.data.frame(apply(mbe_temp[,-1], 2, as.numeric))
  temp1<- cbind(onf, gcc =log_gc[,1], ncp_final, mbe_final)
  temp_col <- c("pos_mbe22", "pos22", "pos20", "CG", "pos_mbe20", "GG", "pos_mbe23",
                "GC", "pos19", "pos12", "pos_mbe12", "CC", "pos18", "pos_mbe18", "TC",
                "pos_mbe19", "pos_mbe13", "pos13", "pos4", "pos3", "pos27", "pos23")
  my_data_temp <- temp1[,temp_col]
  inputData <-as.data.frame(my_data_temp)

  ######################training & Testing data######################



  test_data <- inputData

  ################## Gradient_boost ##############
  predicted_value_gb <- predict(gb_model,newdata = test_data)


  temp_final1<- as.matrix(predicted_value_gb)
  temp_final2<-(ifelse(temp_final1==2,'5mC','Non-5mC'))
  colnames(temp_final2)<- "Status"
  MethStatus<-as.data.frame(cbind(Sequence=res, temp_final2))
  MethStatus$name <- gsub(">","",as.character(MethStatus$name))
  colnames(MethStatus) <- c("IDs", "Sequence", "Status")
  return(MethStatus)
}

