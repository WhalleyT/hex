library("Biostrings")
library("httr")
library("parallel")
library("doParallel")

#data(BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM100, PAM30, PAM40, PAM70, PAM120, PAM250)
utils::data(list=c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"), envir=environment())
utils::globalVariables(names=c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"))

#' Get position specific weights, higher scores for the mid bases compared to tails.
#'
#' Compute the position specific scores for the amino acid bases with heavier weight for the middle bases 
#' and lower score for bases towards the edges. Scores are a function of amino acid bases vector length
#' and the magnitude of score elevation by base distance from edges, bases on the third quartile of edges
#' are simplified to magnitude of 1 beginning from 1.
#'
#' @param x A character vector of Amino Acid bases.
#' @param mag Numeric value to use for seq() function to create distant sequence of weights
#' @return A numeric vector of scores equal to the length of the input vector of Amino Acids
#' @examples
#' \dontrun{
#' get_pos_weights(x=AA.vec, mag=4)
#' }
#' @keywords internal
#' @export
get_pos_weights <- function(x, mag=4){
	posScore <- c()
	xLen <- length(x)
	xMid <- xLen/2
	if(xLen %% 2){
		#print("ODD?")
		xMidCeil <- ceiling(xMid)
		xMidScore <- xMidCeil*mag
		posScore <- seq(1,xMidScore,mag)
		posScore <- c(posScore, rev(posScore[1:length(posScore)-1]))
		xTopFloor <- floor(xLen/3)
		posScore[1:xTopFloor] <- 1:xTopFloor
		xTail <- xLen-xTopFloor+1
		posScore[xTail:xLen] <- xTopFloor:1
		
	}else{
		#print("EVEN?")
		xMidScore <- xMid*mag
		posScore <- seq(1,xMidScore,mag)
		posScore <- c(posScore, rev(posScore))
		xTopFloor <- floor(xLen/3)
		posScore[1:xTopFloor] <- 1:xTopFloor
		xTail <- xLen-xTopFloor+1
		posScore[xTail:xLen] <- xTopFloor:1
	}
	return(posScore)
}

#' Alignment of singular epitope vector against a reference matrix of epitope, rows as epitopes with names and bases in 
#' columns by position.
#'
#' Align the epitope sequence to the reference set of epitope sequences. Sequences are alignment globally
#' from end to end with matches assigned a postitive score as the position specific weight and mis-matches
#' are assigned a negative score as the negative of the position specific weight. Final alignment score is
#' the sum of all the positive and neagative scores from each base. This ensures that the alignment is matches
#' in the middle are scored higher than the alignments with mis-matches in the middle.
#'
#' @param x A character vector of Amino Acid bases.
#' @param ref A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns
#' @param mag Numeric value to use for seq() function to create distant sequence of weights
#' @return A numeric vecdtor of scores for x alignment to the each epitope in ref
#' @examples
#' \dontrun{
#' align_to_ref_epitopes(x, ref, mag=4)
#' }
#' @keywords internal
#' @export
align_to_ref_epitopes <- function(x, ref, mag=4, aln_matrix="BLOSUM62"){
        if(!exists(aln_matrix)){
                #data(aln_matrix)
                data(aln_matrix, envir=environment())
        }
        aln_matrix <- get(aln_matrix)
	xPosWt <- get_pos_weights(x, mag)
	#print(xPosWt)
	xLen <- length(x)
	
	xAlnScoreVec <- c()
	
	if(is.matrix(ref)){
		refSeqCount <- nrow(ref)
	}else if(is.character(ref)){
		refSeqCount <- 1
	}
	else{
		stop("Unsupported object as ref. Must be either matrix or character")
	}
	for(i in 1:refSeqCount)
	{
		xAlnScore <- 0
		#match <- ""
		#misMatch <- ""
		refSeq <- ref[i,]
		refName <- rownames(ref)[i]

                xAlnMat <- diag(aln_matrix[x,refSeq])
                xAlnAmp <- xPosWt * xAlnMat
                xAlnScore <- sum(xAlnAmp)

		#match <- which(x == refSeq)
		##print(match)
		#misMatch <- which(!x == refSeq)
		##print(misMatch)
		#
		#if(length(match)>0){
		#	#print(sum(xPosWt[match]))
		#	xAlnScore <- sum(xPosWt[match])
		#}
		#if(length(misMatch)>0){
		#	#print(sum(0-xPosWt[misMatch]))
		#	xAlnScore <- xAlnScore + (0-sum(xPosWt[misMatch]))
		#}
		xAlnScoreVec[i] <- xAlnScore
		names(xAlnScoreVec)[i] <- refName
	}
	#return(xAlnScoreList)
	return(xAlnScoreVec)
}

#' Read the epitope amino acid sequences from the fasta file
#'
#' Read epitope sequences from a fasta file, provided with the correct path to the file.
#' The read sequences are provided as a matrix with sequence in rows and each base of the
#' sequence as a column.
#'
#' @importFrom Biostrings readAAStringSet as.matrix
#'
#' @param iFile Path to a file containing epitope sequences in fasta format.
#' @return A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns
#' @examples
#' \dontrun{
#' read_epitope_seq(iFile=input.file)
#' }
#' @keywords internal
#' @export
read_epitope_seq <- function(iFile){
	epitopes.AASet <- Biostrings::readAAStringSet(iFile)
	epitopes.mat <- Biostrings::as.matrix(epitopes.AASet, use.names=T)
	return(epitopes.mat)
}

#' Predict the MHC-1 eptiopes from the protein sequence by using the IEDB server API
#'
#' Pass the protein sequence to IEDB server for predicting MHC-1 binding epitopes.
#' The predicted epitope sequences are parsed and formatted as a matrix with 
#' sequence in rows and each base of the sequence as a column.
#'
#' @importFrom Biostrings AAStringSet as.matrix
#' @importFrom httr POST warn_for_status content
#' @importFrom stats setNames
#'
#' @param sequence.name Name of the input protein sequence for predicting epitopes.
#' @param sequence Protein sequence from which the epitopes are predicted.
#' @param method MHC class 1 binding prediction methods.
#' @param allele MHC allele.
#' @param length Length of the peptide to be predicted.
#' @param api.url URL of the IEDB server API, default api is for MHC-1 prediction.
#' @return A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns
#' @examples
#' \dontrun{
#' predict_epitopes(sequence.name="Seq1", 
#'      sequence, method="recommended", 
#'      allele="H-2-Kb", 
#'      length=9, 
#'      api.url"http://tools-cluster-interface.iedb.org/tools_api/mhci/"
#' )
#' }
#' @keywords internal
#' @export
#' 
predict_epitopes <- function(sequence.name, sequence, method="recommended", allele="H-2-Kb", length=9, api.url="http://tools-cluster-interface.iedb.org/tools_api/mhci/"){
	print("in predict epitope")
  res <- httr::POST(api.url,
	    body=list(method = method, 
		    sequence_text = sequence, 
		    allele = allele, 
		    length = length),
	    encode="form"
	)
	httr::warn_for_status(res)

	res.text <- httr::content(res, type="text", encoding = "UTF-8")
	res.text.split <- strsplit(res.text, split="\n")
	res.DF <- t(as.data.frame(lapply(res.text.split[[1]], function(x) strsplit(x, split="\t")), stringsAsFactors=F))
	colnames(res.DF) <- res.DF[1,]
	rownames(res.DF) <- c()
	
	print("matrix")
	print(res.DF)
	if(nrow(res.DF) < 3){
	  res.DF <- as.matrix(t(res.DF[-1,]))
	}else{
	  res.DF <- as.matrix(res.DF[-1,]) 
	}
	rownames(res.DF) <- c(1:nrow(res.DF))
	print(res.DF)
	res.DF <- as.data.frame(res.DF, stringsAsFactors=F)

	if(length(grep("#", sequence.name))>0){
		sequence.name <- gsub("#", "|", sequence.name)
	}
	res.vec <- stats::setNames(res.DF$peptide, paste0(sequence.name, "#", res.DF$start, "#", res.DF$end, "#", res.DF$percentile_rank, "#", res.DF$ann_ic50))
	res.AASet <- Biostrings::AAStringSet(res.vec)
	res.mat <- Biostrings::as.matrix(res.AASet, use.names=T)
	return(res.mat)
}

#' Fetch the MHC-1 eptiopes prediction result from the IEDB server API
#'
#' Pass the protein sequence to IEDB server for predicting MHC-1 binding epitopes.
#' The predicted epitope sequences are pasrsed and formatted as a matrix with 
#' sequence in rows and each base of the sequence as a column.
#'
#' @importFrom httr POST warn_for_status content
#'
#' @param sequence Protein sequence from which the epitopes are predicted.
#' @param method MHC class 1 binding prediction methods.
#' @param allele MHC allele.
#' @param length Length of the peptide to be predicted.
#' @param api.url URL of the IEDB server API, default api is for MHC-1 prediction.
#' @return A data.frame of the IEDB MHC-1 prediction result on the input sequence
#' @examples
#' \dontrun{
#' fetch_epitopes(sequence, 
#'    method="recommended", 
#'    allele="H-2-Kb", 
#'    length=9, 
#'    api.url="http://tools-cluster-interface.iedb.org/tools_api/mhci/")
#' }
#' @keywords internal
#' @export
fetch_epitopes <- function(sequence, method="recommended", allele="H-2-Kb", length=9, api.url="http://tools-cluster-interface.iedb.org/tools_api/mhci/"){
  #httr::POST("https://eo3mlll51sol16x.m.pipedream.net")
	res <- httr::POST(api.url,
	    body=list(method = method, 
		    sequence_text = sequence, 
		    allele = allele, 
		    length = length),
	    encode="form"
	)
	httr::warn_for_status(res)

	res.text <- httr::content(res, type="text")
	res.text.split <- strsplit(res.text, split="\n")
	res.DF <- t(as.data.frame(lapply(res.text.split[[1]], function(x) strsplit(x, split="\t")), stringsAsFactors=F))
	colnames(res.DF) <- res.DF[1,]
	rownames(res.DF) <- c()
	print(res.DF)
	if(nrow(res.DF) < 3){
	  res.DF <- as.matrix(t(res.DF[-1,]))
	}else{
	  res.DF <- as.matrix(res.DF[-1,]) 
	}
	rownames(res.DF) <- c(1:nrow(res.DF))
	res.DF <- as.data.frame(res.DF, stringsAsFactors=F)
	
	return(res.DF)
}

#' Fetch the all MHC alleles for a specified species from the IEDB server API
#'
#' Pass the species information to IEDB server for fetching all the associated MHC alleles.
#'
#' @importFrom httr POST warn_for_status content
#'
#' @param species Common name of the species to query the IEDB server.
#' @param api.url URL of the IEDB server API, default api is for MHC-1 prediction.
#' @return A vector of the MHC allele associated with the specified species
#' @examples
#' \dontrun{
#' fetch_alleles(species="human", api.url="http://tools-cluster-interface.iedb.org/tools_api/mhci/")
#' }
#' @keywords internal
#' @export
fetch_alleles <- function(species="human", api.url="http://tools-cluster-interface.iedb.org/tools_api/mhci/"){
	res <- httr::POST(api.url,
	    body=list(method = "ann", 
		    species = species),
	    encode="form"
	)
	httr::warn_for_status(res)

	res.text <- httr::content(res, type="text")
	res.text.split <- strsplit(res.text, split="\n")
        res.vec <- as.vector(res.text.split[[1]])
        res.len <- length(res.vec)
        toRemove <- c(1, (res.len-2):res.len)
        res.vec <- res.vec[-toRemove]
	
	return(res.vec)
}

#' Filter the epitopes by IC50 score
#'
#' Filter the predicted epitope sequence data frame by the user provided IC50 score threshold.
#' The epitopes with IC50 score <= IC50 threshold are selected and returned.
#'
#' @param epitope.DF Data frame of the epitopes predicted by the IEDB server.
#' @param ic50.threshold A numerical threshold value for the IC50 score of epitope.
#' @return A data.frame of the predicted epitopes filtered by the IC50 threshold score
#' @examples
#' \dontrun{
#' filter_epitopes(epitope.DF, ic50.threshold)
#' }
#' @keywords internal
#' @export
filter_epitopes <- function(epitope.DF, ic50.threshold){
                selIdx <- which(as.numeric(epitope.DF$ann_ic50) <= ic50.threshold)
                if(length(selIdx)==0){
                        print("No epitopes satisfying the IC50 threshold!!!")
                        return(NULL)
                }else{
                        print("Filtering by IC50 threshold....")
	                epitope.DF <- epitope.DF[selIdx,]
                }
                return(epitope.DF)
}

#' Convert the epitopes data frame to amino acid matrix with bases separated into columns
#'
#' Data frame of predict epitopes is converted to a matrix of epitope sequences by rows and their
#' bases in columns.
#'
#' @importFrom stats setNames
#' @importFrom Biostrings AAStringSet as.matrix
#'
#' @param sequence.name Name of the input protein sequence for predicting epitopes.
#' @param epitope.DF Data frame of the epitopes predicted by the IEDB server.
#' @return A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns
#' @examples
#' \dontrun{
#' epitopeDF_to_AAmat(sequence.name, epitope.DF)
#' }
#' @keywords internal
#' @export
epitopeDF_to_AAmat <- function(sequence.name, epitope.DF){
        if(length(grep("#", sequence.name))>0){
		sequence.name <- gsub("#", "|", sequence.name)
	}
	res.vec <- stats::setNames(epitope.DF$peptide, paste0(sequence.name, "#", epitope.DF$start, "#", epitope.DF$end, "#", epitope.DF$percentile_rank, "#", epitope.DF$ann_ic50))
	res.AASet <- Biostrings::AAStringSet(res.vec)
	res.mat <- Biostrings::as.matrix(res.AASet, use.names=T)
	return(res.mat)
}

#' Align set of query epitope sequences against set of reference epitope sequences
#'
#' Aligning set of query epitopes by multi-threaded parallel execution of align_to_red_epitopes()
#' function. The alignment score between all pairs of query and ref epitopes is reported as a matrix.
#'
#' @importFrom parallel makeCluster clusterExport parApply stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @param query.set A matrix of Amino Acid epitope sequences to align against the reference epitope sequences.
#' @param ref.set A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns.
#' @return A matrix of alignment scores for query epitopes against reference epitopes
#' @examples
#' \dontrun{
#' align_sets(query.set, ref.set, env=NULL, num.cores=2)
#' }
#' @keywords internal
#' @export
align_sets <- function(query.set, ref.set, env=NULL, num.cores=2, aln_matrix="BLOSUM62"){
        if(!exists(aln_matrix)){
                print(paste0("Getting data:", aln_matrix))
                #data(aln_matrix)
                data(aln_matrix, envir=environment())
        }
	print("Create Cluster...")
	cl <- parallel::makeCluster(num.cores, outfile="")
	print("Register Cluster...")
	doParallel::registerDoParallel(cl)
	print("Parse Res Set Name...")
	ref.set.name <- deparse(substitute(ref.set))
	ref.set.name <- strsplit(ref.set.name, "\\[")[[1]][1]
	print("Export Functions and Ref Set Name...")
        if(is.null(env)){
            env <- environment()
        }
	#parallel::clusterExport(cl, list("get_pos_weights", "align_to_ref_epitopes", ref.set.name))
	parallel::clusterExport(cl, list("get_pos_weights", "align_to_ref_epitopes", ref.set.name, aln_matrix), envir=env)

	print("par Apply Alignment...")
	t1 <- proc.time()
	alnRes <- parallel::parApply(cl, query.set, 1, function(x) align_to_ref_epitopes(x, ref=ref.set, aln_matrix=aln_matrix))
	print(paste0("Alignment Time : ", (proc.time()-t1)[3]))

	print("Kill Cluster...")
	parallel::stopCluster(cl)
	return(alnRes)
}

#' Get the pairs of query and reference epitopes with the best alignment score or above the 
#' specified alignment score threshold
#'
#' Aligning set of query epitopes by multi-threaded parallel execution of align_to_red_epitopes()
#' function. The alignment score between all pairs of query and ref epitopes is reported as a matrix.
#'
#' @importFrom parallel makeCluster clusterExport parApply stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @param alnRes A matrix of alignment scores for query epitopes against reference epitopes.
#' @param query.set A matrix of Amino Acid epitope sequences to align against the reference epitope sequences.
#' @param ref.set A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns.
#' @param best A boolean flag to switch between selecting the pairs of epitopes with best alignment score or use the 
#' input score value as a threshold.
#' @param score A numerical value to use as alignment score threshold.
#' @return A data frame of alignment result for the selected pairs of aligned epitopes
#' @examples
#' \dontrun{
#' get_best_alignment(alnRes, query.set, ref.set, best=TRUE, score=NULL)
#' }
#' @keywords internal
#' @export
get_best_alignment <- function(alnRes, query.set, ref.set, best=TRUE, score=NULL){
        if(best==TRUE){
                bestAln <- which(alnRes==max(alnRes), arr.ind = TRUE)
        }else if(is.null(score)){
                bestAln <- which(alnRes==max(alnRes), arr.ind = TRUE)
        }else{
               bestAln <- which(alnRes>=score, arr.ind = TRUE) 
        }

        if(length(bestAln)==0){
                return(NULL)
        }

	bestAlnCount <- nrow(bestAln)
	resDF <- data.frame(query=character(), q_start=numeric(), q_end=numeric(), q_epitope_rank=numeric(), q_epitope_ic50=numeric(), q_epitope_seq=character(), ref=character(), r_start=numeric(), r_end=numeric(), r_epitope_rank=numeric(), r_epitope_ic50=numeric(), r_epitope_seq=character(), score=numeric(), stringsAsFactors=F)
	for(i in 1:bestAlnCount){
		bestAlnRowName <- rownames(alnRes)[bestAln[i, "row"]]
		bestAlnColName <- colnames(alnRes)[bestAln[i, "col"]]
		bestAlnScore <- as.numeric(alnRes[bestAln[i,1],bestAln[i,2]])
		qSeq <- paste0(query.set[bestAlnColName,], collapse="")
		rSeq <- paste0(ref.set[bestAlnRowName,], collapse="")

		bestAlnRowName.split <- strsplit(bestAlnRowName, "#")[[1]]
		refName <- bestAlnRowName.split[1]
		rStart <- as.numeric(bestAlnRowName.split[2])
		rEnd <- as.numeric(bestAlnRowName.split[3])
		rRank <- as.numeric(bestAlnRowName.split[4])
		rIC50 <- as.numeric(bestAlnRowName.split[5])

		bestAlnColName.split <- strsplit(bestAlnColName, "#")[[1]]
		queryName <- bestAlnColName.split[1]
		qStart <- as.numeric(bestAlnColName.split[2])
		qEnd <- as.numeric(bestAlnColName.split[3])
		qRank <- as.numeric(bestAlnColName.split[4])
		qIC50 <- as.numeric(bestAlnColName.split[5])

		tmpDF <- data.frame(query=queryName, q_start=qStart, q_end=qEnd, q_epitope_rank=qRank, q_epitope_ic50=qIC50, q_epitope_seq=qSeq, ref=refName, r_start=rStart, r_end=rEnd, r_epitope_rank=rRank, r_epitope_ic50=rIC50, r_epitope_seq=rSeq, score=bestAlnScore, stringsAsFactors=F)
		resDF <- rbind(resDF, tmpDF)
	}
	return(resDF)
}


align_to_cpl <- function(x, ref){
  xlen <- length(x)
  reflen <- nrow(ref)
  
  cpl <- matrix(nrow=xlen, ncol=20)
  aacols <- Biostrings::AAString("ACDEFGHIKLMNPQRSTVWY")
  
  for (i in 1:xlen) {
    for (j in 1:20) {
      if(x[i] == aacols[j]){
        cpl[i,j] <- 100
      }else{
        cpl[i,j] <- 1
      }
    }
  }
  
  colnames(cpl) <- c("A","C","D","E","F","G","H","I","K","L",
                     "M","N","P","Q","R","S","T","V","W","Y")
  
  score_vec <- c()
  
  cpl <- log(cpl / rowSums(cpl))
  
  for (i in 1:reflen) {
    peptide <- ref[i,]
    ref_name <- rownames(ref)[i]
    score <- 0
    
    for (j in 1:20) {
      if(peptide[i] == aacols[j]){
        score <- score + cpl[i,j]
      }
    }
    score_vec[i] <- score
    names(score_vec)[i] <- ref_name
  }
  return(score_vec)
}


align_cpl_sets <- function(query.set, ref.set, env=NULL, num.cores=2, aln_matrix="BLOSUM62"){
  print("Create Cluster...")
  cl <- parallel::makeCluster(num.cores, outfile="")
  print("Register Cluster...")
  doParallel::registerDoParallel(cl)
  print("Parse Res Set Name...")
  ref.set.name <- deparse(substitute(ref.set))
  ref.set.name <- strsplit(ref.set.name, "\\[")[[1]][1]
  print("Export Functions and Ref Set Name...")
  if(is.null(env)){
    env <- environment()
  }
  #parallel::clusterExport(cl, list("get_pos_weights", "align_to_ref_epitopes", ref.set.name))
  parallel::clusterExport(cl, list("align_to_cpl", ref.set.name), envir=env)
  
  print("par Apply Alignment...")
  t1 <- proc.time()
  cplRes <- parallel::parApply(cl, query.set, 1, function(x) align_to_cpl(x, ref.set))
  print(paste0("Alignment Time : ", (proc.time()-t1)[3]))
  
  print("Kill Cluster...")
  parallel::stopCluster(cl)
  return(cplRes)
}

