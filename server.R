#library(epitomeSim)
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(DT)
library(Biostrings)
library(httr)

#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)

source("R/EPItOMe.R")

options(shiny.maxRequestSize=500*1024^2)
utils::data(list=c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"), envir=environment())

shinyServer(
        function(input, output, session){
                myValues <- shiny::reactiveValues(queryAASeq=NULL, refAASeq=NULL, queryEpitopes=NULL, refEpitopes=NULL, cores=2, resultDF=NULL, aln_matrix=NULL)
                
                myValues$queryAASeq <- shiny::reactive({
			if(input$qInputRadio=="file"){
				qFile <- input$qFile
				if (is.null(qFile))
				return(NULL)

				print("Loading Query Protein Sequence File...")
				qAASet <- readAAStringSet(qFile$datapath, format="fasta")
				return(qAASet)
			}else if(input$qInputRadio=="seq"){
				qText <- input$qSeqText
                                print(qText)
                                print(str(qText))
                                print(length(qText))
				if (is.null(qText) || qText=="")
				return(NULL)
				
				print("Loading Query Protein Sequence Text...")
				qTextSplit <- strsplit(qText, ">")
				qSeqList <- list()
				for(i in 2:length(qTextSplit[[1]])){
					qTextSplit2 <- strsplit(qTextSplit[[1]][i], "\n")
					seqName <- qTextSplit2[[1]][1]
					seqAA <- paste0(qTextSplit2[[1]][2:length(qTextSplit2[[1]])], collapse="")
					qSeqList[[seqName]] <- seqAA
				}
				qSeqVec <- unlist(qSeqList)
				qAASet <- AAStringSet(qSeqVec)
				return(qAASet)
			}
                })

                myValues$refAASeq <- shiny::reactive({
			if(input$rInputRadio=="file"){
				rFile <- input$rFile
				if (is.null(rFile))
				return(NULL)

				print("Loading Reference Protein Sequence...")
				rAASet <- readAAStringSet(rFile$datapath, format="fasta")
				return(rAASet)
			}else if(input$rInputRadio=="seq"){
				rText <- input$rSeqText
				if (is.null(rText) || rText=="")
				return(NULL)
				
				print("Loading Reference Protein Sequence Text...")
				rTextSplit <- strsplit(rText, ">")
				rSeqList <- list()
				for(i in 2:length(rTextSplit[[1]])){
					rTextSplit2 <- strsplit(rTextSplit[[1]][i], "\n")
					seqName <- rTextSplit2[[1]][1]
					seqAA <- paste0(rTextSplit2[[1]][2:length(rTextSplit2[[1]])], collapse="")
					rSeqList[[seqName]] <- seqAA
				}
				rSeqVec <- unlist(rSeqList)
				rAASet <- AAStringSet(rSeqVec)
				return(rAASet)
			}
                })

                myValues$cores <- shiny::reactive({
                        if(is.null(input$coresSelect))
                        return(cores<-2)

                        input$coresSelect
                })

                myValues$aln_matrix <- shiny::reactive({
                        data(input$matrixSelect)
                        get(input$matrixSelect)
                })

                observeEvent(input$startButton, {
                        print("Starting Analysis Process...")

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())

			queryAASeqNameList <- names(myValues$queryAASeq())
			refAASeqNameList <- names(myValues$refAASeq())

                        progress$set(message="Predict Epitopes", value=0)

                        updateProgress(value=1/4, detail="Query Sequence(s)")
			print("Predicting Epitopes from Query Protein Sequence...")
			qEpitopesList <- list()
			for(i in 1:length(myValues$queryAASeq())){
				queryAASeq.name <- strsplit(names(myValues$queryAASeq()[i]), " ")[[1]][1]
				query_epitopes <- predict_epitopes(sequence.name=queryAASeq.name, sequence=as.character(myValues$queryAASeq()[[i]]), method=input$methodSelect, allele=input$alleleSelect, length=input$lengthSelect)
                                query_epitopes_df <- fetch_epitopes(sequence=as.character(myValues$queryAASeq()[[i]]), method=input$methodSelect, allele=input$alleleSelect, length=input$lengthSelect)
                                if(input$ic50FilterRadio==TRUE){
                                        query_epitopes_df_filtered <- filter_epitopes(epitope.DF=query_epitopes_df, ic50.threshold=input$ic50Input)
                                        if(!is.null(query_epitopes_df_filtered))
                                        {
                                                query_epitopes <- epitopeDF_to_AAmat(sequence.name=queryAASeq.name, epitope.DF=query_epitopes_df_filtered)
                                                qEpitopesList[[queryAASeq.name]] <- query_epitopes
                                        }
                                }
                                else{
                                        query_epitopes <- epitopeDF_to_AAmat(sequence.name=queryAASeq.name, epitope.DF=query_epitopes_df)
                                        qEpitopesList[[queryAASeq.name]] <- query_epitopes
                                }
			}

                        updateProgress(value=2/4, detail="Refrence Sequence(s)")
			print("Predicting Epitopes from Reference Protein Sequence...")
			rEpitopesList <- list()
			for(i in 1:length(myValues$refAASeq())){
				refAASeq.name <- strsplit(names(myValues$refAASeq()[i]), " ")[[1]][1]
				#ref_epitopes <- predict_epitopes(sequence.name=refAASeq.name, sequence=as.character(myValues$refAASeq()[[i]]), method=input$methodSelect, allele=input$alleleSelect, length=input$lengthSelect, filter=input$ic50FilterRadio, ic50.threshold=input$ic50Input)
				ref_epitopes_df <- fetch_epitopes(sequence=as.character(myValues$refAASeq()[[i]]), method=input$methodSelect, allele=input$alleleSelect, length=input$lengthSelect)
                                if(input$ic50FilterRadio==TRUE){
                                        ref_epitopes_df_filtered <- filter_epitopes(epitope.DF=ref_epitopes_df, ic50.threshold=input$ic50Input)
                                        if(!is.null(ref_epitopes_df_filtered))
                                        {
                                                ref_epitopes <- epitopeDF_to_AAmat(sequence.name=refAASeq.name, epitope.DF=ref_epitopes_df_filtered)
                                                rEpitopesList[[refAASeq.name]] <- ref_epitopes
                                        }
                                }
                                else{
                                        ref_epitopes <- epitopeDF_to_AAmat(sequence.name=refAASeq.name, epitope.DF=ref_epitopes_df)
                                        rEpitopesList[[refAASeq.name]] <- ref_epitopes
                                }
			}

                        if(length(rEpitopesList)==0 || length(qEpitopesList)==0){
                                break
                        }

                        progress$set(message="Align Epitopes", value=3/4, detail="")
			#allBestAlnDF <- data.frame(query=character(), ref=character(), score=integer(), stringsAsFactors=F)
			allBestAlnDF <- data.frame(query=character(), q_start=numeric(), q_end=numeric(), q_epitope_rank=numeric(), q_epitope_ic50=numeric(), q_epitope_seq=character(), ref=character(), r_start=numeric(), r_end=numeric(), r_epitope_rank=numeric(), r_epitope_ic50=numeric(), r_epitope_seq=character(), score=numeric(), stringsAsFactors=F)
			for(i in 1:length(qEpitopesList)){
			    	for(j in 1:length(rEpitopesList)){
                                        updateProgress(value=3/4, detail=paste0("Ref[", j, "] v/s Query[", i, "]"))
					print("Aligning Epitope Sets...")
					alnRes <- align_sets(query.set=qEpitopesList[[i]], ref.set=rEpitopesList[[j]], num.cores=myValues$cores(), env=environment(), aln_matrix=input$matrixSelect)
                                        updateProgress(value=3/4, detail="Get Best Match")

					print("Getting Best Alignment...")
                                        if(input$alignmentFilterRadio=="best"){
                                                bestAlnDF <- get_best_alignment(alnRes, query.set=qEpitopesList[[i]], ref.set=rEpitopesList[[j]])
                                        }else if(input$alignmentFilterRadio=="score"){
                                                bestAlnDF <- get_best_alignment(alnRes, query.set=qEpitopesList[[i]], ref.set=rEpitopesList[[j]], best=F, score=input$alingmentScoreInput)
                                        }

                                        if(is.null(bestAlnDF)){
                                                next
                                        }

					allBestAlnDF <- rbind(allBestAlnDF, bestAlnDF)
				}
			}
                        myValues$resultDF <- allBestAlnDF
                        myValues$queryEpitopes <- qEpitopesList
                        myValues$refEpitopes <- rEpitopesList
                        print(names(myValues$queryEpitopes))
                        print(names(myValues$refEpitopes))

                        #queryAASeq.name <- strsplit(names(myValues$queryAASeq()), " ")[[1]][1]
                        #refAASeq.name <- strsplit(names(myValues$refAASeq()), " ")[[1]][1]

                        #print("Predicting Epitopes from Query Protein Sequence...")
                        #query_epitopes <- predict_epitopes(sequence.name=queryAASeq.name, sequence=as.character(myValues$queryAASeq()[[1]]), method=input$methodSelect, allele=input$alleleSelect, length=input$lengthSelect)
                        #print("Predicting Epitopes from Reference Protein Sequence...")
                        #ref_epitopes <- predict_epitopes(sequence.name=refAASeq.name, sequence=as.character(myValues$refAASeq()[[1]]), method=input$methodSelect, allele=input$alleleSelect, length=input$lengthSelect)

                        #print("Aligning Epitope Sets...")
                        #alnRes <- align_sets(query.set=query_epitopes, ref.set=ref_epitopes, num.cores=myValues$cores(), env=environment())
                        #print("Getting Best Alignment...")
                        #bestAlnDF <- get_best_alignment(alnRes)

                        #myValues$resultDF <- bestAlnDF
                        #myValues$queryEpitopes <- query_epitopes
                        #myValues$refEpitopes <- ref_epitopes
                        progress$set(message="Finished!", value=4/4, detail="")
                })

                output$alnResDT <- DT::renderDataTable({
                        shiny::validate(
                                need(!is.null(myValues$resultDF), "Waiting for alignment results...")
                        )

                        print("In results renderDataTable!")
			#qIDs <- myValues$resultDF[,1]
			#rIDs <- myValues$resultDF[,2]
			#qSeqs <- myValues$queryEpitopes[rownames(myValues$queryEpitopes)==qIDs,]
			#rSeqs <- myValues$refEpitopes[rownames(myValues$refEpitopes)==rIDs,]
			resDF <- myValues$resultDF
                        print(class(resDF$q_start))
                        print(class(resDF$r_start))
                        print(class(resDF$score))

                        #DT::datatable(resDF, filter = list(position='top', clear=FALSE), options = list(search = list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))
                        DT::datatable(resDF, filter = list(position='top', clear=FALSE), options = list(search = list(regex=TRUE, caseInsensitive=FALSE)))
                        #return(resDF)
                })
                #}, filter = list(position='top', clear=FALSE), options = list(search = list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))

                output$alnViewText <- shiny::renderUI({
                        alnMatrix <- myValues$aln_matrix()
                        alnResRow <- input$alnResDT_rows_selected
                        #alnResRow <- input$alnResDT_row_last_clicked

			alnMultiView <- ""
                        if(length(alnResRow)>0){
				print(str(alnResRow))

				for(i in 1:length(alnResRow)){
					qID <- myValues$resultDF[alnResRow[i],1]
					rID <- myValues$resultDF[alnResRow[i],7]
					score <- myValues$resultDF[alnResRow[i],13]
					
					#print(qID)
					#print(rID)
					#print(qEpitopeID)
					#print(rEpitopeID)
					#print(score)
                                        #print(which(rownames(myValues$queryEpitopes[[qID]])==qEpitopeID))
                                        #print(which(rownames(myValues$refEpitopes[[rID]])==rEpitopeID))

					#qEpitopeID <- paste0(myValues$resultDF[alnResRow[i],1:5], collapse="#")
					#rEpitopeID <- paste0(myValues$resultDF[alnResRow[i],7:11], collapse="#")
					#qSeq <- myValues$queryEpitopes[[qID]][rownames(myValues$queryEpitopes[[qID]])==qEpitopeID,]
					#rSeq <- myValues$refEpitopes[[rID]][rownames(myValues$refEpitopes[[rID]])==rEpitopeID,]
					qSeq <- strsplit(myValues$resultDF[alnResRow[i],6], "")[[1]]
					rSeq <- strsplit(myValues$resultDF[alnResRow[i],12], "")[[1]]
					print(qSeq)
					print(rSeq)

					colorVec <- rep("#0000ff", length(qSeq))
					misMatch <- which(!qSeq == rSeq)
                                        if(length(misMatch)>0){
					        colorVec[misMatch] <- "#ff0000"

                                                misMatchMat <- alnMatrix[rSeq[misMatch], qSeq[misMatch]]
                                                if(length(misMatch)>1){
                                                        misMatchMat <- diag(misMatchMat)
                                                        homologs <- which(misMatchMat>0)
                                                        if(length(homologs)>0){
                                                                colorVec[misMatch[homologs]] <- "#006400"
                                                        }
                                                }else{
                                                    if(misMatchMat>0){
                                                            colorVec[misMatch] <- "#006400"
                                                    }
                                                }
                                        }

					#alnView <- "<table style='width:100%'>
					#	<tr style='background-color:beige;'>
					#	    <td><font color='#ff0000'>A</font></td>
					#	    <td><font color='#0000ff'>B</font></td>
					#	    <td><font color='#0000ff'>C</font></td>
					#	    <td><font color='#0000ff'>D</font></td>
					#	    <td><font color='#0000ff'>E</font></td>
					#	    <td><font color='#ff0000'>F</font></td>
					#	    <td><font color='#0000ff'>G</font></td>
					#	    <td><font color='#0000ff'>H</font></td>
					#	    <td><font color='#ff0000'>I</font></td>
					#	</tr>
					#	<tr style='background-color:beige;'>
					#	    <td><font color='#ff0000'>I</font></td>
					#	    <td><font color='#0000ff'>B</font></td>
					#	    <td><font color='#0000ff'>C</font></td>
					#	    <td><font color='#0000ff'>D</font></td>
					#	    <td><font color='#0000ff'>E</font></td>
					#	    <td><font color='#ff0000'>A</font></td>
					#	    <td><font color='#0000ff'>G</font></td>
					#	    <td><font color='#0000ff'>H</font></td>
					#	    <td><font color='#ff0000'>F</font></td>
					#	</tr>
					#</table>"
					tableO <- "<table style='width:100%'>"
					tableC <- "</table>"
					trO <- "<tr style='background-color:beige;'>"
					trC <- "</tr>"
					tableRowR <- paste0("<td><font color='", colorVec, "'>", rSeq, "</font></td>", collapse=" ")
					tableRowQ <- paste0("<td><font color='", colorVec, "'>", qSeq, "</font></td>", collapse=" ")
					alnView <- paste0(tableO,
						trO,
						tableRowR,
						trC,
						trO,
						tableRowQ,
						trC,
						tableC
					)
					#alnView <- paste0(rID, "<BR>", alnView, "<BR>", qID)
					alnView <- paste0(rID, alnView, qID, "<HR>")
					alnMultiView <- paste0(alnMultiView, alnView)
				}
                                #HTML(alnView)
                                HTML(alnMultiView)
                        }
                })

		output$exportEpitopes <- shiny::downloadHandler(
			filename = function(){
				paste(input$exportEpitopesSelect,"_All_Epitopes_", Sys.Date(), '.fasta', sep='')
			},
			content = function(con){
				idx <- which(names(myValues$queryEpitopes)==input$exportEpitopesSelect)
				if(length(idx)>0){
					epitopesMat <- myValues$queryEpitopes[[input$exportEpitopesSelect]]
				}else{
					epitopesMat <- myValues$refEpitopes[[input$exportEpitopesSelect]]
				}
				epitopesAASet <- Biostrings::AAStringSet(apply(epitopesMat, 1, function(x) paste0(x, collapse="")))
				str(epitopesAASet)
				Biostrings::writeXStringSet(epitopesAASet, filepath=con, format="fasta")
			}
		)

		output$exportAlignmentTable <- shiny::downloadHandler(
			filename = function(){
				paste("Alignment_Table_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				write.table(myValues$resultDF, con, quote=FALSE, row.names=FALSE,  sep="\t")
			}
		)

                output$coresUI <- shiny::renderUI({
                        shiny::numericInput(inputId="coresSelect", label="Number of Cores", value=2, min=1, max=detectCores(), step=1)
                })

		output$exportEptiopesUI <- shiny::renderUI({
			if(is.null(myValues$queryEpitopes) || is.null(myValues$refEpitopes)){
				seqNames <- c("NONE")
			}else{
				seqNames <- unique(c(names(myValues$queryEpitopes), names(myValues$refEpitopes)))
			}
			shiny::selectInput(inputId="exportEpitopesSelect", label="Select Input Sequence", choices=seqNames, selected=seqNames[1])
		})

                output$alleleSelectUI <- shiny::renderUI({
			allele_vec <- fetch_alleles(species=input$organismSelect)
			shiny::selectInput(inputId="alleleSelect", label="Select Allele for Prediction", choices=allele_vec, selected=allele_vec[1])
		})

		shinyjs::onclick("showExportLink", shinyjs::toggle(id="exportDiv", anim=TRUE))

		shiny::observe({
			if(is.null(myValues$queryEpitopes) || is.null(myValues$refEpitopes)){
				shinyjs::disable("exportEpitopes")
				shinyjs::disable("exportAlignmentTable")
			}else{
				shinyjs::enable("exportEpitopes")
				shinyjs::enable("exportAlignmentTable")
			}
                        
                        if(is.null(myValues$queryAASeq()) || is.null(myValues$refAASeq())){
				shinyjs::disable("startButton")
			}else{
				shinyjs::enable("startButton")
			}
		})

                observeEvent(input$rInputRadio, {
                        print("Observing Reference Input Type Event!")
                        if(input$rInputRadio=="seq"){
                                shinyjs::hide(id="rFileDiv", anim=TRUE)
                                shinyjs::hide(id="rDatDiv", anim=TRUE)
                                shinyjs::show(id="rSeqDiv", anim=TRUE)
                        }else if(input$rInputRadio=="file"){
                                shinyjs::hide(id="rSeqDiv", anim=TRUE)
                                shinyjs::hide(id="rDatDiv", anim=TRUE)
                                shinyjs::show(id="rFileDiv", anim=TRUE)
                        }else if(input$rInputRadio == "database"){
                                shinyjs::hide(id="rSeqDiv", anim=TRUE)
                                shinyjs::hide(id="rFileDiv", anim=TRUE)
                                shinyjs::show(id="rDatDiv", anim=TRUE)
                        }
                })

		observeEvent(input$qInputRadio, {
                        print("Observing Query Input Type Event!")
                        if(input$qInputRadio=="seq"){
                                shinyjs::hide(id="qFileDiv", anim=TRUE)
                                shinyjs::hide(id="qDatDiv", anim=TRUE)
                                shinyjs::show(id="qSeqDiv", anim=TRUE)
                        }else if(input$qInputRadio=="file"){
                                shinyjs::hide(id="qSeqDiv", anim=TRUE)
                                shinyjs::hide(id="qDatDiv", anim=TRUE)
                                shinyjs::show(id="qFileDiv", anim=TRUE)
                        }else if(input$qInputRadio == "database"){
                          shinyjs::hide(id="qSeqDiv", anim=TRUE)
                          shinyjs::hide(id="qFileDiv", anim=TRUE)
                          shinyjs::show(id="qDatDiv", anim=TRUE)
                        }
                })

                observeEvent(input$ic50FilterRadio, {
                        print("Observing IC50 Filter Input Type Event!")
                        if(input$ic50FilterRadio==TRUE){
				shinyjs::enable("ic50Input")
                        }else if(input$ic50FilterRadio==FALSE){
				shinyjs::disable("ic50Input")
                        }
                })
                
                
                radioButtons(inputId="alignmentFilterRadio", label="Filter Alignment Table", choices=c("Get Best"="best", "By Alignment Score"="score", 
                                                                                                       "By ALS Score" ="cplscore"), selected="best", inline=TRUE)

                observeEvent(input$alignmentFilterRadio, {
                        print("Observing Alignment Filter Input Type Event!")
                        if(input$alignmentFilterRadio=="score"){
        shinyjs::show(id="filtDivAln", anim=TRUE)
        shinyjs::hide(id="filtDivCPL", anim=TRUE)
				shinyjs::enable("alingmentScoreInput")
                        }else if(input$alignmentFilterRadio=="best"){
        shinyjs::show(id="filtDivAln", anim=TRUE)
        shinyjs::hide(id="filtDivCPL", anim=TRUE)
				shinyjs::disable("alingmentScoreInput")
                        }else if(input$alignmentFilterRadio=="cplscore"){
        shinyjs::show(id="filtDivCPL", anim=TRUE)
        shinyjs::hide(id="filtDivAln", anim=TRUE)
        shinyjs::disable("alingmentScoreInput")
                        }
                })
                #shinyBS::addPopover(session, "rankDT_info", "Rank Table Info", content=paste0("<p>Table containing the final rank of the gene from the whole network inferred by INfORM, ",
		#		"HGNC Gene Symbol of the gene, unique identifier from the NCBI ENTREZ Gene database, descriptive name of the gene, user provided differential expression score, ",
		#		"membership of gene in top5, top10 and top20 candidate genes, membership of gene in response modules extracted by using the top5, top10 and top20 candidate genes. </p>",
		#		"<p>Membership information is cumulative, genes of compact membership also belong to the broader membership, eg: gene membership 5+10+20 = total genes in 20. </p>"
		#	),
		#	placement="right",
		#	trigger="focus",
		#	options=list(container="#rt_info")
		#)
        }
)

