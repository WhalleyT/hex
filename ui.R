library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(DT)

dashboardPage(
        dashboardHeader(title="HEX", titleWidth="25%"),
        dashboardSidebar(disable=TRUE),
        dashboardBody(
                useShinyjs(),
                #extendShinyjs(text = jsCode),
                fluidRow(
                        box(
                                title="Upload", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width=12,
                                fluidRow(
                                        column(6,
						wellPanel(
							h3("Reference Protein Sequence"),
							div(id="rFileDiv",
								#fileInput(inputId="rFile", label="Reference Protein Sequence")
								fileInput(inputId="rFile", label="")
							),
							hidden(div(id="rSeqDiv",
								#textAreaInput(inputId="rSeqText", label="Reference Protein Sequence", value="", rows=5)
								tags$textarea(id="rSeqText", placeholder="Paste Protein Sequence", rows=5, cols=120)
							)),
							hidden(div(id="rDatDiv",
							  selectInput(inputId="rDatabase", label="Database", choices=c("Viral"="viral", "Fungal"="fungal", "Category I Bacteria"="catI", 
							                                                               "Category II Bacteria"="catII", "Cancer"="cancer", "Human"="human")))),
							radioButtons(inputId="rInputRadio", label="", choices=c("File"="file", "Sequence"="seq", "Pre-compiled database"="database"), selected="file", inline=TRUE)
						)
                                        ),
                                        column(6,
						wellPanel(
							h3("Query Protein Sequence"),
							div(id="qFileDiv",
								#fileInput(inputId="qFile", label="Query Protein Sequence")
								fileInput(inputId="qFile", label="")
							),
							hidden(div(id="qSeqDiv",
								#textAreaInput(inputId="qSeqText", label="Query Protein Sequence", value="", rows=5)
								tags$textarea(id="qSeqText", placeholder="Paste Protein Sequence", rows=5, cols=120)
							)),
							hidden(div(id="qDatDiv",
							           selectInput(inputId="qDatabase", label="Database", choices=c("Viral"="viral", "Fungal"="fungal", "Category I Bacteria"="catI", 
							                                                                        "Category II Bacteria"="catII", "Cancer"="cancer", "Human"="human")))),
							radioButtons(inputId="qInputRadio", label="", choices=c("File"="file", "Sequence"="seq", "Pre-compiled database"="database"), selected="file", inline=TRUE)
						)
                                        )
                                ),
                                fluidRow(
                                        column(6,
                                                selectInput(inputId="methodSelect", 
                                                        label="Select Method for Prediction", 
                                                        choices=c("recommended", "consensus", "netmhcpan", "ann", "smmpmbec", "smm", "comblib_sidney2008", "netmhccons", "pickpocket"),
                                                        selected="recommended"
                                                ),
                                                selectInput(inputId="organismSelect", 
                                                        label="Select Organism", 
                                                        choices=c("chimpanzee", "cow", "human", "macaque", "mouse", "pig"),
                                                        selected="mouse"
                                                ),
                                                selectInput(inputId="matrixSelect", 
                                                        label="Select Alignment Matrix", 
                                                        choices=c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"),
                                                        selected="BLOSUM62"
                                                )
                                        ),
                                        column(6,
                                                selectInput(inputId="lengthSelect", 
                                                        label="Select Epitope Length for Prediction", 
                                                        choices=c(8:14),
                                                        selected=9
                                                ), 
                                                uiOutput("alleleSelectUI"),
                                                uiOutput("coresUI") 
                                        )
                                ),
                                h3("Filtering Options"),
                                fluidRow(
                                        column(6,
                                                radioButtons(inputId="ic50FilterRadio", label="Filter by IC50", choices=c("YES"=TRUE, "NO"=FALSE), selected=FALSE, inline=TRUE),
                                                numericInput(inputId="ic50Input", label="IC50 Filtering Threshold", value=100, min=0, step=1)
                                        ),
                                        
                                        column(6,
                                                radioButtons(inputId="alignmentFilterRadio", label="Filter Alignment Table", choices=c("Get Best"="best", "By Alignment Score"="score",
                                                                                                                                       "By ALS Score" ="cplscore"), selected="best", inline=TRUE),
                                               hidden(div(id="filtDivAln",
                                                numericInput(inputId="alingmentScoreInput", label="Alignment Score Threshold", value=40, min=0, step=1))),
                                               hidden(div(id="filtDivCPL",
                                                          numericInput(inputId="ALSScoreInput", label="ALS Score Threshold", value=-30, min=0, step=1)))
                                        )
                                ),
                                hr(),
                                actionButton(inputId="startButton", label="Start Analysis")
                        )
                ),
                fluidRow(
                        tabBox(
                                id="viewTabBox", title="Analysis Results", width=12,
                                tabPanel(value="tablePanel", title="Alignment",
					  fluidRow(
						column(12,
							a(id = "showExportLink", "Show/Hide Export Options"),
							hidden(div(id="exportDiv",
								fluidRow(
									column(4,
										wellPanel(
										uiOutput("exportEptiopesUI"),
										downloadButton('exportEpitopes', label='Export Epitopes')
										)
									),
									column(4,
										wellPanel(
										downloadButton('exportAlignmentTable', label='Export Alignment Table')
										)
									)
								)
							))
						)
					 ),
                                         fluidRow(
                                                column(10,
							h3("Alignment Table"), 
                                                    	#div(style = 'overflow-x: scroll', DT::dataTableOutput('alnResDT'))
                                                    	DT::dataTableOutput('alnResDT')
                                                ),
                                                column(2,
							h3("Alignment View"), 
                                                    	uiOutput("alnViewText") 
                                                )
                                         )
                                )
                        )
                )
        )
)

