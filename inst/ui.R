library(shinythemes)
library(shinyBS)
library(threejs)
library(leaflet)

# ------------------------------------------ Custom FileInput function -------------------
# # based on the Shiny fileInput function
# fileInput2 <- function(inputId, label = NULL, labelIcon = NULL, multiple = FALSE, 
#                        accept = NULL, width = NULL, progress = TRUE,...) {
#   
#   # add class fileinput_2 defined in UI to hide the inputTag
#   inputTag <- tags$input(id = inputId, name = inputId, type = "file", 
#                          class = "fileinput_2")
#   if (multiple) 
#     inputTag$attribs$multiple <- "multiple"
#   if (length(accept) > 0) 
#     inputTag$attribs$accept <- paste(accept, collapse = ",")
#   
#   div(..., style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";"), 
#       inputTag,
#       # label customized with an action button
#       tags$label(`for` = inputId, div(icon(labelIcon), label, 
#                                       class = "btn btn-default action-button")),
#       # optionally display a progress bar
#       if(progress)
#         tags$div(id = paste(inputId, "_progress", sep = ""), 
#                  class = "progress shiny-file-input-progress", 
#                  tags$div(class = "progress-bar")
#         )
#   )
# }          


navbarPage(title="Spatial Transcriptomics Downstream App", id="navbar",
           theme = "bootstrap.min.css",

# ------------------------------------------- Data-processing -----------------------
tabPanel("Data Processing",
         
         # # define class fileinput_2 to hide inputTag in fileInput2
         # tags$head(tags$style(HTML(
         #   ".fileinput_2 {
         #   width: 0.1px;
         #   height: 0.1px;
         #   opacity: 0;
         #   overflow: hidden;
         #   position: absolute;
         #   z-index: -1;
         #   }"
         #              ))),
         tags$head(tags$style(HTML('

                                   .modal-lg {
                                   height: 1200px;
                                   width:  1200px;
                                   
                                   }
                                   '))),
       br(),
       br(),
        fluidPage(
          fluidRow(
            column(4),
            column(4,
           actionButton("DataButton", label=img(src="data_icon.png",
                                         height = 200, width = 200)), 
                      align="center"
            ),
           column(4)
           ),
          br(),
          br(),
          
          
         fluidRow(
           column(4,
                  actionButton("FilterButton", label=img(src="filter_icon_with_name.png",
                                                         height = 200, width = 200)), 
                  align="center"
           ),
           column(4,
                  actionButton("ClusterButton", label=img(src="cluster_icon_with_circle_and_name.png",
                                                         height = 200, width = 200)),
                               align="center"
           ),
           column(4,
                  actionButton("NormalizationButton", label=img(src="normalization_icon_circle_and_name.png",
                                                         height = 200, width = 200)),
                               align="center"
           )),
         br(),
         br(),
         fluidRow(
           column(4,
                  actionButton("QCButton", label=img(src="QC_icon.png",
                                                     height = 200, width = 200)),
                  align="center"),
           column(4,
                  actionButton("DEButton", label=img(src="DE_icon.png",
                                                     height = 200, width = 200)),
                  align="center"),
           
           column(4,
                  actionButton("HVGButton", label=img(src="HVG_icon.png",
                                                     height = 200, width = 200)),
                  align="center")
           ),
       br(),
       br(),
       
       fluidRow(
         column(4, 
                actionButton("DIMRedButton", label=img(src="DIM_Reduction_icon.png",
                                                      height = 200, width = 200)),
                align="center"),
         column(4,
                actionButton("HEATMAPButton", label=img(src="heatmap_icon.png",
                                                       height = 200, width = 200)),
                align="center"),
         
         column(4,
                actionButton("batchCorrButton", label=img(src="batchCorr_icon.png",
                                                        height = 200, width = 200)),
                align="center")
       )),
          br(),
          br(),
           
          # Datatab -----
           # popups for all buttons on first page
           bsModal("bsModal_data", title="Data upload", trigger="DataButton",
                   size = "large",
                   sidebarPanel(
                     p("Use the reset button for each new sets of data you work with"),
                     p("Upload ST-output files after running the RNA-seq pipeline"),
                     p("The user can upload multiple data sets which get stored as seperate R-objects
                        (they can later me pooled an manipulated togheter)"),
                     numericInput("nfiles", "number of files", value = 1, min = 1, step = 1),
                     uiOutput("countInputs"),
                     uiOutput("spotInputs"),
                     actionButton("DataCreationButton", "Create Data Sets")

                   ),
                   sidebarPanel(
                     h4("Add annotation to your data set(s)"),
                     p("Add a factor variable to your data set(s), e.g. \"patient\" if you want to control for
                        the patient effect during DE analysis"),
                     textInput("NameDataSetFactorInput", "Factor name:", placeholder="Patient1"),
                     uiOutput("NameDataSetFactorDataSetUI"),
                     actionButton("AddDataSetFactorButton", "Add name"),
                     textOutput("DataSetFactorNameAdded")
                   ),
                    sidebarPanel(
                              htmlOutput("uploadedData")
                    ),
                   sidebarPanel(
                              actionButton("resetDatasets", "Reset Data Set(s)", icon=icon("warning", lib="font-awesome")),
                              htmlOutput("resetDatasetsMsg")
                   )),
          # Filter -----
           bsModal("bsModal_filter", title="Filter options", trigger="FilterButton", 
                   size = "large",
                   fluidRow(
                     column(6,
                            h3("Data Status"),
                            uiOutput("CurrentDataStatus"),
                            br(),
                            h4("Press the update filters button after choosing Data set and/or after you have applied a filter"),
                            actionButton("UpdateFilters", "Update filter status"),
                            br(),
                            h4("Filter away all spots outside of the tissue"),
                            p("In order to use this filter you need to have uploaded spot detection data"),
                            actionButton("insideOutsideTissueFilterButton", "Apply filter"),
                            textOutput("insideOutsideTissueFilterTextOutput"),
                            br(),
                            h4("Filter low transcript abundance genes across all ST-Spots"),
                            p("This filter discard all the genes that have an average count below the given threshold"),
                            numericInput("avgCountFilterInput", "Avg.count:", min=0, value=0.01),
                            actionButton("Apply_avgCount_Filter", "Apply filter"),
                            uiOutput("avgCount_filter_done"),
                            br(),
                            h4("Filter low library size spots"),
                            p("Distribution shown in the graph"),
                            p("Use the slider to choose the threshold level"),
                            uiOutput("librarySizeRemoveText"),
                            actionButton("Apply_librarySize_Filter", "Apply filter"),
                            uiOutput("librarySize_filter_done"),
                            br(),
                            h4("Filter ST-Spots with few identified genes"),
                            p("Use the silder to choose the threshold level"),
                            uiOutput("current_filter_2"),
                            actionButton("Apply_lowgene_spots_filter", "Apply filter"),
                            uiOutput("lowgene_spots_filer_done"),
                            br(),
                            br(),
                            # h4("Filter low abundance genes but keep spatial outliers"),
                            # p("An alternative approach to gene filtering is to keep all genes that have non-zero counts in at least n ST-Spots. 
                            # This provides some more protection against genes with outlier expression patterns, i.e., strong expression in only one or two cells. 
                            # Such outliers are typically uninteresting as they can arise from amplification artifacts that are not replicable across cells. 
                            # (The exception is for studies involving rare cells where the outliers may be biologically relevant.) 
                            # To peform this filtering approach, set the number of n cells below"),
                            #                      
                            #                      sliderInput("n_cells", "Nr of ST-Spots with non.zero counts: ", min=1, max=1000, value=10),
                            #                      textOutput("current_filter2"),
                            #                      actionButton("applyAbundacetypeTwoFilterButton", "Apply filter"),
                            #                      textOutput("abun_filter2"),
                            #                      br(),
                            # 
                  
                           
                            checkboxGroupInput("CommonFilters", "Commonly used quick filters:", choices = c("MALAT1", "Ribosomal genes")),
                            actionButton("Apply_common_filter", "Apply filter"),
                            uiOutput("CommonFiltersOutput"),
                            h4("Filter by gene name"),
                            p("Write the gene names that you want to remove"),
                            textInput("geneID_filter", "GeneIDs: (Obs, no space between IDs)", placeholder = "GeneID1,GeneID2,GendeID3"),
                            actionButton("Apply_geneID_filter", "Apply filter"),
                            uiOutput("geneID_filter_done"),
                            br(),
                            
                            h4("Filter by gene name [Regex]"),
                            p("Use regular expression to filter away gene names"),
                            textInput("RegexGeneID_filter", "Regex:", placeholder = "e.g. (RPS)|(RPL)|(RRP)"),
                            actionButton("Apply_regex_filter", "Apply filter"),
                            uiOutput("regexfilterGenes"),
                            uiOutput("regex_filter_done"),
                            br()
          
                            ),
                            
                     column(6,
                   plotOutput("filterPlot1"),
                   numericInput("filterPlot1_slider", "Library size per ST-Spot", min=0.00001, max=NA, 
                               value=5),
                   plotOutput("filterPlot2"),
                   numericInput("filterPlot2_slider", "Nr of expressed genes / spot", min=1, max=NA, 
                               value=500)
                     )
                   )),
           # ClusterTab -----    
           bsModal("bsModal_cluster", title="Clustering options", trigger="ClusterButton", 
                   size = "large",
                   fluidRow(
                     column(6,
                            h3("Info"),
                            p("Cluster your data prior to normalization of other downstream processess"),
                            p("The output of the clustering are saved in the data set(s)"),
                            p("First, choose your clustering approach and the data set(s) you want to include in the
                              clustering approach. Afterwards you can view the spatial positions of the clusters."),
                            p("OBS. Note that if you merge data sets, genes that are abscent in any of the data sets
                              will be lost for all data sets"),
                            h3("Data set(s) config"),
                            uiOutput("clusterDataSetChooser"),
                            h3("Cluster approach"),
                            selectInput("clusterApproach", "", choices=c("Spearman correlation", "k-means")),
                            p("OBS, k-means not supported atm"),
                            p("min.size = An integer scalar specifying the minimum size of each cluster."),
                            p("If the button \"Apply minSize\" is not pushed, no minimum size is set"),
                            numericInput("clusterMinSizeInput", "min.size: " , min=1, max=100, value=40,   width = "80px"),
                            checkboxInput("clusterMinSizeB", "Use minSize"),
                            actionButton("performClusterButton", "Cluster"),
                            uiOutput("clusterTableOutput")
                            
                      
                     ),
                     
                     column(6,
                              selectInput("typeOfClusterPlot", "What type of cluster do you want to plot", choices = "Spearman"),
                              p("Obs, only spearman avaiable atm"),
                              uiOutput("clusterDataSetPlot"),
                              actionButton("clusterArray_plotB", "Show clusters on array"),
                              plotOutput("clusterArray_plot"),
                              br(),
                              uiOutput("clusterMsgOutput")
                     
                     )
                   )),
       # NormalizationTab -----
           bsModal("bsModal_normalization", title="Normalization options", trigger="NormalizationButton", 
                   size = "large",
                   fluidRow(
                   sidebarPanel(
                            h3("Normalization - \"Scran\" approach"),
                            p("The normalization procedure is adopted from Lun et al (2016) designed for scRNA-seq
                              data. Brefly, expression values are summed across pools of cells and the summed values are then used to calculate normalization size-factors for each individual pool.
                              The size factors are then deconvolved into cell-specific factors that are used for normalization."),
                            p("Here, we provide the pooling alternative presented via the \"Scran\" R-package.
                              if multiple data sets are uploaded, the user can merge these togheter before pooling. (OBS, if genes are missing from any data set relative to another, these will be lost during the merge"),
                            
                            br(),
                            h3("Current selected Data Set"),
                            #uiOutput("CurrentSelectedDataSetNorm"),
                            #p("Obs, this parameter is overriden if \"use merge\" is checked"),
                            
                            h3("Merge data sets prior to normalization"),
                            uiOutput("sceNormUI")
                   ),
                   sidebarPanel(
      
                            h3("Compute size factors"),

                            p("Optional:"),
                            p("An optional factor specifying which cells belong to which cluster, for deconvolution within clusters."),
                            checkboxGroupInput("normclusters", "Clusters:", choices="Spearman", selected=NULL),
                            p("sizes = A numeric vector of pool sizes, i.e., number of cells per pool.
                              The larges size needs to be half the number of the smallest pool"),
                            textInput("NormSizesInput", "Input sizes:", value="NA", placeholder = "ex: 5,10,20,30,40"),
                            
                            actionButton("forcePositiveInfo", "", icon=icon("info-sign", lib="glyphicon")),
                            bsModal("forcePositiveInfoModal", "Info - DE force positive option", trigger="forcePositiveInfo",
                                    htmlOutput("forcePositiveInfoText")),
                            checkboxInput("forcePositive", "Force positive values of size factors"),
                            br(),
                            actionButton("sizefactorB", "Compute size factors"),
                            verbatimTextOutput("sizefactorSummary"),
                            uiOutput("NormalizeMsgOutput")
                            
                   )),
                   fluidRow(
                     sidebarPanel(
                   
                            actionButton("sizeFactorPlotB", "Plot Size Factors"),
                            p("The size factors are based on the clustering and can result in values = 0. This is problematic and interfer with
                              downstream analysis. Check if your size factors looks ok before proceding:"),
                            actionButton("checkNormB", "Check Size Factors"),
                            textOutput("checkNorm"),
                            br(),
                            actionButton("removeSFB", "Remove all ST-features with size factors < 0"),
                            br(),
                            textOutput("removeSF"),
                            br(),
                            actionButton("normalizeB", "Apply Normalization"),
                            textOutput("normalizationDone")
                            ),
                     
                     mainPanel(
                            plotOutput("sizeFactorPlot")
                            
                              )
                     )),
       # QCTab -----
          bsModal("bsModal_qc", title="Quality Control", trigger="QCButton", 
               size = "large",
               fluidRow(
                 sidebarPanel(
                  h4("Inside/Outside tissue diffusion check"),
                  p("In order to run check, spot detected data needs to be uploaded"),
                  
                  actionButton("insideOutsideButton", "Run check"),
                  br(),
                  uiOutput("qcCheckDataSelectionUI"),
                  uiOutput("qcCheckDataSelctionButton")
                )),
                fluidRow(
                  plotOutput("qcInsideOutsidePlots"),
                  plotOutput("qcInsideOutsidePlots2"),
                  uiOutput("qcDataTableOutputUI")
                
               ),
               fluidRow(
                sidebarPanel(
                  h4("Most abundant genes across all features"),
                  uiOutput("QCDataSetRadioButtons"),
                  br(),
                  actionButton("geneQCB", "Plot")
                ),
                mainPanel(
                  plotOutput("geneQC", height = "800px")
                ))
               ),   
       # DETab -----
       bsModal("bsModal_de", title="Differential Expressen (DE) Analysis", trigger="DEButton", 
               size = "large",
               fluidRow(
               sidebarPanel(
                 h3("1. Data set(s) config"),
                 uiOutput("CurrentDataStatusDE"),
                 checkboxGroupInput("DEDataSelectionsParameters", "Only use spots with the following parameters",
                                    choices = c("Spatial Eye Selection", "Spearman Clusters")), 
                 uiOutput("parameterSpatialEyeDE"),
                 uiOutput("parameterClusterDE"),
                 
                 actionButton("dataConfigDE", "Set Data config"),
                 textOutput("dataConfigDEConfirm")
               ),
              sidebarPanel(
                h3("2. Strategy config"),
                p("If the user has performed prior normalization, check that sizeFactors are not equal to zero."),
                p("If normalization have been done prior to DE, normalized counts will be used, otherwise DESeq normalization
                  of counts will be performed. Note that if you are using the pooling strategy, pooling will be performed on raw counts,
                  and they DESeq normalization will be performed on those pooled raw count values"),
                
                checkboxInput("DEPoolStrategy", "I want to use the Spot-Pooling strategy"),
                actionButton("DEPoolStrategyInfoB", "", icon=icon("info-sign", lib="glyphicon")),
                uiOutput("DEPoolStrategyUiOutput1"),
                uiOutput("DEPoolStrategyUiOutput2"),
                uiOutput("DEPoolStrategyUiOutput3"),
                
                bsModal("DEPoolStrategyInfoBModal", "Info - DE pooling strategy", trigger="DEPoolStrategyInfoB", size="small",
                        htmlOutput("DEPoolStrategyInfoText")),
                p("If the pooling strategy is not used, each ST-Spot is treated as an individual sample"),
                actionButton("SetStrategy", "Set Strategy"),
                textOutput("SetDEStrategyConfirm")
              ),
              sidebarPanel(
                h3("Message output box"),
                uiOutput("DEMessegeBoxLostGenes"),
                p("Created contrasts:"),
                uiOutput("DEContrastMessegeBox"),
                br(),
                p("Design used:"),
                uiOutput("DEDesignMessegeBox")
              )
                
                 
               ),
               fluidRow(
               sidebarPanel(
                 h3("3. Design"),
                 uiOutput("DEDesignOutput"),
                 uiOutput("DEVariableOfInterestUI"),
                 uiOutput("SetRefDE"),
                 uiOutput("DEControlTerm"),
                 actionButton("DESeq2RUN", "Run DESeq2")
                
               ),
               sidebarPanel(
                 
                 h3("4. Show results"),
              
                 uiOutput("DEContrastChooser1"),
                 #p("vs"),
                 #uiOutput("DEContrastChooser2"),
                 actionButton("ResultsDE", "Results")

               # p("Choose the factor that you want to perform DE on"),
               # selectInput("DE_Parameter_selection", "Factor:", choices=c("Spearman Clusters", "Spatial Eye selection")),
               # uiOutput("DE_Parameter_selection_part2"),
               # uiOutput("DE_Parameter_selection_part3"),
               # actionButton("RunEdgeR", "Run EdgeR")
               
               ),
               sidebarPanel(
                 h3("5. Downloads"),
                 downloadButton("downloadDETableButton", "DE table output")
               ),
               mainPanel(

               )),
               fluidRow(
                 uiOutput("DETableOutput")
    
                 
      
               ),
               fluidRow(
                 sidebarPanel(
                   h3("6. Plotting options"),
                   selectInput("DE_plot_alt", "Choose plot type:", choices=c("MA", "Volcano")),
                   checkboxInput("labelDEpointsInput", "Label points based on FDR", value=FALSE),
                   uiOutput("labelDEpointsInputFDRselect"),
                   actionButton("DEPlot_execute", "Plot")
                   
                 ),
                 mainPanel(
                   
                   
                   plotOutput("DEPlotOutput",
                              brush= brushOpts(id="MA_brush", delay=300, delayType="debounce", resetOnNew=TRUE))
                  
                 ),
                 
                
                 #verbatimTextOutput("MA_brush_points")
                 uiOutput("MA_brushinfo")
       
               )
       ),
               
       # HVG tab -----
       bsModal("bsModal_hvg", title="Detection of highly variable genes", trigger="HVGButton", 
               size = "large",
               br(),
               h3("Data Selection"),
               uiOutput("CurrentDataStatusHVG"),
               p("Here count data is used to detect highly variable genes (HVGs)
                 that drives heterogeneity across cells in a population. 
                 Keep in mind that the current implementation of the ST pipeline implies above single-cell resolution.
                 This means that cells within an ST-feature is pooled count-wise, which can mask heterogeniety within the ST-feature itself."),
               br(),
               h3("Detection of HVGs"),
               p("In order to detect the HVG genes, we need to estimate the variance in expression for each ST-feature,
                 however, this could derive from both biological and technical variation (e.g. RNA capture noise)."),
               p("One alternative is to fit a trend for the variance estimates of the endogenous genes (i.e. within ST-features.
                 OBS, this assumes that the majority of the genes are NOT DE, such that the technical component dominates the
                 total variance."),
               p("The fitted value is used as an estimate of the technical component, this value is subtracted from the total variance
                  in order to get the biological component for each gene. In this manner, HVGs are defined as the top hits according to highest biological components."),
               
               p("OBS, the following plot is not possible if you have size factors = 0 in your expression object"),
               selectInput("hvgSelectInput", "Choose what to include in your HVG analysis", choices=c("All ST-spots", "Spatial Eye selection")),
               actionButton("RunHVG", "Set Parameters"),
               textOutput("HVGParametersSet"),
               uiOutput("hvgFollowUpQ"),
               numericInput("hvgFDRInput", "Set FDR level", value=0.05, step=0.01),
               numericInput("hvgBioInput", "Set Bio level", value=0.5, step=0.01),
               p("Choose the threshold level for the biological component"),
               actionButton("trendFitB", "Plot"),
               p("(If you get Error msg, it means no gene matched your choosen threshold"),
               plotOutput("trendFit"),
               p("HVGs are defined as genes with biological components that are significantly greater than zero at a false discovery rate (FDR) of 5%."),
               dataTableOutput("HVG_table"),
               br()
               ),
       # DimReductionTab -----
       bsModal("bsModal_dimred", title="Dimensionality reduction", trigger="DIMRedButton", 
               size = "large",
               sidebarPanel(
                 h3("Data Selection"),
                 uiOutput("CurrentDataStatusDimRed"),
                 uiOutput("parameterSpatialEyeDimRed"),
                 checkboxInput("ScreePlotCheckBox", "Include a Scree plot (this will add time to calculation)"),
                 checkboxInput("ArrayPlotCheckBox", "Include spatial array plot(s)"),
                 uiOutput("DimPlotType"),
                 uiOutput("PlotOpt1"),
                 uiOutput("colorList"),
                 selectInput("dimRedInputType", "Count data",
                             choices=c("Raw counts" = "counts",
                                       "Scran normalized counts" = "norm_exprs",
                                       "CPM normalized counts" = "cpm")),
                 uiOutput("dimRedInputTypeNorm"), #depending on input from dimRedInputType
                 p("OBS, if no \"scran\" normalization is saved for the data set(s), CPM will be used instead"),
                 uiOutput("geneID_color_dimRed"),
                 actionButton("dimPlot", "Plot")
               ),
               mainPanel(
                 plotOutput("DimRedPlot",
                            brush= brushOpts(id="DimRed_brush", delay=300, delayType="debounce", resetOnNew=TRUE)),
                 br(),
                 br(),
                 plotOutput("screePlot"),
                 plotOutput("ArrayPlotBelowDimRed")
               ),
               sidebarPanel(
                 h3("Remove points options"),
                 p("This option is available when you choose a single Spatial Eye selection and want to remove ST-spots from the selection"),
                 uiOutput("DimRedSelect2"),
                 uiOutput("RemoveSelectionDimRed"),
                 textOutput("DimRed_brush_points"),
                 verbatimTextOutput("removedOutputDimRed")
               )
       ),
       # HeatmapTab -----
       bsModal("bsModal_heatmap", title="Heatmaps", trigger="HEATMAPButton", 
               size = "large",
               sidebarPanel(
                 h3("Data Selection"),
                 uiOutput("CurrentDataStatusHeatmap"),
                 h4("Only extract Spatial Eye selection"),
                 
                 checkboxInput("HeatDataSelectionsParameters", "Extract Spatial Eye selections", value=FALSE),
                 uiOutput("HeatParameterSpatialEye"),
                 radioButtons("HeatMapTypeRadioButtons", "Type of heatmap", choices=c("Most Variable Genes")),
                 selectInput("HeatCountTypeInput", "Count input type:", c("Raw counts" = "counts",
                                                                                      "Scran normalized counts" = "norm_exprs")),
                 checkboxGroupInput("HeatAnnotationTypeInput", "Display annotations: ", c("Spatial Eye Selection" = "selectionID",
                                                                                "Spearman Clusters" = "spearmanCluster")),
                 numericInput("HeatNrOfGenes", "Number of genes", min=1, max=1000, value=50, step=1),
                 uiOutput("HeatAnnotation"),
                 sliderInput("HeatFontSize", "Fontsize", min=1, max=10, value=7, step=1),
                 checkboxInput("HeatClusterRows", "Cluster Rows", value=FALSE),
                 checkboxInput("HeatClusterCols", "Cluster Columns", value=FALSE),
                 checkboxInput("HeatRowNames", "Row names", value=FALSE),
                 checkboxInput("HeatColNames", "Col names", value=FALSE),
                 uiOutput("HeatGroupBy"), #not added yet
                 actionButton("HeatMapPlotButton", "Plot Heatmap")
               ),
               mainPanel(
                 plotOutput("topVarHeatmap")
               )
       ),
       # BatchCorrect -----
       bsModal("bsModal_batchCorr", title="Batch Correction", trigger="batchCorrButton",
               size="large",
               sidebarPanel(
                 h3("Perform batch correction between data sets"),
                 p("Batch correction is performed with the removeBatchEffects() function from the limma R-packge,
                   Note, it is intended to be used in conjunction with data exploration like a heatmap or PCA plot. 
                   It is not intended to be used in conjunction with a differential expression analysis. 
                   For a DE analysis, you achive the same effect by adding the batches as an additive factor to the linear model."),
                 p("Note that the batch correction can be performed both on raw or normalized counts"),
                 uiOutput("CurrentDataStatusBatchCorr"),
                 selectInput("BatchCorrCountType", "Use count values:", choices = c("Raw counts"="counts",
                                                                                    "Scran normalized counts" = "norm_exprs")),
                 actionButton("ShowBatchCorr", "Perform batch correction"),
                 br(),
                 br(),
                 uiOutput("ApplyBatchCorrButtonUI"),
                 textOutput("ApplyBatchCorrConfirmText")
               ),
               mainPanel(
                 plotOutput("PCABeforeBatchCorr"),
                 plotOutput("PCAAfterBatchCorr")
               )
       # ),
       # 
       # bsModal("bsModal_save", title="Save RData", trigger="SaveButton", 
       #         size = "large",
       #         sidebarPanel(
       #         p("Save your RData for later use"),
       #         downloadButton("SaveData")
       #         
       #         ),
       #         mainPanel(
       #           uiOutput("SavedDataConfirm")
       #         )
       # )
       )),

# --------------------------------------------------------------- Spatial Eye ---------------------------------------------------------------------
tabPanel("Spatial Eye",
         
                    div(class="outer",
                        tags$style(type = "text/css", ".outer {position: fixed; top: 41px; left: 0; right: 0; bottom: 0; overflow: hidden; padding: 0}"),
                        leafletOutput("Map", width="100%", height="100%"),
                   absolutePanel(bottom=50, right=10, width=150,
                                 selectInput("singleClickMode", "Single click mode", c("Selection", "Spot statistics"), selected="")),
                   absolutePanel(bottom=10, right=5, width=200,
                                 actionButton("checkSelection", "Run QC on selection",  icon = icon("eye", lib = "font-awesome"))),
                   
                   absolutePanel(bottom=50, left=10, fixed=TRUE, width="100",
                                 numericInput("SelectionInput", "Selection", value=1, min=1, step=1)),
                   #Change this to numericInput (OBS need to check all server side usage - if chr needed somewere)
                   #downloadButton("imgSave", "Save img")),
                   absolutePanel(bottom=10, width="80%", fixed=TRUE,
                                 actionButton("assignSelect", "Assign Selection"),
                                 actionButton("resetSelect", "Reset Selection")
                                 ),
                                 #actionButton("SetRefPointButton", "Set ref point")),
                  #uiOutput("SpatialEyeDataSetChooser"),
                  absolutePanel(top=50, left=50, width="20%", fixed=TRUE,
                               numericInput(inputId = "SpatialEyeDataSetChooserInput", "Data Set", value = 1, width="30%")),
                                
                  uiOutput("SpotDataMissing"),
                  
                
                
                   
                   
                   bsModal("Plot_and_table", "Plot and Table", "button_plot_and_table", size = "large",
                           plotOutput("TestPlot"),
                           verbatimTextOutput("spotStats"),
                           dataTableOutput("geneStats")),
                   
                   bsModal("testtest", "just a test", "noTrigger", size="large"),
                           
                  bsModal("multiQC", "Plot and Table Multi", "dummie_name", size = "large",
                            plotOutput("multiQCstandard"),
                           
                          
                          
                          
                                     plotOutput("multiQCPCA",  
                                                brush= brushOpts(id="multiQCPCA_brush", delay=300, delayType="debounce", resetOnNew=TRUE)),
                          textOutput("multiQCPCA_brush_points"),
                          actionButton("removeFromMultiQCPCA", "Remove points from selection",  icon = icon("trash", lib = "font-awesome")),
                          verbatimTextOutput("removedOutput"),
                              
                            verbatimTextOutput("multiQCstats"),
                            dataTableOutput("multiQCgeneStats")
                   ))),
          
# -------------------------------------------- INFO Panel ---------------------------------------
          tabPanel("Info", 
                   h4("In this tool you can upload your output after running the ST-pipeline and spot-detection webtool, perform a sanity check on your
                       data, as well as downstream analysis."),
                   h4("Prerequisites"),
                   p("There are a few steps you need to perform prior to useage of this R-tool"),
                   p("1. The output from the ST-bioinformatic pipeline (https://github.com/SpatialTranscriptomicsResearch/st_pipeline).
                         The output is count-values for the detected genes and is in the form of an .tsv file, with genes as columns and spatial coordinates (ST-Spots), as rows."),
                   p("2. The output from the spot-detection tool (https://github.com/SpatialTranscriptomicsResearch/st_aligner),
                         The output is a matrix where the rows consists of each spatial coordinate with the associated pixel cooridnate 
                         corresponding to the correct position of the histology image. This output is also in the form of an .tsv file"),
                   p("3. Tiled histology images. The HE image taken during an ST experiment needs to be tiled for input
                         into this R-tool. Tiling can be performed with the following tool: https://github.com/ludvb/st_tiler ,
                         the output from the tool will be a folder structure of correct format. Place the whole folder into the subfolder
                         \"www\" located in the directory of this R-tool. Inside \"www\" there should be a folder for each Data Set that
                         you want to look at, these folders are named 1,2,3... etc (for dataset 1,2,3 etc...). Exampel: For data set 1, run
                         the tiler-tool and put the output folders in folder \"1\" located inside the \"www\" folder"),
                   h4("Workflow"),
                   img(src='After_Installation_workflow.svg', align = "center"),
                   h4("Dependencies"),
                   p("This tool is build around the ", strong("SCESet object"), "which is the basic data container class adopted from the",
                     a("scater",     href="https://github.com/davismcc/scater"), "package"),
                   p("This object stores your expression values, cell and feature metadata, as well as any normalization and user selections made using the tools provided here.
                     After creation of the SCE object, the user can download it locally, use it for custom made plotting/analysis purposes, or 
                     re-upload it to this tool at a later stage."),
                   
                   
                   
                   h2("Clustering"),
                   h3("Normalization via clustering and deconvolution"),
                   br(),
                   p("Read counts are subject to differences in capture efficiency between ST-features on the array. 
                     Normalization is required to eliminate these feature-specific biases prior to downstream quantitative analyses. 
                     This is often done by assuming that most genes are not differentially expressed (DE) between ST-features
                     Any systematic difference in count size across the non-DE majority of genes between two ST-features is assumed to represent bias and is removed by scaling. 
                     More specifically, size factors are calculated that represent the extent to which counts should be scaled in each library.
                     The characteristics of ST-data can be problematic for these bulk data-based methods due to the dominance of low and zero counts. 
                     To overcome this, we pool counts from many ST-features to increase the count size for accurate size factor estimation.
                     Pool-based size factors are then deconvolved into feature-based factors for feature-specific normalization.
                     This approach is adopted from the singe-cell normalization sugesstion presented by Marioni el.al (2016)."),
                   br(),
                   br(),
                   p("Cluster Method"),
                   p("cluster similar cells together and normalize the cells in each cluster using the deconvolution method"),
                   
                   p("This plot shows the correlation between size factors and library size for all ST-features
                     If we have a tight correlation (i.e. straight line), that would suggesst that systematic differences
                     between ST-features are primarily diven by differences in capture efficiency, 
                     IF one would assume that there is approx. equal total transcriptonal activity in all ST-features.
                     DE between ST-features would however make the correlation trend non-linear, and would result in
                     scatter around the trend line"),
                   br(),
                   h3("Applying the size factors to normalize gene expression"),
                   p("The count data are used to compute normalized log-expression values for use in downstream analyses. 
                     Each value is defined as the log-ratio of each count to the size factor for the corresponding ST-feature, 
                     after adding a prior count of 1 to avoid undefined values at zero counts. 
                     Division by the size factor ensures that any ST-feature-specific biases are removed."),
                   
                   
                   
                   br(),
                   br(),
                  
                   h2("Filtering"),
                   
                   br(),
                   br(),
                   
                   h2("Differential expression"),
                   h4("A word of caution regarding p-values"),
                                            p("It should be stressed that, depending on your clustering/selection method, 
                                             the p-values might be reported as lower than in reality, as the algorithm for DE calculation 
                                              do not take into account the uncertainty and stochasticity introduced during clustering.
                                              Accordingly, these p-values should be used for ranking purposes and not interpretated as 
                                              \"truth\" p-values. However, if groups/clusters are pre-defined, FDR-adjusted p-values 
                                              can be used diretcly to define significat genes after DE comparison")
                   
                   )
           
           )

