library(shiny)
library(shinythemes)
library(shinyBS)
library(DESeq2)
library(ggplot2)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("genefilter")
library(ggrepel)
library(dplyr)
library(tidyr)
library(scater)
library(scran)
library(edgeR)
library(png)
library(grid)
library(maps)
library(threejs)
library(leaflet)
library(leaflet.extras)
#library(mapview) #for mapshot = save image
library(RColorBrewer)
library(pheatmap)
library(sp) #for point in polygon
library(cowplot)
library(Rtsne)
#library(CountClust)

#if(!require("devtools")) install.packages("devtools")
#devtools::install_github("bwlewis/rthreejs")
#increase max upload for images
#devtools::install_github('bhaskarvk/leaflet.extras')
#packrat::init(options = list(ignored.packages = c("gganimate"))
#library(devtools)
#install_github('kkdey/CountClust')
#source("https://bioconductor.org/biocLite.R")
#biocLite()
options(shiny.maxRequestSize=150*1024^2) 

#sceObjs <<- list() #Global variable for saving of all sceObjects
#sceObjsPaths <<- list()
#spotDataObjs <<- list() #Global variable for saving spot data associated with each data set
#spotDataObjsPaths <<- list()
#imgScaleObjs <<- list()


shinyServer(function(input,output, session) {
  
  observe({ print(input$navbar) }) #this function observes when the user changes tabpanels
  
  
# ----------------------------- DATA INPUT -----------------------------
observeEvent(input$resetDatasets,{
    sceObjs <<- list() #Global variable for saving of all sceObjects
    sceObjsPaths <<- list()
    spotDataObjs <<- list() #Global variable for saving spot data associated with each data set
    spotDataObjsPaths <<- list()
    imgScaleObjs <<- list()
    
    output$resetDatasetsMsg <- renderUI({
      HTML(paste("resetted"))
    })
})
  
  output$countInputs=renderUI({
    html_ui = " "
    for (i in 1:input$nfiles){
      html_ui <- paste0(html_ui, fileInput(paste0("count",i), label=paste0("Count file: ",i)))
      
    }
    HTML(html_ui)
  })
  output$spotInputs=renderUI({
    html_ui2 = " "
    for (i in 1:input$nfiles){
      html_ui2 <- paste0(html_ui2, fileInput(paste0("spot",i), label=paste0("Spot file: ",i)))
    }
    HTML(html_ui2)
  })
  
  
  observeEvent(input$DataCreationButton, {
    nrObjs = input$nfiles
    print(nrObjs)
    for (i in 1:nrObjs){
      print(paste("file", i))
      inFilenr = paste0("count",i)
      inFile=input[[inFilenr]]
      print(inFile)
      sceObjsPaths[[i]] <<- inFile[['datapath']]
      
      inFilenr = paste0("spot",i)
      inFile = input[[inFilenr]]
      print(inFile)
      spotDataObjsPaths[[i]] <<- inFile[['datapath']]
    }

      for (n in 1:nrObjs){
          withProgress(message =paste("Creating SCE object nr ", n, "..", sep=""), { #Progress bar START
          st <- read.table(sceObjsPaths[[n]])
    
          st = t(st)
    
          split = strsplit(colnames(st), "x")
          x_cord = sapply(split, "[[", 1)
          y_cord = sapply(split, "[[", 2)
    
          sce <- newSCESet(countData=st)
          sce <- calculateQCMetrics(sce)
          pData(sce)$coord_x = x_cord
          pData(sce)$coord_y = y_cord
          pData(sce)$DataSet = n
          pData(sce)$selectionID = 0 #initierar selection ID = 0
    
          selected<<-NULL
          selected_list <<- NULL
    
    
          pData(sce)$selectionSpotsInput = rep(0, ncol(sce))
          selectionSpotsInputVector <<- rep(0, ncol(sce))
          pData(sce)$AnnotationNameDataSet = "dummy"
    
          symbols = mapIds(org.Hs.eg.db,
                      keys=row.names(sce),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
          rownames(sce) = make.names(symbols, unique=TRUE)
          
          sceObjs[[n]] <<- sce
          
          #load spot data
              spots = read.table(spotDataObjsPaths[[n]], skip=1)
              imgScale = readLines(spotDataObjsPaths[[n]], n=1)
              spots$barcode = paste(spots$V1, spots$V2, sep="x")
              spotDataObjs[[n]] <<- spots #save to global spot data storage
              imgSplit = strsplit(imgScale, "\t")
              imgWidth = sapply(imgSplit, "[[", 3)
              imgHeight = sapply(imgSplit, "[[", 4)
              imgScale = c(as.numeric(imgWidth), as.numeric(imgHeight))
              imgScaleObjs[[n]] <<- imgScale
    
          incProgress(detail = paste("Done!")) #Progress bar END
          Sys.sleep(1)})
          print(paste("DataSet nr", n,  "loaded", sep=": "))
    
      }
      
      outputText = (paste("<p>Data set 1:</p> ", "<p>Nr of ST-Spots: ", dim(sceObjs[[1]])[2],  "</p> <p>Nr of Genes: ", 
                          dim(sceObjs[[1]])[1],"<p>- - - - - - - -</p>", sep=""))
      if (length(sceObjs)>1){
      for (y in 2:length(sceObjs)){
        outputText = append(outputText, paste("<p>Data set ", y,":</p>", "<p>Nr of ST-Spots: ", dim(sceObjs[[y]])[2],  "</p> <p>Nr of Genes: ", 
                            dim(sceObjs[[y]])[1], " <p>- - - - - - - -</p>",sep=""))
      }
      }
      
      output$uploadedData <- renderUI({
        HTML(outputText)
      })
      
      output$NameDataSetFactorDataSetUI = renderUI({
        selectInput("sceFactorNameSelect", "Choose data set", 
                    choices  = seq(1, length(sceObjs), 1))})
    
    
    })
    
observeEvent(input$AddDataSetFactorButton,{
  #Set the sceObject to user selected data set
  sce = sceObjs[[as.numeric(input$sceFactorNameSelect)]]
  pData(sce)$AnnotationNameDataSet = input$NameDataSetFactorInput
  sceObjs[[as.numeric(input$sceFactorNameSelect)]] <<- sce #Save back to global
  
  output$DataSetFactorNameAdded = renderText(print("Added"))
    
  
})

    
# --------------------------Filter ------------------------
  
  observeEvent(input$FilterButton,{
    
    if (length(sceObjs)==0){
      output$CurrentDataStatus = renderUI({
        HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))})
    }else{
      output$CurrentDataStatus = renderUI({
        selectInput("sceFilterSelect", "Which data set do you want to filter", 
                           choices  = seq(1, length(sceObjs), 1))})
    }
  })
  
  # ------- All filter functions under here
observeEvent(input$UpdateFilters, {
  
  #Reset all output texts regarding removed ST-Spots/genes
  output$avgCount_filter_done = renderText({
    print("")
  })
  output$librarySize_filter_done = renderText({
    print("")
  })
  output$insideOutsideTissueFilterTextOutput = renderText({
    print("")
  })
  output$lowgene_spots_filer_done = renderText({
    print("")
  })
  output$CommonFiltersOutput = renderText({
    print("")
  })
  output$geneID_filter_done = renderText({
    print("")
  })
  output$regex_filter_done = renderText({
    print("")
  })
    
    #Set the sceObject to user selected data set
    sce = sceObjs[[as.numeric(input$sceFilterSelect)]]
    
  output$filterPlot1 <- renderPlot({
    
    
    slider1 = input$filterPlot1_slider
    ggtest = ggplot(data=pData(sce), aes(x=total_counts/1e3), colour="#faebd7")+
      geom_histogram(binwidth=0.1, colour="#faebd7", position="dodge")+ylab("Nr of ST-spots")+xlab("Library sizes (tousands)")+
    geom_vline(colour="orange", xintercept=(slider1)) +
    theme(
      legend.background = element_blank(),
      panel.background =element_blank(),
      plot.background = element_rect(fill="#2E3337"),
      axis.title.x = element_text(colour = "#f5925e"),
      axis.text.x = element_text(colour="#f5925e"),
      axis.text.y = element_text(colour="#f5925e"),
      axis.title.y = element_text(colour = "#f5925e"),
      legend.title = element_text(colour = "#f5925e"),
      legend.text = element_text(colour = "#f5925e"),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "bottom")
    ggtest
})

  
  output$filterPlot2 <- renderPlot({
    slider1 = input$filterPlot2_slider

    ggtest = ggplot(data=pData(sce), aes(x=total_features, colour="#faebd7"))+geom_histogram(binwidth=100, colour="#faebd7")
    ggtest = ggtest + geom_vline(colour="orange", xintercept=slider1) +
      xlab("Nr of expressed genes/spot") + ylab("Nr of ST-Spots")+
      theme(
        legend.background = element_blank(),
        panel.background =element_blank(),
        plot.background = element_rect(fill="#2E3337"),
        axis.title.x = element_text(colour = "#f5925e"),
        axis.text.x = element_text(colour="#f5925e"),
        axis.text.y = element_text(colour="#f5925e"),
        axis.title.y = element_text(colour = "#f5925e"),
        legend.title = element_text(colour = "#f5925e"),
        legend.text = element_text(colour = "#f5925e"),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")
    ggtest
  })

  output$librarySizeRemoveText = renderText({
    before = dim(sce)[2]
    libSize <- colSums(counts(sce))
    keep <- libSize >= input$filterPlot1_slider*1000
    sce <- sce[,keep]
    after = dim(sce)[2]
    print(paste("The current filter would remove: ", before-after, " spots",sep=""))
  })
  
  output$current_filter_2 = renderText({
    before = dim(sce)[2]
    x = pData(sce)$total_features
    keep <- x >= input$filterPlot2_slider
    sce <- sce[,keep]
    after = dim(sce)[2]
    print(paste("The current filter would remove: ", before-after, " ST-Spots",sep=""))
  })
})
  
observeEvent(input$Apply_avgCount_Filter, {
    sce = sceObjs[[as.numeric(input$sceFilterSelect)]]
    before = dim(sce)[1]
    ave.counts <- rowMeans(counts(sce))
    keep <- ave.counts >= input$avgCountFilterInput
    sce <- sce[keep,]
    after = dim(sce)[1]
    
    sceObjs[[as.numeric(input$sceFilterSelect)]] <<- sce #save to global
    
    output$avgCount_filter_done = renderText({
      print(paste("Number of genes removed: ", before-after, sep=""))
    })
  })
  
observeEvent(input$Apply_librarySize_Filter, {
    sce = sceObjs[[as.numeric(input$sceFilterSelect)]]
    before = dim(sce)[2]
    libSize <- colSums(counts(sce))
    keep <- libSize >= input$filterPlot1_slider*1000
    sce <- sce[,keep]
    after = dim(sce)[2]
    
    sceObjs[[as.numeric(input$sceFilterSelect)]] <<- sce #save to global
    output$librarySize_filter_done = renderText({
      print(paste("Number of spots removed: ", before-after, sep=""))
    })
  })
  
observeEvent(input$insideOutsideTissueFilterButton, {
    sce = sceObjs[[as.numeric(input$sceFilterSelect)]]
    spot_data = spotDataObjs[[as.numeric(input$sceFilterSelect)]]
    spot_data = spot_data[which(spot_data$barcode %in% colnames(sce)),]
    sce = sce[, which(colnames(sce) %in% spot_data$barcode)]
    #sort to read in spot data into sce object
    sce = sce[, order(colnames(sce))] #string sorted
    spot_data = spot_data[order(spot_data$barcode), ]
    spot_data = spot_data[, 5] #(V5 = inside/outside tissue - yes/true)
    pData(sce)$insideTissue = spot_data
    before = dim(sce)[2]
    sce = sce[, which(pData(sce)$insideTissue == 1)]
    after = dim(sce)[2]
    
    sceObjs[[as.numeric(input$sceFilterSelect)]] <<- sce #save to global
    
    
    output$insideOutsideTissueFilterTextOutput = renderText({
      print(paste("Number of spots removed: ", before-after, sep=""))
    })
  })
  
  # observeEvent(input$applyAbundacetypeTwoFilterButton, {
  #   nrSpotsWithNonZero = input$n_cells
  #   #keep all genes that have <0 zero counts in n Spots
  #   #Count the spots
  #   spotNr = dim(sce)[[2]]
  #   #Count the zeros for each gene
  #   zeros = rowSums(counts(sce))
  # })
  
observeEvent(input$Apply_lowgene_spots_filter, {
    sce = sceObjs[[as.numeric(input$sceFilterSelect)]]
    print(paste("Data set selected: ",input$sceFilterSelect, sep=""))
    
    before = dim(sce)[2]
    x = pData(sce)$total_features
    keep <- x >= input$filterPlot2_slider
    sce <- sce[,keep]
    after = dim(sce)[2]
    
    sceObjs[[as.numeric(input$sceFilterSelect)]] <<- sce #save to global
    
    output$lowgene_spots_filer_done = renderText({
      print(paste("Number of spots removed: ", before-after, sep=""))
    })
  })
  
observeEvent(input$Apply_common_filter, {
    sce = sceObjs[[as.numeric(input$sceFilterSelect)]]
    print(paste("Data set selected: ",input$sceFilterSelect, sep=""))
    before = dim(sce)[1]
    if ("MALAT1" %in% input$CommonFilters){
      remove_idx=which(rownames(sce) == "MALAT1")
      sce <- sce[-remove_idx, ]
    }
    if ("Ribosomal genes" %in% input$CommonFilters){
     #Regex to remove ribosomal genes
     #remove_idx = grep(rownames(sce), pattern="(RPS)|(RPL)|(RRP)")
      ribosomalIDs = c("RPL3",   "RPL4"   ,"RPL5"   ,"RPL6"   ,"RPL7"   ,"RPL7A"  ,"RPL8"   ,"RPL9"   ,"RPL10"  ,"RPL10A" ,"RPL11"  ,"RPL12", 
                       "RPL13",  "RPL13A" ,"RPL14"  ,"RPL15"  ,"RPL17"  ,"RPL18"  ,"RPL18A" ,"RPL19"  ,"RPL21"  ,"RPL22"  ,"RPL23"  ,"RPL23A",
                       "RPL24",  "RPL26"  ,"RPL27"  ,"RPL27A" ,"RPL28"  ,"RPL29"  ,"RPL30"  ,"RPL31"  ,"RPL32"  ,"RPL34"  ,"RPL35"  ,"RPL35A",
                       "RPL36",  "RPL36A" ,"RPL37"  ,"RPL37A" ,"RPL38"  ,"RPL39"  ,"RPL40"  ,"RPL41"  ,"RPLP0"  ,"RPLP1"  ,"RPLP2"  ,"RPSA"  ,
                       "RPS2" ,  "RPS3"   ,"RPS3A"  ,"RPS4X"  ,"RPS4Y"  ,"RPS5"   ,"RPS6"   ,"RPS7"   ,"RPS8"   ,"RPS9"   ,"RPS10"  ,"RPS11" ,
                       "RPS12"  ,"RPS13"  ,"RPS14"  ,"RPS15"  ,"RPS15A" ,"RPS16"  ,"RPS17"  ,"RPS18"  ,"RPS19"  ,"RPS20"  ,"RPS21"  ,"RPS23" ,
                       "RPS24"  ,"RPS25"  ,"RPS26"  ,"RPS27"  ,"RPS27A" ,"RPS28"  ,"RPS29"  ,"RPS30" )
      remove_idx = which(rownames(sce) %in% ribosomalIDs)
      sce = sce[-remove_idx, ]
    }
    after = dim(sce)[1]
    sceObjs[[as.numeric(input$sceFilterSelect)]] <<- sce #save to global
    
    output$CommonFiltersOutput = renderText({
      print(paste("Number of genes removed: ", before-after, sep=""))
    })
  })
  
observeEvent(input$Apply_geneID_filter, {
    sce = sceObjs[[as.numeric(input$sceFilterSelect)]]
    print(paste("Data set selected: ",input$sceFilterSelect, sep=""))
    before = dim(sce)[1]
    print(input$geneID_filter)
    ids = strsplit(input$geneID_filter, ",")
    remove_idx=which(rownames(sce) %in% ids[[1]])
    sce <- sce[-remove_idx, ]
    after = dim(sce)[1]
    
    sceObjs[[as.numeric(input$sceFilterSelect)]] <<- sce #save to global
    
    output$geneID_filter_done = renderText({
      print(paste("Number of genes removed: ", before-after, sep=""))
    })
  })
  
observeEvent(input$Apply_regex_filter, {
    sce = sceObjs[[as.numeric(input$sceFilterSelect)]]
    print(paste("Data set selected: ",input$sceFilterSelect, sep=""))
    before = dim(sce)[1]
    remove_genes_index <<- grep(rownames(sce), pattern=input$RegexGeneID_filter)
    regexTable <<- rownames(sce)[remove_genes_index]
    sce <- sce[-remove_genes_index, ]
    after = dim(sce)[1]
    
    sceObjs[[as.numeric(input$sceFilterSelect)]] <<- sce #save to global
    
    output$regexfilterGenes = renderTable({
      regexTable
    })
   
   
    output$regex_filter_done = renderText({
      print(paste("Number of genes removed: ", before-after, sep=""))
      })
  })

  
  
# ---------------------- Clustering  ----------------

  
observeEvent(input$ClusterButton, {  
  #Initiated when cluster bsModal are pressed
  output$clusterDataSetChooser <- renderUI({
    if (length(sceObjs)==0){
      HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))
    }else{
      checkboxGroupInput("clusterDataSetChooserCheckBox", "Check the data sets you want to include", 
                         choices  = seq(1, length(sceObjs), 1))
    }
  })
})
  
  
observeEvent(input$performClusterButton, {
    
    inputCheck <- tryCatch(
      input$clusterDataSetChooserCheckBox,
      error = function(e) {
        output$clusterMsgOutput = renderText({
          HTML(paste(tags$span(style="color:red", e), sep = ""))
        })
        return(NULL)
      }
    )
    
    if (!is.null(inputCheck)){
    
    withProgress(message =paste("Clustering..."), { #Progress bar START

    #if multiple data sets -> merge
    if (length(input$clusterDataSetChooserCheckBox)>1){
      sce <- mergeSCE(as.numeric(input$clusterDataSetChooserCheckBox)) #function returns merged data sets
    }
    else{
      sce <- sceObjs[[as.numeric(input$clusterDataSetChooserCheckBox)]]
    }
    
    if (dim(sce)[[1]]>2){
    
    #we now have a sce data set, next step choose cluster method and eventual set min.sizes etc
    
    
    if(input$clusterApproach == "Spearman correlation"){
      
      if (input$clusterMinSizeB == TRUE){
          print("min sizes used")
          clusters <- quickCluster(sce, min.size=input$clusterMinSizeInput)
          pData(sce)$spearmanCluster <- clusters  
          
        }else{
          clusters <- quickCluster(sce)
          pData(sce)$spearmanCluster <- clusters  
        }
      
    }
      
    
    #Output the cluster table info
      output$clusterTableOutput = renderTable({
        table(clusters)
      })
      
    #Roll back and save in global sceObjs
      #subset the sceObjs according to pData(sceObjs)$DataSet
      datasets = unique(pData(sce)$DataSet)
      print(paste("datasets:", datasets, sep=""))
      for (i in 1:length(unique(pData(sce)$DataSet))){
        sce_tmp = sce[, which(pData(sce)$DataSet==datasets[i])]
        print(sce_tmp)
        curr_dataset = datasets[i]
        print(paste("curr dataset:", curr_dataset, sep=""))
        #before rolling back we need to "reset" colnames after merge function
        split = strsplit(colnames(sce_tmp), "-")
        colnames(sce_tmp) = sapply(split, "[[", 1)
        sceObjs[[curr_dataset]] <<- sce_tmp #roll back to global sceObjs
      }

      incProgress(detail = paste("Done!")) #Progress bar END
      Sys.sleep(1)}
      else{
        output$clusterMsgOutput = renderUI({
          HTML(paste(tags$span(style="color:red", "Data set to small"), sep = ""))
        })
      }
    })
    }
    })

  

  # -------------- Cluster Plots ---
  
  observeEvent(input$ClusterButton, {
     output$clusterDataSetPlot <- renderUI({
       if (length(sceObjs)==0){
         HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))
       }else{
         checkboxGroupInput("clusterDataSetPlotCheckBox", "Plot the following data set(s)", 
                            choices  = seq(1, length(sceObjs), 1))
       }
     }) 
   })
  
  output$clusterArray_plot = renderPlot({
    if(input$clusterArray_plotB>0){
      
      isolate({

          sceObjs = sceObjs[as.numeric(input$clusterDataSetPlotCheckBox)]
          print(sceObjs)
          nrPlots = length(sceObjs)
          plotList = list()
        
          for (i in 1:nrPlots){
    
          gg <- ggplot(data=pData(sceObjs[[i]]), aes(x=as.numeric(coord_x), 
                                                     y=-as.numeric(coord_y),
                                                     colour=spearmanCluster)) + 
          xlab("X coordinate")+ylab("Y coordinate")+ coord_equal() +
          geom_point()+
            scale_color_discrete(name="Cluster")+
            theme(
              legend.background = element_blank(),
              panel.background =element_blank(),
              plot.background = element_rect(fill="#2E3337"),
              axis.title.x = element_text(colour = "orange"),
              axis.title.y = element_text(colour = "orange"),
              legend.title = element_text(colour = "orange"),
              legend.text = element_text(colour = "orange"),
              panel.grid = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "bottom")
              plotList[[i]] = gg
        }
       
        plot_grid(plotlist=plotList)
      })
  }}, bg = "transparent")
  
  
# ------------------------------- Normalization --------------------
  # OBS need to fix with sceOjbs global
  
  observeEvent(input$NormalizationButton, {
    output$sceNormUI <- renderUI({
      if (length(sceObjs)==0){
        HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))
      }else{
        checkboxGroupInput("sceNormCheckBox", "Check the data sets you want to include", 
                           choices  = seq(1, length(sceObjs), 1))
      }
    }) 
  })
  
  #size factor summary
observeEvent(input$sizefactorB, {
  
  inputCheck <- tryCatch(
    input$sceNormCheckBox,
    error = function(e) {
      output$NormalizeMsgOutput = renderText({
        HTML(paste(tags$span(style="color:red", e), sep = ""))
      })
      return(NULL)
    }
  )
  
  if (!is.null(inputCheck)){
  
  positive = input$forcePositive
  withProgress(message =paste("Computing size factors"), { #Progress bar START
    
  if (input$NormSizesInput=="NA"){ # +++++ If no sizes are set
      
              #if multiple data sets -> merge
              if (length(input$sceNormCheckBox)>1){
                sce <- mergeSCE(as.numeric(input$sceNormCheckBox)) #function returns merged data sets
              }
              else{
               sce <- sceObjs[[as.numeric(input$sceNormCheckBox)]]
              }
              if (input$normclusters=="Spearman"){
                clusters = pData(sce)$spearmanCluster
                print(clusters)
                tryCatch(
                  sce <- computeSumFactors(sce, clusters = clusters, positive=positive),
                  error = function(e) {
                    output$sizefactorSummary = renderUI({
                      HTML(paste(tags$span(style="color:red", e), sep = ""))
                    })
                    return(NULL)
                  }
                )
                
              }else{
                tryCatch(
                  sce <- computeSumFactors(sce, positive=positive),
                  error = function(e) {
                    output$sizefactorSummary = renderUI({
                      HTML(paste(tags$span(style="color:red", e), sep = ""))
                    })
                    return(NULL)
                  }
                )
                
              }
        
            
          len_s <- length(which(sizeFactors(sce)==0))
          output$sizefactorSummary = renderPrint({summary(sizeFactors(sce))})
  }else{# +++++ If sizes are set
    
            sizes=input$NormSizesInput
            sizes = strsplit(sizes, ",")
            sizes = unlist(sizes)
            sizes = as.numeric(sizes)
            print(paste("sizes input: ", sizes))
            
            #if multiple data sets -> merge
              if (length(input$sceNormCheckBox)>1){
              sce <- mergeSCE(as.numeric(input$sceNormCheckBox)) #function returns merged data sets
              }
              else{
              sce <- sceObjs[[as.numeric(input$sceNormCheckBox)]]
              }
              
            
            if (input$normclusters=="Spearman"){
                    clusters = pData(sce)$spearmanCluster
                    print(clusters)
                    print(sizes)
                    tryCatch(
                      sce <- computeSumFactors(sce, clusters=clusters, sizes=sizes, positive=positive),
                      error = function(e) {
                        output$sizefactorSummary = renderUI({
                          HTML(paste(tags$span(style="color:red", e), sep = ""))
                        })
                        return(NULL)
                      }
                    )
                    
              }else{
                tryCatch(
                  sce <- computeSumFactors(sce, sizes=sizes, positive=positive),
                  error = function(e) {
                    output$sizefactorSummary = renderUI({
                      HTML(paste(tags$span(style="color:red", e), sep = ""))
                    })
                    return(NULL)
                  }
                )
              }
              
              len_s <- length(which(sizeFactors(sce)==0))
              output$sizefactorSummary = renderPrint({summary(sizeFactors(sce))})
  }
          
          #We now have a sce object (merged or not) with size factors
          #This sce object are temp saved globaly to sce_norm_tmp to be used in normaliazation function
          #This is done as the user may want to check the result of the size factors etc
          #before applying the normalization values on to global sceObjs object
        
    sce_norm_tmp <<- sce
          
        
        incProgress(detail = paste("Done!")) #Progress bar END
        Sys.sleep(1)})
  }else{
    output$sizefactorSummary = renderUI({
      HTML(paste(tags$span(style="color:red", "Data set to small"), sep = ""))
    })
  }
})
  
  
  
  #Size factor plot
output$sizeFactorPlot <- renderPlot({
  
    if (input$sizeFactorPlotB>0 && exists("sce_norm_tmp")){
      isolate({
    
            fit <- lm(sizeFactors(sce_norm_tmp) ~ sce_norm_tmp$total_counts, data=(sce_norm_tmp))

              df = cbind(sizeFactors(sce_norm_tmp), sce_norm_tmp$total_counts/1e3)
              df = as.data.frame(df)
              ggSize = ggplot(data=df, aes(x=log(V1), y=log(V2)))+geom_point( colour="#faebd7")+
              xlab("Size Factor")+ylab("Library Size (Thousands)")+
              theme(
                legend.background = element_blank(),
                panel.background =element_blank(),
                plot.background = element_rect(fill="#2E3337"),
                axis.title.x = element_text(colour = "orange"),
                axis.text = element_text(colour = "#faebd7"),
                axis.title.y = element_text(colour = "orange"),
                legend.title = element_text(colour = "orange"),
                legend.text = element_text(colour = "orange"),
                panel.grid = element_blank(),
                panel.grid.major = element_blank(),
                legend.position = "bottom")
            ggSize
            })
        }else{
            output$NormalizeMsgOutput = renderUI({
              HTML(paste(tags$span(style="color:red", "No size factors found"), sep = ""))
            })
        }}, bg = "transparent")
      

#Check size factors
observeEvent(input$checkNormB, {
  
  if (exists("sce_norm_tmp")){
  
  sizeFactorsOfZero = length(which(sizeFactors(sce_norm_tmp)==0))
  sizeFactorsNegative = length(which(sizeFactors(sce_norm_tmp)<0))
  
  output$checkNorm = renderPrint({

        #Check if there are size factors = 0 after computing size factors
          if (sizeFactorsOfZero>0 & sizeFactorsNegative==0){
            cat(paste("Detected ", sizeFactorsOfZero),  "ST-Spots with Size Factors = 0")
          }
          else if (sizeFactorsNegative>0){
            cat(paste("Detected ", sizeFactorsNegative, " ST-Spots with negative Size Factors and ", sizeFactorsOfZero, " ST-Spots with Size Factors = 0", sep=""))
          }else{
            cat("Size Factors look ok")
          }
  })
  }else{
    output$NormalizeMsgOutput = renderUI({
      HTML(paste(tags$span(style="color:red", "No data set chosen"), sep = ""))
    })
  }
})
  

#Remove Spots which have non positive size factors

observeEvent(input$removeSFB, {

if(exists("sce_norm_tmp")){
  
  badSF = length(which(sizeFactors(sce_norm_tmp)<=0))
  print(names(which(sizeFactors(sce_norm_tmp)<=0)))
  sce_norm_tmp <<- sce_norm_tmp[,!colnames(sce_norm_tmp) %in% names(which(sizeFactors(sce_norm_tmp)<=0))]
  
  
  output$removeSF = renderPrint({
    cat(paste(badSF, " ST-features removed", sep=""))
      })
}
})
  
 
  
#Apply normalization
observeEvent(input$normalizeB, {
  

        #Normalize
        inputCheck <- tryCatch(
          sce_norm_tmp = normalise(sce_norm_tmp),
        error = function(e) {
          output$NormalizeMsgOutput = renderUI({
            HTML(paste(tags$span(style="color:red", e), sep = ""))
          })
        return(NULL)
        }
        )
        
        if (!is.null(inputCheck)){
        
        
        #Roll back to sceObj global
        datasets = unique(pData(sce_norm_tmp)$DataSet)
        print(paste("datasets:", datasets, sep=""))
          for (i in 1:length(unique(pData(sce_norm_tmp)$DataSet))){
            sce_tmp = sce_norm_tmp[, which(pData(sce_norm_tmp)$DataSet==datasets[i])]
            curr_dataset = datasets[i]
            print(paste("curr dataset:", curr_dataset, sep=""))
            #before rolling back we need to "reset" colnames after merge function
            split = strsplit(colnames(sce_tmp), "-")
            colnames(sce_tmp) = sapply(split, "[[", 1)
            sceObjs[[curr_dataset]] <<- sce_tmp #roll back to global sceObjs
          }
        
        output$normalizationDone = renderPrint({HTML(paste("Done"), sep = "")})
        }
})
        
#Info
output$forcePositiveInfoText = renderUI({
  HTML("<p> Note that the positive=TRUE parameter one can use during the
              scran normalization procedure do not protect against zero size factors.</p>
          <p> It will however avoide negative size factors and keep these ST-Spots from disorting the other estimates.</p> 
          <p> Obs, note that while this \"forced\" non-negative size factors might be of use, the fact that we observe
              size factors that are negative or zero is a major indication low quality ST-Spots, and it's advised that 
              you make sure that your QC and filtering procedures are appropiate. Ideally, positive=TRUE is a measure
              that you should have to use.</p>")
})
 
  
#------------------------------- QC ----------
observeEvent(input$QCButton, {  
    #Initiated when cluster bsModal are pressed
    output$QCDataSetRadioButtons <- renderUI({
      if (length(sceObjs)==0){
        HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))
      }else{
        radioButtons("QCDataSetRadioButtons", "Check the data sets you want to include", 
                           choices  = seq(1, length(sceObjs), 1))
      }
    })
  })

observeEvent(input$insideOutsideButton, {

  withProgress(message ="Running ...", { #Progress bar START
  
  ggInsideList <<- list()
  ggInside2List <<- list()
  qcTable <<- list()
  noGoVar <<- 1
  
  for (i in 1:length(sceObjs)){
    sce = sceObjs[[i]]
    if (dim(sce)[[2]]>1){
      noGoVar <<- 0
    #just consider the spots which we have data from
    spot_data = spotDataObjs[[i]]
    spot_data = spot_data[which(spot_data$barcode %in% colnames(sce)),]
    sce = sce[, which(colnames(sce) %in% spot_data$barcode)]
    #sort to read in spot data into sce object
    sce = sce[, order(colnames(sce))] #string sorted
    spot_data = spot_data[order(spot_data$barcode), ]
    spot_data = spot_data[, 5] #(V5 = inside/outside tissue - yes/true)
    pData(sce)$insideTissue = spot_data
    
    #Total counts inside/outside tissue
    total_inside = sum(pData(sce[,which(pData(sce)$insideTissue == 1)])$total_counts)
    total_outside = sum(pData(sce[,which(pData(sce)$insideTissue == 0)])$total_counts)
    #Number of spots inside/outside tissue
    nr_inside = dim(sce[,which(pData(sce)$insideTissue == 1)])[[2]]
    nr_outside = dim(sce[,which(pData(sce)$insideTissue == 0)])[[2]]
    #Counts/Spot inside/outside tissue
    perSpot_inside = total_inside/nr_inside
    perSpot_outside = total_outside/nr_outside
    
    #create table output
    qcTable[[i]] <<- cbind(total_inside, total_outside, round(perSpot_inside, 0), round(perSpot_outside,0))
   
    #Plots
    datasetName =   pData(sce)$AnnotationNameDataSet
    if (datasetName == "dummy"){
      datasetName = ""
    }
    #-- Inside/Outside Tissue histo
    ggInsideList[[i]] <<- ggplot(data=pData(sce), aes(x=total_features, fill=as.factor(insideTissue)))+
      geom_histogram(binwidth=100, colour="#faebd7", position="fill")+
      ggtitle("Inside Tissue") +
      labs(title=paste("Inside/Outside Tissue DataSet ", i, " | ", datasetName), x="Nr of Expressed genes/spot", y="Percent", fill="Outside/Inside")+
      #xlab("Nr of expressed genes/spot") + ylab("Nr of ST-Spots")+
      theme(
        legend.background = element_blank(),
        panel.background =element_blank(),
        plot.background = element_rect(fill="#2E3337"),
        axis.title.x = element_text(colour = "#f5925e"),
        axis.text.x = element_text(colour="#f5925e"),
        axis.text.y = element_text(colour="#f5925e"),
        title = element_text(colour="#f5925e"),
        axis.title.y = element_text(colour = "#f5925e"),
        legend.title = element_text(colour = "#f5925e"),
        legend.text = element_text(colour = "#f5925e"),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")
    
    
    ggInside2List[[i]] <<- ggplot(data=pData(sce), aes(x=total_features, fill=as.factor(insideTissue)))+
      geom_histogram(binwidth=100, colour="#faebd7")+
      ggtitle("Inside Tissue") +
      labs(title=paste("Inside/Outside Tissue DataSet ",  i, " | ", datasetName), x="Nr of Expressed genes/spot", y="Nr of ST-Spots", fill="Outside/Inside")+
      #xlab("Nr of expressed genes/spot") + ylab("Nr of ST-Spots")+
      theme(
        legend.background = element_blank(),
        panel.background =element_blank(),
        plot.background = element_rect(fill="#2E3337"),
        axis.title.x = element_text(colour = "#f5925e"),
        axis.text.x = element_text(colour="#f5925e"),
        axis.text.y = element_text(colour="#f5925e"),
        title = element_text(colour="#f5925e"),
        axis.title.y = element_text(colour = "#f5925e"),
        legend.title = element_text(colour = "#f5925e"),
        legend.text = element_text(colour = "#f5925e"),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")
    
    }
  }
  
  incProgress(detail = paste("Completed!")) #Progress bar END
  Sys.sleep(1)})
  
  
  #Data table output
  output$qcCheckDataSelectionUI = renderUI({
    checkboxGroupInput("qcCheckSelectionInput", "Show the following data set(s)",  choices  = seq(1, length(sceObjs), 1))
    
  })
  
  #Create action button for show plots and table
  output$qcCheckDataSelctionButton = renderUI({
    actionButton("qcPlotsAndTableButton", "Show plots and table")
  })
})

observeEvent(input$qcPlotsAndTableButton, {
  inputCheck <- tryCatch(
    input$qcCheckSelectionInput,
    error = function(e) {
      output$NormalizeMsgOutput = renderText({
        HTML(paste(tags$span(style="color:red", e), sep = ""))
      })
      return(NULL)
    }
  )
  
  if (!is.null(inputCheck) && noGoVar==0){
  #Create data table output
  outputTable = data.frame("Data set" = numeric(), "Tot. inside"=numeric(), "Tot. outside"=numeric(), "Transcript/spot inside"=numeric(), "Transcript/spot outside"=numeric())
  for (i in 1:length(input$qcCheckSelectionInput)){
    datasetNr = input$qcCheckSelectionInput[i]
    dfRow = qcTable[[i]]
    dfRow = cbind(datasetNr, dfRow)
    colnames(dfRow) = c("Data set", "Tot. inside", "Tot. outside", "Transcript/spot inside", "Transcript/spot outside")
    outputTable = rbind(outputTable, dfRow)
    print(dfRow)
  }
  output$qcDataTableOutputUI = renderUI({
    output$qcDataTable <- renderDataTable(outputTable)
    dataTableOutput("qcDataTable")
  })
  
  #output plots
  ggInside <<- ggInsideList[c(as.numeric(input$qcCheckSelectionInput))]
  ggInside2 <<- ggInside2List[c(as.numeric(input$qcCheckSelectionInput))]
  
  output$qcInsideOutsidePlots = renderPlot({
    plot_grid(plotlist=ggInside)
  }, bg="#2E3337")
  output$qcInsideOutsidePlots2 = renderPlot({
    plot_grid(plotlist=ggInside2)
  }, bg="#2E3337")
  }
})


observeEvent(input$geneQCB, {
      
    output$geneQC = renderPlot({
    isolate({
      sce <- sceObjs[[as.numeric(input$QCDataSetRadioButtons)]]
      
      # ------ Looking at most abundance genes across all features 
      fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
      QC_abund <- plotQC(sce, type = "highest-expression", n=50) + fontsize
      QC_abund +theme(
        legend.background = element_blank(),
        panel.background =element_blank(),
        plot.background = element_blank(),
        axis.title.x = element_text(colour = "orange"),
        axis.title.y = element_text(colour = "orange"),
        legend.title = element_text(colour = "orange"),
        legend.text = element_text(colour = "orange"),
        panel.grid = element_line(linetype="dashed"),
        panel.grid.major = element_line(colour = "gray"))
      
      QC_abund
    })
  
  }, bg = "transparent")
})
    
  
output$downloadQCabundPlot <- downloadHandler(
    filename = "QCabundances_plot.png",
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height,
                       res = 300, units = "in")
      }
      ggsave(file, plot = QC_abund, device = device)
    })
  

  
  
# ------------------------ HVGs ------
      
  observeEvent(input$HVGButton,{
    
    if (length(sceObjs)==0){
      output$CurrentDataStatusHVG = renderUI({
        HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))})
    }else{
      output$CurrentDataStatusHVG = renderUI({
        checkboxGroupInput("sceHVGSelect", "Which data set(s) do you want to include", 
                    choices  = seq(1, length(sceObjs), 1))})
    }
  })
  
  
observeEvent(input$RunHVG, {
  
  inputCheck <- tryCatch(
    input$sceHVGSelect,
    error = function(e) {
      HTML(paste(tags$span(style="color:red", e), sep = ""))
      return(NULL)
    }
  )
  
  if (!is.null(inputCheck)){
  
    
    print(input$sceHVGSelect)
    
    #if multiple data sets -> merge
    if (length(input$sceHVGSelect)>1){
      sceMergedHVG = mergeSCE(as.numeric(input$sceHVGSelect)) #function returns merged data sets
    }
    else{
      sceMergedHVG = sceObjs[[as.numeric(input$sceHVGSelect)]]
    }
    
    if (exists("sceMergedHVG")){
        if (input$hvgSelectInput == "All ST-spots"){
                sce_tmp_hvg_global <<- sceMergedHVG
        }else if(input$hvgSelectInput == "Spatial Eye selection"){
          
          
          output$hvgFollowUpQ = renderUI({
              checkboxGroupInput("hvgFollowUpQCheckbox", "Include selections:",
                                 choices = levels(as.factor(pData(sceMergedHVG)$selectionID)), selected=NULL)
            })
            
              
            
    
        }
    }
  output$HVGParametersSet = renderText({
    print("Parameters Set")
  })
  }
})

observeEvent(input$hvgFollowUpQCheckbox, {
     sce_tmp_hvg_global <<-  sce_tmp_hvg_global[, which(pData(sce_tmp_hvg_global)$selectionID %in% input$hvgFollowUpQCheckbox) ]
    print(sce_tmp_hvg_global)
})
  
    
    
observeEvent(input$trendFitB, {
  
  var.fit <- tryCatch(
    trendVar(sce_tmp_hvg_global, trend="loess", use.spikes=FALSE, span=0.2),
    error = function(e) {
      HTML(paste(tags$span(style="color:red", e), sep = ""))
      return(NULL)
    }
  )
  var.out <- tryCatch(
    decomposeVar(sce_tmp_hvg_global, var.fit),
    error = function(e) {
      HTML(paste(tags$span(style="color:red", e), sep = ""))
      return(NULL)
    }
  )
    
    if (!is.null(var.fit) && !is.null(var.out)){
    
    hvg.out <- var.out[which(var.out$FDR <= input$hvgFDRInput & var.out$bio >= input$hvgBioInput),]
    hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
    nrow(hvg.out)
    
    output$trendFit = renderPlot({
      plotExpression(sce_tmp_hvg_global, rownames(hvg.out)[1:nrow(hvg.out)])+
        theme(axis.title.y = element_text(colour = "orange"),
              axis.text = element_text(colour = "orange"))
         
    }, bg="#2E3337")
    output$HVG_table = renderDataTable({
      hvg.out
    })
    }
})
  
  
  
# ------------------------ Heatmaps -------------
  
observeEvent(input$HEATMAPButton,{
    
    if (length(sceObjs)==0){
      output$CurrentDataStatusHeatmap = renderUI({
        HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))})
    }else{
      output$CurrentDataStatusHeatmap = renderUI({
        checkboxGroupInput("sceHeatSelect", "Which data set do you want to include", 
                    choices  = seq(1, length(sceObjs), 1))})
    }
})
  

observeEvent(input$HeatDataSelectionsParameters, {
  
  print("test")
  
  #If user only want to extract Spatial Eye data
  if (input$HeatDataSelectionsParameters==TRUE){
    #go through the sceObjects and look for different selectionIDs
    selectionID_list = list()
    for (i in 1:length(sceObjs)){
      sce_tmp = sceObjs[[i]]
      selectionID_list[[i]] = unique(pData(sceObjs[[i]])$selectionID)
    }
    extractIDs_SE = unique(unlist(selectionID_list)) 
    
    if (length(extractIDs_SE)==0){
      output$HeatParameterSpatialEye = renderUI({
        HTML(paste(tags$span(style="color:red", "OBS: No selection performed"), sep = ""))})
    }else{
      output$HeatParameterSpatialEye = renderUI({
        checkboxGroupInput("HeatParameterSpatialEyeInput", "Include the following selections from SpatialEye", 
                           choices  = extractIDs_SE[order(extractIDs_SE)])})
        }
    }
  })

 
 
observeEvent(input$HeatMapPlotButton, {
  
  output$topVarHeatmap = renderPlot({
    isolate({
    print(input$sceHeatSelect)
    #if multiple data sets -> merge
    if (length(input$sceHeatSelect)>1){
      sce <- mergeSCE(as.numeric(input$sceHeatSelect)) #function returns merged data sets
    }
    else{
      sce <- sceObjs[[as.numeric(input$sceHeatSelect)]]
    }
      
      
    #Extract only the spots that belong to parameter selections
    selections_SE <- input$HeatParameterSpatialEyeInput
    print(selections_SE)
    
    if (length(selections_SE)>0){
      print("Only Spatial eye set")
      index_extract = pData(sce)$selectionID %in% selections_SE
      sce <- sce[, index_extract] 
    }
   
    
    #Heatmap Type (only most variable genes avaiable right now)
    if (input$HeatMapTypeRadioButtons == "Most Variable Genes"){
    
    #Type of data input for heatmap
    count_type= input$HeatCountTypeInput
    if (count_type=="counts"){
      topVarGenes <- head(order(rowVars(counts(sce)),decreasing=TRUE),input$HeatNrOfGenes)
      mat <- counts(sce)[ topVarGenes, ]
    }else if(count_type=="norm_exprs"){
      topVarGenes <- head(order(rowVars(exprs(sce)),decreasing=TRUE),input$HeatNrOfGenes)
      mat <- exprs(sce)[ topVarGenes, ]
    }
    
    print(length(input$HeatAnnotationTypeInput))
    #Annotation options
    #If user have not selected any annotation -> just use dataset as annotation
    if (length(input$HeatAnnotationTypeInput) < 1){
      annotation = "DataSet"
    }else{
    annotation = input$HeatAnnotationTypeInput
    }
    mat <- mat - rowMeans(mat)
    df <- as.data.frame(as.character(pData(sce)[, annotation]))
    rownames(df)=colnames(sce)
    colnames(df)=annotation
    heat = pheatmap(mat, show_colnames = input$HeatColNames, show_rownames = input$HeatRowNames,
                    cluster_rows = input$HeatClusterRows, 
                    cluster_cols = input$HeatClusterCols, annotation=df,
                    fontsize = input$HeatFontSize)
    }
    })
    },bg = "#2E3337")
  

  
})

# ------------------- DIM REDUCTION SECTION ------
  
observeEvent(input$DIMRedButton,{
    
    if (length(sceObjs)==0){
      output$CurrentDataStatusDimRed = renderUI({
        HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))})
    }else{
      output$CurrentDataStatusDimRed = renderUI({
        checkboxGroupInput("DimRedSelect", "Which data set(s) do you want to include", 
                           choices  = seq(1, length(sceObjs), 1))})
    }
  })

#Let user filter out selections from spatial eye
observeEvent(input$DimRedSelect, {
  sceObjs = sceObjs[as.numeric(input$DimRedSelect)]
  selectionID_list = list()
  for (i in 1:length(sceObjs)){
    if ("selectionID" %in% names(pData(sceObjs[[i]]))){
    #go through the sceObjects and look for different selectionIDs
    sce_tmp = sceObjs[[i]]
    selectionID_list[[i]] = unique(pData(sceObjs[[i]])$selectionID)
    }
  }
    
    extractIDs_SE = unique(unlist(selectionID_list)) 
    
    if (length(extractIDs_SE)==0){
      output$parameterSpatialEyeDimRed = renderUI({
        HTML(paste(tags$span(style="color:red", "OBS: No selection performed"), sep = ""))})
    }else{
      output$parameterSpatialEyeDimRed = renderUI({
        checkboxGroupInput("parameterSpatialEyeDimRedInput", "Optional: Only plot the following selections from SpatialEye", 
                           choices  = extractIDs_SE)})
    }
})
  
  
output$DimPlotType = renderUI({
    selectInput("DimPlotType", "Select dimensionality reduction method", c("PCA", "tSNE"))
  })

#If user wants to do tSNE, also add perplexity option
observeEvent(input$DimPlotType, {
  if (input$DimPlotType == "tSNE"){
    output$PlotOpt1 = renderUI({
      numericInput("perplexity2", "Perplexity: ", value=10)
    })
  }
  if (input$DimPlotType == "PCA"){
    output$PlotOpt1 = renderUI({
      numericInput("ncompInput", "Number of PCs to show:", value=2)
    })
  }
})

output$colorList = renderUI({
    selectInput("colorList", "Color by: ", c("Total counts" = "total_counts", 
                                             "Log total counts" = "log10_total_counts", 
                                             "Spearman cluster" = "spearmanCluster", 
                                             "Spatial Eye Selection" = "selectionID",
                                             "GeneID" = "GeneID", 
                                             "Data set" = "Data set"))
  })
  
#If user wants to color by geneID  
observeEvent(input$colorList, {
  if (input$colorList == "GeneID"){
      output$geneID_color_dimRed = renderUI({
      textInput("geneIDcolorDimRedInput", "GeneID", placeholder = "GeneID")
      })
  }
})
    
observeEvent(input$dimRedInputType, {
  if (input$dimRedInputType == "Normalized counts"){
    output$dimRedInputTypeNorm = renderUI({
      selectInput("dimRedInputTypeNormSelect", "Type of normalization", choices = c("Scran", "CPM"))
    })
  }else{output$dimRedInputTypeNorm = renderUI({HTML(paste(tags$span("", sep = "")))})}
})
  
observeEvent(input$dimPlot, {
  
  inputCheck <- tryCatch(
    input$DimRedSelect,
    error = function(e) {
        HTML(paste(tags$span(style="color:red", e), sep = ""))
      return(NULL)
    }
  )
  
  if (!is.null(inputCheck)){

  #if multiple data sets -> merge
  if (length(input$DimRedSelect)>1){
    sce <- mergeSCE(as.numeric(input$DimRedSelect)) #function returns merged data sets
  }
  else{
    sce <- sceObjs[[as.numeric(input$DimRedSelect)]]
  }
  #If the user have checked any of the Optional "Only plot the following selections from SpatialEye checkboxes"
  selections_SE <- input$parameterSpatialEyeDimRedInput
  if(length(selections_SE)>0){
    print("Only Spatial eye set")
    index_extract = pData(sce)$selectionID %in% selections_SE
    sce <- sce[, index_extract] #This sce object only contains the ST-Spots of interest
  }
 
  output$DimRedPlot <- renderPlot({
   
    count_type = input$dimRedInputType #norm_exprs is created during scran normalization 
   
    if (input$colorList == "GeneID"){
      print(input$geneIDcolorDimRedInput)
      gene_idx = which(rownames(sce) == input$geneIDcolorDimRedInput)
      gene_id_counts = pData(sce[which(rownames(sce) == input$geneIDcolorDimRedInput),])$total_counts
      pData(sce)$GeneIDCounts = gene_id_counts
      colorBy = "GeneIDCounts"
    }else if (input$colorList == "Data set"){
      colorBy = "DataSet" 
    }else{colorBy = input$colorList
    }
    
    
    
    if (input$DimPlotType == "PCA"){

    ncomp = input$ncompInput
    fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
    outPlot <- plotPCA(sce, exprs_values=count_type, colour_by=colorBy, ncomponents=ncomp) + #coord_equal() +
    ggtitle("PCA plot")+
      theme(plot.background = element_rect(fill="#2E3337"),
            axis.title.x = element_text(colour = "orange"),
            axis.ticks = element_blank(),
            axis.title.y = element_text(colour = "orange"),
            legend.title = element_blank(),
            legend.text = element_text(colour="orange", size=5),
            title = element_text(colour="orange"),
            panel.grid = element_blank(),
            panel.grid.major = element_blank())
    
    outPlot <<- outPlot #save for brush function
    
    outPlot
    
    #Also plot a screeplot
    if (input$ScreePlotCheckBox == TRUE){
      if (count_type=="counts"){pc = prcomp(t(counts(sce)))}else{pc = prcomp(t(exprs(sce)))}
    output$screePlot = renderPlot({
      screeplot(pc, type="line", col="orange", main="", npcs=6)
    }, bg = "#2E3337")
    }
    
    }else if (input$DimPlotType == "tSNE"){
      fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
      outPlot <- plotTSNE(sce, exprs_values=count_type, perplexity=input$perplexity2, colour_by=colorBy) + 
        fontsize + ggtitle(paste("Perplexity = ", input$perplexity2, sep="")) +
        theme(plot.background = element_rect(fill="#2E3337"),
              axis.title.x = element_text(colour = "orange"),
              axis.ticks = element_blank(),
              axis.title.y = element_text(colour = "orange"),
              legend.title = element_blank(),
              legend.text = element_text(colour="orange", size=5),
              title = element_text(colour="orange"),
              panel.grid = element_blank(),
              panel.grid.major = element_blank())
      
      outPlot <<- outPlot #save for brush function
    }
    
    outPlot #plot
    
    }, bg="transparent")
  
  #Also plot array plots below
  if (input$ArrayPlotCheckBox == TRUE){
  output$ArrayPlotBelowDimRed <- renderPlot({
    
    if (length(input$DimRedSelect)>1){
      tmp <- sceObjs[as.numeric(input$DimRedSelect)]
    }
    else{
      tmp <- list(sceObjs[[as.numeric(input$DimRedSelect)]])
    }
    
    nrPlots = length(tmp)
    plotListArray = list()
    
    for (i in 1:nrPlots){
      sce = tmp[[i]]
      #If the user have checked any of the Optional "Only plot the following selections from SpatialEye checkboxes"
      selections_SE <- input$parameterSpatialEyeDimRedInput
      if(length(selections_SE)>0){
        index_extract = pData(sce)$selectionID %in% selections_SE
        sce <- sce[, index_extract] #This sce object only contains the ST-Spots of interest
      }
     
      colorBy = input$colorList
      
      if (colorBy == "total_counts"){
        ggColor = pData(sce)$"total_counts"
      }else if(colorBy == "log10_total_counts"){
        ggColor = pData(sce)$"log10_total_counts"
      }else if(colorBy == "Spearman cluster"){
        ggColor = as.factor(pData(sce)$"spearmanCluster")
      }else if(colorBy == "selectionID"){
        ggColor = as.factor(pData(sce)$"selectionID")
      }else if(colorBy == "Data set"){
        ggColor = as.factor(pData(sce)$"DataSet")
      }else if(colorBy == "GeneID"){
        gene_idx = which(rownames(sce) == input$geneIDcolorDimRedInput)
        gene_id_counts = pData(sce[which(rownames(sce) == input$geneIDcolorDimRedInput),])$total_counts
        pData(sce)$GeneIDCounts = gene_id_counts
        ggColor = pData(sce)$"GeneIDCounts"
      } 

      

      gg <- ggplot(data=pData(sce), aes(x=as.numeric(coord_x), y=-as.numeric(coord_y), colour=ggColor)) + 
        xlab("X coordinate")+ylab("Y coordinate")+ coord_equal() +
        geom_point(size = 3, alpha = 0.4 ) +
        theme(
          legend.background = element_blank(),
          panel.background =element_blank(),
          plot.background = element_rect(fill="#2E3337"),
          axis.title.x = element_text(colour = "orange"),
          axis.ticks = element_blank(),
          axis.title.y = element_text(colour = "orange"),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(colour="orange", size=5),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "bottom")
      plotListArray[[i]] = gg
    }
    
    plot_grid(plotlist=plotListArray)
  }, bg = "transparent")
  }
}
})

#Add feature to remove selectionID after plotting DimRed
observeEvent(input$dimPlot, {
  #If the user have checked any of the Optional "Only plot the following selections from SpatialEye checkboxes"
  selections_SE <- input$parameterSpatialEyeDimRedInput
  if(length(selections_SE)>0){
  output$RemoveSelectionDimRed = renderUI({
    actionButton("RemoveSelectionIDfromDimRed", "Remove selected ST-spot(s) from selectionID")
  })
  }
})


output$DimRed_brush_points = renderPrint({
  brush_DimRed <<- tryCatch(
    brushedPoints(outPlot$data,input$DimRed_brush),
    error = function(e) {
        HTML(paste(tags$span(style="color:red", e), sep = ""))
      return(NULL)
    }
  )
    print(rownames(brush_DimRed))
})

# ++++++++++++++ Remove brushed points from current selection group
observeEvent(input$DIMRedButton,{
  
  if (length(sceObjs)==0){
    output$DimRedSelect2 = renderUI({
      HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))})
  }else{
    output$DimRedSelect2 = renderUI({
      radioButtons("DimRedSelect2", "Which data set do you want to include", 
                         choices  = seq(1, length(sceObjs), 1))})
  }
})

observeEvent(input$RemoveSelectionIDfromDimRed, {
  sce = sceObjs[[as.numeric(input$DimRedSelect2)]]
  print(sce)
  remove_idx = which(colnames(sce) %in% rownames(brush_DimRed))
  print(rownames(brush_DimRed))
  print(remove_idx)
  selection_group = pData(sce)$selectionID
  selection_group[remove_idx] = "0"
  pData(sce)$selectionID = selection_group
  sceObjs[[as.numeric(input$DimRedSelect2)]] <<- sce #save to global after removal
  output$removedOutputDimRed = renderPrint({cat("You just removed:", length(remove_idx), " points from selection")})
})



  
# ADD ABLITY TO DOWNLOAD PLOTS !  

  # output$downloadDimPlot <- downloadHandler(
  #   filename = "Dimplot.png",
  #   content = function(file) {
  #     device <- function(..., width, height) {
  #       grDevices::png(..., width = width, height = height,
  #                      res = 300, units = "in")
  #     }
  #     ggsave(file, plot = outPlot, device = device)
  #   })
  
  
  
# ------------------ DE SECTION ----------------------------------------------------------------------
  
  output$DEPoolStrategyInfoText = renderUI({
    HTML("<p> The pooling strategy works in the following manner: </p>
          <p> After identification of two subsets of ST-Spots via the Spatial Eye,
              the spots within each subset are pooled togheter and compared 
              against eachother for DE-analysis. Moreover, the user can define a
              sampling number wihtin each subset, this will result in a random
              sub-sampling of ST-Spots within each pool. These subsets will be treated as
              biological replicates within the selection-group </p>
          <p> For example, if the user have two selections of ST-Spots, and uses a
              \"random-subset\" number of 3, three sub-populations within each selection will be
              created at random, i.e. we will have 6 samples in total and during DE, the analysis
              will compare 3 biological \"psuedo-replicates\".</p>
          <p> By setting \"set.seed()\" to a specific number, the user can redo the same 
              random sampling procedure multiple times. </p>
          <p> </p>
          <p> This could be an alternative if for example your ST-spots are of low quality
              and treating them as individual samples not appropiate </p>")
  })



# 1. Data set(s) config ==================
#uiOutput for available data set(s)
observeEvent(input$DEButton,{
  
  if (length(sceObjs)==0){
    output$CurrentDataStatusDE = renderUI({
      HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))})
  }else{
    output$CurrentDataStatusDE = renderUI({
      checkboxGroupInput("sceDESelect", "Which data set(s) do you want to include", 
                         choices  = seq(1, length(sceObjs), 1))})
  }
})
#uiOutput for available selections
observeEvent(input$DEButton,{
  
  #Depending on the parameter selection, show the alternatives to include
  observeEvent(input$DEDataSelectionsParameters, {
  print(input$DEDataSelectionsParameters)
    
    if ("Spatial Eye Selection" %in% input$DEDataSelectionsParameters){
      #go through the sceObjects and look for different selectionIDs
      selectionID_list = list()
      for (i in 1:length(sceObjs)){
        sce_tmp = sceObjs[[i]]
        selectionID_list[[i]] = unique(pData(sceObjs[[i]])$selectionID)
      }
      extractIDs_SE = unique(unlist(selectionID_list)) 
      
      if (length(extractIDs_SE)==0){
        output$parameterSpatialEyeDE = renderUI({
          HTML(paste(tags$span(style="color:red", "OBS: No selection performed"), sep = ""))})
      }else{
        output$parameterSpatialEyeDE = renderUI({
          checkboxGroupInput("parameterSpatialEyeDEInput", "Include the following selections from SpatialEye", 
                             choices  = extractIDs_SE[order(extractIDs_SE)])})
      }
    }
    if ("Spearman Clusters" %in% input$DEDataSelectionsParameters){
      #go through the sceObjects and look for different clusters
      clustersID_list = list()
      for (i in 1:length(sceObjs)){
        sce_tmp = sceObjs[[i]]
        clustersID_list[[i]] = unique(pData(sceObjs[[i]])$spearmanCluster)
      }
      extractIDs_Cluster = unique(unlist(clustersID_list)) 
      
      if (length(extractIDs_Cluster)==0){
        output$parameterClusterDE = renderUI({
          HTML(paste(tags$span(style="color:red", "OBS: No clustering performed"), sep = ""))})
      }else{
        output$parameterClusterDE = renderUI({
          checkboxGroupInput("parameterClusterDEInput", "Include the following clusters", 
                             choices  = extractIDs_Cluster[order(extractIDs_Cluster)])})
      }
    }
  
  })
})

observeEvent(input$dataConfigDE, {
  inputCheck <- tryCatch(
    input$sceDESelect,
    error = function(e) {
      HTML(paste(tags$span(style="color:red", e), sep = ""))
      return(NULL)
    }
  )
  
  if (!is.null(inputCheck)){
  #if multiple data sets -> merge
   if (length(input$sceDESelect)>1){
     sce_DE <<- mergeSCE(as.numeric(input$sceDESelect)) #function returns merged data sets
   }
   else{
     sce_DE <<- sceObjs[[as.numeric(input$sceDESelect)]]
   }
  
  #Extract only the spots that belong to parameter selections
  selections_SE <- input$parameterSpatialEyeDEInput
  selections_Cluster <- input$parameterClusterDEInput
  
  if (length(selections_SE)>0 & length(selections_Cluster)>0){
    print("SpatialEye + Spearman set")
    index_extract = pData(sce_DE)$selectionID %in% selections_SE & pData(sce_DE)$spearmanCluster %in% selections_Cluster
  }else if(length(selections_SE)>0){
    print("Only Spatial eye set")
    index_extract = pData(sce_DE)$selectionID %in% selections_SE
  }else if(length(selections_Cluster)>0){
    print("Only Spearman set")
    index_extract = pData(sce_DE)$spearmanCluster %in% selections_Cluster
  }else{
    #use all spots
    index_extract = 1:dim(sce_DE)[[2]]
  }
  
  
  sce_DE <<- sce_DE[, index_extract] #This sce_DE object only contains the ST-Spots of interest
  
  output$dataConfigDEConfirm = renderText(print("Data config confirmed"))

  } 
})


# 2. Strategy config =================
output$DEPoolStrategyUiOutput1 = renderUI({ 
    if (input$DEPoolStrategy == TRUE){
    numericInput("DEPoolStrategySetSeed", "Set seed:", value=1)
    }else({return(NULL)})
})

output$DEPoolStrategyUiOutput2 = renderUI({ 
  if (input$DEPoolStrategy == TRUE){
    numericInput("DEPoolStrategyPoolSize", "Pool Size:", value=5)
  }else({return(NULL)})
})

output$DEPoolStrategyUiOutput3 = renderUI({ 
  if (input$DEPoolStrategy == TRUE){
    selectInput("DEPoolStrategyParameter", "Pool based on", choices = c("Spatial Eye Selection", "Spearman Clusters", "Data Sets"))
  }else({return(NULL)})
})


#Set strategy button performs pooling of ST-Spots if wanted
observeEvent(input$SetStrategy, {
  
  #If pooling strategy:
  if (input$DEPoolStrategy == TRUE){
    withProgress(message ="Pooling started ...", { #Progress bar START
      
      poolSize = input$DEPoolStrategyPoolSize #User defined pool size
      
      #Extract individual objects for pooling
      selections_SE <- input$parameterSpatialEyeDEInput
      print(paste("selections_SE:", selections_SE))
      selections_Cluster <- input$parameterClusterDEInput
      dataSet_cluster <- input$sceDESelect
      sceSE_list = list()
      print(input$DEPoolStrategyParameter)
      # if pooling is based on Spatial Eye
      if (input$DEPoolStrategyParameter == "Spatial Eye Selection"){
        print("Pool based on Spatial EYE")
        IDs = selections_SE
        #Create sceSE_list of sce objects so that sce nr1 = sce that belong to first selection group etc.
        for (i in 1:length(selections_SE)){
          sceSE_list[[i]] = sce_DE[, which(pData(sce_DE)$selectionID %in% selections_SE[i])]
        }
        
      }else if(input$DEPoolStrategyParameter == "Spearman Clusters"){
        print("Pool based on Spearman Clusters")
        IDs = selections_Cluster
        print(selections_Cluster)
        for (i in 1:length(selections_Cluster)){
          sceSE_list[[i]] = sce_DE[, which(pData(sce_DE)$spearmanCluster %in% selections_Cluster[i])]
        }
      }else if(input$DEPoolStrategyParameter == "Data Sets"){
        print("Pool based on Data sets")
        IDs = dataSet_cluster
        for (i in 1:length(dataSet_cluster)){
          sceSE_list[[i]] = sce_DE[, which(pData(sce_DE)$DataSet %in% dataSet_cluster[i])]
        }
      }
      
      #We now have a list of SCE object that are input for pooling
      test <<- sceSE_list
      print(paste("sceSE_list: ", sceSE_list))
      print(paste("poolsize: ", poolSize))
      pooled_sce_list = list()
      #Loop through the sce objects and create pools
      for (i in 1:length(sceSE_list)){
        
        split_SE = split(colnames(sceSE_list[[i]]), sample(1:poolSize)) #This is a list of ST-spot names
        print(paste("split_SE:", split_SE))
        
        #pool the counts of the first pool list (collection of ST-Spot names)
        pool_counts = counts(sce_DE[, which(colnames(sce_DE) %in% split_SE[[1]]) ]) #OBS we take raw counts here!
        pool_counts = rowSums(pool_counts)
        #Then go through the rest of the pool list (collection of ST-Spot names) and pool these
        for (j in 2:length(split_SE)){
         pool_counts = cbind(pool_counts, rowSums(counts(sce_DE[, which(colnames(sce_DE) %in% split_SE[[j]]) ]))) #This count matrix is now considered one sample
        }
        
        pool_counts = as.data.frame(pool_counts)
        
        #print(pool_counts)
        colnames(pool_counts) = paste("DE-", IDs[i], "-Sample_", seq(1:ncol(pool_counts)), sep="")
        pooled_sce_list[[i]] = pool_counts
        pooled_sce_list[[i]]$ID = rownames(pooled_sce_list[[i]]) #add column for later joining
        

      }
      
      #Check if user looses rows
      max_rows_before = max(unlist(lapply(pooled_sce_list, FUN  = function(x){max(nrow(x))})))
      
      #Join all into one data frame
      DE_df = plyr::join_all(pooled_sce_list, by="ID", type="full")
      #This will introduce NA for the rows that do not match
      #Remove these for downstream DE analysis
      DE_df = na.omit(DE_df)
      rownames(DE_df) = DE_df$ID
      DE_df$ID = NULL
      
      #Skapa ett SCE objekt fran denna lista 
      
      DE_sce = newSCESet(countData = DE_df)
      
      split = strsplit(colnames(DE_sce), "-")
      split = sapply(split, "[[", 2)
      pData(DE_sce)$Pools = split
      DE_sce <<- DE_sce #Save the DE-sce object global for input to DESeq2/EdgeR functions
      
      #rows after (to check if lost rows)
      rows_after = nrow(DE_sce)
      lost_rows = max_rows_before - rows_after
      if (lost_rows>0){
        output$DEMessegeBoxLostGenes = renderUI({
          verbatimTextOutput(paste("During DE analysis you lost ", lost_rows, " genes that was not present in all data sets"))
        })
      }
      
      incProgress(detail = paste("Pooling completed!")) #Progress bar END
      Sys.sleep(1)})
    print("Pooling Done")
  }else{ #if no pooling should be performed
    DE_sce <<- sce_DE
  }
  
  output$SetDEStrategyConfirm = renderText(print("Strategy Set"))
  
  print("Strategy Set")
  
  output$DEDesignOutput = renderUI({
    namesVector = names(pData(DE_sce))
    namesVector = namesVector[which(namesVector %in% c("selectionID", "spearmanCluster", "DataSet", "Pools", "AnnotationNameDataSet"))]
    checkboxGroupInput("useDEDesign", "Make linear combinations with parameters:", choices=namesVector)
  })
  
  output$DEControlTerm = renderUI({
    namesVector = names(pData(DE_sce))
    namesVector = namesVector[which(namesVector %in% c("selectionID", "spearmanCluster", "DataSet", "Pools", "AnnotationNameDataSet"))]
    checkboxGroupInput("DEInteractionTerm", "Add interaction terms with parameters:", choices=namesVector)
  })
  
  output$DEVariableOfInterestUI = renderUI({
    namesVector = names(pData(DE_sce))
    namesVector = namesVector[which(namesVector %in% c("selectionID", "spearmanCluster", "DataSet", "Pools", "AnnotationNameDataSet"))]
    selectInput("DEVariableOfInterest", "Choose variable of interest:", choices=namesVector)
  })
  
  })

#Can only handle 2 interaction terms atm
observe({
  if (length(input$DEInteractionTerm) > 2){
    
    updateCheckboxGroupInput(session,"DEInteractionTerm", selected = tail(input$DEInteractionTerm, 2) )
  }
})

#Can only handle 3 variables of interest atm
observe({
  if (length(input$useDEDesign) > 3){
    updateCheckboxGroupInput(session,"useDEDesign", selected = tail(input$useDEDesign, 3) )
  }
})

observeEvent(input$useDEDesign, {
  output$SetRefDE = renderUI({
    choice = input$DEVariableOfInterest
    choices = unique(pData(DE_sce)[[choice]])
    choices = choices[order(choices)]
    radioButtons("setDERef", "Set reference", choices=choices)
  })
})

# 3. Running DE =======================

observeEvent(input$DESeq2RUN,{
  
    
  
    withProgress(message ="Running DESeq2 ...", { #Progress bar START
    
      dds <- tryCatch(
        convertTo(DE_sce, type="DESeq2", pData.col = input$useDEDesign),
        error = function(e) {
            HTML(paste(tags$span(style="color:red", e), sep = ""))
          return(NULL)
        }
      )
  
    if(!is.null(dds)){
      
    VariablesToUse = input$useDEDesign
    print(paste("VariableToUse: ", VariablesToUse))
    
    VariableOfInterest = input$DEVariableOfInterest
    VariablesToUse = as.factor(VariablesToUse)
    VariablesToUse <- tryCatch(
      relevel(VariablesToUse, VariableOfInterest),
      error = function(e) {
        HTML(paste(tags$span(style="color:red", e), sep = ""))
        return(NULL)
      }
    )
    if(!is.null(VariablesToUse)){

    #convert to factors for DE design
    for (i in 1:length(input$useDEDesign)){
      colData(dds)[i] = as.factor(colData(dds)[[i]])
    }
    #if(length(input$DEInteractionTerm) < 1){
    #If no interaction terms
    print("No interaction terms used")
      if(length(input$useDEDesign)==2){
          print("2 parameters in design")
          design = paste("~", VariablesToUse[2],  paste("+", VariablesToUse[1], collapse="+"), sep="")
          design(dds) = formula(design)
          print(design(dds))
      }else if(length(input$useDEDesign)==3){
        print("3 parameters in design")
        design = paste("~", VariablesToUse[3],  paste("+", VariablesToUse[1], VariablesToUse[2], collapse="+"), sep="")
        design(dds) = formula(design)
        print(design(dds))
      }else if(length(input$useDEDesign)==1){
        print("1 parameter in design")
        design = paste("~", paste(VariablesToUse, collapse = "+"), sep="")
        
        design(dds) = formula(design)
        print(design(dds))
      }
    
    if(length(input$DEInteractionTerm)>0){
    #If interactionterms is wanted
    print("Adding interaction term")
      design = paste(design, " +", paste(input$DEInteractionTerm[1], input$DEInteractionTerm[2], sep=":"), sep="")
      #dds$group = factor(paste0(dds[[input$DEInteractionTerm[1]]], dds[[input$DEInteractionTerm[2]]]))
      design(dds) = formula(design)
      print(design(dds))
    }
    
    print(paste("First level: ", input$setDERef))
    
    colData(dds)[[as.character(input$DEVariableOfInterest)]] = relevel(colData(dds)[[as.character(input$DEVariableOfInterest)]], input$setDERef)
    dds = tryCatch(
      DESeq(dds),
      error = function(e) {
        output$DEDesignMessegeBox = renderUI({
         HTML(paste(tags$span(style="color:red", e), sep = ""))
        })
        return(NULL)
      }
    )
    print(dds)
    if (!is.null(dds)){
      output$DEContrastMessegeBox = renderUI({
        HTML(paste(resultsNames(dds), sep=""))
      })
    
      output$DEDesignMessegeBox = renderUI({
        HTML(paste(design(dds), sep="" ))
      })
    
      dds_global <<- dds
    
    
      output$DEContrastChooser1 <- renderUI({
        checkboxGroupInput("DEConrastCheckBoxGroupInput1", "", choices=resultsNames(dds))
      })
      observeEvent(input$DEConrastCheckBoxGroupInput1, {
        print(input$DEConrastCheckBoxGroupInput1)
      })
    }
    }
    }
    incProgress(detail = paste("DE completed!")) #Progress bar END
    Sys.sleep(1)})
    print("DE Done")
    }
)


# 4. Show results =================
#let user choose the contrast that are wanted from the DE run
#Contrasts are found in resultsNames(dds)
#these contrasts are output in the messege box in the shiny app

observeEvent(input$ResultsDE, {
  t1 = input$DEConrastCheckBoxGroupInput1
  #t2 = input$DEConstrastCheckBoxGroupInput2
  if (length(t1)<1){
    return(NULL)
  }else{
    res = results(dds_global, contrast=list(t1))
  }
  
  DESeq2_var <<- 1
  res_global <<- res #save global for plotting etc
  
 
  #output a result table
  res_tableOutput <<- as.data.frame(res)
  res_tableOutput$GeneID = rownames(res)
  res_tableOutput = res_tableOutput[c(7,1,2,3,4,5,6)]
  
    output$DETableOutput =  renderUI({
        output$aa <- renderDataTable(res_tableOutput)
        dataTableOutput("aa")
  })
    #Show the download table button handle
    output$downloadDETableButton =
      downloadHandler(
        
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
          paste("DE_result_table", "csv", sep = ".")
        },
        
        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
          res_table_download = as.data.frame(res_tableOutput)
          # Write to a file specified by the 'file' argument
          write.csv(res_table_download, file, sep = ",",
                    row.names = TRUE)
        }
      )
    
 
})

# observeEvent(input$DE_Parameter_selection, {
#     
#       if (input$DE_Parameter_selection=="Spearman Clusters"){
#             
#               output$DE_Parameter_selection_part2 = renderUI({
#               selectInput("ClusterDE1", "Cluster1:",c=levels(pData(sce_DE)$cluster))})
#               output$DE_Parameter_selection_part3 = renderUI({
#                 selectInput("ClusterDE2", "Cluster2:",c=levels(pData(sce_DE)$cluster))})
#       }else if (input$DE_Parameter_selection == "Spatial Eye selection"){
#         
#         #Selection from Spatial Eye stored in pData(sce_DE)$selectionID
#         selectionLevels = levels(as.factor(pData(sce_DE)$selectionID))
#         
#         output$DE_Parameter_selection_part2 = renderUI({
#           numericInput("SelectionDE1", label="Selection1:",choices= selectionLevels)})
#         output$DE_Parameter_selection_part3 = renderUI({
#           numericInput("SelectionDE2", label="Selection2:",choices=selectionLevels)})
#       }
#  
# })


observeEvent(input$RunEdgeR, {
    
    
    #If pooling strategy:
    if (input$DEPoolStrategy == TRUE){
      withProgress(message ="Pooling started ...", { #Progress bar START
      
      c1 = as.character(input$ClusterDE1)
      c2 = as.character(input$ClusterDE2)
      c1_SE = as.character(input$SelectionDE1)
      c2_SE = as.character(input$SelectionDE2)
      
      if (input$DE_Parameter_selection=="Spearman Clusters"){
        
        set.seed(input$DEPoolStrategySetSeed)
        
        sceE1 = sce[, which(pData(sce)$cluster %in% c1)]
        sceE2 = sce[, which(pData(sce)$cluster %in% c2)]
        
      }else if(input$DE_Parameter_selection == "Spatial Eye selection"){
        
        set.seed(input$DEPoolStrategySetSeed)
        
        sceE1 = sce[, which(pData(sce)$selectionID %in% c1_SE)]
        sceE2 = sce[, which(pData(sce)$selectionID %in% c2_SE)]
      }
      
      pool1_1bc = sample(colnames(sceE1), round(length(colnames(sceE1))/2))
      pool1_2bc = colnames(sceE1[, -which(colnames(sceE1) %in% pool1_1bc)])
      
      pool2_1bc = sample(colnames(sceE2), round(length(colnames(sceE2))/2))
      pool2_2bc = colnames(sceE2[, -which(colnames(sceE2) %in% pool2_1bc)]) 
      
      pool1_1_counts=rowSums(counts(sceE1[, which(colnames(sceE1) %in% pool1_1bc)]))
      pool1_2_counts=rowSums(counts(sceE1[, which(colnames(sceE1) %in% pool1_2bc)]))
      
      pool2_1_counts=rowSums(counts(sceE2[, which(colnames(sceE2) %in% pool2_1bc)]))
      pool2_2_counts=rowSums(counts(sceE2[, which(colnames(sceE2) %in% pool2_2bc)]))
      
      pool = as.data.frame(cbind(pool1_1_counts, pool1_2_counts, pool2_1_counts, pool2_2_counts ))
      
      
      sce_pooled <- newSCESet(countData=pool)
      sce_pooled <- calculateQCMetrics(sce_pooled)
      
      pData(sce_pooled)$group = c("1","1","2","2")
      
      de.design = model.matrix(~0 + pData(sce_pooled)$group)
      colnames(de.design) = c("First","Second")
      yE = convertTo(sce_pooled, type="edgeR", group = pData(sce_pooled)$group)
      
      
      incProgress(detail = paste("Pooling completed!")) #Progress bar END
      Sys.sleep(1)})
      print("Pooling Done")
      
    }else if(input$DEPoolStrategy == FALSE){

    #if not pooling stratery:
   
    c1 = as.character(input$ClusterDE1)
    c2 = as.character(input$ClusterDE2)
    c1_SE = as.character(input$SelectionDE1)
    c2_SE = as.character(input$SelectionDE2)
    if (input$DE_Parameter_selection=="Spearman Clusters"){
        sceE=sce[, which(pData(sce)$cluster %in% c(c1,c2))]
        pData(sceE)$cluster = droplevels(pData(sceE)$cluster)
        de.design = model.matrix(~0 + pData(sceE)$cluster)
        colnames(de.design) = c("First","Second")
        yE = convertTo(sceE, type="edgeR", group = pData(sceE)$cluster)
    }else if(input$DE_Parameter_selection == "Spatial Eye selection"){
        sceE=sce[, which(pData(sce)$selectionID %in% c(c1_SE,c2_SE))]
        pData(sceE)$selectionID = droplevels(as.factor(pData(sceE)$selectionID))
        de.design = model.matrix(~0 + pData(sceE)$selectionID)
        colnames(de.design) = c("First", "Second")
        yE = convertTo(sceE, type="edgeR", group = pData(sceE)$selectionID)
    }
    }
    
    
    withProgress(message ="Performing DE analysis ...", { #Progress bar START
   
    # Estimate dispersion parameter for GLM
    yE <- estimateGLMCommonDisp(yE, de.design)
    yE <- estimateGLMTrendedDisp(yE, de.design, method="power")
    yE <- estimateGLMTagwiseDisp(yE,de.design)
    
    # Plot mean-variance
    #plotBCV(yE)
    
    # Model fitting
    fit.edgeR <- glmFit(yE, de.design)
    # Differential expression
    contrasts.edgeR <- makeContrasts(First - Second, levels=de.design)
    lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)
    # Access results tables
    edgeR_results <- lrt.edgeR$table
    sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = 0.05)
    genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
    print(head(edgeR_results[order(edgeR_results$PValue), ], 20))
    
    lrt.edgeR$table = lrt.edgeR$table[order(rownames(lrt.edgeR$table)),]
    FDR = as.data.frame(topTags(lrt.edgeR, n=nrow(lrt.edgeR$table)))[,c(1,5)]
    FDR = FDR[order(rownames(FDR)),]
    FDR = FDR[,2]
    lrt.edgeR$table$FDR = FDR
    sig=ifelse(lrt.edgeR$table$FDR<0.01, "xFDR<0.01", "Not Sig")
    lrt.edgeR$table$sig = sig
    lrt.edgeR$table$GeneID = rownames(lrt.edgeR$table)
    lrt.edgeR <<- lrt.edgeR
    res_tableOutput <<- as.data.frame(lrt.edgeR$table)
    res_tableOutput <<- res_tableOutput[, -6]
    res_tableOutput <<- res_tableOutput[, c(6,1,2,4,5)]
  
    
    
    incProgress(detail = paste("DE completed!")) #Progress bar END
    Sys.sleep(1)})
    print("DE Done")
    
    
     output$DETableOutput =  renderUI({
      
       output$aa <- renderDataTable(res_tableOutput)
       dataTableOutput("aa")
       })

     })
  
  # DE-plots ==========================
  
    observeEvent(input$DEPlot_execute, {
    
      output$DEPlotOutput = renderPlot({
        
        if(DESeq2_var == 1){
          table = as.data.frame(res_global)
          colnames(table) = c("logCPM", "logFC", "lfcSE", "stat", "PValue", "FDR")
          table$GeneID = rownames(table)
        }
        else{
          table = lrt.edgeR$table
        }
        
        if (input$DE_plot_alt == "Volcano"){

          sig = ifelse(table$FDR < 0.01, "Sig", "Non-sig" )
          table$sig = sig
          table$col = ifelse(table$sig=="Sig", "red", "gray")
      
          ggVolcano = ggplot(table, aes(x=table$logFC, y=-log10(table$PValue))) +
            geom_point(size=1, col= table$col) +
            coord_cartesian()+
            xlab("LogFoldChange")+ ylab("-log10(pvalue)")+
            theme(
              legend.background = element_blank(),
              panel.background =element_blank(),
              plot.background = element_rect(fill="#2E3337"),
              axis.title.x = element_text(colour = "orange"),
              axis.text = element_text(colour = "#faebd7"),
              axis.title.y = element_text(colour = "orange"),
              legend.title = element_blank(),
              legend.text = element_blank(),
              legend.key = element_blank(),
              panel.grid = element_blank(),
              panel.grid.major = element_blank())
            
          if (input$labelDEpointsInput == TRUE){
            gene_label_list = table[which(table$FDR < input$labelDEpointsInputFDRselect),]$GeneID
            ggVolcano = ggVolcano +
              geom_text_repel(data=table[which(table$GeneID %in% gene_label_list),],
                              aes(x=table[which(table$GeneID %in% gene_label_list),]$logFC,
                                  y=-log10(table[which(table$GeneID %in% gene_label_list),]$PValue),
                                  label=table[which(table$GeneID %in% gene_label_list),]$GeneID), col="white")
          }
          
          ggVolcano
        }else if (input$DE_plot_alt=="MA"){
          
          sig = ifelse(table$FDR < 0.01, "Sig", "Non-sig" )
          table$sig = sig
          table$col = ifelse(table$sig=="Sig", "red", "gray")
        
          outDEPlot <- ggplot(data=table, aes(x=table$logCPM, y=table$logFC))+
            geom_point(size=1, col= table$col) +
            coord_cartesian()+
            geom_hline(yintercept = 2, col="red", linetype="dashed")+geom_hline(yintercept = -2, col="red", linetype="dashed")+
            xlab("Mean Expression") +
            ylab("Log Fold Change")+
            theme(
              legend.background = element_blank(),
              panel.background = element_blank(),
              plot.background = element_rect(fill="#2E3337"),
              axis.title.x = element_text(colour = "orange"),
              axis.title.y = element_text(colour = "orange"),
              axis.text = element_text(colour = "orange"),
              legend.title = element_text(colour = "orange"),
              legend.text = element_text(colour = "orange"),
              panel.grid = element_blank(),
              panel.grid.major = element_blank()
            )
          
          if (input$labelDEpointsInput == TRUE){
            gene_label_list = table[which(table$FDR < input$labelDEpointsInputFDRselect),]$GeneID
            outDEPlot = outDEPlot +
                              geom_text_repel(data=table[which(table$GeneID %in% gene_label_list),],
                                              aes(x=table[which(table$GeneID %in% gene_label_list),]$logCPM,
                                                  y=table[which(table$GeneID %in% gene_label_list),]$logFC,
                              label=table[which(table$GeneID %in% gene_label_list),]$GeneID), col="white")
          }
          
          outDEPlot
        
        }
          })
      
      #Table below
      output$MA_brushinfo =  renderUI({
        
        if (exists("lrt.edgeR")){ #If edgeR was run
          
          ymin=(input$MA_brush$ymin)
          xmin=(input$MA_brush$xmin)
          ymax=(input$MA_brush$ymax)
          xmax=(input$MA_brush$xmax)
          
          if (input$DE_plot_alt=="MA"){
            
            #subset genes within brush area
            sub_E = lrt.edgeR$table[which(lrt.edgeR$table$logCPM>xmin & lrt.edgeR$table$logCPM<xmax & lrt.edgeR$table$logFC > ymin & lrt.edgeR$table$logFC < ymax),]
            sub_E = sub_E[, c(1,2,4,5,7)]
          }else if (input$DE_plot_alt =="Volcano"){
            #subset genes within brush area<
            sub_E = lrt.edgeR$table[which(lrt.edgeR$table$logFC>xmin & lrt.edgeR$table$logFC<xmax & -log10(lrt.edgeR$table$PValue) > ymin & 
                                            -log10(lrt.edgeR$table$PValue) < ymax),]
            sub_E = sub_E[, c(1,4,5,7)]
            
          }
          
          output$brushInfoTable <- renderDataTable((sub_E))
          dataTableOutput("brushInfoTable")
          
        }else if(exists("res_global")){ #If DESeq2 was run
          ymin=(input$MA_brush$ymin)
          xmin=(input$MA_brush$xmin)
          ymax=(input$MA_brush$ymax)
          xmax=(input$MA_brush$xmax)
          
          if (input$DE_plot_alt=="MA"){
            table = as.data.frame(res_global)
            colnames(table) = c("logCPM", "logFC", "lfcSE", "stat", "PValue", "FDR")
            table$GeneID = rownames(table)
            table = table[c(7,1,2,3,4,5,6)]
            #subset genes within brush area
            sub_E = table[which(table$logCPM>xmin & table$logCPM<xmax & table$logFC > ymin & table$logFC < ymax),]
            print(sub_E)
          }else if (input$DE_plot_alt =="Volcano"){
            #subset genes within brush area<
            sub_E = table[which(table$logFC>xmin & table$logFC<xmax & -log10(table$PValue) > ymin & 
                                  -log10(table$PValue) < ymax),]
            
          }
          
          output$brushInfoTable <- renderDataTable((sub_E))
          dataTableOutput("brushInfoTable")
          
          
        }else{return()}
        
      })
      
      output$MA_brush_points = renderPrint({
        if (input$DE_plot_alt=="MA"){
          brush_MA <<- brushedPoints(outDEPlot$data,input$MA_brush)
          str(input$MA_brush)}else if (input$DE_plot_alt=="Volcano"){
            brush_MA <<- brushedPoints(ggVolcano$data,input$MA_brush)
            str(input$MA_brush)
          }
        
      })
    
    })
    
    output$labelDEpointsInputFDRselect = renderUI({
      if (input$labelDEpointsInput == TRUE){
        numericInput("labelDEpointsInputFDRselect", "FDR threshold", min=0e-10, max=1, value=0.01, step=0.001)
      }else{return(NULL)}
    })
    

   
    


  
    
  #-------------------------- Remove Batch Effect (Limma) --------------------
  
    observeEvent(input$batchCorrButton,{
      
      if (length(sceObjs)==0){
        output$CurrentDataStatusBatchCorr = renderUI({
          HTML(paste(tags$span(style="color:red", "OBS: You have not loaded any data sets yet"), sep = ""))})
      }else{
        output$CurrentDataStatusBatchCorr = renderUI({
          checkboxGroupInput("sceBatchCorrSelect", "Which data set do you want to include", 
                             choices  = seq(1, length(sceObjs), 1))})
      }
    
      
    })
    
observeEvent(input$ShowBatchCorr, {
  
        count_type = input$BatchCorrCountType
  
        output$PCABeforeBatchCorr = renderPlot({
          
          
        
        
        
            print(input$sceBatchCorrSelect)
            #if multiple data sets -> merge
            if (length(input$sceBatchCorrSelect)>1){
              sce <- mergeSCE(as.numeric(input$sceBatchCorrSelect)) #function returns merged data sets
            }
            else{
              sce <- sceObjs[[as.numeric(input$sceBatchCorrSelect)]]
            }
            
            #Plot after-PCA
          output$PCAAfterBatchCorr = renderPlot({
          
            batches = pData(sce)$DataSet
          
            #using raw counts as input
            t.corrected.counts = BatchCorrectedCounts(t(counts(sce)), batches, use_parallel = FALSE)
            corrected.counts = t(t.corrected.counts)
          
            exprs(sce) = corrected.counts
          
            corrected.counts <<- corrected.counts
            sce_batch_corrected <<- sce
          
            outPlot <- plotPCA(sce, exprs_values="exprs", colour_by="DataSet") + #coord_equal() +
              ggtitle("PCA plot after correction")+
              theme(plot.background = element_rect(fill="#2E3337"),
                  axis.title.x = element_text(colour = "orange"),
                  axis.ticks = element_blank(),
                  axis.title.y = element_text(colour = "orange"),
                  legend.title = element_blank(),
                  legend.text = element_text(colour="orange", size=5),
                  title = element_text(colour="orange"),
                  panel.grid = element_blank(),
                  panel.grid.major = element_blank())
          
            outPlot
          
          }, bg = "transparent")
          
          
        #Plot before-PCA
          outPlot <- plotPCA(sce, exprs_values=count_type, colour_by="DataSet") + #coord_equal() +
            ggtitle("PCA plot after correction")+
            theme(plot.background = element_rect(fill="#2E3337"),
                  axis.title.x = element_text(colour = "orange"),
                  axis.ticks = element_blank(),
                  axis.title.y = element_text(colour = "orange"),
                  legend.title = element_blank(),
                  legend.text = element_text(colour="orange", size=5),
                  title = element_text(colour="orange"),
                  panel.grid = element_blank(),
                  panel.grid.major = element_blank())
          
          outPlot
        }, bg = "transparent")
        
        output$ApplyBatchCorrButtonUI = renderUI({
          actionButton("ApplyBatchCorrButton", "Replace with batch corrected counts")
        })
        
       
    })
observeEvent(input$ApplyBatchCorrButton, {
  
          count_type = input$BatchCorrCountType

          if (input$BatchCorrCountType == "counts"){
            counts(sce_batch_corrected) = corrected.counts
          }else if(input$BatchCorrCountType == "norm_exprs"){
            exprs(sce_batch_corrected) = corrected.counts
          }
          
          sce = sce_batch_corrected
          #Roll back from merged data set
          datasets = unique(pData(sce)$DataSet)
          print(paste("datasets:", datasets, sep=""))
          for (i in 1:length(unique(pData(sce)$DataSet))){
            sce_tmp = sce[, which(pData(sce)$DataSet==datasets[i])]
            curr_dataset = datasets[i]
            print(paste("curr dataset:", curr_dataset, sep=""))
            #before rolling back we need to "reset" colnames after merge function
            split = strsplit(colnames(sce_tmp), "-")
            colnames(sce_tmp) = sapply(split, "[[", 1)
            sceObjs[[curr_dataset]] <<- sce_tmp #roll back to global sceObjs
          }
          output$ApplyBatchCorrConfirmText = renderText({
            print("Done")
          })
          
    })
    
    
  # ----------------- SELECTOR PLOT FUNCTION -----------------------------------------------------------
  updatePlot = function(x, y){

    # Depending on the plot the user wants ---
    if (input$selectPlot == "PCA"){
    fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
    pca <- plotPCA(sce, exprs_values="exprs") + fontsize
    x <- pca$data$PC1
    y <- pca$data$PC2
    }
    else if (input$selectPlot == "tSNE"){
      perplex = as.numeric(input$perplexity)
      title = paste("Perplexity = ", perplex)
      fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
      tsne <- plotTSNE(sce, exprs_values="exprs", perplexity=perplex, colour_by="total_features"
                       ) + fontsize
      x = tsne$data$X1
      y= tsne$data$X2
      
    }
    else if (input$selectPlot == "Standard Array"){
      x=sce$coord_x
      y=sce$coord_y
    }
    
    
  N <- length(x)
  if (length(y)!=N) { stop("length of x and y vectors should be equal") }
  collected <- new.env()
  resetValues(collected, data.frame(x=x, y=y))
  
  # Internal functions, to avoid passing many arguments around.
  plotFun1 <- function(output) {
    generatePlot1(output, collected, pch=16)
  }
  plotFun2 <- function(output) {
    generatePlot2(output, collected)
  }
  updateSelect <- function(input, output, setting) {
    brushed <- brushedPoints(collected$coords, xvar="x", yvar="y", input$plot1_brush, allRows=TRUE)
    collected$current.selected[brushed$selected_] <- setting
    plotFun1(output)
  }

  plotFun1(output)
  plotFun2(output)
  
  observeEvent(input$select, { updateSelect(input, output, TRUE) })
  observeEvent(input$unselect, { updateSelect(input, output, FALSE) })
  observeEvent(input$clear, {
    collected$current.selected[] <- FALSE
    plotFun1(output)
  })
  
  observeEvent(input$list_add, {
    n <- length(collected$all.selected)
    collected$all.selected[[n+1]] <- collected$current.selected
    collected$old.selected[collected$current.selected] <- TRUE
    collected$current.selected[] <- FALSE
    plotFun1(output)
    plotFun2(output)
  })
  
  observeEvent(input$reset, {
    resetValues(collected)
    plotFun1(output)
    plotFun2(output)
  })
  
  observeEvent(input$finish, {
    keep.val <<- collected$all.selected
  })
  }
  

  
  # ---------------------------------------------------------------------------------------- DOWNLOAD handlers etc -------------------------------------------
  
    output$SaveData = downloadHandler(
      filename =  paste("SCE-", Sys.Date(), ".RData", sep="")
      ,
      content = function(filename) {
        save(sceObjs, filename)
      }
    )
    
  
  
  
  #=================================================================================== SPATIAL EYE ===================================================================================
  
  #leaflet.extras
  
  # The goal of leaflet.extras is to provide extra functionality to the leaflet R package using various leaflet plugins.
  # Installation
  # 
  # # We need latest leaflet package from Github, as CRAN package is too old.
  # devtools::install_github('rstudio/leaflet')
  # devtools::install_github('bhaskarvk/leaflet.extras')
  
  # -------------------------------------------------------- some info for dev leaflet ------------------------
  
  # 
  # Groups and Layer IDs may appear similar, in that both are used to assign a name to a layer. 
  # However, they differ in that layer IDs are used to provide a unique identifier to individual markers and shapes, etc., while groups are used to give shared labels to many items.
  # You generally provide one group value for the entire addMarkers call, and you can reuse that same group value in future addXXX calls to add to that group’s membership (as in the example above).
  # 
  # layerId arguments are always vectorized: when calling e.g. addMarkers you need to provide one layer ID per marker, and they must all be unique. 
  # If you add a circle with a layer ID of "foo" and later add a different shape with the same layer ID, the original circle will be removed.
  
  # ---------------------------------------------------------------------------------------------------------
  
 
  

    
    #acm_defaults <- function(map, x, y) addCircleMarkers(map, x, y, radius=6, color="black", fillColor="orange", fillOpacity=1, opacity=1, weight=2, stroke=TRUE, layerId="Selected")
    observe({ 
      print(input$navbar) 
      if (input$navbar == "Spatial Eye"){
        
       if (length(spotDataObjs)==0){
         output$SpotDataMissing = renderUI({
           
           modalDialog("SpotDataMissing", "Spot Data is missing")
           
         })
       }else{

        
        if (length(spotDataObjs)>0){
          print(input$SpatialEyeDataSetChooserInput)
          
          if (input$SpatialEyeDataSetChooserInput > length(sceObjs)){
            selectedData = length(sceObjs) #avoid user to scroll to Data set that does not exist
            print(length(sceObjs))
          }else{
            selectedData = input$SpatialEyeDataSetChooserInput
          }
          sce = sceObjs[[as.numeric(selectedData)]]
          #just consider the spots which we have data from
          print(paste("dim before: ", dim(sce)))
          spot_data = spotDataObjs[[as.numeric(selectedData)]]
          #Also read in data for img scaling
          spot_data_dims = imgScaleObjs[[as.numeric(selectedData)]]
          spot_data = spot_data[which(spot_data$barcode %in% colnames(sce)),]
          sce = sce[, which(colnames(sce) %in% spot_data$barcode)]
          #sort to read in spot data into sce object
          sce = sce[, order(colnames(sce))] #string sorted
          spot_data = spot_data[order(spot_data$barcode), ]
          #convert to img coordinates
          spot_data = spot_data[, c(3,4,5)] #just use column V3,V4,v5 (V5 = inside/outside tissue - yes/true)
          #pixel-coords are stored in column V3,V4 in output from spot-detection tool
          scalingFactorImg = 256/max(spot_data_dims)
          
          spot_data$V3 = (spot_data$V3*scalingFactorImg)
          spot_data$V4 = (spot_data_dims[2] - spot_data$V4)*scalingFactorImg
          print(head(spot_data))
          colnames(spot_data) = c("x", "y", "barcode")
          pData(sce)$spot_coord_x = as.numeric(spot_data$x)
          pData(sce)$spot_coord_y = as.numeric(spot_data$y)
          selection_group<<-rep(0, length=length(colnames(sce)))
          print(paste("dim after: ", dim(sce)))
          sceObjs[[as.numeric(selectedData)]] <<- sce #save to global sceObjs
          
          #Run spatial eye
          spatialEye <<- TRUE #not used atm
          
        }else{
          output$SpotDataMissing = renderUI({
            
            modalDialog("SpotDataMissing", "Spot Data is missing")
            
          })
         #just temporary for testing
         x_seq = rep(seq(from=40, to=210, 10), 50)
         x_seq = x_seq[1:ncol(sce)]
         y_seq = rep(seq(from=70, to=220, 10), 50)
         y_seq = x_seq[1:ncol(sce)]
         y_seq=sample(y_seq)
         x_seq=sample(x_seq)
         pData(sce)$coord_x = x_seq
         pData(sce)$coord_y = y_seq
         selection_group<<-rep(0, length=length(colnames(sce)))
         sce <<- sce #save new global var
        } 
      
      
    
    # ------------------------------------------------------ Plot the Default tissue "map" -------------------------------------------------------------------------------------
    
    # if (spatialEye == FALSE){
    #   print("Not Run")
    #   
    # }else{
    
    output$Map <- renderLeaflet({
      
      
      
      #urlTemplate for picture depends on the current selected data set
      #folders with pictures needs to be named in the correct manner
      if (input$SpatialEyeDataSetChooserInput > length(sceObjs)){
        selectedData = length(sceObjs) #avoid user to scroll to Data set that does not exist
        print(length(sceObjs))
      }else{
      selectedData = input$SpatialEyeDataSetChooserInput
      }
      sce = sceObjs[[as.numeric(selectedData)]]
      urlImg = paste("tiles", selectedData, "{z}/{x}/{y}.png", sep="/")
      print(urlImg)
      
      
      m <<- leaflet(options= leafletOptions(
        crs=leafletCRS(crsClass='L.CRS.Simple'))) %>% addTiles(urlTemplate = urlImg,
                                                               option=list(continuousWorld=TRUE, tms=TRUE,
                                                                           tileSize="256",minZoom="1",maxZoom="6")) %>%
        setView(lat=128, lng=128, zoom=1) %>%
        addCircleMarkers(data = pData(sce), lat = sce$spot_coord_y, lng = sce$spot_coord_x, color="black", stroke=FALSE, fillOpacity=0, group="ClickMarkers", layerId = colnames(sce) ) %>%
        addLayersControl(baseGroups = c("Tissue [add experiment name]"), 
                         overlayGroups = c("Spots", "Selected", "draw")) %>%
        addStyleEditor() %>%
        addDrawToolbar(targetGroup="draw",
                       polylineOptions = FALSE,
                       circleOptions = FALSE,
                       markerOptions = FALSE,
                       editOptions=editToolbarOptions(selectedPathOptions=selectedPathOptions()))
      
    })
       }
      }
    })
    
    
    observeEvent(input$Map_zoom, {
      if (input$SpatialEyeDataSetChooserInput > length(sceObjs)){
        selectedData = length(sceObjs) #avoid user to scroll to Data set that does not exist
        print(length(sceObjs))
      }else{
        selectedData = input$SpatialEyeDataSetChooserInput
      }
      sce = sceObjs[[as.numeric(selectedData)]]
      #OK.... this solution is maybe ugly? The radius needs to be adjusted to the zoomlevel, i.e. I will here update everything each time the zoom-button is pressed... 
      proxy <- leafletProxy("Map")
      print(input$Map_zoom)
      radi = 3*2^(input$Map_zoom-1)
      proxy %>%
        addCircleMarkers(data=pData(sce), radius= radi, lat = sce$spot_coord_y, lng = sce$spot_coord_x, layerId =colnames(sce), color="black", group="Spots", fillOpacity=0.4, stroke=FALSE) %>%
        hideGroup("Selected") %>% showGroup("Selected")
      
    })
    
    # ------------------------------------------------- Multiple selection marker - draw tool --------------------------------------------------------------------------------------------
    # ------------------- Draw selection marker
    observeEvent(input$Map_draw_new_feature, {
      if (input$SpatialEyeDataSetChooserInput > length(sceObjs)){
        selectedData = length(sceObjs) #avoid user to scroll to Data set that does not exist
        print(length(sceObjs))
      }else{
        selectedData = input$SpatialEyeDataSetChooserInput
      }
      sce = sceObjs[[as.numeric(selectedData)]]
      d <<- input$Map_draw_new_feature
      
      proxy <- leafletProxy("Map")
      if (d$properties$feature_type=="rectangle"){
        cord_lower_left = unlist(d$geometry$coordinates[[1]][[1]])
        cord_upper_left = unlist(d$geometry$coordinates[[1]][[2]])
        cord_lower_right = unlist(d$geometry$coordinates[[1]][[4]])
        cord_upper_right = unlist(d$geometry$coordinates[[1]][[3]])
        #proxy %>% addCircleMarkers(lng=cord_lower_left[1], lat=cord_lower_left[2], fillColor="red")
        #proxy %>% addCircleMarkers(lng=cord_upper_left[1], lat=cord_upper_left[2], fillColor="green")
        #proxy %>% addCircleMarkers(lng=cord_lower_right[1], lat=cord_lower_right[2], fillColor="yellow")
        #proxy %>% addCircleMarkers(lng=cord_upper_right[1], lat=cord_upper_right[2], fillColor="blue")
        
        
        #Select ST-spots that are within the drawn area
        #Approach idea: Subset the SCE-set for all coord_x, coord_y that are within area (given spot radius?)
        
        # x -bound
        sub_x = sce$spot_coord_x > cord_lower_left[1] & sce$spot_coord_x < cord_lower_right[1]
        
        # y - bound
        sub_y = sce$spot_coord_y > cord_lower_left[2] & sce$spot_coord_y < cord_upper_left[2]
        
        # subset sce object
        sce = sce[, sub_x & sub_y]
        colors <- brewer.pal(4,"Set1")
        color <<- colors[as.numeric(input$SelectionInput)]
        
        selection_group[sub_x & sub_y] <<- input$SelectionInput

  
      }else if(d$properties$feature_type=="polygon"){
        print("test polygon")
        d <- d
        geom <- d$geometry
        
        #Use the "point in polygon" function from the sp-package to determine which spots are inside of poly
        cords = unlist(geom$coordinates)
        pol.x = cords[seq(1, length(cords), 2)]
        pol.y = cords[seq(2, length(cords), 2)]
        spots_inside = point.in.polygon(point.x=pData(sce)$spot_coord_x, point.y=pData(sce)$spot_coord_y, pol.x=pol.x, pol.y=pol.y)
        
        #subset sce object
        sce = sce[, spots_inside==1]
        
        colors <- brewer.pal(4,"Set1")
        color <<- colors[as.numeric(input$SelectionInput)]
        
        selection_group[spots_inside==1] <<- input$SelectionInput
        
      }
     
      
      print(selection_group) #testing
      radi = 3*2^(input$Map_zoom-1)
      proxy %>% addCircleMarkers(radius=radi, data = pData(sce), lat = sce$spot_coord_y, lng = sce$spot_coord_x, color="black", fillColor=color, fillOpacity=1, opacity=1, weight=2, stroke=TRUE, 
                                 group="Selected", layerId = colnames(selected) )
      
      
      
    })
    
    # --------- Single click selectionn marker
    observeEvent(input$Map_marker_click, { 
      if (input$SpatialEyeDataSetChooserInput > length(sceObjs)){
        selectedData = length(sceObjs) #avoid user to scroll to Data set that does not exist
        print(length(sceObjs))
      }else{
        selectedData = input$SpatialEyeDataSetChooserInput
      }
      sce = sceObjs[[as.numeric(selectedData)]]
      
      #Two different modes:
      #(1) "Selection" - When user click on ST-spot, that spots is marked to the current selectionID
      #(2) "Spot statistics" - When user click on ST-spot, RNA-seq data is shown in popup window for that specific spot
      
      # --------- Mode (1)
      if (input$singleClickMode == "Selection"){
        p <- input$Map_marker_click
        print(p)
        proxy <- leafletProxy("Map")
        print(paste("You selected: ", p$id))
        
        colors <- brewer.pal(4,"Set1")
        color = colors[as.numeric(input$SelectionInput)]
        
        # subset sce object
        sce_tmp <- sce[, which(colnames(sce)==p$id)]
        #print(sce)
        selection_group[which(colnames(sce)==p$id)] <<- input$SelectionInput
        print(paste("Selection Group:", selection_group))
        
        radi = 3*2^(input$Map_zoom-1)
        proxy %>% addCircleMarkers(radius=radi, data = pData(sce_tmp), lat = sce_tmp$spot_coord_y, lng = sce_tmp$spot_coord_x, color="black", fillColor=color, fillOpacity=1, opacity=1, weight=2, stroke=TRUE, 
                                   group="Selected", layerId = colnames(selected) )
        
      }
      
      # ----------- Mode (2)
      if (input$singleClickMode == "Spot statistics"){
        print("Statistics mode ON")
        p <- input$Map_marker_click
        print(p$id)
        toggleModal(session, "Plot_and_table", toggle = "open")
        sce_tmp <<- sce[,colnames(sce)==p$id]
        table_row2=counts(sce_tmp[order(-counts(sce_tmp)),])
        table_row1=rownames(sce_tmp[order(-counts(sce_tmp)),])
        table_row3= exprs(sce_tmp[order(-counts(sce_tmp)),])
        table_out = cbind(table_row1, table_row2)
        table_out = cbind(table_out, table_row3)
        colnames(table_out) = c("GeneID", "Raw-count", "Norm. expression")
        
        table2_1 = dim(sce_tmp)[1]
        table2_2 = pData(sce_tmp)$total_counts
        table2_3 = pData(sce_tmp)$total_counts/dim(sce_tmp)[1]
        table2_4 = pData(sce_tmp)$pct_counts_top_50_features
        table2_5 = pData(sce_tmp)$pct_counts_top_100_features
        table2_6 = pData(sce_tmp)$pct_counts_top_200_features
        
        
        h1 = hist(sce$total_counts/1e3, breaks=40)
        cuts1 <- cut(h1$breaks, c(sce_tmp$total_counts/1e3, sce_tmp$total_counts/1e3+(h1$breaks[2]-h1$breaks[1])))
        h2 = hist(sce$total_features, breaks=40)
        cuts2 <- cut(h2$breaks, c(sce_tmp$total_features, sce_tmp$total_features+(h2$breaks[2]-h2$breaks[1])))
        colors <- brewer.pal(4,"Set1")
        
        output$TestPlot <- renderPlot({
          
          par(mfrow=c(1,2))
          plot(h1, col=colors[cuts1], main="Total counts across ST-spots", xlab="Library sizes (thousands)",  ylab="Number of ST-features")
          plot(h2, col=colors[cuts2], xlab="Number of expressed genes", main="Total genes identified across ST-spots",  ylab="Number of ST-features")
          
        })
        
         output$spotStats <- renderPrint({
           
           cat(paste("# Genes\n", table2_1, "\nTotal counts\n", table2_2, "\nMean count per Gene\n", round(table2_3,0), "\nTop 50 genes accont for\n", round(table2_4,0), "%", "\nTop 100 genes accont for\n", round(table2_5,0),"%", "\nTop 200 genes accont for\n", round(table2_6, 0), "%"), sep="")
         })
         
         output$geneStats <- renderDataTable({
           
           table_out
           
           
         }, options = list(pageLength=10))
      }
    })
    
    # observeEvent(input$Map_shape_click, {
    #   v$msg <- paste("Clicked shape", input$Map_shape_click$id)
    #   print(v$msg)
    # })
    
    
    # -------------------------------------------------------------------- Assign selection ----------------------------------------------------------------------------------
    #After all selections are performed, the "Assign Selection" button will finaly put the selectionIDs in the pData of the sce object
    observeEvent(input$assignSelect, {
      if (input$SpatialEyeDataSetChooserInput > length(sceObjs)){
        selectedData = length(sceObjs) #avoid user to scroll to Data set that does not exist
        print(length(sceObjs))
      }else{
        selectedData = input$SpatialEyeDataSetChooserInput
      }
      sce = sceObjs[[as.numeric(selectedData)]]
      pData(sce)$selectionID = selection_group
      sceObjs[[as.numeric(selectedData)]] <<- sce #save to global
      print("Selection assigned")
    })
    
    # -------------------------------------------------------------------- Reset selection ----------------------------------------------------------------------------------
    #The reset selection button will remove all selection (if any were assigned) from the pData of the sce object
    
    observeEvent(input$resetSelect, {
      if (input$SpatialEyeDataSetChooserInput > length(sceObjs)){
        selectedData = length(sceObjs) #avoid user to scroll to Data set that does not exist
        print(length(sceObjs))
      }else{
        selectedData = input$SpatialEyeDataSetChooserInput
      }
      sce = sceObjs[[as.numeric(selectedData)]]
      pData(sce)$selectionID = 0
      proxy <- leafletProxy("Map")
      proxy %>% clearShapes() %>%
                clearMarkers() %>%
                addCircles(radius=5, lat=0, lng=0, stroke=TRUE) %>%
        #Total hidden markers for clicking events?
        addCircleMarkers(data = pData(sce), lat = sce$spot_coord_y, lng = sce$spot_coord_x, color="black", stroke=FALSE, fillOpacity=0, group="ClickMarkers", layerId = colnames(sce) ) 
      
      print(pData(sce)$selectionID)
      print("Selection Reset")
      sceObjs[[as.numeric(selectedData)]] <<- sce
    })
    
    
    # --------------------------------------------------------------------- Selection QC look-up ----------------------------------------------------------------------------
    
    observeEvent(input$checkSelection, {
      observeEvent(input$resetSelect, {
        if (input$SpatialEyeDataSetChooserInput > length(sceObjs)){
          selectedData = length(sceObjs) #avoid user to scroll to Data set that does not exist
          print(length(sceObjs))
        }else{
          selectedData = input$SpatialEyeDataSetChooserInput
        }
        sce = sceObjs[[as.numeric(selectedData)]]
      #Make QC on alla spots that are included in the current selection
      currentSelect = input$SelectionInput
      print(paste("currentSelect: ", currentSelect))
      print(selection_group)
      #subset the SCE object for these ST-spots
      multiQC_idx = which(selection_group==currentSelect)
      print(multiQC_idx)
      sce_multiQC <<- sce[, multiQC_idx]
      
      
      #Toggle on the modal popup-window
      toggleModal(session, "multiQC", toggle = "open")
      
      #Table stats
      # OBS OBS NEED TO CHECK IF THIS IS CORRECT !!!!! ### =¤ =¤?=¤?=¤ ?=?¤=?¤=?¤=?¤=¤?=¤?=
      table_row2=counts(sce_multiQC[order(-rowSums(counts(sce_multiQC))),])
      table_row1=rownames(sce_multiQC[order(-rowSums(counts(sce_multiQC))),])
      table_row3= exprs(sce_multiQC[order(-rowSums(counts(sce_multiQC))),])
      table_out = cbind(table_row1, rowSums(table_row2))
      table_out = cbind(table_out, rowSums(table_row3))
      colnames(table_out) = c("GeneID", "Raw-count", "Norm. expression")
      
      output$multiQCstandard <- renderPlot({
        
        par(mfrow=c(1,2))
        h1 = hist(sce_multiQC$total_counts/1e3, breaks=40)
        h2 = hist(sce_multiQC$total_features, breaks=40)
        plot(h1, main="Total counts across ST-spots", xlab="Library sizes (thousands)",  ylab="Number of ST-features")
        plot(h2, xlab="Number of expressed genes", main="Total genes identified across ST-spots",  ylab="Number of ST-features")
        
      })
      
      output$multiQCPCA <- renderPlot({
        par(mfrow=c(1,2))
        plotPCAmultiQC <<- plotPCA(sce_multiQC, exprs_values="exprs") + coord_fixed()
        plotPCAmultiQC
      })
      
      output$multiQCstats <- renderPrint({
        
        avg_lib_size = sum(colSums(counts(sce_multiQC)))/as.numeric(dim(sce_multiQC)[2])
        avg_geneNr_per_spot = NULL #fix later
        cat(paste("Avg. library size: ", avg_lib_size, "\n", "Avg. nr of genes / ST-Spot: ", sep=""))
      })
      
      output$multiQCgeneStats <- renderDataTable({
        
        table_out
        
        
      }, options = list(pageLength=10))
    })
    })
    
    # ------------------------------------------------------------------------------------ Multi QC PCA - Brush effect -----------------------------------------------------------
    output$multiQCPCA_brushinfo <- renderPrint({
      cat("input$plot_brush:\n")
      str(input$plot_brush)
    })
                                   
    output$multiQCPCA_brush_points = renderPrint({
      brush_PCA <<- brushedPoints(plotPCAmultiQC$data,input$multiQCPCA_brush)
      #print(input$multiQCPCA_brush)
      #print(rownames(brush_PCA))
      
    })
    
    # ++++++++++++++ Remove brushed points from current selection group
    observeEvent(input$removeFromMultiQCPCA, {
      remove_idx = which(colnames(sce) %in% rownames(brush_PCA))
      selection_group[remove_idx] <<- "0"
      print(selection_group)
      output$removedOutput = renderPrint({cat("You just removed:", length(remove_idx), " points from selection")})
    })

      
    
    
    
    # --------------------------------------------------- Img Save -----------------------------------------------------------------------------------------------------------------
    # #http://stackoverflow.com/questions/31336898/how-to-save-leaflet-in-rstudio-map-as-png-or-jpg-file
    #   
    #   output$html_link <- downloadHandler(
    #     filename = 'temp.html',
    #     
    #     content = function(file) {
    #       src <- normalizePath('mymap.Rmd')
    #       
    #       # temporarily switch to the temp dir, in case you do not have write
    #       # permission to the current working directory
    #       owd <- setwd(tempdir())
    #       on.exit(setwd(owd))
    #       file.copy(src, 'mymap.Rmd')
    #       
    #       out <- render('mymap.Rmd',
    #                     html_document()
    #       )
    #       file.rename(out, file)
    #     }
    #   )
    
    
    
  
  
  
  
  
  # -----------------------------------------------------------------------------------------------------------------------------------------------------------
  #END OF SERVER FUNCTION CALL:
  
  })




# ----------------------- EXTRA FUNCTIONS OUTSIDE OF SERVER FUNCTION - FOR SELECTOR PLOT -------------------------------------------------------------

# A battery of internal functions, taken out to reduce the length of the main function.

resetValues <- function(collected, coords=NULL) {
  if (is.null(coords)) { coords <- collected$coords }
  else { collected$coords <- coords }
  N <- nrow(coords)
  collected$old.selected <- logical(N)
  collected$current.selected <- logical(N)
  collected$all.selected <- list()
}

generatePlot1 <- function(output, collected, ...) {
  output$plot1 <- renderPlot({
    cols <- rep("grey", length(collected$current.selected))
    cols[collected$old.selected] <- "orange"
    cols[collected$current.selected] <- "red"
    x <- collected$coords$x
    y <- collected$coords$y
    plot(x, y, col=cols, ...)
  })
}

generatePlot2 <- function(output, collected, ...) {
  output$plot2 <- renderPlot({
    n <- length(collected$all.selected)
    shading <- rev(grey.colors(n))
    x <- collected$coords$x
    y <- collected$coords$y
    plot(x, y, type="n", ...)
    for (i in seq_len(n)) {
      current <- collected$all.selected[[i]]
      text(x[current], y[current], labels=i, col=shading[i])
    }
  })
}

# ------------- SCE merge function -------------

mergeSCE = function(sets){
  withProgress(message =paste("Combining data sets..."), {
  print(paste("SceMerge function sets: ", sets, sep=""))
  list_of_data = sceObjs[sets]
  print(paste("list_of_data: ", list_of_data))
  common_names = Reduce(intersect, lapply(list_of_data, row.names))
  list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })
  #Now we only have the genes that all 3 sets have in common, we can now merge them
  for (i in 1:length(list_of_data)){
    colnames(list_of_data[[i]]) = paste(colnames(list_of_data[[i]]), i, sep="-")
    list_of_data[[i]] = list_of_data[[i]][order(rownames(list_of_data[[i]]))]
  }
  
  sceMerged = mergeSCESet(list_of_data[[1]], list_of_data[[2]], fdata_cols_x = FALSE)
  if (length(sets)>2){
    for (i in 2:(length(list_of_data)-1)){
      fData(sceMerged)$dummy = 1
      sceMerged = mergeSCESet(sceMerged, list_of_data[[i+1]], fdata_cols_x = FALSE, fdata_cols_y = FALSE)
    }
  }
  
  incProgress(detail = paste("Done!")) #Progress bar END
  Sys.sleep(1)})
  
  return(sceMerged)
}


# ----------- Gradient tool from ref points -----------
gradient_matrix = function(sce, refPoint){
  
  #refPoint is given as a selectionID from Spatial EYE
  #create tmp sce with only ref-Spots
  sce_ref = sce[, which(pData(sce)$selectionID == refPoint)]
  
  min_distance_vector = NULL
  x1_vector = pData(sce)$spot_coord_x
  y1_vector = pData(sce)$spot_coord_y
  x2_vector = pData(sce_ref)$spot_coord_x
  y2_vector = pData(sce_ref)$spot_coord_y
  
  #Loop through the ST-Spots
  #pixel coordinates saved in pData(sce)$spot_coord_y and pData(sce)$spot_coord_x
  for (i in 1:ncol(sce)){
    print(i)
    distance=NULL
    
    #loop through each of the ref_points to compare against
    for (j in 1:ncol(sce_ref)){
      
      x_dist = (x2_vector[j] - x1_vector[i])
      y_dist = (y2_vector[j] - y1_vector[i])
      distance[j] = sqrt(x_dist*x_dist + y_dist*y_dist) #fill vector with distances
    }
    #Extract the shortest distance
    min_distance_vector[i] = min(distance)
  }
  #read into pData
  pData(sce)$dist_from_ref = min_distance_vector
  #return the sce object (Obs not saved to global yet)
  return(sce)
  
}
