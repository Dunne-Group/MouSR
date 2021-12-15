
vals <- reactiveValues(count=0) 

## Uploading file's size requirement (<30MB)
options(shiny.maxRequestSize = 30*1024^2)
## Start shiny server  

shinyServer(function(input, output, session) {
  pdf(NULL)
  isolate(vals$count <- vals$count + 1)
######################################################################################################### 
  #undefined columns selected
#######################################################################################################################
#########################################################################################################
  observe({

    chFileFormat <- input$Fileformat
    if (chFileFormat == ""){ 
      output$a1<- renderUI({
          isolate(
            actionButton("stop", label = "Check Input Files"))
      }) 
      
      observeEvent(input$stop,{
          shinyalert("Warning!", "Please Select one Option from above", type = "error")
        })
        }
    
    
    else if (chFileFormat == "1"){
      output$a1<- renderUI({
        isolate(
          actionButton("dataSubmit", label = "Check Input Files"))
      })

      observeEvent(input$dataSubmit,{
        withProgress(message = "Checking your data, please wait",
                     {
                       
                       
##### UpLoding ##################################################################################################################
    
      values <- reactiveValues()
      
          file_to_read <- input$file1
          file_to_read2 <- input$file2
          
          if (is.null(file_to_read ) & is.null(file_to_read2)) { 
            
            queryMagic <- function() {
              print(paste(" Please upload your files first"))
              
              return("Data")
            }
            output$console2 <- renderPrint({
              logText()
              return(print(values[["log"]]))
            })
            
            logText <- reactive({
              values[["log"]] <- capture.output(data <- queryMagic())
            })
            
            return()
          }
      else if (!is.null(file_to_read ) & !is.null(file_to_read2)) { 
      
      UploadData <- reactiveValues(

       mainfile = read.delim(file=file_to_read$datapath,sep=input$Sep1, header=TRUE)) 
  
      Uploadtable <- reactiveValues(
       
        labelfile = read.delim(file=file_to_read2$datapath,sep=input$Sep1, header=TRUE) 
      )
     

     if (length(colnames( UploadData$mainfile)) < 4 || length(colnames(Uploadtable$labelfile))> 3) {
       
         queryMagic <- function() {
           if(input$dataSubmit)
           print(paste("Wrong File Uploaded: Input 1 should be your Main file and Input 2 should be your Labels table"))
           
           return("Data")
         }
         output$console2 <- renderPrint({
           logText()
           return(print(values[["log"]]))
         })
         
         logText <- reactive({
           values[["log"]] <- capture.output(data <- queryMagic())
         })
         
         output$b1<- renderUI({
         if(input$dataSubmit)
           isolate(
           actionButton("wait", label = "Continue to Submit"))
         }) 
      
         output$error <- renderText({
           if(input$dataSubmit)
             isolate({
           paste("If there is any error or Input summary table seems wrong, please make sure: 1)You have uploaded correct files with correct format (CSV (Comma delimited)/ txt (Tab delimited)).
                 2) You have picked correct file Separators.")
                    })
         })
          
     } 
    
      
     else if (length(colnames( UploadData$mainfile)) > 4 && length(colnames(Uploadtable$labelfile))< 4) {
        
        queryMagic <- function() {
         print(paste(" Correct input Uploaded"))
          return("Data")
        }
        output$console2 <- renderPrint({
          logText()
          return(print(values[["log"]]))
        })
        
        logText <- reactive({
          values[["log"]] <- capture.output(data <- queryMagic())
        })
        
        if(input$dataSubmit)
          isolate({
        output$b<- renderUI({
          
         actionButton("dataPCA", label = " Continue to Submit")
        }) 
        output$b1<- renderUI({
   
        }) 
        
        output$error <- renderText({
          if(input$dataSubmit)
            isolate({
              })
        })
        
          }) 
     } 
      
        }
      
        
         
      ######## Summary ##############################################################################################################
      output$Title <- renderText({
        if(input$dataSubmit)
          isolate({
            no.samples <- length(colnames(UploadData$mainfile))-2
            paste("This cohort contains ", as.character(no.samples)," Samples:", sep="")
          })
      })
      ######################################################################################################################
      output$Info <- renderText({
        if(input$dataSubmit)
          isolate({
            paste(as.character( (Uploadtable$labelfile)[,1]), collapse=", " )
          })
      })
      #####################################################################################################################
      output$t1 <- renderText({
        if(input$dataSubmit)
          isolate({
            paste(as.character(colnames(Uploadtable$labelfile)[2]),sep=""," Label:")
          })
      })
      #####################################################################################################################
      output$t2 <- renderText({
        if(input$dataSubmit)
          isolate({
            paste(as.character(levels(as.factor((Uploadtable$labelfile)[,2]))), collapse=", ")
          })
      })
      #######################################################################################################################
      output$Summary <- renderTable({
        if(input$dataSubmit)
          isolate({
            no.samples <- length(colnames(UploadData$mainfile))-2
            no.gene <- dim(UploadData$mainfile)[1]
            res.summary <- rbind(no.samples, no.gene)
            rownames(res.summary) <- c("Samples", "Gene_ID")
            colnames(res.summary) <- "Number"
            res.summary
          })
      },digits=0,rownames = TRUE, align="lc")
      
      #})
      #################################################################################################################  
      
   # }  
     })  
  })    
    
    ##############PCA 2D PLOT  ####################################################
    ##PCA###
    data <- reactiveValues()
    observeEvent(input$dataPCA,{
      withProgress(message = "Preparing PCA & MDS Plots, please wait",
                   {
                     file_to_read=input$file1
                     f6<- read.delim(file_to_read$datapath,sep=input$Sep1, header=T)
                     names(f6)[1] <- "ID"
                     names(f6)[2] <- "symbol"
                   
                     head(f6)
                     dim(f6)
                     # remove rows with zero value
                     f7 <- column_to_rownames(f6, var = "ID")[-1]
                     f7 <- f7[rowSums(f7 != 0)>0,]
                     #### attaching gene symbol to the file
                     f8 <- rownames_to_column(f7, var="ID")
                     exp <- join(f8[,1,drop=F], f6, by="ID", type="inner")
                     data$exp<-exp
                     file_to_read2=input$file2
                     table <- read.delim(file_to_read2$datapath,sep=input$Sep1, header=T)
                     names(table)[1] <- " Types"
                     names(table)[2] <- "Samples"
                     names(table)[3] <- "Labels"
                     rownames(exp) <- exp[,1]
                     exp <- exp[,-c(1,2)]
                     exp.scale <- t(scale(t(exp), scale = F))
                     pc <- prcomp(exp.scale, scale. = F, center = F)
                     percentage <- round(pc$sdev / sum(pc$sdev) * 100, 2)
                     pc$rotation
                     pcr <- data.frame(pc$rotation[,1:3], Group = table[,1] )
                     data$pp <-data.frame(pc$rotation[,1:3], Group = table[,1] )
                     data$percentage <- paste( colnames(pcr), "(", paste( as.character(percentage), "%", ")", sep="") )
                    
                     pcr1 <- data.matrix(pc$rotation[,1:3])
                     pcr2 <- data.frame(pc$rotation[,1:3], Group = table[,1] )
                     pcr$label<- table[,3]
                     data$label<- table[,3]
                     
                     data$gr<-(table[,1])
                     
                     data$plot<-pcr1[,1:3]
                     data$plot3<-pcr2[,1:3]
                     
                     
                     data$table<-table
                     
                     A <- table[,1]
                     B = factor(A)
                     B7=(levels(B)) 
                     data$B7<- length(B7)
                     data$A <-A
                     ##############Group Selection################################                  
                     output$color1 <- renderUI({
                       selectInput('select','Please pick a color for all of your Groups:', as.list(B7), multiple = TRUE) 
                     })
                     
                     output$check1 <- renderUI({
                       selectInput ('Group1comp','Category A', as.list(B7),multiple = TRUE)
                     })
                     
                     output$check2 <- renderUI({
                       selectInput ('Group2comp','Category B', as.list(B7),multiple = TRUE)
                     })
                     
                     output$check3 <- renderUI({
                       selectInput ('Group3comp','Group A', as.list(B7),selected = as.list(B7[1]), multiple = FALSE)
                       #selectInput ('Group3comp','Group A', as.list(B7),selectize=TRUE)
                       
                     })
                     
                     output$check4 <- renderUI({
                       selectInput ('Group4comp','Group B', as.list(B7),selected = as.list(B7[2]), multiple = FALSE)
                     })
                     
                     my_list <- list("both","row","column")    
                     output$checkdan <- renderUI({
                       selectInput ('Group7comp','Choose clustering option', as.list(my_list),multiple = FALSE)
                     }) 
                     
                     output$checkdanselected <- renderUI({
                       selectInput ('Group8comp','Choose clustering option', as.list(my_list),multiple = FALSE)
                     }) 
                     
                     
                     ############## MDS Plot ##########################################
                     distance.matrix <- dist(scale(t(exp), scale = T,center = T ), method="euclidean")
                     
                     mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
                     
                     ## calculate the percentage of variation that each MDS axis accounts for...
                     data$mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
                     #print(data$mds.var.per)
                     mds.values <- mds.stuff$points
                     data$mds.data <- data.frame(Sample=rownames(mds.values),
                                                 X=mds.values[,1],
                                                 Y=mds.values[,2], Group = table[,1] )
                     
                     ############################################################################### 
                     shinyalert("Done!", "Please go to Data analysis section", type = "success")
                   }) #progress bar
    }) # End submit
    
    #############Colour ###########################################################  
    
    gg_fill_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    output$platter <- renderUI({ 
      lev <- sort(unique(input$select)) 
      cols <- gg_fill_hue(length(lev))
      
      lapply(seq_along(lev), function(i) {
        colourInput(inputId = paste0("col", lev[i]),
                    label = paste0("Choose colour for ", lev[i]), 
                    value = cols[i]
        )        
      })
    }) 
    ##############PCA 2D PLOT  ##################################################### 
    observe({
      chp <- input$optionp
      chlab <- input$optionlab
      textpca <- as.numeric(input$label_pca)
      pointpca <- as.numeric(input$point_pca)
      widthpca <- reactive ({ input$plot_widthpca })
      heightpca <- reactive ({ input$plot_heightpca}) 
      if (chp == "1"){
        output$PCAPlot2D <- renderPlot(width = widthpca, height = heightpca,{
          cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
         
          cols <- eval(parse(text = cols))
        
          req(length(cols) == length(input$select))
          if (length(cols) == data$B7){
            if (chlab == "1"){
              sp <- ggplot(data$pp, aes(PC1, PC2,label = data$label,fill = Group))+ geom_point(aes(col = Group),size=pointpca)
              p <- sp+ scale_color_manual(values = cols)+
                geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[2])}
           
            else if (chlab == "2"){
              sp <- ggplot(data$pp, aes(PC1, PC2, fill = Group)) +geom_point(aes(col = Group),size=pointpca) 
              p <-sp + scale_color_manual(values = cols)+
                theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[2])}
          
            vals$p<-p
            print(p)
            vals$cols<-cols
            
          }
          
          else{
            if (chlab == "1"){
              p <-ggplot(data$pp, aes(PC1, PC2, label = data$label, fill = Group))+geom_point(aes(col = Group),size=pointpca)+geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[2])}
            else if (chlab == "2"){
              p <-ggplot(data$pp, aes(PC1, PC2, color = Group))+geom_point(aes(col = Group),size=pointpca)+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[2])}
            vals$p<-p
            print(p)
            
          }
        }) 
      }# Chp=1
      
      else if (chp == "2"){
        output$PCAPlot2D <- renderPlot(width = widthpca, height = heightpca,{
          cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
          cols <- eval(parse(text = cols))
          
          req(length(cols) == length(input$select))
          if (length(cols) == data$B7){ 
            if (chlab == "1"){
              sp<- ggplot(data$pp, aes(PC1, PC3,label = data$label, fill = Group)) + 
               
                geom_point(aes(col = Group),size=pointpca) 
            
              p <- sp+ scale_color_manual(values = vals$cols)+
                geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[3])}
            else if  (chlab == "2"){
              sp<- ggplot(data$pp, aes(PC1, PC3, fill = Group)) + 
              
                geom_point(aes(col = Group),size=pointpca) 
             
              p <- sp+ scale_color_manual(values = vals$cols)+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[3])}
            
            vals$p<-p
            
            print(p)
          }
          
          else{
            if (chlab == "1"){
              p <-ggplot(data$pp, aes(PC1, PC3, label = data$label, fill = Group))+geom_point(aes(col = Group),size=pointpca)+geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[3])}
            else if (chlab == "2"){
              p <-ggplot(data$pp, aes(PC1, PC3, color = Group))+geom_point(aes(col = Group),size=pointpca)+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[3])}
            vals$p<-p
            print(p)
          } 
          
          
        })
      }#Chp=2
      
      else if (chp == "3"){
       
        output$PCAPlot2D <- renderPlot(width = widthpca, height = heightpca,{
          cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
          
          cols <- eval(parse(text = cols))
       
          req(length(cols) == length(input$select))
          if (length(cols) == data$B7){
            if (chlab == "1"){
              sp<- ggplot(data$pp, aes(PC2, PC3,label = data$label, fill = Group)) + 
                geom_point(aes(col = Group),size=pointpca) 
              
              p <- sp +  scale_color_manual(values = vals$cols) +
                geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[2]) + ylab(data$percentage[3])}
            
            else if (chlab == "2"){
              sp<- ggplot(data$pp, aes(PC2, PC3, fill = Group)) + 
                geom_point(aes(col = Group),size=pointpca) 
             
              p <- sp +  scale_color_manual(values = vals$cols) +theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[2]) + ylab(data$percentage[3])
            }
            
            
            vals$p<-p
            print(p)
          }
          else{
            if (chlab == "1"){
              p <-ggplot(data$pp, aes(PC2, PC3, label = data$label, fill = Group))+geom_point(aes(col = Group),size=pointpca)+geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[2]) + ylab(data$percentage[3])}
            
            else if (chlab == "2"){
              p <-ggplot(data$pp, aes(PC2, PC3, color = Group))+geom_point(aes(col = Group),size=pointpca)+theme_bw(base_size = 14)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.background = element_blank())+ xlab(data$percentage[2]) + ylab(data$percentage[3]) }
            vals$p<-p
            print(p)
            
            
          } 
          
        })
      }#Chp=3
    })  # end observe  
    ###########################################################################################################      
    ####### MDS PLOT
    observe({
      chmds <- input$optlab
      textmds<- as.numeric(input$label_mds)
      pointmds<- as.numeric(input$point_mds)
      widthmds <- reactive ({ input$plot_widthmds })
      heightmds <- reactive ({ input$plot_heightmds}) 
      output$MDSPlot <- renderPlot(width = widthmds, height = heightmds,{  
        cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
       
        cols <- eval(parse(text = cols))

        req(length(cols) == length(input$select))
        if (length(cols) == data$B7){
          if (chmds == "1"){
            g<- ggplot(data$mds.data, aes(x=X, y=Y, label = data$label, fill = Group)) + 
              geom_point(aes(col = Group),size=pointmds) 

            s <- g +  scale_color_manual(values = vals$cols) +
              geom_text(size= textmds,position = position_nudge(y = +5))+theme_bw(base_size = 14)+
              theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank())+ xlab(paste("MDS1 - ", data$mds.var.per[1], "%", sep="")) +ylab(paste("MDS2 - ", data$mds.var.per[2], "%", sep="")) }
          
          else if (chmds == "2"){
            g<- ggplot(data$mds.data, aes(x=X, y=Y, fill = Group)) + 
              geom_point(aes(col = Group),size=pointmds) 
        
            s<- g +  scale_color_manual(values = vals$cols) +theme_bw(base_size = 14)+
              theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank())+ xlab(paste("MDS1 - ", data$mds.var.per[1], "%", sep="")) +ylab(paste("MDS2 - ", data$mds.var.per[2], "%", sep="")) }
          
          vals$s<-s
          print(s)
        }
        else{
          if (chmds == "1"){
            s <-ggplot( data$mds.data, aes(x=X, y=Y, label = data$label, fill = Group))+geom_point(aes(col = Group),size=pointmds)+geom_text(size= textmds,nudge_y = +5)+theme_bw(base_size = 14)+
              theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank())+ xlab(paste("MDS1 - ", data$mds.var.per[1], "%", sep="")) +ylab(paste("MDS2 - ", data$mds.var.per[2], "%", sep="")) }
          else if (chmds == "2"){
            s <-ggplot( data$mds.data, aes(x=X, y=Y, fill= Group))+geom_point(aes(col = Group),size=pointmds)+theme_bw(base_size = 14)+
              theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank())+ xlab(paste("MDS1 - ", data$mds.var.per[1], "%", sep="")) +ylab(paste("MDS2 - ", data$mds.var.per[2], "%", sep="")) }
          vals$s<-s
          print(s)
          
        } 
      }) 
    })     
    ################ 3D PLOT #####################################################################  
    output$PCAPlot <- renderPlotly({
      withProgress(message = "Preparing PCA3D, please wait",
                   {
                     get_colors <- function(groups, group.col = palette()){
                       groups <- data$gr
                       ngrps <- length(levels(groups))
                       if(ngrps > length(group.col))
                         group.col <- rep(group.col, ngrps)
                       color <- group.col[as.numeric(groups)]
                       names(color) <- as.vector(groups)
                       return(color)
                     }
                     gr<-data$gr
                     n <- data$plot3
                     s=get_colors(gr)
                     m<-as.factor(s)
                     cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
                     print(cols)
                     cols <- eval(parse(text = cols))
                     req(length(cols) == length(input$select))
                     if (length(cols) == data$B7){
                       fig <- plot_ly(n, x = ~PC1, y = ~PC2, z = ~PC3, color = ~gr, colors = vals$cols,text = rownames(n))
                       fig <- fig %>% add_markers()
                       chs3d <- input$option3d
                       if (chs3d == "1"){
                         
                         fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                            yaxis = list(title = 'PC2'),
                                                            zaxis = list(title = 'PC3'),
                                                            camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
                        
                         fig <- fig %>% layout(autosize = T)
                         fig <- fig %>% config(
                           toImageButtonOptions= list(filename = 'PCA3DPlot'))
                         fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                            yaxis = list(title = 'PC2'),
                                                            zaxis = list(title = 'PC3')
                         ))
                       }
                       else if (chs3d == "2"){
                         fig <- plot_ly(n, x = ~PC1, y = ~PC2, z = ~PC3, color = ~gr, colors = vals$cols,text = rownames(n))
                         fig <- fig %>% add_markers()
                         
                         fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                            yaxis = list(title = 'PC2'),
                                                            zaxis = list(title = 'PC3'),
                                                            camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
                         
                         fig <- fig %>% layout(autosize = T)
                         fig <- fig %>% config(
                           toImageButtonOptions= list(format = "svg",filename = 'PCA3DPlot'))
                         fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                            yaxis = list(title = 'PC2'),
                                                            zaxis = list(title = 'PC3')
                         ))
                       } 
                       
                     }
                     
                     else
                     {
                       fig <- plot_ly(n, x = ~PC1, y = ~PC2, z = ~PC3, color = ~gr, colors= c("#F8766D","#00BFC4","#7CAE00", "#C77CFF", "#009E73",
                                                                                              "#7CAE00", "#C77CFF", "#D55E00", "#00BFC4"),text = rownames(n))
                       fig <- fig %>% add_markers()
                       chs3d <- input$option3d
                       if (chs3d == "1"){
                         fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                            yaxis = list(title = 'PC2'),
                                                            zaxis = list(title = 'PC3'),
                                                            camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
                         
                         fig <- fig %>% layout(autosize = T)
                         fig <- fig %>% config(
                           toImageButtonOptions= list(filename = 'PCA3DPlot'))
                         fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                            yaxis = list(title = 'PC2'),
                                                            zaxis = list(title = 'PC3')
                         ))
                       }
                       else  if (chs3d == "2"){
                         fig <- plot_ly(n, x = ~PC1, y = ~PC2, z = ~PC3, color = ~gr, colors= c("#F8766D","#00BFC4","#7CAE00", "#C77CFF", "#009E73",
                                                                                                "#7CAE00", "#C77CFF", "#D55E00", "#00BFC4"),text = rownames(n))
                         fig <- fig %>% add_markers()
                         
                         fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                            yaxis = list(title = 'PC2'),
                                                            zaxis = list(title = 'PC3'),
                                                            camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
                         
                         fig <- fig %>% layout(autosize = T)
                         fig <- fig %>% config(
                           toImageButtonOptions= list(format = "svg",filename = 'PCA3DPlot'))
                         fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                            yaxis = list(title = 'PC2'),
                                                            zaxis = list(title = 'PC3')
                         ))
                       }
                     }
                     
                   })  #progress bar
    })# ENd Plotly
    
    ###############Downloads PCA2D ###########################################################  
   # print(input$plot_widthpca*4)
    observe({
      chs <- input$optionchs
      if (chs == "1"){
        output$PlotDownloadPCA2D <- downloadHandler(
          filename = function(){paste("PCAPlot",'.png',sep='')},
          content = function(file){
            device <- function(..., width, height) grDevices::png(..., width =input$plot_widthpca*4, height =input$plot_heightpca*4, res = 350, units = "px")
            ggsave(file, plot = (vals$p)  , device = device)
          })
      }
   
      else if (chs == "2"){
        output$PlotDownloadPCA2D <- downloadHandler(
          filename = function(){paste("PCAPlot",'.svg',sep='')},
          content = function(file) {
            device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthpca*4/300), height =round(input$plot_heightpca*4/300))
            ggsave(file, plot = (vals$p)  , device = device)
          
          })
      }
    })
    
    ########### Donwload MDS ###############################################################################   
    observe({
      chd <- input$optchs
      if (chd== "1"){
        output$PlotDownloadMDS2D <- downloadHandler(
          filename = function(){paste("MDSPlot",'.png',sep='')},
          content = function(file){
            device <- function(..., width, height) grDevices::png(..., width =input$plot_widthmds*4, height =input$plot_heightmds*4, res = 350, units = "px")
            ggsave(file, plot = (vals$s)  , device = device)
          })
      }
      else if (chd == "2"){
        output$PlotDownloadMDS2D <- downloadHandler(
          filename = function(){paste("MDSPlot",'.svg',sep='')},
          content = function(file) {
            device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthmds*4/300), height =round(input$plot_heightmds*4/300))
            ggsave(file, plot = (vals$p)  , device = device)
            
          
          })
      }
    })# End MDS
######################################################################################################
    ###############################################################################################################
    ####CRIS classifier
    datacris <- reactiveValues()
    observeEvent(input$CRISAnalysis,{
      withProgress(message = "Preparing CRIS Classifier Table and Plot, please wait",
                   {
                     Collapsed_data <- datar$Coll
                     head(Collapsed_data)
                    
                     
                     ########## Mouse data############################################   
                     chcris <- input$optionCRIS
                     if (chcris  == "Mouse"){   
                       ####### Get the normalized data from the first part
                       Collapsed_cris <- Collapsed_data[,-c(1,2)]
                       col_name8<- paste(datar$n1, sep="")
                       col_name9<-paste(datar$n2, sep="")
                       
                       names(Collapsed_cris)[1:length(datar$n1)]<-col_name8
                       names(Collapsed_cris)[length(datar$n1)+1:length(datar$n2)]<-col_name9
                       
                       
                       expcris <- data.matrix(Collapsed_cris)
                       
                       # adjust matrix before classification
                       mat_exp_adj <- ematAdjust(expcris)
                       datacris$mat_exp_adj<-mat_exp_adj
                       load("CRIS_template_mouse.RData")
                       
                       NTP_CRIS <- ntp(emat=mat_exp_adj,CIRS_template, seed=367707,doPlot =F)
                       NTP_CRIS <- tibble::rownames_to_column(NTP_CRIS,var = "ID")
                       
                       crisl <-(length(colnames(NTP_CRIS)))
                       
                       output$CRIS<- DT::renderDataTable({NTP_CRIS},
                                                         extensions = 'FixedColumns', 
                                                         options = list(pageLength = 20, autoWidth = FALSE,
                                                                        scrollX = TRUE,
                                                                        scrollCollapse = TRUE,
                                                                        lengthMenu = c(20, 50, 100),
                                                                        columnDefs = list(list(targets = c(3:crisl), searchable = FALSE))),
                                                         editable = FALSE)
                       
                       
                       
                       output$TableDownloadCRIS<- downloadHandler(
                         filename = function(){paste("MOUSE_DATA_NTP_CRIS_results",'.csv')},
                         content = function(file){
                           write.csv(NTP_CRIS, file, sep = ",", quote = F)
                         })
                       
                       
                       
                       observe({  
                         widthcris <- reactive ({ input$plot_widthcris})
                         heightcris<- reactive ({ input$plot_heightcris}) 
                         output$CRIS_Plot<- renderPlot(width = widthcris, height = heightcris,{
                           NTP_CRIS <- ntp(emat=mat_exp_adj,CIRS_template, seed=367707,doPlot = TRUE)
                         })  
                         
                       })
                       
                       observe({
                         chdcris <- input$optioncris
                         if ( chdcris  == "1"){ 
                           
                           output$PlotDownloadCRIS<- downloadHandler(
                             file = "CRIS-Plot.png" , # variable with filename
                             content = function(file) {
                               png(file = file)
                               load("CRIS_template_mouse.RData")
                               NTP_CRIS <- ntp(emat=datacris$mat_exp_adj,CIRS_template, seed=367707,doPlot = TRUE)
                               dev.off()
                             })
                         } else if  ( chdcris == "2"){ 
                           
                           output$PlotDownloadCRIS <- downloadHandler(
                             file = "CRIS-Plot.svg" , # variable with filename
                             content = function(file) {
                               svg(file = file)
                               load("CRIS_template_mouse.RData")
                               NTP_CRIS <- ntp(emat=datacris$mat_exp_adj,CIRS_template, seed=367707,doPlot = TRUE)
                               dev.off()
                             })
                         }
                       })
                       
                     }
                     #############################################################################################
                     else if (chcris  == "Human"){  
                       ####### Get the normalized data from the first part
                       Collapsed_cris <- Collapsed_data[,-c(1,2)]
                       col_name8<- paste(datar$n1, sep="")
                       col_name9<-paste(datar$n2, sep="")
                       
                       names(Collapsed_cris)[1:length(datar$n1)]<-col_name8
                       names(Collapsed_cris)[length(datar$n1)+1:length(datar$n2)]<-col_name9
                       
                       
                       expcris <- data.matrix(Collapsed_cris)
                       
                       # adjust matrix before classification
                       exp_adj <- ematAdjust(expcris)
                       datacris$exp_adj<-exp_adj
                       
                       load("Human.CRIS.template.RData")
                       
                       
                       NTP_CRIS <- ntp(emat=exp_adj, Human.CRIS.template, seed=367707,doPlot =F)
                       NTP_CRIS <- tibble::rownames_to_column(NTP_CRIS,var = "ID")
                       
                       crisl <-(length(colnames(NTP_CRIS)))
                       
                       output$CRIS<- DT::renderDataTable({NTP_CRIS},
                                                         extensions = 'FixedColumns', 
                                                         options = list(pageLength = 20, autoWidth = FALSE,
                                                                        scrollX = TRUE,
                                                                        scrollCollapse = TRUE,
                                                                        lengthMenu = c(20, 50, 100),
                                                                        columnDefs = list(list(targets = c(3:crisl), searchable = FALSE))),
                                                         editable = FALSE)
                       
                       
                       
                       output$TableDownloadCRIS<- downloadHandler(
                         filename = function(){paste("human_DATA_NTP_CMS_results",'.csv')},
                         content = function(file){
                           write.csv(NTP_CRIS, file, sep = ",", quote = F)
                         })
                       
                       
                       
                       observe({  
                         widthcris <- reactive ({ input$plot_widthcris})
                         heightcris<- reactive ({ input$plot_heightcris}) 
                         output$CRIS_Plot<- renderPlot(width = widthcris, height = heightcris,{
                           NTP_CRIS <- ntp(emat=exp_adj, Human.CRIS.template, seed=367707,doPlot = TRUE)
                         })  
                         
                       })
                       
                       observe({
                         chdcris <- input$optioncris
                         if ( chdcris  == "1"){ 
                           output$PlotDownloadCRIS<- downloadHandler(
                             file = "CRIS-Plot.png" , # variable with filename
                             content = function(file) {
                               png(file = file)
                               load("Human.CRIS.template.RData")
                               NTP_CRIS <- ntp(emat= datacris$exp_adj, Human.CRIS.template, seed=367707,doPlot = TRUE)
                               dev.off()
                             })
                         } else if  ( chdcris == "2"){ 
                           
                           output$PlotDownloadCRIS <- downloadHandler(
                             file = "CRIS-Plot.svg" , # variable with filename
                             content = function(file) {
                               svg(file = file)
                               load("Human.CRIS.template.RData")
                               NTP_CRIS <- ntp(emat= datacris$exp_adj, Human.CRIS.template, seed=367707,doPlot = TRUE)
                               dev.off()
                             })
                         }
                       })
                       
                     }
                     
                   })   
    }) 
    
    #########################################################################################################################
    #Classifier_CMS
    datacms <- reactiveValues()
    observeEvent(input$CMSAnalysis,{
      withProgress(message = "Preparing CMS Classifier Table and Plot, please wait",
                   {
                     Collapsed_data <- datar$Coll
                     head(Collapsed_data)
                   
                     
                     ########## Mouse data############################################   
                     chcms <- input$optionCMS
                     if (chcms  == "Mouse"){   
                       ####### Get the normalized data from the first part
                       Collapsed_cms <- Collapsed_data[,-c(1,2)]
                       col_name8<- paste(datar$n1, sep="")
                       col_name9<-paste(datar$n2, sep="")
                       
                       names(Collapsed_cms)[1:length(datar$n1)]<-col_name8
                       names(Collapsed_cms)[length(datar$n1)+1:length(datar$n2)]<-col_name9
                       
                       
                       expcms <- data.matrix(Collapsed_cms)
                       
                       # adjust matrix before classification
                       mat_exp_adj <- ematAdjust(expcms)
                       datacms$mat_exp_adj<-mat_exp_adj
                       
                       # load mouse CMS template
                       load("template.CMS.A.RData") # template.CMS.A
                       
                       
                       
                       # CMS NTP and save as text file
                       NTP_CMS <- ntp(emat=mat_exp_adj,template.CMS.A, seed=367707,doPlot = F)
                       NTP_CMS <-tibble:: rownames_to_column(NTP_CMS,var = "ID")
                       # write.table(NTP_CMS,file=CMS_res_nm,sep='\t',row.names = FALSE)
                       
                       
                       
                       cmsl <-(length(colnames(NTP_CMS)))
                       
                       output$CMS<- DT::renderDataTable({NTP_CMS},
                                                        extensions = 'FixedColumns', 
                                                        options = list(pageLength = 20, autoWidth = FALSE,
                                                                       scrollX = TRUE,
                                                                       scrollCollapse = TRUE,
                                                                       lengthMenu = c(20, 50, 100),
                                                                       columnDefs = list(list(targets = c(3:cmsl), searchable = FALSE))),
                                                        editable = FALSE)
                       
                       
                       
                       output$TableDownloadCMS<- downloadHandler(
                         filename = function(){paste("MOUSE_DATA_NTP_CMS_results",'.csv')},
                         content = function(file){
                           write.csv(NTP_CMS, file, sep = ",", quote = F)
                         })
                       
                       
                       
                       observe({  
                         widthcms <- reactive ({ input$plot_widthcms})
                         heightcms<- reactive ({ input$plot_heightcms}) 
                         output$CMS_Plot<- renderPlot(width = widthcms, height = heightcms,{
                           NTP_CMS<-ntp(emat=mat_exp_adj,template.CMS.A, seed=367707,doPlot = T)
                         })  
                         
                       })
                       
                       observe({
                         chdcms <- input$optioncms
                         if ( chdcms  == "1"){ 
                           output$PlotDownloadCMS <- downloadHandler(
                             file = "CMS-Plot.png" , # variable with filename
                             content = function(file) {
                               png(file = file)
                               load("template.CMS.A.RData") # template
                               NTP_CMS <- ntp(emat=datacms$mat_exp_adj,template.CMS.A, seed=367707,doPlot = T)
                               dev.off()
                             })
                         } else if  ( chdcms == "2"){ 
                           
                           output$PlotDownloadCMS <- downloadHandler(
                             file = "CMS-Plot.svg" , # variable with filename
                             content = function(file) {
                               svg(file = file)
                               load("template.CMS.A.RData") # template
                               NTP_CMS <- ntp(emat=datacms$mat_exp_adj,template.CMS.A, seed=367707,doPlot = T)
                               dev.off()
                             })
                         }
                       })
                     }
                     #############################################################################################
                     else if (chcms  == "Human"){  
                       ####### Get the normalized data from the first part
                       Collapsed_cms <- Collapsed_data[,-c(1,2)]
                       col_name8<- paste(datar$n1, sep="")
                       col_name9<-paste(datar$n2, sep="")
                       
                       names(Collapsed_cms)[1:length(datar$n1)]<-col_name8
                       names(Collapsed_cms)[length(datar$n1)+1:length(datar$n2)]<-col_name9
                       
                       
                       expcms <- data.matrix(Collapsed_cms)
                       
                       # adjust matrix before classification
                       exp_adj <- ematAdjust(expcms)
                       datacms$exp_adj<-exp_adj
                       
                       load("Human.CMS.template.RData") # Human.CMS.template
                       
                       
                       # CMS NTP and save as text file
                       NTP_CMS <- ntp(emat=exp_adj, Human.CMS.template, seed=367707,doPlot = F)
                       # datacms$NTP_CMS <- NTP_CMS
                       NTP_CMS <- rownames_to_column(NTP_CMS,var = "ID")
                       # write.table(NTP_CMS,file=CMS_res_nm,sep='\t',row.names = FALSE)
                       
                       cmsl <-(length(colnames(NTP_CMS)))
                       
                       output$CMS<- DT::renderDataTable({NTP_CMS},
                                                        extensions = 'FixedColumns', 
                                                        options = list(pageLength = 20, autoWidth = FALSE,
                                                                       scrollX = TRUE,
                                                                       scrollCollapse = TRUE,
                                                                       lengthMenu = c(20, 50, 100),
                                                                       columnDefs = list(list(targets = c(3:cmsl), searchable = FALSE))),
                                                        editable = FALSE)
                       
                       
                       
                       output$TableDownloadCMS<- downloadHandler(
                         filename = function(){paste("human_DATA_NTP_CMS_results",'.csv')},
                         content = function(file){
                           write.csv(NTP_CMS, file, sep = ",", quote = F)
                         })
                       
                       
                       
                       observe({  
                         widthcms <- reactive ({ input$plot_widthcms})
                         heightcms<- reactive ({ input$plot_heightcms}) 
                         output$CMS_Plot<- renderPlot(width = widthcms, height = heightcms,{
                           NTP_CMS <- ntp(emat=exp_adj, Human.CMS.template, seed=367707,doPlot = T)
                         })  
                         
                       })
                       
                       observe({
                         chdcms <- input$optioncms
                         if ( chdcms  == "1"){ 
                           output$PlotDownloadCMS <- downloadHandler(
                             file = "CMS-Plot.png" , # variable with filename
                             content = function(file) {
                               png(file = file)
                               load("Human.CMS.template.RData")
                               NTP_CMS <- ntp(emat=datacms$exp_adj,Human.CMS.template, seed=367707,doPlot = T)
                               dev.off()
                             })
                         } else if  ( chdcms == "2"){ 
                           
                           output$PlotDownloadCMS <- downloadHandler(
                             file = "CMS-Plot.svg" , # variable with filename
                             content = function(file) {
                               svg(file = file)
                               load("Human.CMS.template.RData")
                               NTP_CMS <- ntp(emat=datacms$exp_adj,Human.CMS.template, seed=367707,doPlot = T)
                               dev.off()
                             })
                         }
                       })
                       
                     }
                     
                   })   
    }) 
    
######################################################################################################################
#######Heat Map##############
          datar <- reactiveValues()
          observeEvent(input$FilterAnalysis,{
            withProgress(message = "Preparing, please wait",
                         {
                         
                           fs<-data$exp
                           head(fs)
                           rownames(fs)<-fs[,1]
                           fs<- fs[order(fs$ID, decreasing = FALSE), ]
                           
                           file_to_read2=input$file2
                           ds <- read.delim(file_to_read2$datapath,sep=input$Sep1, header=T)
                           names(ds)[1] <- "Types"
                           names(ds)[2] <- "Samples"
                           names(ds)[3] <- "Labels"
                           
                           melfs<-melt(fs)
                           
                           colnames(melfs) <- c("ID", "symbol", "Samples","value")
                           
                           fss<-join(melfs, ds, by = "Samples", type = "inner")
                           
                           if (as.numeric(length(as.character(input$Group1comp)))==1){
                             g1 <- subset(fss, Types==as.character(input$Group1comp)[1])
                             g1 = subset(g1, select = -c(Types,Labels) )
                             colnames(g1) <- c("ID", "symbol", "variable","value")
                             g1<-cast(g1)
                             g1 = subset(g1, select = -c(ID,symbol) )
                            
                             Big1<-g1
                             datar$lg1<-length(g1)
                           }
                           
                           else if (as.numeric(length(as.character(input$Group1comp)))==2){
                             g1 <- subset(fss, Types==as.character(input$Group1comp)[1])
                             g1 = subset(g1, select = -c(Types,Labels) )
                             colnames(g1) <- c("ID", "symbol", "variable","value")
                             g1<-cast(g1)
                             g1 = subset(g1, select = -c(ID,symbol) )
                             datar$lg1<-length(g1)
                             
                             g2 <- subset(fss, Types==as.character(input$Group1comp)[2])
                             g2 = subset(g2, select = -c(Types,Labels) )
                             colnames(g2) <- c("ID", "symbol", "variable","value")
                             g2<-cast(g2)
                             g2 = subset(g2, select = -c(ID,symbol) )
                             datar$lg2<-length(g2)
                            
                             Big1<-cbind(g1,g2)
                           }
                           
                           
                           else if(as.numeric(length(as.character(input$Group1comp)))==3){
                             g1 <- subset(fss, Types==as.character(input$Group1comp)[1])
                             g1 = subset(g1, select = -c(Types,Labels) )
                             colnames(g1) <- c("ID", "symbol", "variable","value")
                             g1<-cast(g1)
                             g1 = subset(g1, select = -c(ID,symbol) )
                             datar$lg1<-length(g1)
                             g2 <- subset(fss, Types==as.character(input$Group1comp)[2])
                             g2 = subset(g2, select = -c(Types,Labels) )
                             colnames(g2) <- c("ID", "symbol", "variable","value")
                             g2<-cast(g2)
                             g2 = subset(g2, select = -c(ID,symbol) )
                             datar$lg2<-length(g2)
                             
                             g3 <- subset(fss, Types==as.character(input$Group1comp)[3])
                             g3 = subset(g3, select = -c(Types,Labels) )
                             colnames(g3) <- c("ID", "symbol", "variable","value")
                             g3<-cast(g3)
                             g3 = subset(g3, select = -c(ID,symbol) )
                             datar$lg3<-length(g3)
                             
                          
                             Big1<-cbind(g1,g2,g3)
                           }
                           
                           else if (as.numeric(length(as.character(input$Group1comp)))==4){
                             g1 <- subset(fss, Types==as.character(input$Group1comp)[1])
                             g1 = subset(g1, select = -c(Types,Labels) )
                             colnames(g1) <- c("ID", "symbol", "variable","value")
                             g1<-cast(g1)
                             g1 = subset(g1, select = -c(ID,symbol) )
                             datar$lg1<-length(g1)
                             g2 <- subset(fss, Types==as.character(input$Group1comp)[2])
                             g2 = subset(g2, select = -c(Types,Labels) )
                             colnames(g2) <- c("ID", "symbol", "variable","value")
                             g2<-cast(g2)
                             g2 = subset(g2, select = -c(ID,symbol) )
                             datar$lg2<-length(g2)
                             g3 <- subset(fss, Types==as.character(input$Group1comp)[3])
                             g3 = subset(g3, select = -c(Types,Labels) )
                             colnames(g3) <- c("ID", "symbol", "variable","value")
                             g3<-cast(g3)
                             g3 = subset(g3, select = -c(ID,symbol) )
                             datar$lg3<-length(g3)
                             g4 <- subset(fss, Types==as.character(input$Group1comp)[4])
                             g4 = subset(g4, select = -c(Types,Labels) )
                             colnames(g4) <- c("ID", "symbol", "variable","value")
                             g4<-cast(g4)
                             g4 = subset(g4, select = -c(ID,symbol) )
                             datar$lg4<-length(g4)
                             
                             Big1<-cbind(g1,g2,g3,g4)
                             
                           }
                           else if (as.numeric(length(as.character(input$Group1comp)))==5){
                             g1 <- subset(fss, Types==as.character(input$Group1comp)[1])
                             g1 = subset(g1, select = -c(Types,Labels) )
                             colnames(g1) <- c("ID", "symbol", "variable","value")
                             g1<-cast(g1)
                             g1 = subset(g1, select = -c(ID,symbol) )
                             datar$lg1<-length(g1)
                             g2 <- subset(fss, Types==as.character(input$Group1comp)[2])
                             g2 = subset(g2, select = -c(Types,Labels) )
                             colnames(g2) <- c("ID", "symbol", "variable","value")
                             g2<-cast(g2)
                             g2 = subset(g2, select = -c(ID,symbol) )
                             datar$lg2<-length(g2)
                             g3 <- subset(fss, Types==as.character(input$Group1comp)[3])
                             g3 = subset(g3, select = -c(Types,Labels) )
                             colnames(g3) <- c("ID", "symbol", "variable","value")
                             g3<-cast(g3)
                             g3 = subset(g3, select = -c(ID,symbol) )
                             datar$lg3<-length(g3)
                             g4 <- subset(fss, Types==as.character(input$Group1comp)[4])
                             g4 = subset(g4, select = -c(Types,Labels) )
                             colnames(g4) <- c("ID", "symbol", "variable","value")
                             g4<-cast(g4)
                             g4 = subset(g4, select = -c(ID,symbol) )
                             datar$lg4<-length(g4)
                             g5 <- subset(fss, Types==as.character(input$Group1comp)[5])
                             g5 = subset(g5, select = -c(Types,Labels) )
                             colnames(g5) <- c("ID", "symbol", "variable","value")
                             g5<-cast(g5)
                             g5 = subset(g5, select = -c(ID,symbol) )
                             datar$lg5<-length(g5)
                             Big1<-cbind(g1,g2,g3,g4,g5)
                             
                           }
                           
                          
                           if (as.numeric(length(as.character(input$Group2comp)))==1){
                             
                             k1 <- subset(fss, Types==as.character(input$Group2comp)[1])
                             k1 = subset(k1, select = -c(Types,Labels) )
                             colnames(k1) <- c("ID", "symbol", "variable","value")
                             k1<-cast(k1)
                             k1 = subset(k1, select = -c(ID,symbol) )
                             datar$lk1<-length(k1)
                             
                             
                             Big2<-k1
                           }
                           
                           else if(as.numeric(length(as.character(input$Group2comp)))==2){
                             
                             k1 <- subset(fss, Types==as.character(input$Group2comp)[1])
                             k1 = subset(k1, select = -c(Types,Labels) )
                             colnames(k1) <- c("ID", "symbol", "variable","value")
                             k1<-cast(k1)
                             k1 = subset(k1, select = -c(ID,symbol) )
                             datar$lk1<-length(k1)
                             k2 <- subset(fss, Types==as.character(input$Group2comp)[2])
                             k2 = subset(k2, select = -c(Types,Labels) )
                             colnames(k2) <- c("ID", "symbol", "variable","value")
                             k2<-cast(k2)
                             k2 = subset(k2, select = -c(ID,symbol) )
                             datar$lk2<-length(k2)
                             
                             
                             Big2<-cbind(k1,k2)
                           }
                           
                           
                           else if (as.numeric(length(as.character(input$Group2comp)))==3){
                             
                             k1 <- subset(fss, Types==as.character(input$Group2comp)[1])
                             k1 = subset(k1, select = -c(Types,Labels) )
                             colnames(k1) <- c("ID", "symbol", "variable","value")
                             k1<-cast(k1)
                             k1 = subset(k1, select = -c(ID,symbol) )
                             datar$lk1<-length(k1)
                             k2 <- subset(fss, Types==as.character(input$Group2comp)[2])
                             k2 = subset(k2, select = -c(Types,Labels) )
                             colnames(k2) <- c("ID", "symbol", "variable","value")
                             k2<-cast(k2)
                             k2 = subset(k2, select = -c(ID,symbol) )
                             datar$lk2<-length(k2)
                             
                             k3 <- subset(fss, Types==as.character(input$Group2comp)[3])
                             k3 = subset( k3, select = -c(Types,Labels) )
                             colnames( k3) <- c("ID", "symbol", "variable","value")
                             k3<-cast( k3)
                             k3 = subset(k3, select = -c(ID,symbol) )
                             datar$lk3<-length( k3)
                             
                             Big2<-cbind(k1,k2,k3)
                           }
                           
                           
                           else if (as.numeric(length(as.character(input$Group2comp)))==4){
                             k1 <- subset(fss, Types==as.character(input$Group2comp)[1])
                             k1 = subset(k1, select = -c(Types,Labels) )
                             colnames(k1) <- c("ID", "symbol", "variable","value")
                             k1<-cast(k1)
                             k1 = subset(k1, select = -c(ID,symbol) )
                             datar$lk1<-length(k1)
                             k2 <- subset(fss, Types==as.character(input$Group2comp)[2])
                             k2 = subset(k2, select = -c(Types,Labels) )
                             colnames(k2) <- c("ID", "symbol", "variable","value")
                             k2<-cast(k2)
                             k2 = subset(k2, select = -c(ID,symbol) )
                             datar$lk2<-length(k2)
                             
                             k3 <- subset(fss, Types==as.character(input$Group2comp)[3])
                             k3 = subset( k3, select = -c(Types,Labels) )
                             colnames( k3) <- c("ID", "symbol", "variable","value")
                             k3<-cast( k3)
                             k3 = subset(k3, select = -c(ID,symbol) )
                             datar$lk3<-length( k3)
                             
                            
                             
                             k4<- subset(fss, Types==as.character(input$Group2comp)[4])
                             k4 = subset(k4, select = -c(Types,Labels) )
                             colnames(k4) <- c("ID", "symbol", "variable","value")
                             k4<-cast(k4)
                             k4 = subset(k4, select = -c(ID,symbol) )
                             datar$lk4<-length(k4)
                             Big2<-cbind(k1,k2,k3,k4)
                           }
                           
                           else if (as.numeric(length(as.character(input$Group2comp)))==5){
                             
                             k1 <- subset(fss, Types==as.character(input$Group2comp)[1])
                             k1 = subset(k1, select = -c(Types,Labels) )
                             colnames(k1) <- c("ID", "symbol", "variable","value")
                             k1<-cast(k1)
                             k1 = subset(k1, select = -c(ID,symbol) )
                             datar$lk1<-length(k1)
                             k2 <- subset(fss, Types==as.character(input$Group2comp)[2])
                             k2 = subset(k2, select = -c(Types,Labels) )
                             colnames(k2) <- c("ID", "symbol", "variable","value")
                             k2<-cast(k2)
                             k2 = subset(k2, select = -c(ID,symbol) )
                             datar$lk2<-length(k2)
                             
                             k3 <- subset(fss, Types==as.character(input$Group2comp)[3])
                             k3 = subset( k3, select = -c(Types,Labels) )
                             colnames( k3) <- c("ID", "symbol", "variable","value")
                             k3<-cast( k3)
                             k3 = subset(k3, select = -c(ID,symbol) )
                             datar$lk3<-length( k3)
                             
                             k4<- subset(fss, Types==as.character(input$Group2comp)[4])
                             k4 = subset(k4, select = -c(Types,Labels) )
                             colnames(k4) <- c("ID", "symbol", "variable","value")
                             k4<-cast(k4)
                             k4 = subset(k4, select = -c(ID,symbol) )
                             datar$lk4<-length(k4)
                             
                             k5<- subset(fss, Types==as.character(input$Group2comp)[5])
                             k5 = subset(k5, select = -c(Types,Labels) )
                             colnames(k5) <- c("ID", "symbol", "variable","value")
                             k5<-cast(k5)
                             k5 = subset(k5, select = -c(ID,symbol) )
                             datar$lk5<-length(k5)
                             Big2<-cbind(k1,k2,k3,k4,k5)
                           }
                           
                           
                           datar$n1<-colnames(Big1)
                           
                           col_name <- paste("A",1:length(Big1), sep="")
                           names(Big1) <-col_name
                           head(Big1)
                           
                           datar$n2<-colnames(Big2)
                           
                           col_name2 <- paste("B",1:length(Big2), sep="")
                           names(Big2) <-col_name2
                           head(Big2)
                           
                           
                           
                           main1<-cbind(Big1,Big2)
                           main<-cbind(fs[c(1,2)],main1)
                           
                           # print(main)
                           
                           
                           
                         })# end a progress
            
            gr <- c(rep("A",length(Big1) ), rep("B",length(Big2)))
            grss = factor(gr, levels = unique(gr))
         

           
                           
                          num <- sapply( main1, is.integer)
                        if(all(num)!= TRUE){
                            shinyalert("Warning!", "Input Data includes Decimal number", type = "error") 
                            Sys.sleep(10)
                          }
                       
                 if (all(num)== TRUE){
                         withProgress(message = "DESeq2 Calculations, please wait",
                        {
                           colData <- data.frame(group=grss)
                           rownames(colData) <- colnames(main[,-c(1,2)])
                           cds <- DESeqDataSetFromMatrix(main[,-c(1,2)], colData, design = ~group)
                           class(cds)
                           cds <- DESeq(cds)
                          
                           DEGs<- results(cds, c("group", levels(grss)))
                           
                           DEGs$symbol<- fs$symbol
                           datar$DEGs<-DEGs
                           datar$head<- head(DEGs)
                           ## add symbol column to DEGs file
                           datar$symbol <- DEGs$symbol
                           datar$cds<-cds
                           datar$ds <- ds
                           ###DEGS Table  #######################################
                           pdat <- data.frame(datar$DEGs)
                           minus_log10_pvalue <- -log(pdat$pvalue, base = 10)
                          
                           pdat<- cbind(pdat,minus_log10_pvalue)
                           pdat<- pdat[order(pdat$padj, decreasing = TRUE), ]
                          
                           
                           datar$Dm <- pdat[, c("symbol","log2FoldChange","minus_log10_pvalue","pvalue", "padj","baseMean","lfcSE","stat")]
                           ##################################################
                           
                           ######################################################
                           dim(datar$head)
                           DEGs <- na.omit(datar$DEGs)
                           
                           DEGs$selectedRowID <- rownames(DEGs)
                           normalized_df <- counts(datar$cds, normalized=TRUE)
                           #### log transform
                           normalized_df <- log2(normalized_df+1)
                           head(normalized_df)
                           rowGroup <- as.vector(datar$symbol)
                           rowID <- as.vector(rownames(normalized_df))
                           collapse.object <- collapseRows(datET=normalized_df, rowGroup=rowGroup, rowID=rowID, method = "MaxMean")
                           class(collapse.object)
                           names(collapse.object)
                           Collapsed_data <- data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
                           head(Collapsed_data)
                           
                           datar$Collapsed_data<-Collapsed_data
                           
                           s=factor(colnames(head(normalized_df)))
                           # print(s)
                           #### Mathing
                           match <- join(Collapsed_data, as.data.frame(DEGs), by = "selectedRowID", type = "inner")
                           dim(match)
                           rownames(match) <- match$symbol
                           datar$match <- match
                           datar$s<-s
                           datar$Coll <- Collapsed_data
                           
                           number1<-match %>% filter(log2FoldChange> 2)
                           number1<-count(number1%>% filter(padj < 0.05))
                           number1<-as.numeric(number1$n)
                           
                           number2<-match %>% filter(log2FoldChange< -2)
                           number2<-count(number2%>% filter(padj < 0.05))
                           number2<-as.numeric(number2$n)
                                 
                           output$upregulated <- renderText({
                             isolate({
                               upregulated<-match %>% filter(log2FoldChange> 2)
                               upregulated<-count(upregulated%>% filter(padj < 0.05))
                               paste("Number of Up-regulated genes:", as.numeric(upregulated$n), sep="")
                             })
                           })
                          
                           output$Downregulated <- renderText({
                             isolate({
                               Downregulated<-match %>% filter(log2FoldChange< -2)
                               Downregulated<-count(Downregulated%>% filter(padj < 0.05))
                               paste("Number of Down-regulated genes:", as.numeric(Downregulated$n),sep="")
                             })
                           })
                           
                        if(number1==0 && number2==0){
                          shinyalert("Warning!", "No gene found for this Threshold, Change the values and click on Create plot.", type = "error")
                          
                        }
                           
                       else if(number1!=0 && number2!=0){
                        shinyalert("All good!", "Please click on create Plot", type = "success")
                             
                        }  
                              
                     })# with progresss 
                      
                 }
          
          })  
                      
            t<-reactiveValues()
            observeEvent(input$findheat,{
                
                           DEGs_top <-subset(datar$match ,log2FoldChange> as.numeric(input$MSGfc1) | log2FoldChange< as.numeric(input$MSGfc))
                           DEGs_top <- subset(DEGs_top, padj <as.numeric (input$MSGpval))  
                           datar$Gene<-DEGs_top$group
                           DEGs_top <-subset( DEGs_top, select = -c(group,selectedRowID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj, symbol))
                           colnames(DEGs_top) <-datar$s
                           
                           output$upregulated <- renderText({
                             isolate({
                               upregulated<-datar$match%>% filter(log2FoldChange> as.numeric(input$MSGfc1))
                               upregulated<-count(upregulated%>% filter(padj < as.numeric (input$MSGpval)))
                               paste("Number of Up-regulated genes:", as.numeric(upregulated$n), sep="")
                             })
                           })
                           
                           output$Downregulated <- renderText({
                             isolate({
                               Downregulated<-datar$match%>% filter(log2FoldChange< as.numeric(input$MSGfc))
                               Downregulated<-count(Downregulated%>% filter(padj <as.numeric (input$MSGpval)))
                               paste("Number of Down-regulated genes:", as.numeric(Downregulated$n),sep="")
                             })
                           })   
                           
                           
                           
                           DEGs_top_scale <-  DEGs_top-rowMeans(DEGs_top) 
                           datar$DEGs_top_scale5<-DEGs_top_scale
                           col_nameA<- paste(datar$n1, sep="")
                           col_nameB<-paste(datar$n2, sep="")
                           
                           names(DEGs_top_scale)[1:length(datar$n1)]<-col_nameA
                           names(DEGs_top_scale)[length(datar$n1)+1:length(datar$n2)]<-col_nameB
                           
                           datar$DEGs_top_scale<-DEGs_top_scale
                           
                           
                           ma <-round(max(DEGs_top_scale, na.rm=TRUE))
                           mi <-round(min(DEGs_top_scale, na.rm=TRUE))
                           
                           gradient_col <- ggplot2::scale_fill_gradient2(
                             low = "blue", high = "red", 
                             midpoint = 0.0, limits = c(mi, ma)
                           )
                           
             
            
            
            observe({
              chl <- input$optionl
              chheat <- input$optionheat
              
              output$heatmapMSGPlot <- renderPlotly({
                withProgress(message = "Drawing Heatmap plot, please wait",
                             {
                               
                               if (chl == "1"){
                                 phe<-heatmaply(DEGs_top_scale, scale_fill_gradient_fun = gradient_col,showticklabels = T, dendrogram =as.character(input$Group7comp),main =paste(as.character(input$Name)), margins =c(0,0,50,0),
                                                fontsize_col = 14)
                                 phe<-phe %>% layout(height = input$plot_heightheat, width = input$plot_widthheat)
                               
                                 if  (chheat == "1"){
                                   phe <- phe %>% config(
                                     toImageButtonOptions= list(filename = as.character(input$Name),res =350, units = "px"))
                                   
                                 }
                                 else if  (chheat == "2"){
                                   phe <- phe %>% config(
                                     toImageButtonOptions= list(format = "svg",filename = as.character(input$Name)))
                                  
                                 }
                               }
                               else if (chl == "2"){
                                 phe<-heatmaply(DEGs_top_scale, scale_fill_gradient_fun = gradient_col,showticklabels=c(TRUE, FALSE),dendrogram =as.character(input$Group7comp),main =paste(as.character(input$Name)),margins =c(0,0,50,0),
                                                fontsize_col = 14)
                                 phe<-phe %>% layout(height = input$plot_heightheat, width = input$plot_widthheat)
                                
                                 
                                 if  (chheat == "1"){
                                   phe <- phe %>% config(
                                     toImageButtonOptions= list(filename = as.character(input$Name)))
                                   
                                 }
                                 else if (chheat == "2"){
                                   phe <- phe %>% config(
                                     toImageButtonOptions= list(format ="svg",filename = as.character(input$Name)))
                                   
                                 }
                                 
                               }
                               
                             })#progress bar
              })
            })# End observe
            
            })
        
            
         # })#END OF HEatmap
          
          
          ###################Heatmap table
          observeEvent(input$FilterAnalysis,{
            DGStable <- na.omit(datar$Dm)
            DGStable <-subset(DGStable, select = -c(minus_log10_pvalue))
            siz <-(length(colnames(DGStable)))
          
            output$DGStable <- DT::renderDataTable({DGStable},
                                                   extensions = 'FixedColumns', 
                                                   options = list(pageLength = 20, autoWidth = FALSE,
                                                                  scrollX = TRUE,
                                                                  scrollCollapse = TRUE,
                                                                  lengthMenu = c(20, 30, 50),
                                                                  columnDefs = list(list(targets = c(2:siz), searchable = FALSE))),
                                                   editable = FALSE)
            
          })
          output$DGEtable <- downloadHandler(
            filename = function(){paste("DGE",'.txt')},
            content = function(file){
              write.csv(na.omit(datar$Dm), file, row.names = TRUE)
            }) 
          
          
          
          ###################    
          
 ################################################################################################################# 
          #heatmap selected list
          observe({
            chzscore <-input$ optiongroupsample
            chselectplot <- input$ optionsheatplot
            if  (chzscore== "1"){
              if  (chselectplot== "1"){
                datselecheat <- reactiveValues()
                observeEvent(input$FilterAnalysis,{
                  match1 <- datar$match 
                 
                  Gene<-match1$group
                  s1<-datar$s
                  DEGs_top<-subset(match1, select = -c(group,selectedRowID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj, symbol))
                  colnames(DEGs_top) <- s1
                  DEGs_top_scale2 <-  DEGs_top-rowMeans(DEGs_top) 
                  col_nameA<- paste(datar$n1, sep="")
                  col_nameB<-paste(datar$n2, sep="")
                  names(DEGs_top_scale2)[1:length(datar$n1)]<-col_nameA
                  names(DEGs_top_scale2)[length(datar$n1)+1:length(datar$n2)]<-col_nameB
                  DEGs_top_scale3<-cbind(Gene,DEGs_top_scale2)
                  names(DEGs_top_scale3)[names(DEGs_top_scale3) == "Gene"] <- "Genes"
                  
                
                  output$Heatmapselect = DT::renderDataTable(DEGs_top_scale3,
                                                             #server = FALSE,
                                                             extensions = 'FixedColumns', 
                                                             options = list(pageLength = 10, autoWidth = FALSE,
                                                                            scrollX = TRUE,
                                                                            scrollCollapse = TRUE,
                                                                            lengthMenu = c(10,50, 100),
                                                                            columnDefs = list(list(targets = c(1:length( DEGs_top_scale3)), 
                                                                                                   searchable = FALSE))
                                                             ))
                  ma <-round(max(DEGs_top_scale2, na.rm=TRUE))
                  # print(ma)
                  mi <-round(min(DEGs_top_scale2, na.rm=TRUE))
                  # print(mi)
                  gradient_col <- ggplot2::scale_fill_gradient2(
                    low = "blue", high = "red", 
                    midpoint = 0.0, limits = c(mi, ma)
                  )
                  
                  # highlight selected rows in the scatterplot
                  observe({
                    chselectl <- input$optionlselectheat
                    chselectheat <- input$optionselecheat
                    output$dotplot = renderPlotly({
                      if ( chselectl  == "1"){
                        s = input$Heatmapselect_rows_selected
                        vals$s<-s
                        par(mar = c(4, 4, 1, .1))
                        heat<- heatmaply( DEGs_top_scale2[vals$s, , drop = FALSE], scale_fill_gradient_fun = gradient_col,showticklabels = T,dendrogram =as.character(input$Group8comp), main =paste(as.character(input$Nameheat)),fontsize_col = 14,fontsize_row = 14) 
                        heat<-heat %>% layout(height = input$plot_heightselectheat, width = input$plot_widthselectheat) 
                        if  (chselectheat == "1"){
                          heat<- heat %>% config(
                            toImageButtonOptions= list(filename = as.character(input$Nameheat)))
                          #p <- data8$plot 
                        }
                        else if  (chselectheat == "2"){
                          heat3 <- heat %>% config(
                            toImageButtonOptions= list(format = "svg",filename = as.character(input$Nameheat)))
                          #p <- data8$plot 
                        }
                      }
                      else if (chselectl == "2"){
                        heat2<- heatmaply( DEGs_top_scale2[vals$s, , drop = FALSE], scale_fill_gradient_fun = gradient_col,showticklabels=c(TRUE, FALSE),dendrogram =as.character(input$Group8comp), main =paste(as.character(input$Nameheat)),fontsize_col = 14,fontsize_row = 14) 
                        heat2<-heat2 %>% layout(height = input$plot_heightselectheat, width = input$plot_widthselectheat)  
                        vals$heat2<-heat2
                        if  (chselectheat== "1"){
                          heat2 <- heat2 %>% config(
                            toImageButtonOptions= list(filename = as.character(input$Nameheat)))
                          
                        }
                        else if (chselectheat == "2"){
                          heat2 <- vals$heat2 %>% config(
                            toImageButtonOptions= list(format ="svg",filename = as.character(input$Nameheat)))
                          
                        }
                        
                      }                                                                                                                    #,margins =c(0,0,50,0)))
                      
                    })
                  })
                  
                })
              }
              
              else if (chselectplot== "2"){
                datselecheat <- reactiveValues()
                observeEvent(  input$FilterAnalysis,{
                  DEGs_top_scale <- datar$DEGs_top_scale
                  DEGs_top_scale1<-cbind(datar$Gene,DEGs_top_scale)
                  names(DEGs_top_scale1)[names(DEGs_top_scale1) == "datar$Gene"] <- "Genes"
                  # print(DEGs_top_scale1)
                  output$Heatmapselect = DT::renderDataTable(DEGs_top_scale1,
                                                             #server = FALSE,
                                                             extensions = 'FixedColumns', 
                                                             options = list(pageLength = 10, autoWidth = FALSE,
                                                                            scrollX = TRUE,
                                                                            scrollCollapse = TRUE,
                                                                            lengthMenu = c(10,50, 100),
                                                                            columnDefs = list(list(targets = c(1:length( DEGs_top_scale1)), 
                                                                                                   searchable = FALSE))
                                                             ))
                  ma <-round(max(DEGs_top_scale, na.rm=TRUE))
                  mi <-round(min(DEGs_top_scale, na.rm=TRUE))
                  
                  gradient_col <- ggplot2::scale_fill_gradient2(
                    low = "blue", high = "red", 
                    midpoint = 0.0, limits = c(mi, ma)
                  )
                  
                  # highlight selected rows in the scatterplot
                  observe({
                    chselectl <- input$optionlselectheat
                    chselectheat <- input$optionselecheat
                    output$dotplot = renderPlotly({
                      if ( chselectl  == "1"){
                        s = input$Heatmapselect_rows_selected
                        vals$s<-s
                        par(mar = c(4, 4, 1, .1))
                        heat<- heatmaply( DEGs_top_scale[vals$s, , drop = FALSE], scale_fill_gradient_fun = gradient_col,showticklabels = T,dendrogram =as.character(input$Group8comp), main =paste(as.character(input$Nameheat)),fontsize_col = 14,fontsize_row = 14) 
                        heat<-heat %>% layout(height = input$plot_heightselectheat, width = input$plot_widthselectheat) 
                        if  (chselectheat == "1"){
                          heat<- heat %>% config(
                            toImageButtonOptions= list(filename = as.character(input$Nameheat)))
                          #p <- data8$plot 
                        }
                        else if  (chselectheat == "2"){
                          heat3 <- heat %>% config(
                            toImageButtonOptions= list(format = "svg",filename = as.character(input$Nameheat)))
                          #p <- data8$plot 
                        }
                      }
                      else if (chselectl == "2"){
                        heat2<- heatmaply( DEGs_top_scale[vals$s, , drop = FALSE], scale_fill_gradient_fun = gradient_col,showticklabels=c(TRUE, FALSE),dendrogram =as.character(input$Group8comp), main =paste(as.character(input$Nameheat)),fontsize_col = 14,fontsize_row = 14) 
                        heat2<-heat2 %>% layout(height = input$plot_heightselectheat, width = input$plot_widthselectheat)  
                        vals$heat2<-heat2
                        if  (chselectheat== "1"){
                          heat2 <- heat2 %>% config(
                            toImageButtonOptions= list(filename = as.character(input$Nameheat)))
                          
                        }
                        else if (chselectheat == "2"){
                          heat2 <- vals$heat2 %>% config(
                            toImageButtonOptions= list(format ="svg",filename = as.character(input$Nameheat)))
                          
                        }
                        
                      }                                                                                                                    #,margins =c(0,0,50,0)))
                      
                    })
                  })
                  
                })
                
              }#chselectplot== "1"
              
            }#endof sampleplot
            if  (chzscore== "2"){
              if  (chselectplot== "1"){
                datselecheat <- reactiveValues()
                observeEvent( input$FilterAnalysis,{
                  match1 <- datar$match 
                  # print(match1)
                  Gene<-match1$group
                  s1<-datar$s
                  DEGs_top<-subset(match1, select = -c(group,selectedRowID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj, symbol))
                  colnames(DEGs_top) <- s1
                  DEGs_top_scale2 <-  DEGs_top-rowMeans(DEGs_top) 
                  
                  if (as.numeric(length(as.character(input$Group1comp)))==1){
                    names(DEGs_top_scale2)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scale2[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    Big3<-c1
                  }
                  
                  else if (as.numeric(length(as.character(input$Group1comp)))==2){
                    names(DEGs_top_scale2)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scale2[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    
                    names(DEGs_top_scale2)[(datar$lg1):datar$lg2-1]<-as.character(input$Group1comp)[2]
                    gname2<-(as.character(input$Group1comp)[2])
                    c2<-data.frame(gname2 =rowMeans(DEGs_top_scale2[(datar$lg1):datar$lg2-1]))
                    names(c2)[names(c2) == "gname2"] <- as.character(input$Group1comp)[2]
                    Big3<-cbind(c1,c2)
                  }
                  
                  
                  else if(as.numeric(length(as.character(input$Group1comp)))==3){
                    
                    names(DEGs_top_scale2)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scale2[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    
                    names(DEGs_top_scale2)[(datar$lg1):datar$lg2-1]<-as.character(input$Group1comp)[2]
                    gname2<-(as.character(input$Group1comp)[2])
                    c2<-data.frame(gname2 =rowMeans(DEGs_top_scale2[(datar$lg1):datar$lg2-1]))
                    names(c2)[names(c2) == "gname2"] <- as.character(input$Group1comp)[2]
                    
                    names(DEGs_top_scale2)[ (datar$lg1+datar$lg2):datar$lg3-1]<-as.character(input$Group1comp)[3]
                    gname3<-(as.character(input$Group1comp)[3])
                    c3<-data.frame(gname3 =rowMeans(DEGs_top_scale2[(datar$lg1+datar$lg2):datar$lg3-1]))
                    names(c3)[names(c3) == "gname3"] <- as.character(input$Group1comp)[3]
                    Big3<-cbind(c1,c2,c3)
                  }
                  
                  
                  else if (as.numeric(length(as.character(input$Group1comp)))==4){
                    
                    names(DEGs_top_scale2)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scale2[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    
                    names(DEGs_top_scale2)[(datar$lg1):datar$lg2-1]<-as.character(input$Group1comp)[2]
                    gname2<-(as.character(input$Group1comp)[2])
                    c2<-data.frame(gname2 =rowMeans(DEGs_top_scale2[(datar$lg1):datar$lg2-1]))
                    names(c2)[names(c2) == "gname2"] <- as.character(input$Group1comp)[2]
                    
                    names(DEGs_top_scale2)[ (datar$lg1+datar$lg2):datar$lg3-1]<-as.character(input$Group1comp)[3]
                    gname3<-(as.character(input$Group1comp)[3])
                    c3<-data.frame(gname3 =rowMeans(DEGs_top_scale2[(datar$lg1+datar$lg2):datar$lg3-1]))
                    names(c3)[names(c3) == "gname3"] <- as.character(input$Group1comp)[3]
                    
                    names(DEGs_top_scale2)[ (datar$lg1+datar$lg2+datar$lg3):datar$lg4-1]<-as.character(input$Group1comp)[4]
                    gname4<-(as.character(input$Group1comp)[4])
                    c4<-data.frame(gname4 =rowMeans(DEGs_top_scale2[ (datar$lg1+datar$lg2+datar$lg3):datar$lg4-1]))
                    names(c4)[names(c4) == "gname4"] <- as.character(input$Group1comp)[4]
                    Big3<-cbind(c1,c2,c3,c4)
                    
                  }
                  
                  else if (as.numeric(length(as.character(input$Group1comp)))==5){
                    names(DEGs_top_scale2)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scale2[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    
                    names(DEGs_top_scale2)[(datar$lg1):datar$lg2-1]<-as.character(input$Group1comp)[2]
                    gname2<-(as.character(input$Group1comp)[2])
                    c2<-data.frame(gname2 =rowMeans(DEGs_top_scale2[(datar$lg1):datar$lg2-1]))
                    names(c2)[names(c2) == "gname2"] <- as.character(input$Group1comp)[2]
                    
                    names(DEGs_top_scale2)[ (datar$lg1+datar$lg2):datar$lg3-1]<-as.character(input$Group1comp)[3]
                    gname3<-(as.character(input$Group1comp)[3])
                    c3<-data.frame(gname3 =rowMeans(DEGs_top_scale2[(datar$lg1+datar$lg2):datar$lg3-1]))
                    names(c3)[names(c3) == "gname3"] <- as.character(input$Group1comp)[3]
                    
                    names(DEGs_top_scale2)[ (datar$lg1+datar$lg2+datar$lg3):datar$lg4-1]<-as.character(input$Group1comp)[4]
                    gname4<-(as.character(input$Group1comp)[4])
                    c4<-data.frame(gname4 =rowMeans(DEGs_top_scale2[ (datar$lg1+datar$lg2+datar$lg3):datar$lg4-1]))
                    names(c4)[names(c4) == "gname4"] <- as.character(input$Group1comp)[4]
                    
                    names(DEGs_top_scale2)[(datar$lg1+datar$lg2+datar$lg3+datar$lg4):datar$lg5-1]<-as.character(input$Group1comp)[5]
                    gname5<-(as.character(input$Group1comp)[5])
                    c5<-data.frame(gname5 =rowMeans(DEGs_top_scale2[(datar$lg1+datar$lg2+datar$lg3+datar$lg4):datar$lg5-1]))  
                    names(c5)[names(c5) == "gname5"] <- as.character(input$Group1comp)[5]
                    Big3<-cbind(c1,c2,c3,c4,c5)
                    
                  }
                  
                  
                  if (as.numeric(length(as.character(input$Group2comp)))==1){
                    names(DEGs_top_scale2)[names(DEGs_top_scale2) == "B1"] <- as.character(input$Group2comp)[1]
                    # colnames(DEGs_top_scale2)[which(names(DEGs_top_scale2) == "B1")] <-as.character(input$Group2comp)[1]
                    #number<-as.numeric(which( colnames(names(DEGs_top_scale2)== as.character(input$Group2comp)[1])))
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scale2)))
                    names(DEGs_top_scale2)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scale2[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]               
                    Big4<-d1
                    
                  }
                  else if(as.numeric(length(as.character(input$Group2comp)))==2){
                    
                    names(DEGs_top_scale2)[names(DEGs_top_scale2) == "B1"] <- as.character(input$Group2comp)[1]
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scale2)))
                    names(DEGs_top_scale2)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scale2[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]    
                    
                    names(DEGs_top_scale2)[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]<-as.character(input$Group2comp)[2]
                    kname2<-(as.character(input$Group2comp)[2])
                    d2<-data.frame(kname2 =rowMeans(DEGs_top_scale2[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]))
                    names(d2)[names(d2) == "kname2"] <- as.character(input$Group2comp)[2]  
                    Big4<-cbind(d1,d2) 
                    
                  }
                  
                  
                  else if (as.numeric(length(as.character(input$Group2comp)))==3){
                    names(DEGs_top_scale2)[names(DEGs_top_scale2) == "B1"] <- as.character(input$Group2comp)[1]
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scale2)))
                    names(DEGs_top_scale2)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scale2[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]      
                    
                    names(DEGs_top_scale2)[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]<-as.character(input$Group2comp)[2]
                    kname2<-(as.character(input$Group2comp)[2])
                    d2<-data.frame(kname2 =rowMeans(DEGs_top_scale2[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]))
                    names(d2)[names(d2) == "kname2"] <- as.character(input$Group2comp)[2]  
                    
                    
                    names(DEGs_top_scale2)[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]<-as.character(input$Group2comp)[3]
                    kname3<-(as.character(input$Group2comp)[3])
                    d3<-data.frame(kname3 =rowMeans(DEGs_top_scale2[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]))
                    names(d3)[names(d3) == "kname3"] <- as.character(input$Group2comp)[3]  
                    Big4<-cbind(d1,d2,d3) 
                    
                  }
                  
                  else if (as.numeric(length(as.character(input$Group2comp)))==4){
                    names(DEGs_top_scale2)[names(DEGs_top_scale2) == "B1"] <- as.character(input$Group2comp)[1]
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scale2)))
                    names(DEGs_top_scale2)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scale2[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]      
                    
                    names(DEGs_top_scale2)[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]<-as.character(input$Group2comp)[2]
                    kname2<-(as.character(input$Group2comp)[2])
                    d2<-data.frame(kname2 =rowMeans(DEGs_top_scale2[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]))
                    names(d2)[names(d2) == "kname2"] <- as.character(input$Group2comp)[2]  
                    
                    
                    names(DEGs_top_scale2)[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]<-as.character(input$Group2comp)[3]
                    kname3<-(as.character(input$Group2comp)[3])
                    d3<-data.frame(kname3 =rowMeans(DEGs_top_scale2[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]))
                    names(d3)[names(d3) == "kname3"] <- as.character(input$Group2comp)[3]  
                    
                    names(DEGs_top_scale2)[(number+datar$lk1+datar$lk2+datar$lk3):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4-1)]<-as.character(input$Group2comp)[4]
                    kname4<-(as.character(input$Group2comp)[4])
                    d4<-data.frame(kname4 =rowMeans(DEGs_top_scale2[(number+datar$lk1+datar$lk2+datar$lk3):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4-1)]))
                    names(d4)[names(d4) == "kname4"] <- as.character(input$Group2comp)[4]
                    Big4<-cbind(d1,d2,d3,d4) 
                  }
                  
                  else if (as.numeric(length(as.character(input$Group2comp)))==5){
                    names(DEGs_top_scale2)[names(DEGs_top_scale2) == "B1"] <- as.character(input$Group2comp)[1]
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scale2)))
                    names(DEGs_top_scale2)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scale2[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]       
                    
                    names(DEGs_top_scale2)[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]<-as.character(input$Group2comp)[2]
                    kname2<-(as.character(input$Group2comp)[2])
                    d2<-data.frame(kname2 =rowMeans(DEGs_top_scale2[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]))
                    names(d2)[names(d2) == "kname2"] <- as.character(input$Group2comp)[2]  
                    
                    
                    names(DEGs_top_scale2)[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]<-as.character(input$Group2comp)[3]
                    kname3<-(as.character(input$Group2comp)[3])
                    d3<-data.frame(kname3 =rowMeans(DEGs_top_scale2[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]))
                    names(d3)[names(d3) == "kname3"] <- as.character(input$Group2comp)[3]  
                    
                    names(DEGs_top_scale2)[(number+datar$lk1+datar$lk2+datar$lk3):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4-1)]<-as.character(input$Group2comp)[4]
                    kname4<-(as.character(input$Group2comp)[4])
                    d4<-data.frame(kname4 =rowMeans(DEGs_top_scale2[(number+datar$lk1+datar$lk2+datar$lk3):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4-1)]))
                    names(d4)[names(d4) == "kname4"] <- as.character(input$Group2comp)[4]
                    
                    
                    names(DEGs_top_scale2)[(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4+datar$lk5-1)]<-as.character(input$Group2comp)[5]
                    kname5<-(as.character(input$Group2comp)[5])
                    d5<-data.frame(kname5 =rowMeans(DEGs_top_scale2[(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4+datar$lk5-1)]))
                    names(d5)[names(d5) == "kname5"] <- as.character(input$Group2comp)[5]
                    Big4<-cbind(d1,d2,d3,d4,d5)
                  }
                  
                  dd3<-cbind(Big3,Big4)
                  
                  rownames(dd3)  <- rownames(DEGs_top_scale2)
                  
                  z_score <- function(x){
                    (x) / sd(x)
                    # (x - mean(x)) / sd(x)
                  }
                  
                  data_subset_norm <- t(apply(dd3, 1, z_score))
                  data_subset_norm2 <-data.frame( t(apply(dd3, 1, z_score)))
                 
                  
                  DEGs_top_scale3<-cbind(Gene,data_subset_norm2)
                  names(DEGs_top_scale3)[names(DEGs_top_scale3) == "Gene"] <- "Genes"
                  
                  output$Heatmapselect = DT::renderDataTable(DEGs_top_scale3,
                                                             #server = FALSE,
                                                             extensions = 'FixedColumns', 
                                                             options = list(pageLength = 10, autoWidth = FALSE,
                                                                            scrollX = TRUE,
                                                                            scrollCollapse = TRUE,
                                                                            lengthMenu = c(10,50, 100),
                                                                            columnDefs = list(list(targets = c(1:length( DEGs_top_scale3)), 
                                                                                                   searchable = FALSE))
                                                             ))
                  
                  ma <-(1.25)
                  mi <-(-1.25)
                 
                  
                  gradient_col <- ggplot2::scale_fill_gradient2(
                    low = "#1f1fd8",  mid = "#fffcfc",high = "#fb4141",
                    # midpoint = 0.0,
                    limits = c(mi, ma) 
                  )
                  
                  
                  # highlight selected rows in the scatterplot
                  observe({
                    chselectl <- input$optionlselectheat
                    chselectheat <- input$optionselecheat
                    output$dotplot = renderPlotly({
                      if ( chselectl  == "1"){
                        s = input$Heatmapselect_rows_selected
                        vals$s<-s
                        par(mar = c(4, 4, 1, .1))
                        
                        heat<- heatmaply(data_subset_norm[vals$s, , drop = FALSE],scale_fill_gradient_fun = gradient_col,
                                      
                                         showticklabels = T,dendrogram =as.character(input$Group8comp), main =paste(as.character(input$Nameheat)),fontsize_row = 14,
                                         fontsize_col = 14) 
                        heat<-heat %>% layout(height = input$plot_heightselectheat, width = input$plot_widthselectheat) 
                        if  (chselectheat == "1"){
                          heat<- heat %>% config(
                            toImageButtonOptions= list(filename = as.character(input$Nameheat)))
                          #p <- data8$plot 
                        }
                        else if  (chselectheat == "2"){
                          heat3 <- heat %>% config(
                            toImageButtonOptions= list(format = "svg",filename = as.character(input$Nameheat)))
                          #p <- data8$plot 
                        }
                      }
                      else if (chselectl == "2"){
                        heat2<- heatmaply(data_subset_norm[vals$s, , drop = FALSE],scale_fill_gradient_fun = gradient_col,
                                          # scale = "none", col = bluered(100), trace = "none",
                                          # density.info = "none",
                                          showticklabels=c(TRUE, FALSE),dendrogram =as.character(input$Group8comp), main =paste(as.character(input$Nameheat)),fontsize_col = 14,fontsize_row = 14) 
                        heat2<-heat2 %>% layout(height = input$plot_heightselectheat, width = input$plot_widthselectheat)  
                        vals$heat2<-heat2
                        if  (chselectheat== "1"){
                          heat2 <- heat2 %>% config(
                            toImageButtonOptions= list(filename = as.character(input$Nameheat)))
                          
                        }
                        else if (chselectheat == "2"){
                          heat2 <- vals$heat2 %>% config(
                            toImageButtonOptions= list(format ="svg",filename = as.character(input$Nameheat)))
                          
                        }
                        
                      }                                                                                                                    #,margins =c(0,0,50,0)))
                      
                    })
                  })
                  
                })
              }                        
              
              else if (chselectplot== "2"){
                datselecheat <- reactiveValues()
                observeEvent(input$FilterAnalysis,{
                  DEGs_top_scaleg <- datar$DEGs_top_scale5
                  
                  
                  if (as.numeric(length(as.character(input$Group1comp)))==1){
                    names(DEGs_top_scaleg)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scaleg[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    Big3<-c1
                  }
                  
                  else if (as.numeric(length(as.character(input$Group1comp)))==2){
                    names(DEGs_top_scaleg)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scaleg[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    
                    names(DEGs_top_scaleg)[(datar$lg1):datar$lg2-1]<-as.character(input$Group1comp)[2]
                    gname2<-(as.character(input$Group1comp)[2])
                    c2<-data.frame(gname2 =rowMeans(DEGs_top_scaleg[(datar$lg1):datar$lg2-1]))
                    names(c2)[names(c2) == "gname2"] <- as.character(input$Group1comp)[2]
                    Big3<-cbind(c1,c2)
                  }
                  
                  
                  else if(as.numeric(length(as.character(input$Group1comp)))==3){
                    
                    names(DEGs_top_scaleg)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scaleg[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    
                    names(DEGs_top_scaleg)[(datar$lg1):datar$lg2-1]<-as.character(input$Group1comp)[2]
                    gname2<-(as.character(input$Group1comp)[2])
                    c2<-data.frame(gname2 =rowMeans(DEGs_top_scaleg[(datar$lg1):datar$lg2-1]))
                    names(c2)[names(c2) == "gname2"] <- as.character(input$Group1comp)[2]
                    
                    names(DEGs_top_scaleg)[ (datar$lg1+datar$lg2):datar$lg3-1]<-as.character(input$Group1comp)[3]
                    gname3<-(as.character(input$Group1comp)[3])
                    c3<-data.frame(gname3 =rowMeans(DEGs_top_scaleg[(datar$lg1+datar$lg2):datar$lg3-1]))
                    names(c3)[names(c3) == "gname3"] <- as.character(input$Group1comp)[3]
                    Big3<-cbind(c1,c2,c3)
                  }
                  
                  
                  else if (as.numeric(length(as.character(input$Group1comp)))==4){
                    
                    names(DEGs_top_scaleg)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scaleg[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    
                    names(DEGs_top_scaleg)[(datar$lg1):datar$lg2-1]<-as.character(input$Group1comp)[2]
                    gname2<-(as.character(input$Group1comp)[2])
                    c2<-data.frame(gname2 =rowMeans(DEGs_top_scaleg[(datar$lg1):datar$lg2-1]))
                    names(c2)[names(c2) == "gname2"] <- as.character(input$Group1comp)[2]
                    
                    names(DEGs_top_scaleg)[ (datar$lg1+datar$lg2):datar$lg3-1]<-as.character(input$Group1comp)[3]
                    gname3<-(as.character(input$Group1comp)[3])
                    c3<-data.frame(gname3 =rowMeans(DEGs_top_scaleg[(datar$lg1+datar$lg2):datar$lg3-1]))
                    names(c3)[names(c3) == "gname3"] <- as.character(input$Group1comp)[3]
                    
                    names(DEGs_top_scaleg)[ (datar$lg1+datar$lg2+datar$lg3):datar$lg4-1]<-as.character(input$Group1comp)[4]
                    gname4<-(as.character(input$Group1comp)[4])
                    c4<-data.frame(gname4 =rowMeans(DEGs_top_scaleg[ (datar$lg1+datar$lg2+datar$lg3):datar$lg4-1]))
                    names(c4)[names(c4) == "gname4"] <- as.character(input$Group1comp)[4]
                    Big3<-cbind(c1,c2,c3,c4)
                    
                  }
                  
                  else if (as.numeric(length(as.character(input$Group1comp)))==5){
                    names(DEGs_top_scaleg)[1:datar$lg1-1]<-as.character(input$Group1comp)[1]
                    gname<-(as.character(input$Group1comp)[1])
                    c1<-data.frame( gname=rowMeans(DEGs_top_scaleg[1:datar$lg1-1]))
                    names(c1)[names(c1) == "gname"] <- as.character(input$Group1comp)[1]
                    
                    names(DEGs_top_scaleg)[(datar$lg1):datar$lg2-1]<-as.character(input$Group1comp)[2]
                    gname2<-(as.character(input$Group1comp)[2])
                    c2<-data.frame(gname2 =rowMeans(DEGs_top_scaleg[(datar$lg1):datar$lg2-1]))
                    names(c2)[names(c2) == "gname2"] <- as.character(input$Group1comp)[2]
                    
                    names(DEGs_top_scaleg)[ (datar$lg1+datar$lg2):datar$lg3-1]<-as.character(input$Group1comp)[3]
                    gname3<-(as.character(input$Group1comp)[3])
                    c3<-data.frame(gname3 =rowMeans(DEGs_top_scaleg[(datar$lg1+datar$lg2):datar$lg3-1]))
                    names(c3)[names(c3) == "gname3"] <- as.character(input$Group1comp)[3]
                    
                    names(DEGs_top_scaleg)[ (datar$lg1+datar$lg2+datar$lg3):datar$lg4-1]<-as.character(input$Group1comp)[4]
                    gname4<-(as.character(input$Group1comp)[4])
                    c4<-data.frame(gname4 =rowMeans(DEGs_top_scaleg[ (datar$lg1+datar$lg2+datar$lg3):datar$lg4-1]))
                    names(c4)[names(c4) == "gname4"] <- as.character(input$Group1comp)[4]
                    
                    names(DEGs_top_scaleg)[(datar$lg1+datar$lg2+datar$lg3+datar$lg4):datar$lg5-1]<-as.character(input$Group1comp)[5]
                    gname5<-(as.character(input$Group1comp)[5])
                    c5<-data.frame(gname5 =rowMeans(DEGs_top_scaleg[(datar$lg1+datar$lg2+datar$lg3+datar$lg4):datar$lg5-1]))  
                    names(c5)[names(c5) == "gname5"] <- as.character(input$Group1comp)[5]
                    Big3<-cbind(c1,c2,c3,c4,c5)
                    
                  }
                  
                  
                  if (as.numeric(length(as.character(input$Group2comp)))==1){
                    names(DEGs_top_scaleg)[names(DEGs_top_scaleg) == "B1"] <- as.character(input$Group2comp)[1]
                    # colnames(DEGs_top_scale2)[which(names(DEGs_top_scale2) == "B1")] <-as.character(input$Group2comp)[1]
                    #number<-as.numeric(which( colnames(names(DEGs_top_scale2)== as.character(input$Group2comp)[1])))
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scaleg)))
                    names(DEGs_top_scaleg)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scaleg[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]               
                    Big4<-d1
                    
                  }
                  else if(as.numeric(length(as.character(input$Group2comp)))==2){
                    
                    names(DEGs_top_scaleg)[names(DEGs_top_scaleg) == "B1"] <- as.character(input$Group2comp)[1]
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scaleg)))
                    names(DEGs_top_scaleg)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scaleg[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]     
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]<-as.character(input$Group2comp)[2]
                    kname2<-(as.character(input$Group2comp)[2])
                    d2<-data.frame(kname2 =rowMeans(DEGs_top_scaleg[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]))
                    names(d2)[names(d2) == "kname2"] <- as.character(input$Group2comp)[2]  
                    Big4<-cbind(d1,d2) 
                    
                  }
                  
                  
                  else if (as.numeric(length(as.character(input$Group2comp)))==3){
                    names(DEGs_top_scaleg)[names(DEGs_top_scaleg) == "B1"] <- as.character(input$Group2comp)[1]
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scaleg)))
                    names(DEGs_top_scaleg)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scaleg[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]     
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]<-as.character(input$Group2comp)[2]
                    kname2<-(as.character(input$Group2comp)[2])
                    d2<-data.frame(kname2 =rowMeans(DEGs_top_scaleg[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]))
                    names(d2)[names(d2) == "kname2"] <- as.character(input$Group2comp)[2]  
                    
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]<-as.character(input$Group2comp)[3]
                    kname3<-(as.character(input$Group2comp)[3])
                    d3<-data.frame(kname3 =rowMeans(DEGs_top_scaleg[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]))
                    names(d3)[names(d3) == "kname3"] <- as.character(input$Group2comp)[3]  
                    Big4<-cbind(d1,d2,d3) 
                    
                  }
                  
                  else if (as.numeric(length(as.character(input$Group2comp)))==4){
                    names(DEGs_top_scaleg)[names(DEGs_top_scaleg) == "B1"] <- as.character(input$Group2comp)[1]
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scaleg)))
                    names(DEGs_top_scaleg)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scaleg[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]     
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]<-as.character(input$Group2comp)[2]
                    kname2<-(as.character(input$Group2comp)[2])
                    d2<-data.frame(kname2 =rowMeans(DEGs_top_scaleg[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]))
                    names(d2)[names(d2) == "kname2"] <- as.character(input$Group2comp)[2]  
                    
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]<-as.character(input$Group2comp)[3]
                    kname3<-(as.character(input$Group2comp)[3])
                    d3<-data.frame(kname3 =rowMeans(DEGs_top_scaleg[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]))
                    names(d3)[names(d3) == "kname3"] <- as.character(input$Group2comp)[3]  
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1+datar$lk2+datar$lk3):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4-1)]<-as.character(input$Group2comp)[4]
                    kname4<-(as.character(input$Group2comp)[4])
                    d4<-data.frame(kname4 =rowMeans(DEGs_top_scaleg[(number+datar$lk1+datar$lk2+datar$lk3):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4-1)]))
                    names(d4)[names(d4) == "kname4"] <- as.character(input$Group2comp)[4]
                    Big4<-cbind(d1,d2,d3,d4) 
                  }
                  
                  else if (as.numeric(length(as.character(input$Group2comp)))==5){
                    names(DEGs_top_scaleg)[names(DEGs_top_scaleg) == "B1"] <- as.character(input$Group2comp)[1]
                    number<-as.numeric(match(as.character(input$Group2comp)[1],names(DEGs_top_scaleg)))
                    names(DEGs_top_scaleg)[number:(number+datar$lk1-1)]<-as.character(input$Group2comp)[1]
                    kname<-(as.character(input$Group2comp)[1])
                    d1<-data.frame(kname =rowMeans(DEGs_top_scaleg[number:(number+datar$lk1-1)]))
                    names(d1)[names(d1) == "kname"] <- as.character(input$Group2comp)[1]     
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]<-as.character(input$Group2comp)[2]
                    kname2<-(as.character(input$Group2comp)[2])
                    d2<-data.frame(kname2 =rowMeans(DEGs_top_scaleg[(number+datar$lk1):(number+datar$lk1+datar$lk2-1)]))
                    names(d2)[names(d2) == "kname2"] <- as.character(input$Group2comp)[2]  
                    
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]<-as.character(input$Group2comp)[3]
                    kname3<-(as.character(input$Group2comp)[3])
                    d3<-data.frame(kname3 =rowMeans(DEGs_top_scaleg[(number+datar$lk1+datar$lk2):(number+datar$lk1+datar$lk2+datar$lk3-1)]))
                    names(d3)[names(d3) == "kname3"] <- as.character(input$Group2comp)[3]  
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1+datar$lk2+datar$lk3):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4-1)]<-as.character(input$Group2comp)[4]
                    kname4<-(as.character(input$Group2comp)[4])
                    d4<-data.frame(kname4 =rowMeans(DEGs_top_scaleg[(number+datar$lk1+datar$lk2+datar$lk3):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4-1)]))
                    names(d4)[names(d4) == "kname4"] <- as.character(input$Group2comp)[4]
                    
                    
                    names(DEGs_top_scaleg)[(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4+datar$lk5-1)]<-as.character(input$Group2comp)[5]
                    kname5<-(as.character(input$Group2comp)[5])
                    d5<-data.frame(kname5 =rowMeans(DEGs_top_scaleg[(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4):(number+datar$lk1+datar$lk2+datar$lk3+datar$lk4+datar$lk5-1)]))
                    names(d5)[names(d5) == "kname5"] <- as.character(input$Group2comp)[5]
                    Big4<-cbind(d1,d2,d3,d4,d5)
                  }
                  
                  dd3<-cbind(Big3,Big4)
                  
                  rownames(dd3)  <- rownames(DEGs_top_scaleg)
                  
                  z_score <- function(x){
                    (x) / sd(x)
                    # (x - mean(x)) / sd(x)
                  }
                  
                  data_subset_norm <- t(apply(dd3, 1, z_score))
                  data_subset_norm2 <-data.frame( t(apply(dd3, 1, z_score)))
                  # ma <-max( data_subset_norm , na.rm=TRUE)
                  # mi <-min( data_subset_norm , na.rm=TRUE)
                  
                  DEGs_top_scale3<-cbind(datar$Gene,data_subset_norm2)
                  names(DEGs_top_scale3)[names(DEGs_top_scale3) == "datar$Gene"] <- "Genes"
                  
                  
                  output$Heatmapselect = DT::renderDataTable(DEGs_top_scale3,
                                                             #server = FALSE,
                                                             extensions = 'FixedColumns', 
                                                             options = list(pageLength = 10, autoWidth = FALSE,
                                                                            scrollX = TRUE,
                                                                            scrollCollapse = TRUE,
                                                                            lengthMenu = c(10,50, 100),
                                                                            columnDefs = list(list(targets = c(1:length( DEGs_top_scale3)), 
                                                                                                   searchable = FALSE))
                                                             ))
                  
                  ma <-(1.25)
                  mi <-(-1.25)
               
                  
                  gradient_col <- ggplot2::scale_fill_gradient2(
                    low = "#1f1fd8",  mid = "#fffcfc",high = "#fb4141",
                    # midpoint = 0.0,
                    limits = c(mi, ma) 
                  )
                  
                  
                  # highlight selected rows in the scatterplot
                  observe({
                    chselectl <- input$optionlselectheat
                    chselectheat <- input$optionselecheat
                    output$dotplot = renderPlotly({
                      if ( chselectl  == "1"){
                        s = input$Heatmapselect_rows_selected
                        vals$s<-s
                        par(mar = c(4, 4, 1, .1))
                        
                        heat<- heatmaply(data_subset_norm[vals$s, , drop = FALSE],scale_fill_gradient_fun = gradient_col,
                                       
                                         showticklabels = T,dendrogram =as.character(input$Group8comp), main =paste(as.character(input$Nameheat)),fontsize_col = 14,fontsize_row = 14) 
                        heat<-heat %>% layout(height = input$plot_heightselectheat, width = input$plot_widthselectheat) 
                        if  (chselectheat == "1"){
                          heat<- heat %>% config(
                            toImageButtonOptions= list(filename = as.character(input$Nameheat)))
                          #p <- data8$plot 
                        }
                        else if  (chselectheat == "2"){
                          heat3 <- heat %>% config(
                            toImageButtonOptions= list(format = "svg",filename = as.character(input$Nameheat)))
                          #p <- data8$plot 
                        }
                      }
                      else if (chselectl == "2"){
                        heat2<- heatmaply(data_subset_norm[vals$s, , drop = FALSE],scale_fill_gradient_fun = gradient_col,
                                        
                                          showticklabels=c(TRUE, FALSE),dendrogram =as.character(input$Group8comp), main =paste(as.character(input$Nameheat)),fontsize_col = 14,fontsize_row = 14) 
                        heat2<-heat2 %>% layout(height = input$plot_heightselectheat, width = input$plot_widthselectheat)  
                        vals$heat2<-heat2
                        if  (chselectheat== "1"){
                          heat2 <- heat2 %>% config(
                            toImageButtonOptions= list(filename = as.character(input$Nameheat)))
                          
                        }
                        else if (chselectheat == "2"){
                          heat2 <- vals$heat2 %>% config(
                            toImageButtonOptions= list(format ="svg",filename = as.character(input$Nameheat)))
                          
                        }
                        
                      }                                                                                                                    #,margins =c(0,0,50,0)))
                      
                    })
                  })
                  
                })
                
              }#chselectplot== "1"
              
              
            }#endof group and zscore
            
            
          })# End Observe
################################################################################################################
          #Gene Level Expression 
          dat <- reactiveValues()
          observeEvent(input$FilterAnalysis,{
            withProgress(message = "Preparing list of Gene ID please wait",
                         {
                           
                           observe({
                             chgene <- input$optiongene
                             
                             if (chgene == "1"){             
                               
                               bb<-datar$Collapsed_data  
                               m<- bb[,-c(1,2)]
                               col_namem1<- paste(datar$n1, sep="")
                               col_namem2<-paste(datar$n2, sep="")
                               
                               names(m)[1:length(datar$n1)]<-col_namem1
                               names(m)[length(datar$n1)+1:length(datar$n2)]<-col_namem2
                               
                               m<-cbind(bb[c(1,2)],m)
                               names(m)[names(m) == "group"] <- "Genes"
                               names(m)[names(m) == "selectedRowID"] <- "Ensembl_ID"
                               melted <- melt(m)
                               # # bID<-datar$K
                               colnames(melted) <- c("gene", "ID","Samples", "normalized_counts")                 
                               
                               output$ID = renderDT(m,
                                                    #server = FALSE,
                                                    
                                                    selection = 'single',
                                                    extensions = 'FixedColumns', 
                                                    options = list(pageLength = 10, autoWidth = FALSE,
                                                                   scrollX = TRUE,
                                                                   scrollCollapse = TRUE,
                                                                   lengthMenu = c(10,50, 100),
                                                                   columnDefs = list(list(targets = c(2:length(m)), 
                                                                                          searchable = FALSE))
                                                    ))
                               observe({                  
                                 req(input$ID_rows_selected)
                                 selRow <- m[input$ID_rows_selected,]
                                 # print(selRow[[1]])
                                 # print(as.character(selRow[[1]]))
                                 mo <- subset(melted, gene == as.character(selRow[[1]]))
                                 mm <- join(mo, datar$ds, by = "Samples", type = "inner")
                                 
                                 widthgene <- reactive ({ input$plot_widthgene })
                                 heightgene <- reactive ({ input$plot_heightgene}) 
                                 output$genePlot <- renderPlot(width =  widthgene, height = heightgene,{  
                                   boxplot<-ggplot(mm, aes(x=Types,y=normalized_counts))+
                                     geom_boxplot(width = 0.5) + 
                                     geom_dotplot(
                                       #aes(color=Types),
                                       aes(fill = Types),
                                       binaxis = "y", stackdir = "center",dotsize=0.5
                                     ) +geom_text_repel(aes(label = Samples),position=position_jitter(w=0.15,h=0)) +
                                     ggtitle(selRow[[1]]) +
                                     theme_bw(base_size = 14)+
                                     theme(plot.title=element_text(hjust=0.5))
                                   #,shape=18
                                   boxplot<-  boxplot + theme(axis.title.x = element_blank(), axis.text = element_text(size = 12))
                                   boxplot<-boxplot + stat_summary(fun=mean,  aes(shape="mean"), geom="point", colour = "red", size=4)+scale_shape_manual("", values=c("mean"="x"))
                                   vals$boxplot <- boxplot
                                   print(boxplot)
                                 
                                   
                                 })
                                 
                                 
                                 choptiongene<- input$optionchgene
                                 if (choptiongene == "1"){
                                   output$PlotDownloadgene <- downloadHandler(
                                     filename = function(){paste("CountPlot",'.png',sep='')},
                                    
                                     content = function(file){
                                       device <- function(..., width, height) grDevices::png(..., width =input$plot_widthgene*4, height =input$plot_heightgene*4, res = 300, units = "px")
                                      
                                       ggsave(file, plot = (vals$boxplot)  , device = device)
                                     })
                                 }
                                 else if (choptiongene == "2"){
                                   output$PlotDownloadgene <- downloadHandler(
                                     filename = function(){paste("CountPlot",'.svg',sep='')},
                                     content = function(file) {
                                       device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthgene*4/300), height =round(input$plot_heightgene*4/300))
                                       ggsave(file, plot = (vals$boxplot)  , device = device)
                                      
                                       
                                     })
                                 }  
                                 
                               })#end of Observe 
                               
                             } #end of Count Gene  
                             
                             
                             
                           }) # end of Observe
                         }) # end pf Progresss bar                
          }) # End of Gene Level
          #################################################################################################################    
          ################################################################################################################
          #Volcano Plot
          observeEvent(input$FilterAnalysis,{
         
            output$volc = DT::renderDataTable(na.omit(datar$Dm),
                                             
                                              rownames = FALSE,extensions = 'FixedColumns',
                                              options = list(pageLength = 10, autoWidth = FALSE,
                                                             scrollX = TRUE,
                                                             scrollCollapse = TRUE,
                                                             lengthMenu = c(10,30,60),
                                                             columnDefs = list(list(targets = c(1:7), 
                                                                                    searchable = FALSE))))
            #editable = FALSE
          })                      
          observe({     
            chselectvolc <- input$volc   
            if (chselectvolc   == "1"){     
              datavolc<-reactiveValues()  
              observeEvent(input$FilterAnalysis,{
                observe({
                  chvolc <- input$optionvolc
                  pointvolc<- as.numeric(input$point_volc)
                  hivolc <-input$plot_heightvolc
                  wvolc <-input$plot_widthvolc
            
                  
                  
                  fc1<- -5
                  fc2<- 5
                  p_cutoff<- -log10(0.01)
                  
                  pcut <-(0.01)
                  de <- na.omit(datar$Dm)
                  
                  
                  tooltip = ifelse(is.na(de$symbol), de$pvalue, paste(de$symbol, sep = ""))
                  
                  
                  de$Group <- "Unchanged"
                  
                 
                  de$Group[(de$log2FoldChange) > (fc2) & (de$pvalue) < pcut] <- "Increased"
                  #minus_log10_padj
                  de$Group[(de$log2FoldChange) < (fc1)& (de$pvalue) < pcut] <- "Decreased"
                 
                  
                  output$volcanoPlot <- renderPlotly({
                    withProgress(message = "Drawing Volacno plot, please wait",
                                 {
                                   plot <-  de %>%
                                     ggplot(aes(x = log2FoldChange,
                                                y = minus_log10_pvalue,
                                                text = tooltip,
                                                key = row.names( de))) +
                                     geom_point(aes(col = Group),size=pointvolc) +
                                     xlab("log2FoldChange") +
                                     ylab("-log10(P-value)")+labs(title= as.character(input$volcname))+theme(plot.title = element_text(hjust = 0.5)) 
                                   
                                   p2 <- plot + geom_vline(xintercept=c( fc1,fc2),linetype="dashed", col="black", size=0.1)+ geom_hline(yintercept=p_cutoff,linetype="dashed", col="black", size=0.1) 
                                 
                                   
                                   mycolors <- c("blue", "red", "gray")
                                   names(mycolors) <- c("Decreased", "Increased", "Unchanged")
                                   p3 <- p2 + scale_colour_manual(values = mycolors)
                                   
                                   if ( chvolc  == "1"){
                                     #p3 %>%
                                     ggplotly(p3, width = wvolc, height = hivolc,tooltip = "tooltip") %>%
                                       layout(dragmode = "select")%>%
                                       config(
                                         toImageButtonOptions= list(filename = as.character(input$volcname),res =350, units = "px")) 
                                     
                                     
                                   }
                                   
                                   else if (chvolc == "2"){
                                     p3 %>%
                                       ggplotly(tooltip = "tooltip") %>%
                                       layout(dragmode = "select",autosize=TRUE)%>%
                                       config(
                                         toImageButtonOptions= list(format ="svg",filename = as.character(input$volcname))) 
                                   }
                                 })             
                  })
                  
                })             
                
              })
              
              volcupdate<-reactiveValues()
              observeEvent(input$volcAnalysis,{
                # mylist<- NULL
                observe({
                  chvolc <- input$optionvolc
                  pointvolc<- as.numeric(input$point_volc)
                  hivolc <-input$plot_heightvolc
                  wvolc <-input$plot_widthvolc
                  
                  # pcut<- as.numeric(input$point_adj)
                  fc1<- as.numeric(input$fc_cutoff)
                  fc2<- as.numeric(input$fc_cutoff1)
                  p_cutoff<- -log10(as.numeric(input$ p_cutoff))
                  
                  pcut <-(as.numeric(input$ p_cutoff))
                  de <- na.omit(datar$Dm)
                  
                  
                  tooltip = ifelse(is.na(de$symbol), de$pvalue, paste(de$symbol, sep = ""))
                  
                  
                  de$Group <- "Unchanged"
                  
                
                  
                  de$Group[(de$log2FoldChange) > (fc2) & (de$pvalue) < pcut] <- "Increased"
                  
                  de$Group[(de$log2FoldChange) < (fc1)& (de$pvalue) < pcut] <- "Decreased"
                 
                  
                  output$volcanoPlot <- renderPlotly({
                    withProgress(message = "Drawing Volacno plot, please wait",
                                 {
                                   plot <-  de %>%
                                     ggplot(aes(x = log2FoldChange,
                                                y = minus_log10_pvalue,
                                                # y = minus_log10_padj,
                                                #col=Group,
                                                text = tooltip,
                                                key = row.names( de))) +
                                     geom_point(aes(col = Group),size=pointvolc) +
                                     xlab("log2FoldChange") +
                                     ylab("-log10(P-value)")+labs(title= as.character(input$volcname))+theme(plot.title = element_text(hjust = 0.5)) 
                                   
                                   p2 <- plot + geom_vline(xintercept=c( fc1,fc2),linetype="dashed", col="black", size=0.1)+ geom_hline(yintercept=p_cutoff,linetype="dashed", col="black", size=0.1) 
                                   # scale_color_continuous(guide ='none')
                                   
                                   mycolors <- c("blue", "red", "gray")
                                   names(mycolors) <- c("Decreased", "Increased", "Unchanged")
                                   p3 <- p2 + scale_colour_manual(values = mycolors)
                                   
                                   if ( chvolc  == "1"){
                                     #p3 %>%
                                     ggplotly(p3, width = wvolc, height = hivolc,tooltip = "tooltip") %>%
                                       layout(dragmode = "select")%>%
                                       config(
                                         toImageButtonOptions= list(filename = as.character(input$volcname),res =350, units = "px")) 
                                  
                                   }
                                   
                                   else if (chvolc == "2"){
                                     p3 %>%
                                       ggplotly(tooltip = "tooltip") %>%
                                       layout(dragmode = "select",autosize=TRUE)%>%
                                       config(
                                         toImageButtonOptions= list(format ="svg",filename = as.character(input$volcname))) 
                                   }
                                   
                                 })
                  })
                  
                })
              })  
              
            }
            
            
            else if (chselectvolc   == "2"){
              datavolcp1<-reactiveValues()  
              observeEvent((input$FilterAnalysis),{
                observe({
                  # chvolc <- input$optionvolc
                  pointvolc<- as.numeric(input$point_volc)
                  hivolc <-input$plot_heightvolc
                  wvolc <-input$plot_widthvolc
                  # pcut<- as.numeric(input$point_adj)
                  
                  
                  fc1<- -5
                  fc2<- 5
                  p_cutoff<- -log10(0.01)
                  
                  pcut <-(0.01)
                  de <- na.omit(datar$Dm)
                  
                  
                  tooltip = ifelse(is.na(de$symbol), de$pvalue, paste(de$symbol, sep = ""))
                  
                  
                  de$Group <- "Unchanged"
                  
                  # print(fc1)
                  #print( p_cutoff)
                  
                  de$Group[(de$log2FoldChange) > (fc2) & (de$pvalue) < pcut] <- "Increased"
                  #minus_log10_padj
                  de$Group[(de$log2FoldChange) < (fc1)& (de$pvalue) < pcut] <- "Decreased"
                  
                  output$volcanoPlot2 <- renderPlot(width=wvolc,height=hivolc,{
                    withProgress(message = "Drawing Volacno plot, please wait",
                                 {
                                   mylist<- (input$user_gene_list)
                               
                                   
                                   de$genelabels <- ""
                                   
                                   de$genelabels <- ifelse(de$symbol == mylist[1]
                                                           |de$symbo == mylist[2]
                                                           |de$symbo == mylist[3]
                                                           |de$symbo == mylist[4]
                                                           |de$symbo == mylist[5]
                                                           |de$symbo == mylist[6]
                                                           |de$symbo == mylist[7]
                                                           |de$symbo == mylist[8]
                                                           |de$symbo == mylist[9]
                                                           |de$symbo == mylist[10], T, F)
                                   
                                   plot <- ggplot(de) + 
                                     geom_point(aes(x = log2FoldChange,y = minus_log10_pvalue,col = Group),size=pointvolc) +
                                     geom_text_repel(aes(log2FoldChange, minus_log10_pvalue),size =6 ,label = ifelse(de$genelabels == TRUE, as.character(de$symbol),""), box.padding = unit(0.9, "lines"),hjust= 0.3) + 
                                     xlab("log2FoldChange") +
                                     ylab("-log10(P-value)")+
                                     labs(title= as.character(input$volcname))+theme(plot.title = element_text(hjust = 0.5))+theme_bw(base_size = 16)
                                   
                                   p2 <- plot + geom_vline(xintercept=c( fc1,fc2),linetype="dashed", col="black", size=0.6)+ geom_hline(yintercept=p_cutoff,linetype="dashed", col="black", size=0.6) 
                                  
                                   
                                   mycolors <- c("blue", "red", "gray")
                                   names(mycolors) <- c("Decreased", "Increased", "Unchanged")
                                   p3 <- p2 + scale_colour_manual(values = mycolors)
                                   vals$p3 <-p3
                                   print(p3)
                                   
                                   
                                 })
                    
                  })
                  
                  ########### Donwload volcplot ###############################################################################   
                  observe({
                    chvolc <- input$optionvolc
                    if (chvolc== "1"){
                      output$volcdownload<- downloadHandler(
                        filename = function(){paste("VolcanoPlot",'.png',sep='')},
                        content = function(file){
                          device <- function(..., width, height) grDevices::png(..., width =input$plot_widthvolc*4, height =input$plot_heightvolc*4, res = 350, units = "px")
                          ggsave(file, plot = (vals$p3)  , device = device)
                        })
                    }
                    else if (chvolc == "2"){
                      output$volcdownload <- downloadHandler(
                        filename = function(){paste("VolcanoPlot",'.svg',sep='')},
                        content = function(file) {
                          device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthvolc*4/300), height =round(input$plot_heightvolc*4/300))
                          ggsave(file, plot = (vals$p3)  , device = device)
                          
                          
                        })
                    }
                  })# End plot
                })
              })
              
              datavolcp2<-reactiveValues()  
              observeEvent(input$volcAnalysis,{
                observe({
                  # chvolc <- input$optionvolc
                  pointvolc<- as.numeric(input$point_volc)
                  hivolc <-input$plot_heightvolc
                  wvolc <-input$plot_widthvolc
                  fc1<- as.numeric(input$fc_cutoff)
                  fc2<- as.numeric(input$fc_cutoff1)
                  p_cutoff<- -log10(as.numeric(input$ p_cutoff))
                  
                  pcut <-(as.numeric(input$ p_cutoff))
                  de <- na.omit(datar$Dm)
                  
                  
                  tooltip = ifelse(is.na(de$symbol), de$pvalue, paste(de$symbol, sep = ""))
                  
                  
                  de$Group <- "Unchanged"
                  
                
                  
                  de$Group[(de$log2FoldChange) > (fc2) & (de$pvalue) < pcut] <- "Increased"
                  #minus_log10_padj
                  de$Group[(de$log2FoldChange) < (fc1)& (de$pvalue) < pcut] <- "Decreased"
                  
                  
                  
                  output$volcanoPlot2 <- renderPlot(width=wvolc,height=hivolc,{
                    withProgress(message = "Drawing Volacno plot, please wait",
                                 {
                                   mylist<- (input$user_gene_list)
                                   
                                   de$genelabels <- ""
                                  
                                   de$genelabels <- ifelse(de$symbol == mylist[1]
                                                           |de$symbo == mylist[2]
                                                           |de$symbo == mylist[3]
                                                           |de$symbo == mylist[4]
                                                           |de$symbo == mylist[5]
                                                           |de$symbo == mylist[6]
                                                           |de$symbo == mylist[7]
                                                           |de$symbo == mylist[8]
                                                           |de$symbo == mylist[9]
                                                           |de$symbo == mylist[10], T, F)
                                   
                                   plot <- ggplot(de) + 
                                     geom_point(aes(x = log2FoldChange,y = minus_log10_pvalue,col = Group),size=pointvolc) +
                                     geom_text_repel(aes(log2FoldChange, minus_log10_pvalue),label = ifelse(de$genelabels == TRUE, as.character(de$symbol),""), box.padding = unit(0.7, "lines"),hjust= 0.3) + 
                                     theme(legend.title=element_blank(),text = element_text(size= 14))+  xlab("log2FoldChange") +
                                     ylab("-log10(P-value)")+labs(title= as.character(input$volcname))+theme(plot.title = element_text(hjust = 0.5)) +theme_bw(base_size = 14)
                                   
                                   p2 <- plot + geom_vline(xintercept=c( fc1,fc2),linetype="dashed", col="black", size=0.3)+ geom_hline(yintercept=p_cutoff,linetype="dashed", col="black", size=0.1) 
                                   
                                   
                                   mycolors <- c("blue", "red", "gray")
                                   names(mycolors) <- c("Decreased", "Increased", "Unchanged")
                                   p3 <- p2 + scale_colour_manual(values = mycolors)
                                   vals$p3 <-p3
                                   print(p3)
                                   
                                   
                                 })
                    
                  })
########### Donwload volcplot ###############################################################################   
                  observe({
                    chvolc <- input$optionvolc
                    if (chvolc== "1"){
                      output$volcdownload<- downloadHandler(
                        filename = function(){paste("VolcanoPlot",'.png',sep='')},
                        content = function(file){
                          device <- function(..., width, height) grDevices::png(..., width =input$plot_widthvolc*4, height =input$plot_heightvolc*4, res = 350, units = "px")
                          ggsave(file, plot = (vals$p3)  , device = device)
                        })
                    }
                    else if (chvolc == "2"){
                      output$volcdownload <- downloadHandler(
                        filename = function(){paste("VolcanoPlot",'.svg',sep='')},
                        content = function(file) {
                          device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthvolc*4/300), height =round(input$plot_heightvolc*4/300))
                          ggsave(file, plot = (vals$p3)  , device = device)
                        
                        })
                    }
                  })# End plot
                  
                })
              })  
            }
            
          })      
          
          
          
          
          output$TableDownloadMSG <- downloadHandler(
            filename = function(){paste("MSGRes",'.txt')},
            content = function(file){
              write.csv(na.omit(datar$Dm), file, row.names = TRUE)
            }) 
          
          
          
          #   
          
          #######ENd Volcano
          ################################################################################################################
          ###################################################################################################################
          ## Analysis tab panel content for mMCP-counter
          datamcp <- reactiveValues()
          observeEvent(input$bb,{
            Collapsed_data <- datar$Coll
            head(Collapsed_data)
            chconter<- input$conter       
            if (chconter == "Mouse"){     
              testdata<- Collapsed_data[,-c(1,2)]
              
              col_name10<- paste(datar$n1, sep="")
              col_name11<-paste(datar$n2, sep="")
              
              names(testdata)[1:length(datar$n1)]<-col_name10
              names(testdata)[length(datar$n1)+1:length(datar$n2)]<-col_name11
              
              #######################################################################################################################      
              cit.dfAggregate <- function (data, partition, MARGIN = 1, fAggreg = median.na) {
                cMARGIN <- setdiff(c(1, 2), MARGIN)
                n <- length(partition)
                N <- dim(data)[MARGIN]
                p <- dim(data)[cMARGIN]
                if (n != N)
                  stop("ERROR - cit.dfAggregate : size of partition doesn't correspond to data dimension")
                l <- split(1:N, partition)
                d <- data
                if (MARGIN == 2)d <- t(data)
                
                d <- matrix(sapply(l, function(i) if (length(i) == 1) {
                  unlist(d[i, ])
                }else {
                  apply(d[i, ], 2, fAggreg)
                }), ncol = p, byrow = TRUE)
                
                
                
                d <- as.data.frame(d)
                rownames(d) <- names(l)
                names(d) <- dimnames(data)[[cMARGIN]]
                if (MARGIN == 2)
                  d <- as.data.frame(t(d))
                d
              }
              
              median.na <- function (x) {
                return(median(x,na.rm=TRUE))
              }
              
              #######      
              mMCPcounter.estimate <- function(exp, features = c("Gene.Symbol","ENSEMBL.ID","Probes")[1]){
                data("mMCPcounter_signatures",  envir=sys.frame(sys.nframe()))
                foundGenes <- intersect(mMCPcounter_signatures[,features],rownames(exp))
                if(length(foundGenes)==0){stop("No signature found in input row names. Please ensure the features argument is accurately set.")}
                absentSignatures <- setdiff(unique(mMCPcounter_signatures$Denomination),unique(mMCPcounter_signatures[mMCPcounter_signatures[,features]%in%rownames(exp),"Denomination"]))
                values <- reactiveValues()
                if(length(absentSignatures)==0)
                  
                  queryMagic <- function() {
                    print(paste("All genes were found for population(s)."))
                    
                    return("Data")
                  }
                
                if(length(absentSignatures)>0)
                  # {warning(paste("No genes were found for population(s): ",paste(absentSignatures,collapse = ", "),".",sep=""))}
                  
                  queryMagic <- function() {
                    print(paste("No genes were found for population(s): ",paste(absentSignatures,collapse = ", "),".",sep=""))
                    
                    return("Data")
                  }
                output$console <- renderPrint({
                  logText()
                  return(print(values[["log"]]))
                })
                
                logText <- reactive({
                  values[["log"]] <- capture.output(data <- queryMagic())
                })
                localSig <- mMCPcounter_signatures[mMCPcounter_signatures[,features] %in% foundGenes,]
                expAg <- exp[localSig[,features],]
                expAg <- cit.dfAggregate(expAg,localSig$Denomination,fAggreg = median.na)
                expAg <- expAg[c("T cells", "CD8 T cells", "NK cells", "B derived", "Memory B cells", "Monocytes / macrophages", "Monocytes", "Granulocytes", "Mast cells", "Eosinophils", "Neutrophils", "Basophils", "Vessels", "Lymphatics", "Endothelial cells", "Fibroblasts"),]
                expAg <- expAg[apply(expAg,1,function(x){sum(is.na(x))})<ncol(expAg),]
                return(expAg)
              }      
              ###################################################################      
              #####  Mouse MCP-counter
              m_mcp <- mMCPcounter.estimate(testdata)
              
              
              ############################################################            
              ### scale MCP result before plotting heatmap
              m_mcp_scaled <- t(scale(t(m_mcp)))
              
              ## Remove row with NA value
              m_mcp_scaled <- na.omit(m_mcp_scaled)
              
              #print(m_mcp_scaled)
              #datamcp$m_mcp <- m_mcp
              size <-(length(colnames(m_mcp)))
              
              my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
              
              mcpplot<- pheatmap(m_mcp_scaled, show_rownames = T, angle_col = "90", cluster_cols = F, color = my_palette, cluster_rows = T, fontsize_row = 14, fontsize_col = 14)
              
              
              output$mMCPRes <- DT::renderDataTable({m_mcp},extensions = 'FixedColumns', options = list(pageLength = 10, autoWidth =FALSE,scrollX = TRUE,
                                                                                                        scrollCollapse = TRUE,
                                                                                                        lengthMenu = c(10, 20, 50),
                                                                                                        columnDefs = list(list(targets = c(1:size), searchable = FALSE))),
                                                    
                                                    editable = FALSE)
              #End of MCP
              
              output$TableDownloadmMCP <- downloadHandler(  
                filename = function(){paste("mMCPRes",'.csv')},
                content = function(file){
                  write.csv(m_mcp, file, row.names = TRUE)
                }) 
              
              observe({  
                widthmcp <- reactive ({ input$plot_widthmcp})
                heightmcp <- reactive ({ input$plot_heightmcp}) 
                output$mMCPlot <- renderPlot(width = widthmcp, height = heightmcp,{
                  mcpplot
                })  
                
              })
              
              observe({
                chmcp <- input$optionmcp
                if (chmcp == "1"){ 
                  output$PlotDownloadmMCP <- downloadHandler(
                    filename <- function() {
                      paste("mMCP-Plot",".png", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::png(..., width = input$plot_widthmcp*4, height =input$plot_heightmcp*4 , res = 350, units = "px")
                      ggsave(file, plot =  mcpplot, device = device)
                      # dev.off()
                    })
                } else if  (chmcp == "2"){ 
                  
                  output$PlotDownloadmMCP <- downloadHandler(
                    filename <- function() {
                      paste("mMCP-Plot",".svg", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthmcp*4/300), height =round(input$plot_heightmcp*4/300))
                      ggsave(file, plot =  mcpplot , device = device)
                      #device <- function(..., width, height) grDevices::svg(..., width = 16, height = 10)
                     # ggsave(file, plot =  mcpplot, device = device)
                     # dev.off()
                    })
                }
              }) #End of Observe download
            } else if (chconter== "Human"){
              testdata<- Collapsed_data[,-c(1,2)]
              
              col_name10<- paste(datar$n1, sep="")
              col_name11<-paste(datar$n2, sep="")
              
              names(testdata)[1:length(datar$n1)]<-col_name10
              names(testdata)[length(datar$n1)+1:length(datar$n2)]<-col_name11
              
              #####  Human MCP-counter
              mcp <- MCPcounter.estimate(testdata, featuresType = 'HUGO_symbols')
              
              
              ### scale MCP result before plotting heatmap
              mcp_scaled <- t(scale(t(mcp)))
              
              ## Remove row with NA value
              mcp_scaled <- na.omit(mcp_scaled)
              #datamcp$mcp <- mcp
              size <-(length(colnames(mcp)))
              
              my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
              
              mcpplot<- pheatmap(mcp_scaled, show_rownames = T, angle_col = "90", cluster_cols = F, color = my_palette, cluster_rows = T, fontsize_row = 14, fontsize_col = 14)  
              
              output$mMCPRes <- DT::renderDataTable({mcp},extensions = 'FixedColumns', options = list(pageLength = 10, autoWidth =FALSE,scrollX = TRUE,
                                                                                                      scrollCollapse = TRUE,
                                                                                                      lengthMenu = c(10, 20, 50),
                                                                                                      columnDefs = list(list(targets = c(1:size), searchable = FALSE))),
                                                    
                                                    editable = FALSE)
              
              output$TableDownloadmMCP <- downloadHandler(  
                filename = function(){paste("MCPRes",'.csv')},
                content = function(file){
                  write.csv(mcp, file, row.names = TRUE)
                }) 
              
              observe({  
                widthmcp <- reactive ({ input$plot_widthmcp})
                heightmcp <- reactive ({ input$plot_heightmcp}) 
                output$mMCPlot <- renderPlot(width = widthmcp, height = heightmcp,{
                  mcpplot
                })  
                
              })
              
              observe({
                chmcp <- input$optionmcp
                if (chmcp == "1"){ 
                  output$PlotDownloadmMCP <- downloadHandler(
                    filename <- function() {
                      paste("MCP-Plot",".png", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::png(..., width = input$plot_widthmcp*4, height =input$plot_heightmcp*4 , res = 350, units = "px")
                      ggsave(file, plot =  mcpplot, device = device)
                      # dev.off()
                    })
                } else if  (chmcp == "2"){ 
                  
                  output$PlotDownloadmMCP <- downloadHandler(
                    filename <- function() {
                      paste("MCP-Plot",".svg", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthmcp*4/300), height =round(input$plot_heightmcp*4/300))
                      ggsave(file, plot =  mcpplot , device = device)
                     # device <- function(..., width, height) grDevices::svg(..., width = 16, height = 10)
                      #ggsave(file, plot =  mcpplot, device = device)
                     # dev.off()
                    })
                }
              }) #End of Observe download
              
            } 
            
            
          })#End of MCP
          
#################################################################################################################  
######## Gene Enrichment ########################################################################################################      
          datass <- reactiveValues()
          observeEvent(input$ssGSEArAnalysis,{
            Collapsed_data <- datar$Coll
            head(Collapsed_data)
           
            ########## Mouse data############################################   
            chspecies <- input$optionspecies
            if (chspecies == "Mouse"){   
              ####### Get the normalized data from the first part
              Mouse.Data <- Collapsed_data[,-c(1,2)]
              col_name8<- paste(datar$n1, sep="")
              col_name9<-paste(datar$n2, sep="")
              
              names(Mouse.Data)[1:length(datar$n1)]<-col_name8
              names(Mouse.Data)[length(datar$n1)+1:length(datar$n2)]<-col_name9
              

              
              ##########################################ssGSEA for Mouse data
              ch <- input$option
              if (ch == "Hallmark"){
                withProgress(message = "Preparing Hallmark Table and Plot, please wait",
                             {
                               load("HallMark_Mouse.rdata")
                               datass$ssgsea.mouse <- gsva(as.matrix(Mouse.Data), x,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                               ssgsea.mouse.scaled <- t(scale(t(datass$ssgsea.mouse)))
                               my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
                               datass$hallp<- pheatmap(ssgsea.mouse.scaled, angle_col = "90", color=my_palette, fontsize_row = 8, fontsize_col = 10, 
                                                       border_color = NA, family = "Helvetica",treeheight_row = 8, legend = T, cluster_cols = F)
                             })
              } 
              
              else if (ch == "GeneOntology") {
                withProgress(message = "Preparing GeneOntology Table, please wait",
                             {
                               load("GO_Mouse.rdata")
                               datass$ssgsea.mouse <- gsva(as.matrix(Mouse.Data), x,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                             })
              }
              
              else if (ch=="USER_Geneset"){
                withProgress(message = "Preparing USER_Geneset Table, please wait",
                             {
                               file_to_read5 <- input$file5
                               if(is.null(file_to_read5)){return("")}
                               file6 <-  file_to_read5$datapath
                               e = new.env()
                               name <- load(file6, envir = e)
                               data <- e[[name]]
                               datass$ssgsea.mouse <- gsva(as.matrix(Mouse.Data), data,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                             })
              }
              
              
              datass$s <-(length(colnames(datass$ssgsea.mouse)))
              
              output$ssGSEA <- DT::renderDataTable({datass$ssgsea.mouse},
                                                   extensions = 'FixedColumns', 
                                                   options = list(pageLength = 20, autoWidth = FALSE,
                                                                  scrollX = TRUE,
                                                                  scrollCollapse = TRUE,
                                                                  lengthMenu = c(20, 50, 100),
                                                                  columnDefs = list(list(targets = c(1:datass$s), searchable = FALSE))),
                                                   editable = FALSE)
              
              output$TableDownloadssGSEA <- downloadHandler(
                filename = function(){paste("ssGSEA.Mouse",'.csv')},
                content = function(file){
                  write.csv(datass$ssgsea.mouse, file, sep = ",", quote = F)
                })
              
              
              
              observe({  
                widthhall <- reactive ({ input$plot_widthhall})
                heighthall <- reactive ({ input$plot_heighthall}) 
                output$hallPlot <- renderPlot(width = widthhall, height = heighthall,{
                  datass$hallp
                })  
                
              })
              
              observe({
                chhall <- input$optionhall
                if (chhall == "1"){ 
                  output$PlotDownloadhall <- downloadHandler(
                    filename <- function() {
                      paste("Hallmark-Plot",".png", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::png(..., width = input$plot_widthhall*4, height =input$plot_heighthall*4 , res = 350, units = "px")
                      ggsave(file, plot = datass$hallp, device = device)
                      # dev.off()
                    })
                } else if  (chhall== "2"){ 
                  
                  output$PlotDownloadhall <- downloadHandler(
                    filename <- function() {
                      paste("Hallmark-Plot",".svg", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthhall*4/300), height =round(input$plot_heighthall*4/300))
                      ggsave(file, plot = datass$hallp , device = device)
                    
                    })
                }
              })
              
            }
            
            #############################################################################################
            else if (chspecies == "Human"){  
              Human.Data <- Collapsed_data[,-c(1,2)]
              col_name8<- paste(datar$n1, sep="")
              col_name9<-paste(datar$n2, sep="")
              
              names(Human.Data)[1:length(datar$n1)]<-col_name8
              names(Human.Data)[length(datar$n1)+1:length(datar$n2)]<-col_name9
              
              ch <- input$option
              if (ch == "Hallmark"){
                withProgress(message = "Preparing Hallmark Table and Plot, please wait",
                             {
                               load("HallMark_Human.rdata")
                               datass$ssgsea.human <- gsva(as.matrix(Human.Data), x,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                               ssgsea.human.scaled <- t(scale(t(datass$ssgsea.human)))
                               my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
                               datass$hallp<- pheatmap(ssgsea.human.scaled, angle_col = "90", color=my_palette, fontsize_row = 8, fontsize_col = 10, 
                                                       border_color = NA, family = "Helvetica",treeheight_row = 8, legend = T, cluster_cols = F)
                             })
              } 
              
              
              
              else if (ch == "GeneOntology") {
                withProgress(message = "Preparing GeneOntology Table, please wait",
                             {
                               load("GO_Human.rdata")
                               datass$ssgsea.human <- gsva(as.matrix(Human.Data), x,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                             })
              }
              
              else if (ch=="USER_Geneset"){
                withProgress(message = "Preparing USER_Geneset Table, please wait",
                             {
                               file_to_read5 <- input$file5
                               if(is.null(file_to_read5)){return("")}
                               file6 <-  file_to_read5$datapath
                               e = new.env()
                               name <- load(file6, envir = e)
                               data <- e[[name]]
                               datass$ssgsea.human<- gsva(as.matrix(Human.Data), data,
                                                          min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                             })
              }
              
              
              datass$s <-(length(colnames(datass$ssgsea.human)))
              
              output$ssGSEA <- DT::renderDataTable({datass$ssgsea.human},
                                                   extensions = 'FixedColumns', 
                                                   options = list(pageLength = 20, autoWidth = FALSE,
                                                                  scrollX = TRUE,
                                                                  scrollCollapse = TRUE,
                                                                  lengthMenu = c(20, 50, 100),
                                                                  columnDefs = list(list(targets = c(1:datass$s), searchable = FALSE))),
                                                   editable = FALSE)
              
              
              output$TableDownloadssGSEA <- downloadHandler(
                filename = function(){paste("ssGSEA.Human",'.csv')},
                content = function(file){
                  write.csv(datass$ssgsea.human, file, sep = ",", quote = F)
                })
              
              
              
              observe({  
                widthhall <- reactive ({ input$plot_widthhall})
                heighthall <- reactive ({ input$plot_heighthall}) 
                output$hallPlot <- renderPlot(width = widthhall, height = heighthall,{
                  datass$hallp
                })  
                
              })
              
              
              
              observe({
                chhall <- input$optionhall
                if (chhall == "1"){ 
                  output$PlotDownloadhall <- downloadHandler(
                    filename <- function() {
                      paste("Hallmark-Plot",".png", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::png(..., width = input$plot_widthhall*4, height =input$plot_heighthall*4 , res = 350, units = "px")
                      ggsave(file, plot = datass$hallp, device = device)
                      # dev.off()
                    })
                } else if  (chhall== "2"){ 
                  
                  output$PlotDownloadhall <- downloadHandler(
                    filename <- function() {
                      paste("Hallmark-Plot",".svg", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthhall*4/300), height =round(input$plot_heighthall*4/300))
                      ggsave(file, plot = datass$hallp , device = device)
                     # device <- function(..., width, height) grDevices::svg(..., width = 16, height = 10)
                      #ggsave(file, plot = datass$hallp, device = device)
                     # dev.off()
                    })
                }
              })
            }  
          })#End Gene Enrichment 
          ##################################################  
          ###GESA PLOT######################## 
          datagesa <- reactiveValues()
          observeEvent(input$gseaplot,{
            withProgress(message = "Preparing list of Hallmark-Pathways, please wait...",
                         {
                           fs1<-data$exp
                           head(fs1)
                           rownames(fs1)<-fs1[,1]
                           fs1<- fs1[order(fs1$ID, decreasing = FALSE), ]
                           
                           file_to_read2=input$file2
                           ds1 <- read.delim(file_to_read2$datapath,sep=input$Sep1, header=T)
                           names(ds1)[1] <- "Types"
                           names(ds1)[2] <- "Samples"
                           names(ds1)[3] <- "Labels"
                           melfs1<-melt(fs1)
                           
                           colnames(melfs1) <- c("ID", "symbol", "Samples","value")
                           
                           fss1<-join(melfs1, ds1, by = "Samples", type = "inner")
                           
                           g11 <- subset(fss1, Types==as.character(input$Group3comp)[1])
                           g11 = subset(g11, select = -c(Types,Labels) )
                           colnames(g11) <- c("ID", "symbol", "variable","value")
                           g11<-cast(g11)
                           g11 = subset(g11, select = -c(ID,symbol) )
                           Big11<-g11
                           
                           k11 <- subset(fss1, Types==as.character(input$Group4comp)[1])
                           k11 = subset(k11, select = -c(Types,Labels) )
                           colnames(k11) <- c("ID", "symbol", "variable","value")
                           k11<-cast(k11)
                           k11 = subset(k11, select = -c(ID,symbol) )
                           Big12<-k11
                           main11<-cbind(Big11,Big12)
                           
                           num2 <- sapply( main11, is.integer)
                           if(all(num2)!= TRUE){
                             shinyalert("Warning!", "Input Data includes Decimal number", type = "error") 
                             Sys.sleep(10)
                           }
                           if (all(num2)== TRUE){
                             
                           main1<-cbind(fs1[c(1,2)],main11)
                           }
                           gr1 <- c(rep(as.character(input$Group3comp)[1],length(Big11) ), rep(as.character(input$Group4comp)[1],length(Big12)))
                           grss1 = factor(gr1, levels = unique(gr1))
                           
                           
                           #gr1 = factor(ds1[,1])
                           colData1 <- data.frame(group=grss1)
                           # grss1 = factor(ds1[,1], levels = unique(ds1[,1]))
                           rownames(colData1) <- colnames(main1[,-c(1,2)])
                           cds1 <- DESeqDataSetFromMatrix(main1[,-c(1,2)], colData1, design = ~group)
                           class(cds1)
                           cds1 <- DESeq(cds1)
                           DEGs1<- results(cds1, c("group", levels(grss1)))
                           DEGs1$symbol<- fs1$symbol
                           DEGs1 <- na.omit(DEGs1)
                           pdat <- data.frame(DEGs1)
                           # pdat<- subset(pdat,pdat$pvalue<0.05)
                           names(pdat)[names(pdat) == "selectedRowID"] <- "row"
                           rownames(pdat) <- c()
                           res<-pdat
                           res2 <- res %>%
                             dplyr::select(symbol, log2FoldChange) %>%
                             #na.omit() %>%
                             distinct() %>%
                             group_by(symbol) %>%
                             summarize(log2FoldChange=mean(log2FoldChange))
                           #arrange(des(stat))
                           res2
                           
                           ranks <- deframe(res2)
                           ranks.2 <- sort(ranks, decreasing = T) 
                           
                           chspeciesplot <- input$optionspeciesplot 
                           if (chspeciesplot == "Mouse"){                  
                             load("HallMark_Mouse.rdata")
                             pat <- x
                             set.seed(127)
                             fgseaRes <- fgsea(pathways=pat, stats=ranks.2, nperm=10000,minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) %>%
                               as_tibble() %>%
                               arrange(padj)
                             
                             filtered_pathway <- fgseaRes
                             
                             filt_p <- as.vector(filtered_pathway$pathway)
                             
                             output$check5 <- renderUI({
                               selectInput ('Group5comp','List of Pathways', as.list(filtered_pathway$pathway),multiple = FALSE)
                             })
                             
                             
                            
                           }
                           else if  (chspeciesplot == "Human"){
                             load("HallMark_Human.rdata")
                             pat <- x
                             set.seed(127)
                             fgseaRes <- fgsea(pathways=pat, stats=ranks.2, nperm=10000,minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) %>%
                               as_tibble() %>%
                               arrange(padj)
                             
                             filtered_pathway <- fgseaRes
                             
                             
                             filt_p <- as.vector(filtered_pathway$pathway)
                             
                             output$check5 <- renderUI({
                               selectInput ('Group5comp','List of Pathways', as.list(filtered_pathway$pathway),multiple = FALSE)
                             })
                             
                           }
                           
                           
                           m<-filt_p[1]
                           
                           pannel<-filtered_pathway %>% filter_all(any_vars(. %in% c(m)))
                           pannel<-data.frame(rbind(pannel)) 
                           
                           
                           Pval<-pannel$pval
                           padj <-pannel$padj
                           ES<-pannel$ES
                           NES <-pannel$NES
                           
                           pp = formatC(Pval, digits = 4, format = "f")
                           dd = formatC( padj, digits = 4, format = "f")
                           ES = formatC(ES, digits = 4, format = "f")
                           NES = formatC( NES, digits = 4, format = "f")
                           
                           
                           pathway=pat[[m]]
                           stats=ranks.2 
                           gseaParam = 1
                           ticksSize = 0.3
                           
                           rnk <- rank(-stats)
                           ord <- order(rnk)
                           statsAdj <- stats[ord]
                           statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
                           statsAdj <- statsAdj/max(abs(statsAdj))
                           pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
                           pathway <- sort(pathway)
                           gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                                                   returnAllExtremes = TRUE)
                           bottoms <- gseaRes$bottoms
                           tops <- gseaRes$tops
                           n <- length(statsAdj)
                           xs <- as.vector(rbind(pathway - 1, pathway))
                           ys <- as.vector(rbind(bottoms, tops))
                           toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
                           diff <- (max(tops) - min(bottoms))/8
                           x = y = NULL
                           
                           pen<-ggplot(toPlot, aes(x = x, y = y))+geom_point(color = "green",  size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                                                                        linetype = "dashed") + geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
                             geom_hline(yintercept = 0,colour = "black") + geom_line(color = "green")+ theme_bw() + geom_segment(data = data.frame(x = pathway), 
                                                                                                                                 mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + theme(panel.border = element_blank(), 
                                                                                                                                                                                                                       panel.grid.minor = element_blank()) + labs(x = "Rank", y = "Enrichment Score (ES)",tag = paste(as.character(input$Group3comp)[1]))+
                             labs(title =m,subtitle = paste("P-value:",pp,", FDR:",dd,", ES:",ES,", NES:",NES))+ theme(plot.title = element_text(hjust = 0.5, face="bold")
                                                                                                                       ,plot.subtitle  = element_text(hjust = 0.5, face="bold", color="red"))+ theme(plot.tag.position = c(0.15, 0.02))                                                                             
                           
                           datagesa$gseaplot<-pen+labs(x=paste(as.character(input$Group4comp)[1])) + theme(axis.title.x = element_text(size=13,hjust=0.9))
                           
                           lead<-pannel$leadingEdge
                           test<- lead[[1]]
                           numgene<-length(test)
                           
                           datagesa$lead <- lead
                           
                           output$leading <- renderText({
                             isolate({
                               paste("Number of leadingEdge Genes in this pathway:", as.numeric(numgene), sep="")
                             })
                             
                           })       
                           
                               
                           observe({
                             i<-as.character(input$Group5comp)[1]
                             pannel<-filtered_pathway %>% filter_all(any_vars(. %in% c(i)))
                             pannel<-data.frame(rbind(pannel)) 
                             
                             # observe({ 
                             if (!is.na(i)){
                               lead<-pannel$leadingEdge
                               test<- lead[[1]]
                               numgene<-length(test)
                               
                               output$leading <- renderText({
                                 isolate({
                                   paste("Number of leadingEdge Genes in this pathway:", as.numeric(numgene), sep="")
                                 })
                               })
                               datagesa$lead <- lead
                             }
                             # })
                             
                             Pval<-pannel$pval
                             padj <-pannel$padj
                             ES<-pannel$ES
                             NES <-pannel$NES
                             
                             pp = formatC(Pval, digits = 4, format = "f")
                             dd = formatC( padj, digits = 4, format = "f")
                             ES = formatC(ES, digits = 4, format = "f")
                             NES = formatC( NES, digits = 4, format = "f")
                             
                             pathway=pat[[i]]
                             stats=ranks.2 
                             gseaParam = 1
                             ticksSize = 0.3
                             
                             rnk <- rank(-stats)
                             ord <- order(rnk)
                             statsAdj <- stats[ord]
                             statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
                             statsAdj <- statsAdj/max(abs(statsAdj))
                             pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
                             pathway <- sort(pathway)
                             gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                                                     returnAllExtremes = TRUE)
                             bottoms <- gseaRes$bottoms
                             tops <- gseaRes$tops
                             n <- length(statsAdj)
                             xs <- as.vector(rbind(pathway - 1, pathway))
                             ys <- as.vector(rbind(bottoms, tops))
                             toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
                             diff <- (max(tops) - min(bottoms))/8
                             x = y = NULL
                             
                             
                             
                             pen<-ggplot(toPlot, aes(x = x, y = y))+geom_point(color = "green",  size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                                                                          linetype = "dashed") + geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
                               geom_hline(yintercept = 0,colour = "black") + geom_line(color = "green")+ theme_bw() + geom_segment(data = data.frame(x = pathway), 
                                                                                                                                   mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + theme(panel.border = element_blank(), 
                                                                                                                                                                                                                         panel.grid.minor = element_blank()) + labs(x = "Rank", y = "Enrichment Score (ES)",tag = paste(as.character(input$Group3comp)[1]))+
                               labs(title =i,subtitle = paste("P-value:",pp,", FDR:",dd,", ES:",ES,", NES:",NES))+ theme(plot.title = element_text(hjust = 0.5, face="bold")
                                                                                                                         ,plot.subtitle  = element_text(hjust = 0.5, face="bold", color="red"))+ theme(plot.tag.position = c(0.15, 0.02))                                                                             
                             
                             datagesa$gseaplot<-pen+labs(x=paste(as.character(input$Group4comp)[1])) + theme(axis.title.x = element_text(size=13,hjust=0.9))
                             
                             
                           })
                           
                         }) 
          })
          observeEvent(input$gegoplot,{
            withProgress(message = "Preparing list of GeneOntology-Pathways, please wait...",
                         {
                           fs1<-data$exp
                           head(fs1)
                           rownames(fs1)<-fs1[,1]
                           fs1<- fs1[order(fs1$ID, decreasing = FALSE), ]
                           
                           file_to_read2=input$file2
                           ds1 <- read.delim(file_to_read2$datapath,sep=input$Sep1, header=T)
                           names(ds1)[1] <- "Types"
                           names(ds1)[2] <- "Samples"
                           names(ds1)[3] <- "Labels"
                           
                           melfs1<-melt(fs1)
                           
                           colnames(melfs1) <- c("ID", "symbol", "Samples","value")
                           
                           fss1<-join(melfs1, ds1, by = "Samples", type = "inner")
                           
                           g11 <- subset(fss1, Types==as.character(input$Group3comp)[1])
                           g11 = subset(g11, select = -c(Types,Labels) )
                           colnames(g11) <- c("ID", "symbol", "variable","value")
                           g11<-cast(g11)
                           g11 = subset(g11, select = -c(ID,symbol) )
                           Big11<-g11
                           
                           k11 <- subset(fss1, Types==as.character(input$Group4comp)[1])
                           k11 = subset(k11, select = -c(Types,Labels) )
                           colnames(k11) <- c("ID", "symbol", "variable","value")
                           k11<-cast(k11)
                           k11 = subset(k11, select = -c(ID,symbol) )
                           Big12<-k11
                           main11<-cbind(Big11,Big12)
                           
                           num2 <- sapply( main11, is.integer)
                           if(all(num2)!= TRUE){
                             shinyalert("Warning!", "Input Data includes Decimal number", type = "error") 
                             Sys.sleep(10)
                           }
                           if (all(num2)== TRUE){
                             
                             main1<-cbind(fs1[c(1,2)],main11)
                           }
                           
                          
                           
                           gr1 <- c(rep(as.character(input$Group3comp)[1],length(Big11) ), rep(as.character(input$Group4comp)[1],length(Big12)))
                           grss1 = factor(gr1, levels = unique(gr1))
                           
                           
                         
                           colData1 <- data.frame(group=grss1)
                         
                           rownames(colData1) <- colnames(main1[,-c(1,2)])
                           cds1 <- DESeqDataSetFromMatrix(main1[,-c(1,2)], colData1, design = ~group)
                           class(cds1)
                           cds1 <- DESeq(cds1)
                           DEGs1<- results(cds1, c("group", levels(grss1)))
                           DEGs1$symbol<- fs1$symbol
                           DEGs1 <- na.omit(DEGs1)
                           pdat <- data.frame(DEGs1)
                         
                           names(pdat)[names(pdat) == "selectedRowID"] <- "row"
                           rownames(pdat) <- c()
                           res<-pdat
                           res2 <- res %>%
                             dplyr::select(symbol, log2FoldChange) %>%
                             #na.omit() %>%
                             distinct() %>%
                             group_by(symbol) %>%
                             summarize(log2FoldChange=mean(log2FoldChange))
                           #arrange(des(stat))
                           res2
                           
                           ranks <- deframe(res2)
                           ranks.2 <- sort(ranks, decreasing = T) 
                           
                           chspeciesplot <- input$optionspeciesplot 
                           if (chspeciesplot == "Mouse"){
                             load("GO_Mouse.rdata")
                             pat <- x
                             set.seed(127)
                             fgseaRes <- fgsea(pathways=pat, stats=ranks.2, nperm=10000,minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) %>%
                               as_tibble() %>%
                               arrange(padj)
                             
                             filtered_pathway<- fgseaRes
                             
                             filtered_pathway<- filtered_pathway[order( filtered_pathway$ES, decreasing = TRUE), ]
                             filtered_pathway<- filtered_pathway[1:50,]
                             
                             
                             
                             filt_p <- as.vector(filtered_pathway$pathway)
                             
                             output$check5 <- renderUI({
                               selectInput ('Group6comp','List of Pathways', as.list(filtered_pathway$pathway),multiple = FALSE)
                             })
                             
                           }
                           else if  (chspeciesplot == "Human"){
                             load("GO_Human.rdata")
                             pat <- x
                             set.seed(127)
                             fgseaRes <- fgsea(pathways=pat, stats=ranks.2, nperm=10000,minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) %>%
                               as_tibble() %>%
                               arrange(padj)
                             
                             filtered_pathway<- fgseaRes
                             
                             filtered_pathway<- filtered_pathway[order( filtered_pathway$ES, decreasing = TRUE), ]
                             filtered_pathway<- filtered_pathway[1:50,]
                             
                             
                             
                             filt_p <- as.vector(filtered_pathway$pathway)
                             
                             output$check5 <- renderUI({
                               selectInput ('Group6comp','List of Pathways', as.list(filtered_pathway$pathway),multiple = FALSE)
                             })
                             
                           }
                           m<-filt_p[1]
                           
                           pannel2<-filtered_pathway %>% filter_all(any_vars(. %in% c(m)))
                           pannel2<-data.frame(rbind(pannel2)) 
                           
                           
                           Pval<-pannel2$pval
                           padj <-pannel2$padj
                           ES<-pannel2$ES
                           NES <-pannel2$NES
                           
                           pp = formatC(Pval, digits = 4, format = "f")
                           dd = formatC( padj, digits = 4, format = "f")
                           ES = formatC(ES, digits = 4, format = "f")
                           NES = formatC( NES, digits = 4, format = "f")
                           
                           
                           pathway=pat[[m]]
                           stats=ranks.2 
                           gseaParam = 1
                           ticksSize = 0.3
                           
                           rnk <- rank(-stats)
                           ord <- order(rnk)
                           statsAdj <- stats[ord]
                           statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
                           statsAdj <- statsAdj/max(abs(statsAdj))
                           pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
                           pathway <- sort(pathway)
                           gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                                                   returnAllExtremes = TRUE)
                           bottoms <- gseaRes$bottoms
                           tops <- gseaRes$tops
                           n <- length(statsAdj)
                           xs <- as.vector(rbind(pathway - 1, pathway))
                           ys <- as.vector(rbind(bottoms, tops))
                           toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
                           diff <- (max(tops) - min(bottoms))/8
                           x = y = NULL
                           
                           pen<-ggplot(toPlot, aes(x = x, y = y))+geom_point(color = "green",  size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                                                                        linetype = "dashed") + geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
                             geom_hline(yintercept = 0,colour = "black") + geom_line(color = "green")+ theme_bw() + geom_segment(data = data.frame(x = pathway), 
                                                                                                                                 mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + theme(panel.border = element_blank(), 
                                                                                                                                                                                                                       panel.grid.minor = element_blank()) + labs(x = "Rank", y = "Enrichment Score (ES)",tag = paste(as.character(input$Group3comp)[1]))+
                             labs(title =m,subtitle = paste("P-value:",pp,", FDR:",dd,", ES:",ES,", NES:",NES))+ theme(plot.title = element_text(hjust = 0.5, face="bold")
                                                                                                                       ,plot.subtitle  = element_text(hjust = 0.5, face="bold", color="red"))+ theme(plot.tag.position = c(0.15, 0.02))                                                                             
                           
                           datagesa$gseaplot<-pen+labs(x=paste(as.character(input$Group4comp)[1])) + theme(axis.title.x = element_text(size=13,hjust=0.9))
                           
                           lead2<-pannel2$leadingEdge
                           test2<- lead2[[1]]
                           numgene2<-length(test2)
                           
                           datagesa$lead <- lead2
                           
                           output$leading <- renderText({
                             isolate({
                               paste("Number of leadingEdge Genes in this pathway:", as.numeric(numgene2), sep="")
                             })
                             
                           })       
                           
                           observe({
                             i<-as.character(input$Group6comp)[1]
                             pannel2<-filtered_pathway %>% filter_all(any_vars(. %in% c(i)))
                             pannel2<-data.frame(rbind(pannel2)) 
                             
                             if (!is.na(i)){
                               lead2<-pannel2$leadingEdge
                               test2<- lead2[[1]]
                               numgene2<-length(test2)
                               
                               output$leading <- renderText({
                                 isolate({
                                   paste("Number of leadingEdge Genes in this pathway:", as.numeric(numgene2), sep="")
                                 })
                               })
                               datagesa$lead <- lead2
                             }
                             # })
                             
                             Pval<-pannel2$pval
                             padj <-pannel2$padj
                             ES<-pannel2$ES
                             NES <-pannel2$NES
                             
                             pp = formatC(Pval, digits = 4, format = "f")
                             dd = formatC( padj, digits = 4, format = "f")
                             ES = formatC(ES, digits = 4, format = "f")
                             NES = formatC( NES, digits = 4, format = "f")
                             
                             pathway=pat[[i]]
                             stats=ranks.2 
                             gseaParam = 1
                             ticksSize = 0.3
                             
                             rnk <- rank(-stats)
                             ord <- order(rnk)
                             statsAdj <- stats[ord]
                             statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
                             statsAdj <- statsAdj/max(abs(statsAdj))
                             pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
                             pathway <- sort(pathway)
                             gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                                                     returnAllExtremes = TRUE)
                             bottoms <- gseaRes$bottoms
                             tops <- gseaRes$tops
                             n <- length(statsAdj)
                             xs <- as.vector(rbind(pathway - 1, pathway))
                             ys <- as.vector(rbind(bottoms, tops))
                             toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
                             diff <- (max(tops) - min(bottoms))/8
                             x = y = NULL
                             
                             
                             
                             pen<-ggplot(toPlot, aes(x = x, y = y))+geom_point(color = "green",  size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                                                                          linetype = "dashed") + geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
                               geom_hline(yintercept = 0,colour = "black") + geom_line(color = "green")+ theme_bw() + geom_segment(data = data.frame(x = pathway), 
                                                                                                                                   mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + theme(panel.border = element_blank(), 
                                                                                                                                                                                                                         panel.grid.minor = element_blank()) + labs(x = "Rank", y = "Enrichment Score (ES)",tag = paste(as.character(input$Group3comp)[1]))+
                               labs(title =i,subtitle = paste("P-value:",pp,", FDR:",dd,", ES:",ES,", NES:",NES))+ theme(plot.title = element_text(hjust = 0.5, face="bold")
                                                                                                                         ,plot.subtitle  = element_text(hjust = 0.5, face="bold", color="red"))+ theme(plot.tag.position = c(0.15, 0.02))                                                                             
                             
                             datagesa$gseaplot<-pen+labs(x=paste(as.character(input$Group4comp)[1])) + theme(axis.title.x = element_text(size=13,hjust=0.9))
                             
                             
                           })     
                           
                           
                         })
          })
          
          
          output$listofgenes<- downloadHandler(
            filename = function(){paste("LeadingEdgeGenes",'.csv')},
            content = function(file){
              write.table(as.data.frame(datagesa$lead), file, quote=F ,row.names=FALSE,col.names=FALSE)
            })
          
          
          observe({  
            widthgsea <- reactive ({ input$plot_widthgsea })
            heightgsea <- reactive ({ input$plot_heightgsea}) 
            output$gseaPlot <- renderPlot(width = widthgsea, height =  heightgsea,{  
              gseaplot <-datagesa$gseaplot
              vals$gseaplot <-gseaplot
              print(gseaplot)
            })
          }) 
          observe({  
            choptiongsea<- input$optiongsea
            
            if (choptiongsea == "1"){
              output$PlotDownloadgsea <- downloadHandler(
                filename = function(){paste("GSEAPlot",'.png',sep='')},
                content = function(file){
                  device <- function(..., width, height) grDevices::png(..., width =input$plot_widthgsea*4, height =input$plot_heightgsea*4, res = 350, units = "px")
                  ggsave(file, plot = ( vals$gseaplot) , device = device)
                })
            }
            else if (choptiongsea == "2"){
              output$PlotDownloadgsea <- downloadHandler(
                filename = function(){paste("GSEAPlot",'.svg',sep='')},
                content = function(file) {
                  device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthgsea*4/300), height =round(input$plot_heightgsea*4/300))
                  ggsave(file, plot = vals$gseaplot , device = device)
                 # svg(file, width =16, height =10)
                 # print( vals$gseaplot)
                 # dev.off()
                  
                })
            }  
            
          })#observe
          
          
          #output$gseaPlot <- renderPlot({
          #  datagesa$gseaplot
          # })   
          
          
          #################################################################################################################  
          
        } ##End of Raw Count choice
        
        
        #} #End of Normalize File
        
     else if (chFileFormat == "2"){ 
       output$a1<- renderUI({
         isolate(
           actionButton("dataSubmit2", label = "Check Input Files"))
       })
       # if(input$dataSubmit)({ 
       
       observeEvent(input$dataSubmit2,{
            withProgress(message = "Checking your data, please wait",
                         {
##### UpLoding ##################################################################################################################
##### UpLoding ##################################################################################################################
                                                
                                                
                           
values <- reactiveValues()

file_to_read <- input$file1
file_to_read2 <- input$file2

if (is.null(file_to_read ) & is.null(file_to_read2)) { 
  
  queryMagic <- function() {
    print(paste(" Please upload your files first"))
    
    return("Data")
  }
  output$console2 <- renderPrint({
    logText()
    return(print(values[["log"]]))
  })
  
  logText <- reactive({
    values[["log"]] <- capture.output(data <- queryMagic())
  })
  
  return()
}
else if (!is.null(file_to_read ) & !is.null(file_to_read2)) { 
  
  UploadData <- reactiveValues(
    
    mainfile = read.delim(file=file_to_read$datapath,sep=input$Sep1, header=TRUE)) 
  
  Uploadtable <- reactiveValues(
    
    labelfile = read.delim(file=file_to_read2$datapath,sep=input$Sep1, header=TRUE) 
  )
  
  
  if (length(colnames( UploadData$mainfile)) < 4 || length(colnames(Uploadtable$labelfile))> 3) {
    
    queryMagic <- function() {
      if(input$dataSubmit2)
        print(paste("Wrong File Uploaded: Input 1 should be your Main file and Input 2 should be your Labels table"))
      
      return("Data")
    }
    output$console2 <- renderPrint({
      logText()
      return(print(values[["log"]]))
    })
    
    logText <- reactive({
      values[["log"]] <- capture.output(data <- queryMagic())
    })
    
    output$b1<- renderUI({
      if(input$dataSubmit2)
        isolate(
          actionButton("wait", label = "Continue to Submit"))
    }) 
    
    output$error <- renderText({
      if(input$dataSubmit2)
        isolate({
          paste("If there is any error or Input summary table seems wrong, please make sure: 1) You have uploaded correct files with correct format (CSV (Comma delimited)/ txt (Tab delimited)). 2) You have picked correct file Seperators.")
        })
    })
    
  } 
  
  
  else if (length(colnames( UploadData$mainfile)) > 4 && length(colnames(Uploadtable$labelfile))< 4) {
    
    queryMagic <- function() {
      print(paste(" Correct input Uploaded"))
      
      return("Data")
    }
    output$console2 <- renderPrint({
      logText()
      return(print(values[["log"]]))
    })
    
    logText <- reactive({
      values[["log"]] <- capture.output(data <- queryMagic())
    })
    
    if(input$dataSubmit2)
      isolate({
        output$b<- renderUI({
          
          actionButton("dataPCA2", label = " Continue to Submit")
        }) 
        output$b1<- renderUI({
          
        }) 
        
        output$error <- renderText({
          if(input$dataSubmit2)
            isolate({
            })
        })
        
      }) 
  } 
  
}

######## Summary ##############################################################################################################
                             output$Title <- renderText({
                               if(input$dataSubmit2)
                                 isolate({
                                   no.samples <- length(colnames(UploadData$mainfile))-2
                                   paste("This cohort contains ", as.character(no.samples)," Samples:", sep="")
                                 })
                             })
                             ######################################################################################################################
                             output$Info <- renderText({
                               if(input$dataSubmit2)
                                 isolate({
                                   paste(as.character( (Uploadtable$labelfile)[,1]), collapse=", " )
                                 })
                             })
                             #####################################################################################################################
                             output$t1 <- renderText({
                               if(input$dataSubmit2)
                                 isolate({
                                   paste(as.character(colnames(Uploadtable$labelfile)[2]),sep=""," Label:")
                                 })
                             })
                             #####################################################################################################################
                             output$t2 <- renderText({
                               if(input$dataSubmit2)
                                 isolate({
                                   paste(as.character(levels(as.factor((Uploadtable$labelfile)[,2]))), collapse=", ")
                                 })
                             })
                             #######################################################################################################################
                             output$Summary <- renderTable({
                               if(input$dataSubmit2)
                                 isolate({
                                   no.samples <- length(colnames(UploadData$mainfile))-2
                                   no.gene <- dim(UploadData$mainfile)[1]
                                   res.summary <- rbind(no.samples, no.gene)
                                   rownames(res.summary) <- c("Samples", "Gene_ID")
                                   colnames(res.summary) <- "Number"
                                   res.summary
                                 })
     
                 },digits=0,rownames = TRUE, align="lc")
    
    #})
    #################################################################################################################  
    
    # }  
  })  
  
  ##############PCA 2D PLOT  ####################################################
  ##PCA###
  data <- reactiveValues()
  observeEvent(input$dataPCA2,{
    withProgress(message = "Preparing PCA & MDS Plots, please wait",
                 {                  
                           file_to_read=input$file1
                           f6<- read.delim(file_to_read$datapath,sep=input$Sep1, header=T)
                          
                           head(f6)
                           dim(f6)
                           names(f6)[1] <- "ID"
                           names(f6)[2] <- "symbol"
                           # remove rows with zero value
                           f7 <- column_to_rownames(f6, var = "ID")[-1]
                           f7 <- f7[rowSums(f7 != 0)>0,]
                           #### attaching gene symbol to the file
                           f8 <- rownames_to_column(f7, var="ID")
                           exp <- join(f8[,1,drop=F], f6, by="ID", type="inner")
                           data$exp<-exp
                           file_to_read2=input$file2
                           table <- read.delim(file_to_read2$datapath,sep=input$Sep1, header=T)
                           names(table)[1] <- " Types"
                           names(table)[2] <- "Samples"
                           names(table)[3] <- "Labels"
                           
                           rownames(exp) <- exp[,1]
                           exp <- exp[,-c(1,2)]
                           exp.scale <- t(scale(t(exp), scale = F))
                           pc <- prcomp(exp.scale, scale. = F, center = F)
                           percentage <- round(pc$sdev / sum(pc$sdev) * 100, 2)
                           pc$rotation
                           pcr <- data.frame(pc$rotation[,1:3], Group = table[,1] )
                           data$pp <-data.frame(pc$rotation[,1:3], Group = table[,1] )
                           data$percentage <- paste( colnames(pcr), "(", paste( as.character(percentage), "%", ")", sep="") )
                       
                           pcr1 <- data.matrix(pc$rotation[,1:3])
                           pcr2 <- data.frame(pc$rotation[,1:3], Group = table[,1] )
                           pcr$label<- table[,3]
                           data$label<- table[,3]
                           
                           data$gr<-(table[,1])
                           
                           data$plot<-pcr1[,1:3]
                           data$plot3<-pcr2[,1:3]
                           
                          
                           data$table<-table
                           
                           A <- table[,1]
                           B = factor(A)
                           B7=(levels(B)) 
                           data$B7<- length(B7)
                           data$A <-A
                           ##############Group Selection################################                  
                           output$color1 <- renderUI({
                             selectInput('select','Please pick a color for all of your Groups:', as.list(B7), multiple = TRUE) 
                           })
                           
                           output$check1 <- renderUI({
                             selectInput ('Group1comp','Group A', as.list(B7),multiple = TRUE)
                           })
                           
                           output$check2 <- renderUI({
                             selectInput ('Group2comp','Group B', as.list(B7),multiple = TRUE)
                           })
                           
                           output$check3 <- renderUI({
                             selectInput ('Group3comp','Group A', as.list(B7),multiple = TRUE)
                           })
                           
                           output$check4 <- renderUI({
                             selectInput ('Group4comp','Group B', as.list(B7),multiple = TRUE)
                           })
                           
                           norm_list <- list("none","quantile")
                           output$norm <- renderUI({
                             selectInput ('Groupnorm','Choose one option to Normalize your decimal numbers (if needed)', as.list(norm_list),multiple = FALSE)
                           })
                           ############## MDS Plot ##########################################
                           distance.matrix <- dist(scale(t(exp), scale = T,center = T ), method="euclidean")
                           
                           mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
                           
                           ## calculate the percentage of variation that each MDS axis accounts for...
                           data$mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
                           
                           mds.values <- mds.stuff$points
                           data$mds.data <- data.frame(Sample=rownames(mds.values),
                                                       X=mds.values[,1],
                                                       Y=mds.values[,2], Group = table[,1] )
                           
                           ###############################################################################                   
                         }) #progress bar
          }) # End submit
          
          #############Colour ###########################################################  
          
          gg_fill_hue <- function(n) {
            hues = seq(15, 375, length = n + 1)
            hcl(h = hues, l = 65, c = 100)[1:n]
          }
          output$platter <- renderUI({ 
            lev <- sort(unique(input$select)) 
            cols <- gg_fill_hue(length(lev))
            
            lapply(seq_along(lev), function(i) {
              colourInput(inputId = paste0("col", lev[i]),
                          label = paste0("Choose colour for ", lev[i]), 
                          value = cols[i]
              )        
            })
          }) 
          ##############PCA 2D PLOT  ##################################################### 
          observe({
            chp <- input$optionp
            chlab <- input$optionlab
            textpca <- as.numeric(input$label_pca)
            pointpca <- as.numeric(input$point_pca)
            widthpca <- reactive ({ input$plot_widthpca })
            heightpca <- reactive ({ input$plot_heightpca}) 
            if (chp == "1"){
              output$PCAPlot2D <- renderPlot(width = widthpca, height = heightpca,{
                cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
               
                cols <- eval(parse(text = cols))
               
                req(length(cols) == length(input$select))
                if (length(cols) == data$B7){
                  if (chlab == "1"){
                    sp <- ggplot(data$pp, aes(PC1, PC2,label = data$label,fill = Group))+ geom_point(aes(col = Group),size=pointpca)
                    p <- sp+ scale_color_manual(values = cols)+
                      geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[2])}
                  
                  else if (chlab == "2"){
                    sp <- ggplot(data$pp, aes(PC1, PC2, fill = Group)) +geom_point(aes(col = Group),size=pointpca) 
                    p <-sp + scale_color_manual(values = cols)+
                      theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[2])}
               
                  vals$p<-p
                  print(p)
                  vals$cols<-cols
                  
                }
                
                else{
                  if (chlab == "1"){
                    p <-ggplot(data$pp, aes(PC1, PC2, label = data$label, fill = Group))+geom_point(aes(col = Group),size=pointpca)+geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[2])}
                  else if (chlab == "2"){
                    p <-ggplot(data$pp, aes(PC1, PC2, color = Group))+geom_point(aes(col = Group),size=pointpca)+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[2])}
                  vals$p<-p
                  print(p)
                }
              }) 
            }# Chp=1
            
            else if (chp == "2"){
              output$PCAPlot2D <- renderPlot(width = widthpca, height = heightpca,{
                cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
                cols <- eval(parse(text = cols))
              
                req(length(cols) == length(input$select))
                if (length(cols) == data$B7){ 
                  if (chlab == "1"){
                    sp<- ggplot(data$pp, aes(PC1, PC3,label = data$label, fill = Group)) + 
                      
                      geom_point(aes(col = Group),size=pointpca) 
                    
                    p <- sp+ scale_color_manual(values = vals$cols)+
                      geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[3])}
                  else if  (chlab == "2"){
                    sp<- ggplot(data$pp, aes(PC1, PC3, fill = Group)) + 
                      
                      geom_point(aes(col = Group),size=pointpca) 
                  
                    p <- sp+ scale_color_manual(values = vals$cols)+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[3])}
                  
                  vals$p<-p
                
                  print(p)
                }
                
                else{
                  if (chlab == "1"){
                    p <-ggplot(data$pp, aes(PC1, PC3, label = data$label, fill = Group))+geom_point(aes(col = Group),size=pointpca)+geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[3])}
                  else if (chlab == "2"){
                    p <-ggplot(data$pp, aes(PC1, PC3, color = Group))+geom_point(aes(col = Group),size=pointpca)+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[1]) + ylab(data$percentage[3])}
                  vals$p<-p
                  print(p)
                } 
                
                
              })
            }#Chp=2
            
            else if (chp == "3"){
        
              output$PCAPlot2D <- renderPlot(width = widthpca, height = heightpca,{
                cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
               
                cols <- eval(parse(text = cols))
               
                
                # To prevent errors
                req(length(cols) == length(input$select))
                if (length(cols) == data$B7){
                  if (chlab == "1"){
                    sp<- ggplot(data$pp, aes(PC2, PC3,label = data$label, fill = Group)) + 
                      geom_point(aes(col = Group),size=pointpca) 
                    
                    p <- sp +  scale_color_manual(values = vals$cols) +
                      geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[2]) + ylab(data$percentage[3])}
                  
                  else if (chlab == "2"){
                    sp<- ggplot(data$pp, aes(PC2, PC3, fill = Group)) + 
                      geom_point(aes(col = Group),size=pointpca) 
                  
                    p <- sp +  scale_color_manual(values = vals$cols) +theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[2]) + ylab(data$percentage[3])}
                   
                  
                  vals$p<-p
                  print(p)
                }
                else{
                  if (chlab == "1"){
                    p <-ggplot(data$pp, aes(PC2, PC3, label = data$label, fill = Group))+geom_point(aes(col = Group),size=pointpca)+geom_text(size= textpca,position = position_nudge(y = -0.05))+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[2]) + ylab(data$percentage[3])}
                  else if (chlab == "2"){
                    p <-ggplot(data$pp, aes(PC2, PC3, color = Group))+geom_point(aes(col = Group),size=pointpca)+theme_bw(base_size = 14)+
                      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.background = element_blank())+ xlab(data$percentage[2]) + ylab(data$percentage[3])}
                  vals$p<-p
                  print(p)
                  
                } 
                
              })
            }#Chp=3
          })  # end observe  
          ###########################################################################################################      
          ####### MDS PLOT
          observe({
            chmds <- input$optlab
            textmds<- as.numeric(input$label_mds)
            pointmds<- as.numeric(input$point_mds)
            widthmds <- reactive ({ input$plot_widthmds })
            heightmds <- reactive ({ input$plot_heightmds}) 
            output$MDSPlot <- renderPlot(width = widthmds, height = heightmds,{  
              cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
             
              cols <- eval(parse(text = cols))
          
              
            
              req(length(cols) == length(input$select))
              if (length(cols) == data$B7){
                if (chmds == "1"){
                  g<- ggplot(data$mds.data, aes(x=X, y=Y, label = data$label, fill = Group)) + 
                    geom_point(aes(col = Group),size=pointmds) 
                  
                  s <- g +  scale_color_manual(values = vals$cols) +
                    geom_text(size= textmds,position = position_nudge(y = +5))+theme_bw(base_size = 14)+
                    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          panel.background = element_blank())+ xlab(paste("MDS1 - ", data$mds.var.per[1], "%", sep="")) +ylab(paste("MDS2 - ", data$mds.var.per[2], "%", sep="")) }
                
                else if (chmds == "2"){
                  g<- ggplot(data$mds.data, aes(x=X, y=Y, fill = Group)) + 
                    geom_point(aes(col = Group),size=pointmds) 
                
                  s<- g +  scale_color_manual(values = vals$cols) +theme_bw(base_size = 14)+
                    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          panel.background = element_blank())+ xlab(paste("MDS1 - ", data$mds.var.per[1], "%", sep="")) +ylab(paste("MDS2 - ", data$mds.var.per[2], "%", sep="")) }
                
                vals$s<-s
                print(s)
              }
              else{
                if (chmds == "1"){
                  s <-ggplot( data$mds.data, aes(x=X, y=Y, label = data$label, fill = Group))+geom_point(aes(col = Group),size=pointmds)+geom_text(size= textmds,nudge_y = +5)+theme_bw(base_size = 14)+
                    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          panel.background = element_blank())+ xlab(paste("MDS1 - ", data$mds.var.per[1], "%", sep="")) +ylab(paste("MDS2 - ", data$mds.var.per[2], "%", sep="")) }
                else if (chmds == "2"){
                  s <-ggplot( data$mds.data, aes(x=X, y=Y, fill= Group))+geom_point(aes(col = Group),size=pointmds)+theme_bw(base_size = 14)+
                    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          panel.background = element_blank())+ xlab(paste("MDS1 - ", data$mds.var.per[1], "%", sep="")) +ylab(paste("MDS2 - ", data$mds.var.per[2], "%", sep="")) }
                vals$s<-s
                print(s)
                
              } 
            }) 
          })     
          ################ 3D PLOT #####################################################################  
          output$PCAPlot <- renderPlotly({
            withProgress(message = "Preparing PCA3D, please wait",
                         {
                           get_colors <- function(groups, group.col = palette()){
                             groups <- data$gr
                             ngrps <- length(levels(groups))
                             if(ngrps > length(group.col))
                               group.col <- rep(group.col, ngrps)
                             color <- group.col[as.numeric(groups)]
                             names(color) <- as.vector(groups)
                             return(color)
                           }
                           gr<-data$gr
                           n <- data$plot3
                           s=get_colors(gr)
                           m<-as.factor(s)
                           cols <- paste0("c(", paste0("input$col", sort(input$select), collapse = ", "), ")")
                           print(cols)
                           cols <- eval(parse(text = cols))
                           req(length(cols) == length(input$select))
                           if (length(cols) == data$B7){
                             fig <- plot_ly(n, x = ~PC1, y = ~PC2, z = ~PC3, color = ~gr, colors = vals$cols,text = rownames(n))
                             fig <- fig %>% add_markers()
                             chs3d <- input$option3d
                             if (chs3d == "1"){
                               
                               fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                                  yaxis = list(title = 'PC2'),
                                                                  zaxis = list(title = 'PC3'),
                                                                  camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
                               
                               fig <- fig %>% layout(autosize = T)
                               fig <- fig %>% config(
                                 toImageButtonOptions= list(filename = 'PCA3DPlot'))
                               fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                                  yaxis = list(title = 'PC2'),
                                                                  zaxis = list(title = 'PC3')
                               ))
                             }
                             else if (chs3d == "2"){
                               fig <- plot_ly(n, x = ~PC1, y = ~PC2, z = ~PC3, color = ~gr, colors = vals$cols,text = rownames(n))
                               fig <- fig %>% add_markers()
                               
                               fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                                  yaxis = list(title = 'PC2'),
                                                                  zaxis = list(title = 'PC3'),
                                                                  camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
                               
                               fig <- fig %>% layout(autosize = T)
                               fig <- fig %>% config(
                                 toImageButtonOptions= list(format = "svg",filename = 'PCA3DPlot'))
                               fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                                  yaxis = list(title = 'PC2'),
                                                                  zaxis = list(title = 'PC3')
                               ))
                             } 
                             
                           }
                           
                           else
                           {
                             fig <- plot_ly(n, x = ~PC1, y = ~PC2, z = ~PC3, color = ~gr, colors= c("#F8766D","#00BFC4","#7CAE00", "#C77CFF", "#009E73",
                                                                                                    "#7CAE00", "#C77CFF", "#D55E00", "#00BFC4"),text = rownames(n))
                             fig <- fig %>% add_markers()
                             chs3d <- input$option3d
                             if (chs3d == "1"){
                               fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                                  yaxis = list(title = 'PC2'),
                                                                  zaxis = list(title = 'PC3'),
                                                                  camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
                               
                               fig <- fig %>% layout(autosize = T)
                               fig <- fig %>% config(
                                 toImageButtonOptions= list(filename = 'PCA3DPlot'))
                               fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                                  yaxis = list(title = 'PC2'),
                                                                  zaxis = list(title = 'PC3')
                               ))
                             }
                             else  if (chs3d == "2"){
                               fig <- plot_ly(n, x = ~PC1, y = ~PC2, z = ~PC3, color = ~gr, colors= c("#F8766D","#00BFC4","#7CAE00", "#C77CFF", "#009E73",
                                                                                                      "#7CAE00", "#C77CFF", "#D55E00", "#00BFC4"),text = rownames(n))
                               fig <- fig %>% add_markers()
                               
                               fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                                  yaxis = list(title = 'PC2'),
                                                                  zaxis = list(title = 'PC3'),
                                                                  camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
                               
                               fig <- fig %>% layout(autosize = T)
                               fig <- fig %>% config(
                                 toImageButtonOptions= list(format = "svg",filename = 'PCA3DPlot'))
                               fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                                  yaxis = list(title = 'PC2'),
                                                                  zaxis = list(title = 'PC3')
                               ))
                             }
                           }
                           
                         })  #progress bar
          })# ENd Plotly
          
          ###############Downloads PCA2D ###########################################################  
          observe({
            chs <- input$optionchs
            if (chs == "1"){
              output$PlotDownloadPCA2D <- downloadHandler(
                filename = function(){paste("PCAPlot",'.png',sep='')},
                content = function(file){
                  device <- function(..., width, height) grDevices::png(..., width =input$plot_widthpca*4, height =input$plot_heightpca*4, res = 350, units = "px")
                  ggsave(file, plot = (vals$p)  , device = device)
                })
            }
            else if (chs == "2"){
              output$PlotDownloadPCA2D <- downloadHandler(
                filename = function(){paste("PCAPlot",'.svg',sep='')},
                content = function(file) {
                  device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthpca*4/300), height =round(input$plot_heightpca*4/300))
                  ggsave(file, plot = vals$p , device = device)
             
                })
            }
          })
          
          ########### Donwload MDS ###############################################################################   
          observe({
            chd <- input$optchs
            if (chd== "1"){
              output$PlotDownloadMDS2D <- downloadHandler(
                filename = function(){paste("MDSPlot",'.png',sep='')},
                content = function(file){
                  device <- function(..., width, height) grDevices::png(..., width =input$plot_widthmds*4, height =input$plot_heightmds*4, res = 350, units = "px")
                  ggsave(file, plot = (vals$s)  , device = device)
                })
            }
            else if (chd == "2"){
              output$PlotDownloadMDS2D <- downloadHandler(
                filename = function(){paste("MDSPlot",'.svg',sep='')},
                content = function(file) {
                  device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthmds*4/300), height =round(input$plot_heightmds*4/300))
                  ggsave(file, plot = vals$s , device = device)
                
                })
            }
          })# End MDS
          ######################################
          ####CRIS classifier
          datacris <- reactiveValues()
          observeEvent(input$CRISAnalysis,{
            withProgress(message = "Preparing CRIS Classifier Table and Plot, please wait",
                         {
                           fcris<-data$exp
                           head(fcris)
                           rownames(fcris)<-fcris[,1]
                           fcris<-fcris[order( fcris$ID, decreasing = FALSE), ]
                           Collcris <- fcris[,-c(1,2)]
                           collcris<-data.matrix(Collcris)
                           Collcris1<-collcris+1
                           normalized_df <- log2(Collcris1)
                           head(normalized_df)
                           rowGroup <- as.vector(fcris$symbol)
                           rowID <- as.vector(rownames(normalized_df))
                           collapse.object <- collapseRows(datET=normalized_df, rowGroup=rowGroup, rowID=rowID, method = "MaxMean")
                           class(collapse.object)
                           names(collapse.object)
                           Collapsed_data <- data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
                           head(Collapsed_data) 
                           
                           
                           ########## Mouse data############################################   
                           chcris <- input$optionCRIS
                           if (chcris  == "Mouse"){   
                             ####### Get the normalized data from the first part
                             Collapsed_cris <- Collapsed_data[,-c(1,2)]
                             
                             expcris <- data.matrix(Collapsed_cris)
                             
                             # adjust matrix before classification
                             mat_exp_adj <- ematAdjust(expcris)
                             datacris$mat_exp_adj<-mat_exp_adj
                             load("CRIS_template_mouse.RData")
                             
                             NTP_CRIS <- ntp(emat=mat_exp_adj,CIRS_template, seed=367707,doPlot =F)
                             NTP_CRIS <- tibble::rownames_to_column(NTP_CRIS,var = "ID")
                             
                             crisl <-(length(colnames(NTP_CRIS)))
                             
                             output$CRIS<- DT::renderDataTable({NTP_CRIS},
                                                               extensions = 'FixedColumns', 
                                                               options = list(pageLength = 20, autoWidth = FALSE,
                                                                              scrollX = TRUE,
                                                                              scrollCollapse = TRUE,
                                                                              lengthMenu = c(20, 50, 100),
                                                                              columnDefs = list(list(targets = c(3:crisl), searchable = FALSE))),
                                                               editable = FALSE)
                             
                             
                             
                             output$TableDownloadCRIS<- downloadHandler(
                               filename = function(){paste("MOUSE_DATA_NTP_CRIS_results",'.csv')},
                               content = function(file){
                                 write.csv(NTP_CRIS, file, sep = ",", quote = F)
                               })
                             
                             
                             
                             observe({  
                               widthcris <- reactive ({ input$plot_widthcris})
                               heightcris<- reactive ({ input$plot_heightcris}) 
                               output$CRIS_Plot<- renderPlot(width = widthcris, height = heightcris,{
                                 NTP_CRIS <- ntp(emat=mat_exp_adj,CIRS_template, seed=367707,doPlot = TRUE)
                               })  
                               
                             })
                             
                             observe({
                               chdcris <- input$optioncris
                               if ( chdcris  == "1"){ 
                                 
                                 output$PlotDownloadCRIS<- downloadHandler(
                                   file = "CRIS-Plot.png" , # variable with filename
                                   content = function(file) {
                                     png(file = file)
                                     load("CRIS_template_mouse.RData")
                                     NTP_CRIS <- ntp(emat=datacris$mat_exp_adj,CIRS_template, seed=367707,doPlot = TRUE)
                                     dev.off()
                                   })
                               } else if  ( chdcris == "2"){ 
                                 
                                 output$PlotDownloadCRIS <- downloadHandler(
                                   file = "CRIS-Plot.svg" , # variable with filename
                                   content = function(file) {
                                     svg(file = file)
                                     load("CRIS_template_mouse.RData")
                                     NTP_CRIS <- ntp(emat=datacris$mat_exp_adj,CIRS_template, seed=367707,doPlot = TRUE)
                                     dev.off()
                                   })
                               }
                             })
                             
                           }
                           #############################################################################################
                           else if (chcris  == "Human"){  
                             ####### Get the normalized data from the first part
                             Collapsed_cris <- Collapsed_data[,-c(1,2)]
                             
                             
                             expcris <- data.matrix(Collapsed_cris)
                             
                             # adjust matrix before classification
                             exp_adj <- ematAdjust(expcris)
                             datacris$exp_adj<-exp_adj
                             
                             load("Human.CRIS.template.RData")
                             
                             
                             NTP_CRIS <- ntp(emat=exp_adj, Human.CRIS.template, seed=367707,doPlot =F)
                             NTP_CRIS <- tibble::rownames_to_column(NTP_CRIS,var = "ID")
                             
                             crisl <-(length(colnames(NTP_CRIS)))
                             
                             output$CRIS<- DT::renderDataTable({NTP_CRIS},
                                                               extensions = 'FixedColumns', 
                                                               options = list(pageLength = 20, autoWidth = FALSE,
                                                                              scrollX = TRUE,
                                                                              scrollCollapse = TRUE,
                                                                              lengthMenu = c(20, 50, 100),
                                                                              columnDefs = list(list(targets = c(3:crisl), searchable = FALSE))),
                                                               editable = FALSE)
                             
                             
                             
                             output$TableDownloadCRIS<- downloadHandler(
                               filename = function(){paste("human_DATA_NTP_CMS_results",'.csv')},
                               content = function(file){
                                 write.csv(NTP_CRIS, file, sep = ",", quote = F)
                               })
                             
                             
                             
                             observe({  
                               widthcris <- reactive ({ input$plot_widthcris})
                               heightcris<- reactive ({ input$plot_heightcris}) 
                               output$CRIS_Plot<- renderPlot(width = widthcris, height = heightcris,{
                                 NTP_CRIS <- ntp(emat=exp_adj, Human.CRIS.template, seed=367707,doPlot = TRUE)
                               })  
                               
                             })
                             
                             observe({
                               chdcris <- input$optioncris
                               if ( chdcris  == "1"){ 
                                 output$PlotDownloadCRIS<- downloadHandler(
                                   file = "CRIS-Plot.png" , # variable with filename
                                   content = function(file) {
                                     png(file = file)
                                     load("Human.CRIS.template.RData")
                                     NTP_CRIS <- ntp(emat= datacris$exp_adj, Human.CRIS.template, seed=367707,doPlot = TRUE)
                                     dev.off()
                                   })
                               } else if  ( chdcris == "2"){ 
                                 
                                 output$PlotDownloadCRIS <- downloadHandler(
                                   file = "CRIS-Plot.svg" , # variable with filename
                                   content = function(file) {
                                     svg(file = file)
                                     load("Human.CRIS.template.RData")
                                     NTP_CRIS <- ntp(emat= datacris$exp_adj, Human.CRIS.template, seed=367707,doPlot = TRUE)
                                     dev.off()
                                   })
                               }
                             })
                             
                           }
                           
                         })   
          }) 
          
          #########################################################################################################################
          #Classifier_CMS
          datacms <- reactiveValues()
          observeEvent(input$CMSAnalysis,{
            withProgress(message = "Preparing CMS Classifier Table and Plot, please wait",
                         {
                           fcms<-data$exp
                           head(fcms)
                           rownames(fcms)<-fcms[,1]
                           fcms<-fcms[order( fcms$ID, decreasing = FALSE), ]
                           Collcms <- fcms[,-c(1,2)]
                           collcms<-data.matrix(Collcms)
                           Collcms1<-collcms+1
                           normalized_df <- log2(Collcms1)
                           head(normalized_df)
                           rowGroup <- as.vector(fcms$symbol)
                           rowID <- as.vector(rownames(normalized_df))
                           collapse.object <- collapseRows(datET=normalized_df, rowGroup=rowGroup, rowID=rowID, method = "MaxMean")
                           class(collapse.object)
                           names(collapse.object)
                           Collapsed_data <- data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
                           head(Collapsed_data) 
                           
                           # print( head(Collapsed_data))
                           ########## Mouse data############################################   
                           chcms <- input$optionCMS
                           if (chcms  == "Mouse"){   
                             ####### Get the normalized data from the first part
                             Collapsed_cms <- Collapsed_data[,-c(1,2)]
                             
                             expcms <- data.matrix(Collapsed_cms)
                             
                             # adjust matrix before classification
                             mat_exp_adj <- ematAdjust(expcms)
                             datacms$mat_exp_adj<-mat_exp_adj
                             
                             # load mouse CMS template
                             load("template.CMS.A.RData") # template.CMS.A
                             
                             
                             
                             # CMS NTP and save as text file
                             NTP_CMS <- ntp(emat=mat_exp_adj,template.CMS.A, seed=367707,doPlot = F)
                             NTP_CMS <-tibble:: rownames_to_column(NTP_CMS,var = "ID")
                             # write.table(NTP_CMS,file=CMS_res_nm,sep='\t',row.names = FALSE)
                             
                             
                             
                             cmsl <-(length(colnames(NTP_CMS)))
                             
                             output$CMS<- DT::renderDataTable({NTP_CMS},
                                                              extensions = 'FixedColumns', 
                                                              options = list(pageLength = 20, autoWidth = FALSE,
                                                                             scrollX = TRUE,
                                                                             scrollCollapse = TRUE,
                                                                             lengthMenu = c(20, 50, 100),
                                                                             columnDefs = list(list(targets = c(3:cmsl), searchable = FALSE))),
                                                              editable = FALSE)
                             
                             
                             
                             output$TableDownloadCMS<- downloadHandler(
                               filename = function(){paste("MOUSE_DATA_NTP_CMS_results",'.csv')},
                               content = function(file){
                                 write.csv(NTP_CMS, file, sep = ",", quote = F)
                               })
                             
                             
                             
                             observe({  
                               widthcms <- reactive ({ input$plot_widthcms})
                               heightcms<- reactive ({ input$plot_heightcms}) 
                               output$CMS_Plot<- renderPlot(width = widthcms, height = heightcms,{
                                 NTP_CMS<-ntp(emat=mat_exp_adj,template.CMS.A, seed=367707,doPlot = T)
                               })  
                               
                             })
                             
                             observe({
                               chdcms <- input$optioncms
                               if ( chdcms  == "1"){ 
                                 output$PlotDownloadCMS <- downloadHandler(
                                   file = "CMS-Plot.png" , # variable with filename
                                   content = function(file) {
                                     png(file = file)
                                     load("template.CMS.A.RData") # template
                                     NTP_CMS <- ntp(emat=datacms$mat_exp_adj,template.CMS.A, seed=367707,doPlot = T)
                                     dev.off()
                                   })
                               } else if  ( chdcms == "2"){ 
                                 
                                 output$PlotDownloadCMS <- downloadHandler(
                                   file = "CMS-Plot.svg" , # variable with filename
                                   content = function(file) {
                                     svg(file = file)
                                     load("template.CMS.A.RData") # template
                                     NTP_CMS <- ntp(emat=datacms$mat_exp_adj,template.CMS.A, seed=367707,doPlot = T)
                                     dev.off()
                                   })
                               }
                             })
                           }
                           #############################################################################################
                           else if (chcms  == "Human"){  
                             ####### Get the normalized data from the first part
                             Collapsed_cms <- Collapsed_data[,-c(1,2)]
                             
                             expcms <- data.matrix(Collapsed_cms)
                             
                             # adjust matrix before classification
                             exp_adj <- ematAdjust(expcms)
                             datacms$exp_adj<-exp_adj
                             
                             load("Human.CMS.template.RData") # Human.CMS.template
                             
                             
                             # CMS NTP and save as text file
                             NTP_CMS <- ntp(emat=exp_adj, Human.CMS.template, seed=367707,doPlot = F)
                             # datacms$NTP_CMS <- NTP_CMS
                             NTP_CMS <- rownames_to_column(NTP_CMS,var = "ID")
                             # write.table(NTP_CMS,file=CMS_res_nm,sep='\t',row.names = FALSE)
                             
                             cmsl <-(length(colnames(NTP_CMS)))
                             
                             output$CMS<- DT::renderDataTable({NTP_CMS},
                                                              extensions = 'FixedColumns', 
                                                              options = list(pageLength = 20, autoWidth = FALSE,
                                                                             scrollX = TRUE,
                                                                             scrollCollapse = TRUE,
                                                                             lengthMenu = c(20, 50, 100),
                                                                             columnDefs = list(list(targets = c(3:cmsl), searchable = FALSE))),
                                                              editable = FALSE)
                             
                             
                             
                             output$TableDownloadCMS<- downloadHandler(
                               filename = function(){paste("human_DATA_NTP_CMS_results",'.csv')},
                               content = function(file){
                                 write.csv(NTP_CMS, file, sep = ",", quote = F)
                               })
                             
                             
                             
                             observe({  
                               widthcms <- reactive ({ input$plot_widthcms})
                               heightcms<- reactive ({ input$plot_heightcms}) 
                               output$CMS_Plot<- renderPlot(width = widthcms, height = heightcms,{
                                 NTP_CMS <- ntp(emat=exp_adj, Human.CMS.template, seed=367707,doPlot = T)
                               })  
                               
                             })
                             
                             observe({
                               chdcms <- input$optioncms
                               if ( chdcms  == "1"){ 
                                 output$PlotDownloadCMS <- downloadHandler(
                                   file = "CMS-Plot.png" , # variable with filename
                                   content = function(file) {
                                     png(file = file)
                                     load("Human.CMS.template.RData")
                                     NTP_CMS <- ntp(emat=datacms$exp_adj,Human.CMS.template, seed=367707,doPlot = T)
                                     dev.off()
                                   })
                               } else if  ( chdcms == "2"){ 
                                 
                                 output$PlotDownloadCMS <- downloadHandler(
                                   file = "CMS-Plot.svg" , # variable with filename
                                   content = function(file) {
                                     svg(file = file)
                                     load("Human.CMS.template.RData")
                                     NTP_CMS <- ntp(emat=datacms$exp_adj,Human.CMS.template, seed=367707,doPlot = T)
                                     dev.off()
                                   })
                               }
                             })
                             
                           }
                           
                         })   
          }) 
          
          ######################################################################################################### 
          #heatmap 
          #######Heat Map##############
          #datar <- reactiveValues()
          observeEvent(input$FilterAnalysis,{
            showModal( modalDialog("These Options including, Heatmap Plot, Selected Gene List Heatmap Plot,Gene Level Expression and Volacno Plot are not available for The Normalized Read count uploaded Data Matrix",
                                   easyClose = TRUE,
                                   footer = NULL
            ))
            Sys.sleep(10)
            fs<-data$exp
            head(fs)
            
            removeModal()
          })
          #}
          #})  
          ##############################################################################################################  
          #################################################### 
          #Gene Level Expression 
          dat <- reactiveValues()
          observeEvent(input$FilterAnalysis,{
            observe({
              chgene <- input$optiongene
              if (chgene == "1"){
                
                showModal( modalDialog("The Selected Gene List Heatmap Plot Option is not available for The Normalized Read count uploaded Data Matrix",
                                       easyClose = TRUE,
                                       footer = NULL
                ))
                Sys.sleep(5)
                removeModal()
                
              }
              
              else if (chgene == "2"){  
                showModal( modalDialog("The Gene Level Expression Option is not available for The Normalized Read count uploaded Data Matrix",
                                       easyClose = TRUE,
                                       footer = NULL
                ))
                Sys.sleep(5)
                
                removeModal()
              }
              
              
              
            }) # end of Observe
            # }) # end pf Progresss bar  
            #dat$Collapsed_data <- Collapsed_data
          }) # End of Gene Level  
          
          ################################################################# 
          ##########################################################################################################
          ###########MCP COUNTER ############################################################################################# Analysis tab panel content for mMCP-counter
          datamcp <- reactiveValues()
          observeEvent(input$bb,{
            withProgress(message = "Preparing data..., please wait",
                         {
                           fmcp<-data$exp
                           head(fmcp)
                           rownames(fmcp)<-fmcp[,1]
                           fmcp<-fmcp[order( fmcp$ID, decreasing = FALSE), ]
                           Coll2 <- fmcp[,-c(1,2)]
                           coll2<-data.matrix(Coll2)
                           Coll1<-coll2+1
                           normalized_df <- log2(Coll1)
                           head(normalized_df)
                           rowGroup <- as.vector(fmcp$symbol)
                           rowID <- as.vector(rownames(normalized_df))
                           collapse.object <- collapseRows(datET=normalized_df, rowGroup=rowGroup, rowID=rowID, method = "MaxMean")
                           class(collapse.object)
                           names(collapse.object)
                           Collapsed_data <- data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
                           head(Collapsed_data)                
                           
                         })      
            
            chconter<- input$conter       
            if (chconter == "Mouse"){     
              testdata<- Collapsed_data[,-c(1,2)]
              #######################################################################################################################      
              cit.dfAggregate <- function (data, partition, MARGIN = 1, fAggreg = median.na) {
                cMARGIN <- setdiff(c(1, 2), MARGIN)
                n <- length(partition)
                N <- dim(data)[MARGIN]
                p <- dim(data)[cMARGIN]
                if (n != N)
                  stop("ERROR - cit.dfAggregate : size of partition doesn't correspond to data dimension")
                l <- split(1:N, partition)
                d <- data
                if (MARGIN == 2)d <- t(data)
                
                d <- matrix(sapply(l, function(i) if (length(i) == 1) {
                  unlist(d[i, ])
                }else {
                  apply(d[i, ], 2, fAggreg)
                }), ncol = p, byrow = TRUE)
                
                
                
                d <- as.data.frame(d)
                rownames(d) <- names(l)
                names(d) <- dimnames(data)[[cMARGIN]]
                if (MARGIN == 2)
                  d <- as.data.frame(t(d))
                d
              }
              
              median.na <- function (x) {
                return(median(x,na.rm=TRUE))
              }
              
              #######      
              mMCPcounter.estimate <- function(exp, features = c("Gene.Symbol","ENSEMBL.ID","Probes")[1]){
                data("mMCPcounter_signatures",  envir=sys.frame(sys.nframe()))
                foundGenes <- intersect(mMCPcounter_signatures[,features],rownames(exp))
                if(length(foundGenes)==0){stop("No signature found in input row names. Please ensure the features argument is accurately set.")}
                absentSignatures <- setdiff(unique(mMCPcounter_signatures$Denomination),unique(mMCPcounter_signatures[mMCPcounter_signatures[,features]%in%rownames(exp),"Denomination"]))
                values <- reactiveValues()
                if(length(absentSignatures)==0)
                  
                  queryMagic <- function() {
                    print(paste("All genes were found for population(s)."))
                    
                    return("Data")
                  }
                
                if(length(absentSignatures)>0)
                  # {warning(paste("No genes were found for population(s): ",paste(absentSignatures,collapse = ", "),".",sep=""))}
                  
                  queryMagic <- function() {
                    print(paste("No genes were found for population(s): ",paste(absentSignatures,collapse = ", "),".",sep=""))
                    
                    return("Data")
                  }
                output$console <- renderPrint({
                  logText()
                  return(print(values[["log"]]))
                })
                
                logText <- reactive({
                  values[["log"]] <- capture.output(data <- queryMagic())
                })
                localSig <- mMCPcounter_signatures[mMCPcounter_signatures[,features] %in% foundGenes,]
                expAg <- exp[localSig[,features],]
                expAg <- cit.dfAggregate(expAg,localSig$Denomination,fAggreg = median.na)
                expAg <- expAg[c("T cells", "CD8 T cells", "NK cells", "B derived", "Memory B cells", "Monocytes / macrophages", "Monocytes", "Granulocytes", "Mast cells", "Eosinophils", "Neutrophils", "Basophils", "Vessels", "Lymphatics", "Endothelial cells", "Fibroblasts"),]
                expAg <- expAg[apply(expAg,1,function(x){sum(is.na(x))})<ncol(expAg),]
                return(expAg)
              }      
              ###################################################################      
              #####  Mouse MCP-counter
              m_mcp <- mMCPcounter.estimate(testdata)
              
              
              ############################################################     
              ### scale MCP result before plotting heatmap
              m_mcp_scaled <- t(scale(t(m_mcp)))
              
              ## Remove row with NA value
              m_mcp_scaled <- na.omit(m_mcp_scaled)
              #datamcp$m_mcp <- m_mcp
              size <-(length(colnames(m_mcp)))
              
              my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
              
              mcpplot<- pheatmap(m_mcp_scaled, show_rownames = T, angle_col = "90", cluster_cols = F, color = my_palette, cluster_rows = T,  fontsize_row = 14, fontsize_col = 14)
              
              
              output$mMCPRes <- DT::renderDataTable({m_mcp},extensions = 'FixedColumns', options = list(pageLength = 10, autoWidth =FALSE,scrollX = TRUE,
                                                                                                        scrollCollapse = TRUE,
                                                                                                        lengthMenu = c(10, 20, 50),
                                                                                                        columnDefs = list(list(targets = c(1:size), searchable = FALSE))),
                                                    
                                                    editable = FALSE)
              #End of MCP
              
              output$TableDownloadmMCP <- downloadHandler(  
                filename = function(){paste("mMCPRes",'.csv')},
                content = function(file){
                  write.csv(m_mcp, file, row.names = TRUE)
                }) 
              
              observe({  
                widthmcp <- reactive ({ input$plot_widthmcp})
                heightmcp <- reactive ({ input$plot_heightmcp}) 
                output$mMCPlot <- renderPlot(width = widthmcp, height = heightmcp,{
                  mcpplot
                })  
                
              })
              
              observe({
                chmcp <- input$optionmcp
                if (chmcp == "1"){ 
                  output$PlotDownloadmMCP <- downloadHandler(
                    filename <- function() {
                      paste("mMCP-Plot",".png", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::png(..., width = input$plot_widthmcp*4, height =input$plot_heightmcp*4 , res = 350, units = "px")
                      ggsave(file, plot =  mcpplot, device = device)
                       dev.off()
                    })
                } else if  (chmcp == "2"){ 
                  
                  output$PlotDownloadmMCP <- downloadHandler(
                    filename <- function() {
                      paste("mMCP-Plot",".svg", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthmcp*4/300), height =round(input$plot_heightmcp*4/300))
                      ggsave(file, plot = mcpplot , device = device)
                      #dev.off()
                    })
                }
              }) #End of Observe download
            } else if (chconter== "Human"){
              
              testdata<-   Collapsed_data[,-c(1,2)]
              
              #print(testdata)
              #####  Human MCP-counter
              mcp <- MCPcounter::MCPcounter.estimate(testdata, featuresType = 'HUGO_symbols')
              
              
              ### scale MCP result before plotting heatmap
              mcp_scaled <- t(scale(t(mcp)))
              
              ## Remove row with NA value
              mcp_scaled <- na.omit(mcp_scaled)
              #datamcp$mcp <- mcp
              size <-(length(colnames(mcp)))
              
              my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
              
              mcpplot<- pheatmap(mcp_scaled, show_rownames = T, angle_col = "90",cluster_cols = F, color = my_palette, cluster_rows = T, fontsize_row = 14, fontsize_col = 14)  
              
              output$mMCPRes <- DT::renderDataTable({mcp},extensions = 'FixedColumns', options = list(pageLength = 10, autoWidth =FALSE,scrollX = TRUE,
                                                                                                      scrollCollapse = TRUE,
                                                                                                      lengthMenu = c(10, 20, 50),
                                                                                                      columnDefs = list(list(targets = c(1:size), searchable = FALSE))),
                                                    
                                                    editable = FALSE)
              
              output$TableDownloadmMCP <- downloadHandler(  
                filename = function(){paste("MCPRes",'.csv')},
                content = function(file){
                  write.csv(mcp, file, row.names = TRUE)
                }) 
              
              observe({  
                widthmcp <- reactive ({ input$plot_widthmcp})
                heightmcp <- reactive ({ input$plot_heightmcp}) 
                output$mMCPlot <- renderPlot(width = widthmcp, height = heightmcp,{
                  mcpplot
                })  
                
              })
              
              observe({
                chmcp <- input$optionmcp
                if (chmcp == "1"){ 
                  output$PlotDownloadmMCP <- downloadHandler(
                    filename <- function() {
                      paste("MCP-Plot",".png", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::png(..., width = input$plot_widthmcp*4, height =input$plot_heightmcp*4 , res = 350, units = "px")
                      ggsave(file, plot =  mcpplot, device = device)
                      # dev.off()
                    })
                } else if  (chmcp == "2"){ 
                  
                  output$PlotDownloadmMCP <- downloadHandler(
                    filename <- function() {
                      paste("MCP-Plot",".svg", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthmcp*4/300), height =round(input$plot_heightmcp*4/300))
                      ggsave(file, plot = mcpplot , device = device)
                     # device <- function(..., width, height) grDevices::svg(..., width = 16, height = 10)
                     # ggsave(file, plot =  mcpplot, device = device)
                      #dev.off()
                    })
                }
              }) #End of Observe download
              
            } 
            
            
          })#End of MCP
          
          ################################################ 
          ######## Gene Enrichment ########################################################################################################      
          datass <- reactiveValues()
          observeEvent(input$ssGSEArAnalysis,{
            
            # Collapsed_data <- dat$Collapsed_data
            #  print( head(Collapsed_data))
            ########## Mouse data############################################   
            chspecies <- input$optionspecies
            if (chspecies == "Mouse"){   
              withProgress(message = "Preparing data..., please wait",
                           {
                             fge<-data$exp
                             head(fge)
                             rownames(fge)<-fge[,1]
                             fge<-fge[order( fge$ID, decreasing = FALSE), ]
                             Coll3 <- fge[,-c(1,2)]
                             coll3<-data.matrix(Coll3)
                             Coll<-coll3+1
                             normalized_df <- log2(Coll)
                             head(normalized_df)
                             rowGroup <- as.vector(fge$symbol)
                             rowID <- as.vector(rownames(normalized_df))
                             collapse.object <- collapseRows(datET=normalized_df, rowGroup=rowGroup, rowID=rowID, method = "MaxMean")
                             class(collapse.object)
                             names(collapse.object)
                             Collapsed_data <- data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
                             head(Collapsed_data)  
                             
                           })
              Mouse.Data <- Collapsed_data[,-c(1,2)]
              ####### Get the normalized data from the first part
              # Mouse.Data <- Collapsed_data
              # Mouse.Data <- Collapsed_data
              
              ##########################################ssGSEA for Mouse data
              ch <- input$option
              if (ch == "Hallmark"){
                withProgress(message = "Preparing Hallmark Table and Plot, please wait",
                             {
                               load("HallMark_Mouse.rdata")
                               datass$ssgsea.mouse <- gsva(as.matrix(Mouse.Data), x,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                               ssgsea.mouse.scaled <- t(scale(t(datass$ssgsea.mouse)))
                               my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
                               datass$hallp<- pheatmap(ssgsea.mouse.scaled, angle_col = "90", color=my_palette, fontsize_row = 8, fontsize_col = 10, 
                                                       border_color = NA, family = "Helvetica",treeheight_row = 8, legend = T, cluster_cols = F)
                             })
              } 
              
              else if (ch == "GeneOntology") {
                withProgress(message = "Preparing GeneOntology Table, please wait",
                             {
                               load("GO_Mouse.rdata")
                               datass$ssgsea.mouse <- gsva(as.matrix(Mouse.Data), x,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                             })
              }
              
              else if (ch=="USER_Geneset"){
                withProgress(message = "Preparing USER_Geneset Table, please wait",
                             {
                               file_to_read5 <- input$file5
                               if(is.null(file_to_read5)){return("")}
                               file6 <-  file_to_read5$datapath
                               e = new.env()
                               name <- load(file6, envir = e)
                               data <- e[[name]]
                               datass$ssgsea.mouse <- gsva(as.matrix(Mouse.Data), data,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                             })
              }
              
              
              datass$s <-(length(colnames(datass$ssgsea.mouse)))
              
              output$ssGSEA <- DT::renderDataTable({datass$ssgsea.mouse},
                                                   extensions = 'FixedColumns', 
                                                   options = list(pageLength = 20, autoWidth = FALSE,
                                                                  scrollX = TRUE,
                                                                  scrollCollapse = TRUE,
                                                                  lengthMenu = c(20, 50, 100),
                                                                  columnDefs = list(list(targets = c(1:datass$s), searchable = FALSE))),
                                                   editable = FALSE)
              
              output$TableDownloadssGSEA <- downloadHandler(
                filename = function(){paste("ssGSEA.Mouse",'.csv')},
                content = function(file){
                  write.csv(datass$ssgsea.mouse, file, sep = ",", quote = F)
                })
              
              
              
              observe({  
                widthhall <- reactive ({ input$plot_widthhall})
                heighthall <- reactive ({ input$plot_heighthall}) 
                output$hallPlot <- renderPlot(width = widthhall, height = heighthall,{
                  datass$hallp
                })  
                
              })
              
              observe({
                chhall <- input$optionhall
                if (chhall == "1"){ 
                  output$PlotDownloadhall <- downloadHandler(
                    filename <- function() {
                      paste("Hallmark-Plot",".png", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::png(..., width = input$plot_widthhall*4, height =input$plot_heighthall*4 , res = 350, units = "px")
                      ggsave(file, plot = datass$hallp, device = device)
                      # dev.off()
                    })
                } else if  (chhall== "2"){ 
                  
                  output$PlotDownloadhall <- downloadHandler(
                    filename <- function() {
                      paste("Hallmark-Plot",".svg", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthhall*4/300), height =round(input$plot_heighthall*4/300))
                      ggsave(file, plot = datass$hallp , device = device)
                     # device <- function(..., width, height) grDevices::svg(..., width = 16, height = 10)
                     # ggsave(file, plot = datass$hallp, device = device)
                     # dev.off()
                    })
                }
              })
              
            }
            
            
            #############################################################################################
            #############################################################################################
            else if (chspecies == "Human"){ 
              withProgress(message = "Preparing data..., please wait",
                           {
                             fge<-data$exp
                             head(fge)
                             rownames(fge)<-fge[,1]
                             fge<-fge[order( fge$ID, decreasing = FALSE), ]
                             Coll3 <- fge[,-c(1,2)]
                             coll3<-data.matrix(Coll3)
                             Coll<-coll3+1
                             normalized_df <- log2(Coll)
                             head(normalized_df)
                             rowGroup <- as.vector(fge$symbol)
                             rowID <- as.vector(rownames(normalized_df))
                             collapse.object <- collapseRows(datET=normalized_df, rowGroup=rowGroup, rowID=rowID, method = "MaxMean")
                             class(collapse.object)
                             names(collapse.object)
                             Collapsed_data <- data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
                             head(Collapsed_data)  
                           })
              
              Human.Data <- Collapsed_data[,-c(1,2)]
              
              ch <- input$option
              if (ch == "Hallmark"){
                withProgress(message = "Preparing Hallmark Table and Plot, please wait",
                             {
                               load("HallMark_Human.rdata")
                               datass$ssgsea.human <- gsva(as.matrix(Human.Data), x,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                               ssgsea.human.scaled <- t(scale(t(datass$ssgsea.human)))
                               my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
                               datass$hallp<- pheatmap(ssgsea.human.scaled, angle_col = "90", color=my_palette, fontsize_row = 8, fontsize_col = 10, 
                                                       border_color = NA, family = "Helvetica",treeheight_row = 8, legend = T, cluster_cols = F)
                             })
              } 
              
              
              
              else if (ch == "GeneOntology") {
                withProgress(message = "Preparing GeneOntology Table, please wait",
                             {
                               load("GO_Human.rdata")
                               datass$ssgsea.human <- gsva(as.matrix(Human.Data), x,
                                                           min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                             })
              }
              
              else if (ch=="USER_Geneset"){
                withProgress(message = "Preparing USER_Geneset Table, please wait",
                             {
                               file_to_read5 <- input$file5
                               if(is.null(file_to_read5)){return("")}
                               file6 <-  file_to_read5$datapath
                               e = new.env()
                               name <- load(file6, envir = e)
                               data <- e[[name]]
                               datass$ssgsea.human<- gsva(as.matrix(Human.Data), data,
                                                          min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea',parallel.sz=0)
                             })
              }
              
              
              datass$s <-(length(colnames(datass$ssgsea.human)))
              
              output$ssGSEA <- DT::renderDataTable({datass$ssgsea.human},
                                                   extensions = 'FixedColumns', 
                                                   options = list(pageLength = 20, autoWidth = FALSE,
                                                                  scrollX = TRUE,
                                                                  scrollCollapse = TRUE,
                                                                  lengthMenu = c(20, 50, 100),
                                                                  columnDefs = list(list(targets = c(1:datass$s), searchable = FALSE))),
                                                   editable = FALSE)
              
              
              output$TableDownloadssGSEA <- downloadHandler(
                filename = function(){paste("ssGSEA.Human",'.csv')},
                content = function(file){
                  write.csv(datass$ssgsea.human, file, sep = ",", quote = F)
                })
              
              
              
              observe({  
                widthhall <- reactive ({ input$plot_widthhall})
                heighthall <- reactive ({ input$plot_heighthall}) 
                output$hallPlot <- renderPlot(width = widthhall, height = heighthall,{
                  datass$hallp
                })  
                
              })
              
              
              
              observe({
                chhall <- input$optionhall
                if (chhall == "1"){ 
                  output$PlotDownloadhall <- downloadHandler(
                    filename <- function() {
                      paste("Hallmark-Plot",".png", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::png(..., width = input$plot_widthhall*4, height =input$plot_heighthall*4 , res = 350, units = "px")
                      ggsave(file, plot = datass$hallp, device = device)
                      # dev.off()
                    })
                } else if  (chhall== "2"){ 
                  
                  output$PlotDownloadhall <- downloadHandler(
                    filename <- function() {
                      paste("Hallmark-Plot",".svg", sep = "")
                    },
                    content <- function(file) {
                      device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthhall*4/300), height =round(input$plot_heighthall*4/300))
                      ggsave(file, plot = datass$hallp , device = device)
                     # device <- function(..., width, height) grDevices::svg(..., width = 16, height = 10)
                     # ggsave(file, plot = datass$hallp, device = device)
                     # dev.off()
                    })
                }
              })
            }  
          })#End Gene Enrichment 
          ######################################################################## 
          ##########################  
          ####GESA PLOT################################################################################################### 
          ###GESA PLOT######################## 
          datagesa <- reactiveValues()
          observeEvent(input$gseaplot,{
            withProgress(message = "Preparing list of Hallmark-Pathways, please wait...",
                         {
                           fgsea<-data$exp
                           head(fgsea)
                           rownames(fgsea)<-fgsea[,1]
                           fgsea<- fgsea[order(fgsea$ID, decreasing = FALSE), ]
                           
                           file_to_read2=input$file2
                           dsgsea <- read.delim(file_to_read2$datapath,sep=input$Sep1, header=T)
                           names(dsgsea)[1] <- "Types"
                           names(dsgsea)[2] <- "Samples"
                           names(dsgsea)[3] <- "Labels"
                           melfs1<-melt(fgsea)
                           
                           colnames(melfs1) <- c("ID", "symbol", "Samples","value")
                           
                           fss1<-join(melfs1, dsgsea, by = "Samples", type = "inner")
                           
                           g11 <- subset(fss1, Types==as.character(input$Group3comp)[1])
                           g11 = subset(g11, select = -c(Types,Labels) )
                           colnames(g11) <- c("ID", "symbol", "variable","value")
                           g11<-cast(g11)
                           g11 = subset(g11, select = -c(ID,symbol) )
                           Big11<-g11
                           
                           k11 <- subset(fss1, Types==as.character(input$Group4comp)[1])
                           k11 = subset(k11, select = -c(Types,Labels) )
                           colnames(k11) <- c("ID", "symbol", "variable","value")
                           k11<-cast(k11)
                           k11 = subset(k11, select = -c(ID,symbol) )
                           Big12<-k11
                           main1<-cbind(Big11,Big12)
                           main1<-cbind(fgsea[c(1,2)],main1)
                           
                           gr1 <- c(rep(as.character(input$Group3comp)[1],length(Big11) ), rep(as.character(input$Group4comp)[1],length(Big12)))
                           grss1 = factor(gr1, levels = unique(gr1))
                           
                           design <- model.matrix(~0+grss1)  
                           rownames(design) <- colnames(main1[,-c(1,2)])
                           
                           colnames(design) <- levels(grss1)
                           par(mar=c(0,0,0,0))
                           v <- voom(main1[,-c(1,2)], design=design, plot=F, normalize=as.character(input$Groupnorm))
                           fit <- lmFit(v, design)
                           
                           comp <- makeContrasts(contrasts= paste(as.character(input$Group3comp), as.character(input$Group4comp), sep="-"), levels=design )
                           
                           contrast.fit <- contrasts.fit(fit, contrasts=comp)
                           contrast.fit <- eBayes(contrast.fit)  
                           
                           
                           
                           contrast.fit <- na.omit(contrast.fit)
                           pdat2 <- data.frame(contrast.fit)
                           
                           pdat2 <- topTable(contrast.fit, number=Inf, adjust="BH", sort.by="none") 
                           
                           
                           pdat2$symbol<- fgsea$symbol
                           
                           
                           res<-pdat2
                           res2 <- res %>%
                             dplyr::select(symbol, logFC) %>%
                             #na.omit() %>%
                             distinct() %>%
                             group_by(symbol) %>%
                             summarize(logFC=mean(logFC))
                           #arrange(des(stat))
                           res2
                           
                           ranks <- deframe(res2)
                           ranks.2 <- sort(ranks, decreasing = T) 
                           
                           
                           chspeciesplot <- input$optionspeciesplot 
                           if (chspeciesplot == "Mouse"){                  
                             load("HallMark_Mouse.rdata")
                             pat <- x
                             set.seed(127)
                             fgseaRes <- fgsea(pathways=pat, stats=ranks.2, nperm=10000,minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) %>%
                               as_tibble() %>%
                               arrange(padj)
                             
                             filtered_pathway <- fgseaRes
                             
                             # print(filtered_pathway)
                             
                             # filtered_pathway <- subset(fgseaRes, pval < 0.05)
                             
                             filt_p <- as.vector(filtered_pathway$pathway)
                             
                             output$check5 <- renderUI({
                               selectInput ('Group5comp','List of Pathways', as.list(filtered_pathway$pathway),multiple = FALSE)
                             })
                             
                             
                             #print(filt_p[2])
                           }
                           else if  (chspeciesplot == "Human"){
                             load("HallMark_Human.rdata")
                             pat <- x
                             set.seed(127)
                             fgseaRes <- fgsea(pathways=pat, stats=ranks.2, nperm=10000,minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) %>%
                               as_tibble() %>%
                               arrange(padj)
                             
                             filtered_pathway <- fgseaRes
                             
                             
                             filt_p <- as.vector(filtered_pathway$pathway)
                             
                             output$check5 <- renderUI({
                               selectInput ('Group5comp','List of Pathways', as.list(filtered_pathway$pathway),multiple = FALSE)
                             })
                             
                           }
                           
                           
                           m<-filt_p[1]
                           
                           pannel<-filtered_pathway %>% filter_all(any_vars(. %in% c(m)))
                           pannel<-data.frame(rbind(pannel)) 
                           
                           
                           Pval<-pannel$pval
                           padj <-pannel$padj
                           ES<-pannel$ES
                           NES <-pannel$NES
                           
                           pp = formatC(Pval, digits = 4, format = "f")
                           dd = formatC( padj, digits = 4, format = "f")
                           ES = formatC(ES, digits = 4, format = "f")
                           NES = formatC( NES, digits = 4, format = "f")
                           
                           
                           pathway=pat[[m]]
                           stats=ranks.2 
                           gseaParam = 1
                           ticksSize = 0.3
                           
                           rnk <- rank(-stats)
                           ord <- order(rnk)
                           statsAdj <- stats[ord]
                           statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
                           statsAdj <- statsAdj/max(abs(statsAdj))
                           pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
                           pathway <- sort(pathway)
                           gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                                                   returnAllExtremes = TRUE)
                           bottoms <- gseaRes$bottoms
                           tops <- gseaRes$tops
                           n <- length(statsAdj)
                           xs <- as.vector(rbind(pathway - 1, pathway))
                           ys <- as.vector(rbind(bottoms, tops))
                           toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
                           diff <- (max(tops) - min(bottoms))/8
                           x = y = NULL
                           
                           pen<-ggplot(toPlot, aes(x = x, y = y))+geom_point(color = "green",  size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                                                                        linetype = "dashed") + geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
                             geom_hline(yintercept = 0,colour = "black") + geom_line(color = "green")+ theme_bw() + geom_segment(data = data.frame(x = pathway), 
                                                                                                                                 mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + theme(panel.border = element_blank(), 
                                                                                                                                                                                                                       panel.grid.minor = element_blank()) + labs(x = "Rank", y = "Enrichment Score (ES)",tag = paste(as.character(input$Group3comp)[1]))+
                             labs(title =m,subtitle = paste("P-value:",pp,", FDR:",dd,", ES:",ES,", NES:",NES))+ theme(plot.title = element_text(hjust = 0.5, face="bold")
                                                                                                                       ,plot.subtitle  = element_text(hjust = 0.5, face="bold", color="red"))+ theme(plot.tag.position = c(0.15, 0.02))                                                                             
                           
                           datagesa$gseaplot<-pen+labs(x=paste(as.character(input$Group4comp)[1])) + theme(axis.title.x = element_text(size=13,hjust=0.9))
                           
                           lead<-pannel$leadingEdge
                           test<- lead[[1]]
                           numgene<-length(test)
                           
                           datagesa$lead <- lead
                           
                           output$leading <- renderText({
                             isolate({
                               paste("Number of leadingEdge Genes in this pathway:", as.numeric(numgene), sep="")
                             })
                             
                           })       
                           
                           observe({
                             i<-as.character(input$Group5comp)[1]
                             pannel<-filtered_pathway %>% filter_all(any_vars(. %in% c(i)))
                             pannel<-data.frame(rbind(pannel)) 
                             
                             # observe({ 
                             if (!is.na(i)){
                               lead<-pannel$leadingEdge
                               test<- lead[[1]]
                               numgene<-length(test)
                               
                               output$leading <- renderText({
                                 isolate({
                                   paste("Number of leadingEdge Genes in this pathway:", as.numeric(numgene), sep="")
                                 })
                               })
                               datagesa$lead <- lead
                             }
                             # })
                             
                             Pval<-pannel$pval
                             padj <-pannel$padj
                             ES<-pannel$ES
                             NES <-pannel$NES
                             
                             pp = formatC(Pval, digits = 4, format = "f")
                             dd = formatC( padj, digits = 4, format = "f")
                             ES = formatC(ES, digits = 4, format = "f")
                             NES = formatC( NES, digits = 4, format = "f")
                             
                             pathway=pat[[i]]
                             stats=ranks.2 
                             gseaParam = 1
                             ticksSize = 0.3
                             
                             rnk <- rank(-stats)
                             ord <- order(rnk)
                             statsAdj <- stats[ord]
                             statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
                             statsAdj <- statsAdj/max(abs(statsAdj))
                             pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
                             pathway <- sort(pathway)
                             gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                                                     returnAllExtremes = TRUE)
                             bottoms <- gseaRes$bottoms
                             tops <- gseaRes$tops
                             n <- length(statsAdj)
                             xs <- as.vector(rbind(pathway - 1, pathway))
                             ys <- as.vector(rbind(bottoms, tops))
                             toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
                             diff <- (max(tops) - min(bottoms))/8
                             x = y = NULL
                             
                             
                             
                             pen<-ggplot(toPlot, aes(x = x, y = y))+geom_point(color = "green",  size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                                                                          linetype = "dashed") + geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
                               geom_hline(yintercept = 0,colour = "black") + geom_line(color = "green")+ theme_bw() + geom_segment(data = data.frame(x = pathway), 
                                                                                                                                   mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + theme(panel.border = element_blank(), 
                                                                                                                                                                                                                         panel.grid.minor = element_blank()) + labs(x = "Rank", y = "Enrichment Score (ES)",tag = paste(as.character(input$Group3comp)[1]))+
                               labs(title =i,subtitle = paste("P-value:",pp,", FDR:",dd,", ES:",ES,", NES:",NES))+ theme(plot.title = element_text(hjust = 0.5, face="bold")
                                                                                                                         ,plot.subtitle  = element_text(hjust = 0.5, face="bold", color="red"))+ theme(plot.tag.position = c(0.15, 0.02))                                                                             
                             
                             datagesa$gseaplot<-pen+labs(x=paste(as.character(input$Group4comp)[1])) + theme(axis.title.x = element_text(size=13,hjust=0.9))
                             
                             
                           })
                           
                         }) 
          })
          observeEvent(input$gegoplot,{
            withProgress(message = "Preparing list of GeneOntology-Pathways, please wait...",
                         {
                           fgsea<-data$exp
                           head(fgsea)
                           rownames(fgsea)<-fgsea[,1]
                           fgsea<- fgsea[order(fgsea$ID, decreasing = FALSE), ]
                           
                           file_to_read2=input$file2
                           dsgsea <- read.delim(file_to_read2$datapath,sep=input$Sep1, header=T)
                           names(dsgsea)[1] <- "Types"
                           names(dsgsea)[2] <- "Samples"
                           names(dsgsea)[3] <- "Labels"
                           
                           melfs1<-melt(fgsea)
                           
                           colnames(melfs1) <- c("ID", "symbol", "Samples","value")
                           
                           fss1<-join(melfs1, dsgsea, by = "Samples", type = "inner")
                           
                           g11 <- subset(fss1, Types==as.character(input$Group3comp)[1])
                           g11 = subset(g11, select = -c(Types,Labels) )
                           colnames(g11) <- c("ID", "symbol", "variable","value")
                           g11<-cast(g11)
                           g11 = subset(g11, select = -c(ID,symbol) )
                           Big11<-g11
                           
                           k11 <- subset(fss1, Types==as.character(input$Group4comp)[1])
                           k11 = subset(k11, select = -c(Types,Labels) )
                           colnames(k11) <- c("ID", "symbol", "variable","value")
                           k11<-cast(k11)
                           k11 = subset(k11, select = -c(ID,symbol) )
                           Big12<-k11
                           main1<-cbind(Big11,Big12)
                           main1<-cbind(fgsea[c(1,2)],main1)
                           
                           gr1 <- c(rep(as.character(input$Group3comp)[1],length(Big11) ), rep(as.character(input$Group4comp)[1],length(Big12)))
                           grss1 = factor(gr1, levels = unique(gr1))
                           
                           design <- model.matrix(~0+grss1)  
                           rownames(design) <- colnames(main1[,-c(1,2)])
                           
                           colnames(design) <- levels(grss1)
                           par(mar=c(0,0,0,0))
                           v <- voom(main1[,-c(1,2)], design=design, plot=F, normalize=as.character(input$Groupnorm))
                           fit <- lmFit(v, design)
                           
                           comp <- makeContrasts(contrasts= paste(as.character(input$Group3comp), as.character(input$Group4comp), sep="-"), levels=design )
                           
                           contrast.fit <- contrasts.fit(fit, contrasts=comp)
                           contrast.fit <- eBayes(contrast.fit)  
                           
                           
                           
                           contrast.fit <- na.omit(contrast.fit)
                           pdat2 <- data.frame(contrast.fit)
                           
                           pdat2 <- topTable(contrast.fit, number=Inf, adjust="BH", sort.by="none") 
                           
                           
                           pdat2$symbol<- fgsea$symbol
                           
                           
                           res<-pdat2
                           res2 <- res %>%
                             dplyr::select(symbol, logFC) %>%
                             #na.omit() %>%
                             distinct() %>%
                             group_by(symbol) %>%
                             summarize(logFC=mean(logFC))
                           #arrange(des(stat))
                           res2
                           
                           ranks <- deframe(res2)
                           ranks.2 <- sort(ranks, decreasing = T) 
                           
                           chspeciesplot <- input$optionspeciesplot 
                           if (chspeciesplot == "Mouse"){
                             load("GO_Mouse.rdata")
                             pat <- x
                             set.seed(127)
                             fgseaRes <- fgsea(pathways=pat, stats=ranks.2, nperm=10000,minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) %>%
                               as_tibble() %>%
                               arrange(padj)
                             
                             filtered_pathway<- fgseaRes
                             
                             filtered_pathway<- filtered_pathway[order( filtered_pathway$ES, decreasing = TRUE), ]
                             filtered_pathway<- filtered_pathway[1:50,]
                             
                             
                             
                             filt_p <- as.vector(filtered_pathway$pathway)
                             
                             output$check5 <- renderUI({
                               selectInput ('Group6comp','List of Pathways', as.list(filtered_pathway$pathway),multiple = FALSE)
                             })
                             
                           }
                           else if  (chspeciesplot == "Human"){
                             load("GO_Human.rdata")
                             pat <- x
                             set.seed(127)
                             fgseaRes <- fgsea(pathways=pat, stats=ranks.2, nperm=10000,minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) %>%
                               as_tibble() %>%
                               arrange(padj)
                             
                             filtered_pathway<- fgseaRes
                             
                             filtered_pathway<- filtered_pathway[order( filtered_pathway$ES, decreasing = TRUE), ]
                             filtered_pathway<- filtered_pathway[1:50,]
                             
                             
                             
                             filt_p <- as.vector(filtered_pathway$pathway)
                             
                             output$check5 <- renderUI({
                               selectInput ('Group6comp','List of Pathways', as.list(filtered_pathway$pathway),multiple = FALSE)
                             })
                             
                           }
                           m<-filt_p[1]
                           
                           pannel2<-filtered_pathway %>% filter_all(any_vars(. %in% c(m)))
                           pannel2<-data.frame(rbind(pannel2)) 
                           
                           
                           Pval<-pannel2$pval
                           padj <-pannel2$padj
                           ES<-pannel2$ES
                           NES <-pannel2$NES
                           
                           pp = formatC(Pval, digits = 4, format = "f")
                           dd = formatC( padj, digits = 4, format = "f")
                           ES = formatC(ES, digits = 4, format = "f")
                           NES = formatC( NES, digits = 4, format = "f")
                           
                           
                           pathway=pat[[m]]
                           stats=ranks.2 
                           gseaParam = 1
                           ticksSize = 0.3
                           
                           rnk <- rank(-stats)
                           ord <- order(rnk)
                           statsAdj <- stats[ord]
                           statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
                           statsAdj <- statsAdj/max(abs(statsAdj))
                           pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
                           pathway <- sort(pathway)
                           gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                                                   returnAllExtremes = TRUE)
                           bottoms <- gseaRes$bottoms
                           tops <- gseaRes$tops
                           n <- length(statsAdj)
                           xs <- as.vector(rbind(pathway - 1, pathway))
                           ys <- as.vector(rbind(bottoms, tops))
                           toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
                           diff <- (max(tops) - min(bottoms))/8
                           x = y = NULL
                           
                           pen<-ggplot(toPlot, aes(x = x, y = y))+geom_point(color = "green",  size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                                                                        linetype = "dashed") + geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
                             geom_hline(yintercept = 0,colour = "black") + geom_line(color = "green")+ theme_bw() + geom_segment(data = data.frame(x = pathway), 
                                                                                                                                 mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + theme(panel.border = element_blank(), 
                                                                                                                                                                                                                       panel.grid.minor = element_blank()) + labs(x = "Rank", y = "Enrichment Score (ES)",tag = paste(as.character(input$Group3comp)[1]))+
                             labs(title =m,subtitle = paste("P-value:",pp,", FDR:",dd,", ES:",ES,", NES:",NES))+ theme(plot.title = element_text(hjust = 0.5, face="bold")
                                                                                                                       ,plot.subtitle  = element_text(hjust = 0.5, face="bold", color="red"))+ theme(plot.tag.position = c(0.15, 0.02))                                                                             
                           
                           datagesa$gseaplot<-pen+labs(x=paste(as.character(input$Group4comp)[1])) + theme(axis.title.x = element_text(size=13,hjust=0.9))
                           
                           lead2<-pannel2$leadingEdge
                           test2<- lead2[[1]]
                           numgene2<-length(test2)
                           
                           datagesa$lead <- lead2
                           
                           output$leading <- renderText({
                             isolate({
                               paste("Number of leadingEdge Genes in this pathway:", as.numeric(numgene2), sep="")
                             })
                             
                           })       
                           
                           observe({
                             i<-as.character(input$Group6comp)[1]
                             pannel2<-filtered_pathway %>% filter_all(any_vars(. %in% c(i)))
                             pannel2<-data.frame(rbind(pannel2)) 
                             # print(i)
                             # observe({ 
                             if (!is.na(i)){
                               lead2<-pannel2$leadingEdge
                               test2<- lead2[[1]]
                               numgene2<-length(test2)
                               
                               output$leading <- renderText({
                                 isolate({
                                   paste("Number of leadingEdge Genes in this pathway:", as.numeric(numgene2), sep="")
                                 })
                               })
                               datagesa$lead <- lead2
                             }
                             # })
                             
                             Pval<-pannel2$pval
                             padj <-pannel2$padj
                             ES<-pannel2$ES
                             NES <-pannel2$NES
                             
                             pp = formatC(Pval, digits = 4, format = "f")
                             dd = formatC( padj, digits = 4, format = "f")
                             ES = formatC(ES, digits = 4, format = "f")
                             NES = formatC( NES, digits = 4, format = "f")
                             
                             pathway=pat[[i]]
                             stats=ranks.2 
                             gseaParam = 1
                             ticksSize = 0.3
                             
                             rnk <- rank(-stats)
                             ord <- order(rnk)
                             statsAdj <- stats[ord]
                             statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
                             statsAdj <- statsAdj/max(abs(statsAdj))
                             pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
                             pathway <- sort(pathway)
                             gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                                                     returnAllExtremes = TRUE)
                             bottoms <- gseaRes$bottoms
                             tops <- gseaRes$tops
                             n <- length(statsAdj)
                             xs <- as.vector(rbind(pathway - 1, pathway))
                             ys <- as.vector(rbind(bottoms, tops))
                             toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
                             diff <- (max(tops) - min(bottoms))/8
                             x = y = NULL
                             
                             
                             
                             pen<-ggplot(toPlot, aes(x = x, y = y))+geom_point(color = "green",  size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                                                                          linetype = "dashed") + geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
                               geom_hline(yintercept = 0,colour = "black") + geom_line(color = "green")+ theme_bw() + geom_segment(data = data.frame(x = pathway), 
                                                                                                                                   mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + theme(panel.border = element_blank(), 
                                                                                                                                                                                                                         panel.grid.minor = element_blank()) + labs(x = "Rank", y = "Enrichment Score (ES)",tag = paste(as.character(input$Group3comp)[1]))+
                               labs(title =i,subtitle = paste("P-value:",pp,", FDR:",dd,", ES:",ES,", NES:",NES))+ theme(plot.title = element_text(hjust = 0.5, face="bold")
                                                                                                                         ,plot.subtitle  = element_text(hjust = 0.5, face="bold", color="red"))+ theme(plot.tag.position = c(0.15, 0.02))                                                                             
                             
                             datagesa$gseaplot<-pen+labs(x=paste(as.character(input$Group4comp)[1])) + theme(axis.title.x = element_text(size=13,hjust=0.9))
                             
                             
                           })     
                           
                           
                         })
          })
          
          
          output$listofgenes<- downloadHandler(
            filename = function(){paste("LeadingEdgeGenes",'.csv')},
            content = function(file){
              write.table(as.data.frame(datagesa$lead), file, quote=F ,row.names=FALSE,col.names=FALSE)
            })
          
          
          observe({  
            widthgsea <- reactive ({ input$plot_widthgsea })
            heightgsea <- reactive ({ input$plot_heightgsea}) 
            output$gseaPlot <- renderPlot(width = widthgsea, height =  heightgsea,{  
              gseaplot <-datagesa$gseaplot
              vals$gseaplot <-gseaplot
              print(gseaplot)
            })
          }) 
            observe({  
            choptiongsea<- input$optiongsea
            
            if (choptiongsea == "1"){
              output$PlotDownloadgsea <- downloadHandler(
                filename = function(){paste("GSEAPlot",'.png',sep='')},
                content = function(file){
                  device <- function(..., width, height) grDevices::png(..., width =input$plot_widthgsea*4, height =input$plot_heightgsea*4, res = 350, units = "px")
                  ggsave(file, plot = ( vals$gseaplot) , device = device)
                })
            }
            else if (choptiongsea == "2"){
              output$PlotDownloadgsea <- downloadHandler(
                filename = function(){paste("GSEAPlot",'.svg',sep='')},
                content = function(file) {
                  device <- function(..., width, height) grDevices::svg(..., width =round(input$plot_widthgsea*4/300), height =round(input$plot_heightgsea*4/300))
                  ggsave(file, plot = vals$gseaplot , device = device)
               
                  
                })
            }  
            
          })#observe
          
          #output$gseaPlot <- renderPlot({
          #  datagesa$gseaplot
          # })   
          
          
          
          
       })  #end of observe
  #############################################################################################################  
        }#End of Normalize File
        
      })#End of Observe for datafile choice  
      
      
      ##############################################################################################################
      ########### Update count #########
      # Reactively update the client.
      output$count <- renderText({
        vals$count
      })
      
      # When a session ends, decrement the counter.
      session$onSessionEnded(function(){
        isolate(vals$count <- vals$count - 1)
      })  
      
      ###########################################################################################################################
    }
    )##END shiny server
    
###########################################################################################################################
    
    
    
    
    
    
  