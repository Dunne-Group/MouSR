## This App is used for interactivaly visualizing RNA-seq analysis 
#getwd();
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinycustomloader)
library(ggplot2)
library(CMScaller)
library(tibble)
library(DESeq2)
library(plyr)
library(pheatmap)
library(biomaRt)
library(heatmaply)
library(reshape)
library(plotly)
library(WGCNA)
library(lattice)
library(RColorBrewer)
library(GSVA)
library(rlist)
library(msigdbr)
library(tidyverse)
library(mMCPcounter)
library(magrittr)
library(dplyr)
library(ggrepel)
library(shinycssloaders)
library(readxl)
library(DT)
library(colourpicker)
library(fgsea)
library(enrichplot)
library(MCPcounter)
library(limma)
library(shinyalert)

header <- dashboardHeader(
  title = "Molecular Subtyping Resource ",tags$li(a(href='http://www.beatson.gla.ac.uk/',target="_blank",img(src = "logo9.png",align="right", height = 102, width = 420),img(src = "logo23.png",align="left", height = 102, width = 620),
                                                  style = "padding-top:1px; padding-bottom:1px;padding-right:0px"),
                                                  class ="dropdown")
)
#########################################################################################################################################700####################################################700####
sidebar <- dashboardSidebar(
                             sidebarMenu(id="menu1",
                                         menuItem("Menu Items", icon = icon("fas fa-feather-alt")
                                         ),
                                         menuItem("Introduction", icon = icon("user"),tabName = "intro"
                                         ),
                                         menuItem("Data Input", icon = icon("file"),
                                                  menuSubItem("Upload Your Data", tabName = "dataupload")
                                         ),
                                         menuItem("Data Analysis", icon = icon("area-chart"), 
                                                  menuItem("PCA & MDS ", icon = icon("fas fa-save"),
                                                           menuSubItem("PCA 2D Plot", tabName="PCA"),
                                                           menuSubItem("PCA 3D Plot", tabName="3DPCA"),
                                                           menuSubItem("MDS 2D Plot", tabName="MDS")),
                                                  
                                                  menuItem("Differential Gene Expressions",icon = icon("fas fa-microscope"),
                                                           menuSubItem("Heatmap", tabName="MSG-heatmap"),
                                                           menuSubItem("Heatmap Selected Genes", tabName="heatmap"),
                                                           menuSubItem("Gene Expression Levels", tabName="Gene-plot"),
                                                           menuSubItem("Volcano", tabName="Volcano")
                                                           
                                                  ),
                                                  menuItem("GSEA",icon = icon("fas fa-dna"),
                                                           menuSubItem("GSEA Plot",tabName="GSEA-plot"),
                                                           menuSubItem("ssGSEA Table",tabName="ssGSEA"),
                                                           menuSubItem("Hallmark Plot ", tabName="Hallmarkplot")
                             
                                                  ),
                                                  
                                                  menuItem("mMCP/MCP-counter",icon = icon("fas fa-paste"),
                                                           menuSubItem("mMCP/MCP Table", tabName="mMCP"),
                                                           menuSubItem("mMCP/MCP Plot ", tabName="mMCPplot")
                                                  )
                                         ),
                                         
                                         menuItem("Molecular Subtyping(v1.0)", icon = icon("fad fa-layer-group", lib = "font-awesome"), 
                                                  menuSubItem("CMS-Classifier(Beta version)", tabName="CMS_T"),
                                                  menuSubItem("CRIS-Classifier(Beta version)", tabName="CRIS_T")
                                         ),
                                         
                                        
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         div(style="width:100px;padding-left:50px;", h5("Version 1.0")),
                                         div(h6(" ShinyApp created by Baharak Ahmaderaghi"))
                                        
                             )
) #endof sidebar <- dashboardSidebar
#########################################################################################################################################
body <- dashboardBody(
  tabItems(
    ## Introduction tab panel################################################################################################################
    tabItem(tabName="intro",
            br(),
            h2(strong(span("INTRODUCTION:"))),
            br(),
            p(strong(span("Welcome to the Molecular Subtyping Resource (",HTML(paste("MouS", tags$sup("R"), sep = "")),"):",style = "color:green"), style = "font-family: 'times'; font-si16pt")),
            br(),
            p("Number of user(s) currently connected to this App",div(style="width:200px;padding-right:100px;", verbatimTextOutput("count"))),
            p(strong("Please Note:"),"Chrome, Edge and then Firefox are the recommended browsers for the best App experience."),
            br(),
            p(strong(span("Introduction-Video"),a(strong("Watch it here"),href="https://youtu.be/_fFXYMp1o7g", target="_blank"),style = "color:red")),
            br(),
            p(" This Molecular Subtyping Resource has been developed within the Molecular Pathology Research Group",
              a("Dunne-Lab", href="https://dunne-lab.com", target="_blank"),
              "as part of the CRUK-funded",strong("ACRCelerate"),"colorectal cancer pre-clinical accelerator programme."),
            p("The ACRCelerate project will bring together an european-wide consortium of basic and clinical scientists at the forefront of CRC research to interrogate a suite of state-of-the-art preclinical models. The overarching aim of ACRCelerate is to generate robust and reproducible preclinical data to de-risk future clinical trials via patient stratification. Specifically, the various models will be categorised into subtypes based on their gene activity so that treatments can be aimed at particular subgroups of patients. In doing so, the network will be able to accelerate the next generation of stratified trials for CRC through accurate disease subtype positioning. Using this preclinical testing platform, the consortium aims to accelerate the identification and development of new drugs for CRC."),
            br(), 
            p("This portal accepts uploads of either; text, csv files. The format of your Transcriptional Data Matrix",strong("RNA-seq-Read Count(Both Decimal and Integer numbers)"),"should be similar to this table, with a size limit of up to 30 MB, length no restriction:
            If your Transcriptional Data Matrix file is not in the same orientation to the table below, you will need to transpose it. A transposing tool can be accessed via this link:",
              a("Click here", href="https://www.convertcsv.com/transpose-csv.htm", target="_blank")
              
            ),
            br(), 
            br(),
            HTML("
          <html>
            <head>
            <style>
            table, th, td {
              border: 1px solid black;
              border-collapse: collapse;
            }
          th, td {
            padding: 5px;
            text-align: left;
          }
          </style>
            </head>
            <body>
            <TABLE>
            <CAPTION><EM>Transcriptional Data Matrix Table</EM></CAPTION>
             <tr>
                   <th>ID</th>
                   <th>symbol</th>
                    <th>KPNt1</th>
                    <th>KPNt2</th>
                    <th>KPNt3</th>
                    <th>KPt1</th>
                    <th>KPt2</th>
                     <th>KPt3</th>
                    <th>APNt1</th>
                    <th>APNt2</th>
                    <th>APNt3</th>
                    <th>......</th>
                   
                    
               </tr>
               <tr>
    <td>ENSMUSG00000000001</td>
    <td>Gnai3</td>
    <td>2792</td>
    <td>867</td>
    <td>3074</td>
    <td>2096</td>
    <td>4537</td>
    <td>5269</td>
    <td>5635</td>
    <td>6359</td>
    <td>9770</td>
    <td>......</td>
  </tr>
  
                 <tr>
    <td>ENSMUSG00000000028</td>
    <td>Cdc45</td>
    <td>510</td>
    <td>143</td>
    <td>131</td>
    <td>178</td>
    <td>289</td>
    <td>316</td>
    <td>479</td>
    <td>608</td>
    <td>1065</td>
    <td>......</td>
  </tr>
                 <tr>
    <td>ENSMUSG00000000031</td>
    <td>H19</td>
    <td>212</td>
    <td>128</td>
    <td>796</td>
    <td>300</td>
    <td>634</td>
    <td>598</td>
    <td>100495</td>
    <td>2082</td>
    <td>18134</td>
    <td>......</td>
  </tr>
  <tr>
    <td>ENSMUSG00000000037</td>
    <td>Scml2</td>
    <td>8</td>
    <td>0</td>
    <td>0</td>
    <td>11</td>
    <td>1</td>
    <td>2</td>
    <td>27</td>
    <td>122</td>
    <td>135</td>
    <td>......</td>
  </tr>
    <tr>
    <td>ENSMUSG00000000049</td>
    <td>Apoh</td>
    <td>10</td>
    <td>26</td>
    <td>1</td>
    <td>25</td>
    <td>27</td>
    <td>15</td>
    <td>0</td>
    <td>23</td>
    <td>40</td>
    <td>......</td>
  </tr>
  <tr>
    <td>ENSMUSG00000000056</td>
    <td>Narf</td>
    <td>925</td>
    <td>176</td>
    <td>730</td>
    <td>569</td>
    <td>850</td>
    <td>830</td>
    <td>1790</td>
    <td>1632</td>
    <td>2115</td>
    <td>......</td>
  </tr>
                 <tr>
    <td>ENSMUSG00000000058</td>
    <td>Cav2</td>
    <td>2195</td>
    <td>568</td>
    <td>460</td>
    <td>1480</td>
    <td>2586</td>
    <td>1211</td>
    <td>303</td>
    <td>457</td>
    <td>605</td>
    <td>...</td>
  </tr>
  <tr>
    <td>ENSMUSG00000000078</td>
    <td>Klf6</td>
    <td>3651</td>
    <td>838</td>
    <td>2349</td>
    <td>4263</td>
    <td>5373</td>
    <td>5483</td>
    <td>3156</td>
    <td>5324</td>
    <td>6919</td>
    <td>......</td>
  </tr>
  <tr>
    <td>ENSMUSG00000000085</td>
    <td>Scmh1</td>
    <td>764</td>
    <td>301</td>
    <td>491</td>
    <td>1175</td>
    <td>1099</td>
    <td>601</td>
    <td>1420</td>
    <td>814</td>
    <td>1634</td>
    <td>......</td>
  </tr>
  <tr>
    <td>ENSMUSG00000000088</td>
    <td>Cox5a</td>
    <td>2127</td>
    <td>693</td>
    <td>1716</td>
    <td>3127</td>
    <td>2506</td>
    <td>3221</td>
    <td>1333</td>
    <td>3077</td>
    <td>4507</td>
    <td>......</td>
  </tr>
  
  <td>ENSMUSG00000000093</td>
    <td>Tbx2</td>
    <td>127</td>
    <td>31</td>
    <td>86</td>
    <td>250</td>
    <td>269</td>
    <td>192</td>
    <td>106</td>
    <td>206</td>
    <td>344</td>
    <td>......</td>
  </tr>
  <td>ENSMUSG00000000094</td>
    <td>Tbx4</td>
    <td>5</td>
    <td>0</td>
    <td>24</td>
    <td>20</td>
    <td>42</td>
    <td>14</td>
    <td>2028</td>
    <td>143</td>
    <td>209</td>
    <td>......</td>
  </tr>
<td>ENSMUSG00000000120</td>
    <td>Ngfr</td>
    <td>28</td>
    <td>9</td>
    <td>108</td>
    <td>47</td>
    <td>96</td>
    <td>45</td>
    <td>571</td>
    <td>71</td>
    <td>76</td>
    <td>......</td>
  </tr>
<td>ENSMUSG00000000125</td>
    <td>Wnt3</td>
    <td>0</td>
    <td>0</td>
    <td>12</td>
    <td>1</td>
    <td>2</td>
    <td>6</td>
    <td>20</td>
    <td>99</td>
    <td>100</td>
    <td>......</td>
  </tr>
<td>ENSMUSG00000000126</td>
    <td>Wnt9a</td>
    <td>89</td>
    <td>5</td>
    <td>18</td>
    <td>190</td>
    <td>171</td>
    <td>50</td>
    <td>24</td>
    <td>25</td>
    <td>22</td>
    <td>......</td>
  </tr>
  

  
 
   <tr>
    <td>......</td>
    <td>......</td>
    <td>......</td>
    <td>......</td>
    <td>......</td>
    <td>.......</td>
     <td>......</td>
    <td>......</td>
    <td>......</td>
    <td>......</td>
    <td>......</td>
    <td>......</td>
  </tr>
            </TABLE>
            </body>
            "),
            br(),
            p("The format of your Sample Information and Labels Table should be similar to this,
            If your types are two words, such as Early tumor, you should write them as ",strong("Early_tumor."),"Also please note the order of your Sample input should be exact same match with your Transcriptional Data Matrix columns."),
            HTML("
          <html>
            <head>
            <style>
            table, th, td {
              border: 1px solid black;
              border-collapse: collapse;
            }
          th, td {
            padding: 5px;
            text-align: left;
          }
          </style>
            </head>
            <body>
            <TABLE>
            <CAPTION><EM>Sample Information and Labels Table</EM></CAPTION>
             <tr>
                   <th>Group</th>
                    <th>Samples	</th>
                    <th>Group-ID</th>
                    
               </tr>
               <tr>
    <td>KPN_tumour</td>
    <td>KPNt1</td>
    <td>KPNt1 </td>

  </tr>
                <tr>
    <td>KPN_tumour</td>
    <td>KPNt2</td>
    <td> KPNt2</td>

  </tr>
                <tr>
    <td>KPN_tumour</td>
    <td>KPNt3</td>
    <td>KPNt3 </td>

  </tr>
  <tr>
    <td>KP_tumour</td>
    <td>KPt1</td>
    <td>KPt1 </td>

  </tr>
  <tr>
    <td>KP_tumour</td>
    <td>KPt2</td>
    <td>KPt2 </td>

  </tr>
  <tr>
    <td>KP_tumour</td>
    <td>KPt3</td>
    <td>KPt3 </td>

  </tr>
   <tr>
     <td>APN_tumour</td>
    <td>APNt1</td>
    <td>APNt1 </td>

  </tr>
 <td>APN_tumour</td>
    <td>APNt2</td>
    <td>APNt2</td>

  </tr>
 <td>APN_tumour</td>
    <td>APNt3</td>
    <td>APNt3 </td>

  </tr>
 <tr>
  <td>......</td>
    <td>......</td>
    <td>...... </td>

  </tr>
   
            </TABLE>
            </body>
    
                 "),
            br(),
            br(),
            
            h5 ('Source of inspirations:'),
            br(),
            a("DEApp",href="https://yanli.shinyapps.io/DEApp/", target="_blank"),
            
            
            br()
            
    ), # End introduction tab panel
    #########################################################################################################################################      
    ## Upload data
    tabItem(tabName="dataupload",
            fluidRow(
              box(title = "Input Data: Please Upload your Data",
                  solidHeader = T, status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,
                  p(strong(span("Data Input-Video"),a(strong("Watch it here"),href="https://youtu.be/Y05M5CE75I4", target="_blank"),style = "color:red")),
                  radioButtons(inputId="Fileformat",label="Select one Option based on the Format of your Transcriptional Data Matrix ",choices=c("None selected" = "", " Read count (Integer Numbers)"="1"," Read count (Decimal Numbers) "="2")),
                   p(strong("Please note")," In Read count (Decimal Numbers) option, the Entire Differential Gene Expressions section is not available in the current version of the App."
                  ,style="color:red; text-align:left ;"),
               
                  
                   fluidRow(
                    ###########################################################                
                    box(title = "Input 1: Transcriptional Data Matrix",
                        solidHeader = T, status = "info",
                        width = 6,
                        fileInput("file1","Upload Your Transcriptional Data Matrix File here:",
                                  multiple = TRUE,
                                  accept=c('text/tab-separated-values',
                                           'text/csv',
                                           'text/comma-separated-values',
                                           'text/tab-separated-values',
                                           '.txt',
                                           '.csv',
                                           '.tsv')),
                        radioButtons(inputId="Sep1",label="Separator(file)",choices=c(Comma=",",Semicolon=";", Tab="\t"), selected='\t')
                        ,helpText("The Examplar file of 'Read count (Integer Numbers)' is available", a("here", href="https://drive.google.com/file/d/1XsSJVfUixYuBdcFdL9gxo4ODeDXqO0-R/view?usp=sharing", target="_blank")
                                 ,style="color:green; text-align:center; font-size: 16px;")
                    ),

                    #####################################################               
                    box(title = "Input 2: Sample Information and Labels",
                        solidHeader = T, status = "info",
                        width = 6,
                        fileInput("file2","Upload Your Sample Information and Labels File here:",
                                  multiple = TRUE,
                                  accept=c('text/tab-separated-values',
                                           'text/csv',
                                           'text/comma-separated-values',
                                           'text/tab-separated-values',
                                           '.txt',
                                           '.csv',
                                           '.tsv')),
                       
                       radioButtons(inputId="Sep",label="Separator(file)",choices=c(Comma=",",Semicolon=";", Tab="\t"), selected='\t')
                       ,helpText("The Examplar file of 'Sample Information and Labels'is available", a("here", href=" https://drive.google.com/file/d/11psEllJudRFBSybJOX3pB43VMboEL92z/view?usp=sharing", target="_blank")
                                 ,style="color:green; text-align:center;font-size: 16px;")
                    ),
                    
                    ######################################################
                    uiOutput("a1"), 
                    tags$style("button#stop{margin-left:auto;margin-right:auto;display:block; padding: 10px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                    tags$style("button#dataSubmit {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                    tags$style("button#dataSubmit2 {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                    helpText("Clicking Check Input Files again before the upload is complete will reset the process"
                             ,style="color:red; text-align:center ;"),
                    
                  ),
                 
                verbatimTextOutput("console2"),
                 tags$head(tags$style(HTML("
                 #console2 { 
                  color:red;
                 margin-left:auto;
                 margin-right:auto;
                 font-size: 16px;
                 font-weight:bold;
                 text-align:center;
                 }
                "))),
                
                 
                  # end fluidRow(
                  # shinycssloaders::withSpinner(
                  ####################################################
                  fluidRow(
                    box(title = "Input Details",
                        solidHeader = T, status = "primary",
                        width = 12, 
                        fluidRow(
                          column(6,
                                 h4(textOutput("Title")),
                                 verbatimTextOutput("Info")
                          ),
                          column(6,
                                 h4(textOutput("t1")),
                                 verbatimTextOutput("t2")
                          )
                        )
                    )
                  ),
                  ##################################################
                  fluidRow(
                    box(title = "Input Information Summary",
                        solidHeader = T,status = "primary",
                        width = 6,
                        tableOutput("Summary"),
                        tags$style("#Summary table {border: 1px solid black; align: center; margin:auto;}","#Summary th {border: 1px solid black;}","#Summary td {border: 1px solid black;}")
                    ),
                    
                    box(title = "Everything looks correct?",
                            solidHeader = T,status = "primary",
                        width = 6,
                       tags$style("#error{font-size:14px;
                                   color:red;
                                  font-weight:bold; }"),
                      textOutput("error"),
                      br(),
                       uiOutput("b1"), 
                       tags$style("button#wait{margin-left:auto;margin-right:auto;display:block; padding: 10px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                      uiOutput("b"),
                       tags$style("button#dataPCA {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 10px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                       tags$style("button#dataPCA2 {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 10px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                    )
                    
                    
                  ) #end of fluidRow(
                  ################################################
              )#end of  box(title 
            )#end of fluidRow
            
    ), ##End of First tab content for data input
    #########################################################################################################################################
    #################################################################################################################
    # PCA Tab
    tabItem(tabName="3DPCA",
            fluidRow(          
              box(title = "PCA 3D Analysis ",
                  solidHeader = T,status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,height = 800,
                  column(3,
                         p(strong(span("PCA 3D-Video"),a(strong("Watch it here"),href="https://youtu.be/GI-KpZEi74M", target="_blank"),style = "color:red")),
                         radioButtons(inputId="option3d",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                         p(strong("Plotly Modebar control:")," From left to right"),
                         p( "1)Download plot as a selected format."),
                         p( "2) Zoom on selection. Click and hold your mouse button to select a region to zoom on."),
                         p( "3)Image panning."),
                         p( "4)Rotates the plot around its middle point in three-dimensional space.") ,
                         p( "5)Rotates the plot around its middle point while constraining one axis slightly."),
                         p( "6)Reset axis to the original axis ranges."),
                         p( "7)Reset to the last saved axis ranges."),
                         p( "8)Show closest data on hover -- this function is always on."),
                         br(),
                         
                  ),
                  column(8,                   
                         helpText("Group color can be changed using PCA2D Pick a color option ",style="color:red; text-align:center ;"),
                         plotlyOutput("PCAPlot")
                  )
              ))),
    ###########################2D PLOT
    tabItem(tabName="PCA",
            fluidRow(  
              box(title = "PCA 2D Analysis ",
                  solidHeader = T, status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,height = "2800",
                  box(title = " Different Filtering Options",
                      width = 12,
                      solidHeader = T,
                      status = "info",
                      fluidRow(
                        column(4,
                               p(strong(span("PCA 2D-Video"),a(strong("Watch it here"),href="https://youtu.be/gLJZ80muYr0", target="_blank"),style = "color:red")),
                               radioButtons(inputId="optionp",label="Select one option",choices=c("PC1&PC2"="1","PC1&PC3"="2","PC2&PC3"="3"),selected = "1"),
                               radioButtons(inputId="optionlab",label="Select label option",choices=c(ON="1",OFF="2"), selected='1'),
                               p(strong("PCA")," :A dimensionality-reduction method to assess quality and clustering characteristics of data sets. We used the prcomp function in the R stats package to perform PCA analysis."), 
                               a(" PCA References",href="https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html", target="_blank")
                        ),
                        column(3,
                               radioButtons(inputId="optionchs",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                               sliderInput("label_pca", "Size of the labels", 0, 10, 4),  
                               sliderInput("point_pca", "Size of the points", 0, 10, 3),  
                               numericInput("plot_heightpca", "Plot_Height (# pixels): ",min = 300, max = 600, value = 500),
                               numericInput("plot_widthpca", "Plot_Width (# pixels):", min = 600, max = 1000,value = 900)
                               
                        ),
                        
                        column(4,
                               p(strong("Put plot option on PC1&PC2 first:"),style="color:red;"),
                               uiOutput("color1"),
                               
                               uiOutput('platter')
                               
                        )
                      )
                  ),
                  box(title = " PCA Analysis 2D Plot",
                      width = 12,height = "1000",
                      solidHeader = T, status = "info",
                      plotOutput("PCAPlot2D"),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      downloadButton("PlotDownloadPCA2D", 
                                     label = "Download"),
                      tags$style("#PlotDownload{float:left; }")
                  )
                  
              )
            )    
    ),
    #########################################################################################################
    ######Multidimensional Scaling 
    tabItem(tabName="MDS",
            fluidRow(  
              box(title = "MDS 2D Analysis ",
                  solidHeader = T, status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,height = "1400",
                  box(title = " Different Filtering Options",
                      width = 12,
                      solidHeader = T,
                      status = "info",
                      fluidRow(
                        column(4,
                               p(strong(span("MDS 2D-Video"),a(strong("Watch it here"),href="https://youtu.be/azpMiC0HonQ", target="_blank"),style = "color:red")),
                               radioButtons(inputId="optlab",label="Select label option",choices=c(ON="1",OFF="2"), selected='1'),
                               p(strong("MDS:"),"Multidimensional Scaling is a dimension-reduction technique designed to project high dimensional data down to 2 dimensions while preserving relative distances between observations. We used the cmdscale function in the R stats package to perform MDS analysis."), 
                               a(" MDS References",href="https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html", target="_blank")
                        ),
                        column(3,
                               radioButtons(inputId="optchs",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                               sliderInput("label_mds", "Size of the labels", 0, 10, 4),  
                               sliderInput("point_mds", "Size of the points", 0, 10, 3),  
                               numericInput("plot_heightmds", "Plot_Height (# pixels): ",min = 300, max = 600, value = 500),
                               numericInput("plot_widthmds", "Plot_Width (# pixels):", min = 600, max = 1000,value = 900)
                               
                        ))
                  ),
                  box(title = " MDS Analysis 2D Plot",
                      width = 12,height = "800",
                      solidHeader = T, status = "info",
                      plotOutput("MDSPlot"),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      downloadButton("PlotDownloadMDS2D", 
                                     label = "Download"),
                      tags$style("#PlotDownload{float:left; }")
                  )
                  
              )
            )    
    ),
    #####################################################################################################
    # ssGSEA Tab  gene set enrichment analysis
    tabItem(tabName="ssGSEA",
            fluidRow(          
              box(title = "ssGSEA Analysis ",
                  solidHeader = T, status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,
                  box(title = " Analysis Filtering Criteria",
                      width = 12,
                      solidHeader = T,
                      status = "info",
                      fluidRow(
                        column(4,
                               p(strong(span("ssGSEA-Video"),a(strong("Watch it here"),href="https://youtu.be/FQcy4nR1ZjE", target="_blank"),style = "color:red")),
                               radioButtons(inputId="optionspecies",label="Select one Species",choices=c(Mouse="Mouse",Human="Human"), selected='Mouse'),
                               radioButtons(inputId="option",label="Select one option",choices=c(Hallmark="Hallmark",GeneOntology="GeneOntology",USER_Geneset="USER_Geneset"), selected='Hallmark'),
                               fileInput("file5","Upload Your list of Gene's File here (.rdata only):"),
                               a("Convert yout .txt file to .rdata here", href="//mouseclassifier.shinyapps.io/SSGSEA/", target="_blank")
                        ),
                        
                        column(6,
                               p(strong("Please pick a related Species according to your already Uploaded Transcriptional Data Matrix."),style="color:red;"),
                               p(strong("If your Uploaded data is Read count (Integer Numbers) and you didn't run the Differential Gene Expression Analysis yet, first you need to run Differential Gene Expression Analysis (Heatmap section) to get your data normalized then come back to run this analysis."),style="color:red;"),
                               p(strong("ssGSEA:"),"The modification of standard gene set enrichment analysis (GSEA) specifically for single sample classification (ssGSEA) is performed using the GSVA package (version 1.32.0)"),a("GSVA ",href="https://doi.org/10.1186/1471-2105-14-7", target="_blank"),p("The R package msigdbr (version 7.1.1) is used to retrieve mouse Hallmark and biological processes (GO_BP) gene sets and applied to samples."),
                               a("msigdbr",href="https://doi.org/10.1016/j.cels.2015.12.004", target="_blank"), p("Enrichment score from ssGSEA demonstrate the degree to which the genes in a particular gene set are co-ordinately up- or down-regulated within a sample"), a("Enrichment score",href="https://doi.org/10.1038/nature08460", target="_blank")
                               
                        )
                      ),
                      actionButton(inputId="ssGSEArAnalysis", label="Submit"),
                      tags$style("button#ssGSEArAnalysis {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                  ),
                  
                  box(title = " Analysis Results",
                      width = 12,
                      solidHeader = T, status = "info",
                      br(),
                      downloadButton("TableDownloadssGSEA", 
                                     label = "Download"),
                      tags$style("#TableDownload {float:left; }"),
                      
                      DT::dataTableOutput("ssGSEA")
                  )
              )
            )    
    ),
    ##############Hall Mark
    tabItem(tabName="Hallmarkplot",
            fluidRow(          
              box(title = "ssGSEA Hallmark Plot ",
                  solidHeader = T,status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,height = "1000px",
                  column(3,
                         radioButtons(inputId="optionhall",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                         br(),
                         numericInput("plot_heighthall", "Plot height (# pixels): ", min = 300, max = 800, value = 544),
                         numericInput("plot_widthhall", "Plot width (# pixels):", min = 300, max = 750, value = 600),
                         downloadButton("PlotDownloadhall", 
                                        label = "Download"),
                         br(),
                         tags$style("#PlotDownload {float:left; }")),
                  column(8,  
                         (plotOutput("hallPlot")))           
                  
                  
              ))),
#############################################
    #GSEA PLots
    #Differential gene expressionS Tab Heatmap
    tabItem(tabName="GSEA-plot",
            fluidRow(
              box(title = " Gene Set Enrichment Analysis (GSEA) Filtering",
                  box(title = "Please Choose Two Groups",
                      width = 6,
                      solidHeader = T,
                      status = "info",
                      fluidRow(
                        column(3,
                               radioButtons(inputId="optionspeciesplot",label="Select one Species",choices=c(Mouse="Mouse",Human="Human"), selected='Mouse'),
                              
                        ),
                        column(4,
                               uiOutput("check3"),
                               
                        ),
                        column(4,
                               uiOutput("check4")
                        ),
                        
                        column(6,
                               uiOutput("norm")
                        ),
                      )
                      
                  ),  
                  width = 12,
                  solidHeader = T,
                  status = "info",
                  fluidRow(
                    column(6,
                           p(strong("Please Note for GeneOntology option, Only the first 50 pathways having highest ES will be displayed here."),style="color:red;"),
                           uiOutput("check5"),
                           textOutput("leading"),
                           tags$head(tags$style("#leading{color: blue;
                                 font-size: 15px;
                                 font-style: italic;
                                  font-weight: bold;
                                 }")),
                          br(),
                          br(),
                          
                           downloadButton("listofgenes", 
                                          label = "Download-Genelists"),
                           br(),
                           tags$style("#TableDownload {float:right; }")

                           
                    )
                    
                  ),
                  fluidRow(
                    column(6,
                          actionButton(inputId="gseaplot", label="Hallmark"),
                          tags$style("button#gseaplot {float:left;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.0em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                          actionButton(inputId="gegoplot", label="GeneOntology"),
                          tags$style("button#gegoplot {float:right;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.0em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                           br(),
                           br(),
                           p(strong("GSEA:"),"Gene set enrichment analysis (GSEA) is a computational method, used to determine whether a priori defined set of genes display significant differences between two biological phenotypes. This is performed using the fgsea package (version 1.16.0) and enrichplot(version 1.10.0)"),a("GSEA",href="https://www.pnas.org/content/102/43/15545.short/", target="_blank"),
                           br(),
                           br(),
                           p(strong(span("GSEA-Video"),a(strong("Watch it here"),href="https://youtu.be/tjpPELDqLfI", target="_blank"),style = "color:red")),
                    )  
                  )
              )
              
            ),
            
            box(title = "Gene Set Enrichment Analysis Plot",
                width = 16, 
                height = "1000px",
                solidHeader = T,
                status = "info",
                column(3,
                       radioButtons(inputId="optiongsea",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                       br(),
                       numericInput("plot_heightgsea", "Plot height (# pixels): ", min = 300, max = 800, value = 500),
                       numericInput("plot_widthgsea", "Plot width (# pixels):", min = 300, max = 750, value = 600),
                       verbatimTextOutput("info"),
                       downloadButton("PlotDownloadgsea", 
                                      label = "Download"),
                       br(),
                       tags$style("#PlotDownload {float:left; }")),
                column(8,  
                       (plotOutput("gseaPlot"))
                       
                )
                
            )), #end heatmap
    
    
    #####################################################################################################
    #####################################################################################################
    tabItem(tabName="CRIS_T",
            fluidRow(          
              box(title = "CRIS-Classifier ",
                  solidHeader = T, status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,
                  box(title = " Analysis Filtering Criteria",
                      width = 12,
                      solidHeader = T,
                      status = "info",
                      fluidRow(
                        column(4,
                               radioButtons(inputId="optionCRIS",label="Select one Species",choices=c(Mouse="Mouse",Human="Human"), selected='Mouse'),
                        ),
                        
                        column(6 ,
                               p(strong("If your Uploaded data is a Read count (Integer Numbers) and you didn't run the Differential Gene Expression Analysis yet, first you need to run Differential Gene Expression Analysis (Heatmap section) to get your data normalized then come back to run this analysis."),style="color:red;"),
                               p(strong("Warning, you should have at least 20 samples to get reliable result using current classifer method."),style="color:blue;"),
                               p(strong("CRIS-Classifier for Human option:"),"This method performs CRC intrinsic subtypes (CRIS) classification based on nearest template prediction (NTP) method and CRIS template used in CMScaller R package."),
                               p(strong("CRIS-Classifier for Mouse option:"),"In order to do CRC intrinsic subtypes (CRIS) classification in mouse, the human CRIS template, embedded in CMScalller, was converted to mouse orthologues using biomaRt, R package, and NTP method was applied to call CRIS subtypes in mouse tissues."),
                               a("biomaRt",href="https://www.nature.com/articles/nprot.2009.97", target="_blank"),"/",a("CMScaller",href="https://www.nature.com/articles/s41598-017-16747-x", target="_blank"),"/",a("NTP",href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015543", target="_blank"),
                               "/",a("CRIS",href="https://www.nature.com/articles/ncomms15107", target="_blank")
                               
                        )
                      ),
                      br(),
                      br(),
                      actionButton(inputId="CRISAnalysis", label="Submit"),
                      tags$style("button#CRISAnalysis {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                  ),
                  
                  box(title = "CRIS-Classifier Table and Plot ",
                      solidHeader = T,status = "info",
                      collapsible = F, collapsed = F,
                      width = 12,height = "1000px",
                      tabsetPanel(type = "tabs",
                                  tabPanel("CRIS_Table",
                                           br(),
                                           downloadButton("TableDownloadCRIS", 
                                                          label = "Download"),
                                           tags$style("#TableDownload {float:left; }"),
                                           
                                           DT::dataTableOutput("CRIS")    
                                           
                                           
                                  ),
                                  tabPanel("CRIS_Plot",
                                           column(3,
                                                  radioButtons(inputId="optioncris",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                                                  br(),
                                                  numericInput("plot_heightcris", "Plot height (# pixels): ", min = 300, max = 800, value = 500),
                                                  numericInput("plot_widthcris", "Plot width (# pixels):", min = 300, max = 750, value = 600),
                                                  downloadButton("PlotDownloadCRIS", 
                                                                 label = "Download"),
                                                  br(),
                                                  tags$style("#PlotDownload {float:left; }")),
                                           column(8,  
                                                  (plotOutput("CRIS_Plot")))
                                           
                                  )
                      )
                      
                  )
              )
            )
    ),
########################## CMS classiffier##############################################################
#####################################################################################################
    tabItem(tabName="CMS_T",
            fluidRow(          
              box(title = "CMS-Classifier ",
                  solidHeader = T, status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,
                  box(title = " Analysis Filtering Criteria",
                      width = 12,
                      solidHeader = T,
                      status = "info",
                      fluidRow(
                        column(4,
                               radioButtons(inputId="optionCMS",label="Select one Species",choices=c(Mouse="Mouse",Human="Human"), selected='Mouse'),
                        ),
                        
                        column(6 ,
                               p(strong("If your Uploaded data is Read count (Integer Numbers) and you didn't run the Differential Gene Expression Analysis yet, first you need to run Differential Gene Expression Analysis (Heatmap section) to get your data normalized then come back to run this analysis."),style="color:red;"),
                               p(strong("Warning, you should have at least 20 samples to get reliable result using current classifer method."),style="color:blue;"),
                                p(strong("CMS-Classifier for Human option:"),"This method performs Consensus Molecular Subtypes (CMS) classification in human based on nearest template prediction (NTP) method and CMS template used in CMScaller R package."),
                               p(strong("CMS-Classifier for Mouse option:"),"In order to do CMS classification in mouse, the human CMS template, embedded in CMScalller, is converted to mouse orthologues using biomaRt, R package, and intersected mouse genes across CMS subtypes were removed. Then, NTP method was applied to call CMS in mouse tissues."),
                               a("biomaRt",href="https://www.nature.com/articles/nprot.2009.97", target="_blank"),"/",a("CMScaller",href="https://www.nature.com/articles/s41598-017-16747-x", target="_blank"),"/",a("NTP",href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015543", target="_blank")
                               
                        )
                      ),
                      br(),
                      br(),
                      actionButton(inputId="CMSAnalysis", label="Submit"),
                      tags$style("button#CMSAnalysis {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                  ),
                  box(title = "CMS-Classifier Table and Plot ",
                      solidHeader = T,status = "info",
                      collapsible = F, collapsed = F,
                      width = 12,height = "1000px",
                      tabsetPanel(type = "tabs",
                                  tabPanel("CMS_Table",
                                           br(),
                                           downloadButton("TableDownloadCMS", 
                                                          label = "Download"),
                                           tags$style("#TableDownload {float:left; }"),
                                           
                                           DT::dataTableOutput("CMS")
                                  ),
                                  tabPanel("CMS_Plot",
                                           column(3,
                                                  radioButtons(inputId="optioncms",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                                                  br(),
                                                  numericInput("plot_heightcms", "Plot height (# pixels): ", min = 300, max = 800, value = 500),
                                                  numericInput("plot_widthcms", "Plot width (# pixels):", min = 300, max = 750, value = 600),
                                                  downloadButton("PlotDownloadCMS", 
                                                                 label = "Download"),
                                                  br(),
                                                  tags$style("#PlotDownload {float:left; }")),
                                           column(8,  
                                                  (plotOutput("CMS_Plot")))           
                                           
                                           
                                  )
                      ))
              ))),         
    
##################################################################################################################  
#####################################################################################################################################
    #Gene Expression Levels
    tabItem(tabName="Gene-plot",
            fluidRow(
              box(title = " Gene Expression Levels",
                  width = 12,
                  solidHeader = T,
                  status = "info",
                  box(title = "Gene Expression Tables",
                      width = 12, 
                      solidHeader = T,
                      status = "info",
                      fluidRow(
                        column(12,
                               p(strong(span("Gene Expression Levels-Video"),a(strong("Watch it here"),href="https://youtu.be/Ya-nnvgcL_w ", target="_blank"),style = "color:red")),
                                br(),
                               div(dataTableOutput("ID")))
                      )
                  ),
                  
                  box(title = "Gene Expression Plot",
                      width = 12, 
                      height = "800px",
                      solidHeader = T,
                      status = "info",
                      fluidRow(
                        column(3,
                               radioButtons(inputId="optiongene",label=" Pick on Option for Plot ",choices=c(CountPlot="1"), selected='1'),
                                p(strong("CountPlot:")," click on a gene of interest from the table above. (The Table is only for this option)",style="color:red;"),
                               br(),
                               radioButtons(inputId="optionchgene",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                               numericInput("plot_heightgene", "Plot_Height (# pixels): ",min = 300, max = 600, value = 500),
                               numericInput("plot_widthgene", "Plot_Width (# pixels):", min = 600, max = 750,value = 700)
                        ),
                        column(8,
                               plotOutput("genePlot"),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               downloadButton("PlotDownloadgene", 
                                              label = "Download"),
                               tags$style("#PlotDownload {float:left; }")
                               
                        ))
                  ))
            )),
    #####################################################
    #Differential gene expressionS Tab Heatmap
    tabItem(tabName="MSG-heatmap",
            fluidRow(
              box(title = " Differential Gene Expression Analysis",
                  box(title = "Differential Categories comparison Analysis",
                      width = 6,
                      solidHeader = T,
                      status = "info",
                      textOutput("Groupcomparison"),
                      tags$style("#Groupcomparison{ font-weight: bold; color: #0033FF; padding-top: .3cm; padding-bottom: .3cm;}"),
                      p("Please select Maximum 5  groups for analysis in each Differential Category A and  B "
                        , style="font-weight: bold"),
                      p("Please DO Not compare 2 same groups",style="color:red;"),
                      fluidRow(
                        column(6,
                               uiOutput("check1")),
                        column(6,
                               uiOutput("check2")
                        )
                      ),
                      p(strong("Differential gene expression analysis:"),"DESeq2 package (version 1.24.0) was used to identify genes which are differentially expressed between the two groups selected by user."),
                      a("DESeq2",href="https://doi.org/10.1186/s13059-014-0550-8", target="_blank"),
                      br(),
                      p(strong("Heatmap:"),"The heatmaply package (version 1.1.1) is used to visualize the differentially expressed genes."),
                      a("heatmaply",href=" https://cran.r-project.org/web/packages/heatmaply/heatmaply.pdf", target="_blank"),
                      br(),
                      br(),
                      p(strong(span("Heatmap-Video"),a(strong("Watch it here"),href="https://youtu.be/2WcGevvoItU", target="_blank"),style = "color:red")),
                  ),                     
                  
                  width = 12,
                  solidHeader = T,
                  status = "info",
                  fluidRow(
                    column(2,
                           p("By clicking Create Plot button, the Heatmap plot will be created using default values as below. Then you can change these values and press it again to Update plot."
                             ,style="color:red;"),
                           textInput(inputId="MSGpval", 
                                     label="Padj-value:", 
                                     value="0.05"),
                           
                           textInput(inputId="MSGfc", 
                                     label="log2FoldChange(min):", 
                                     value="-2"),
                           
                           textInput(inputId="MSGfc1", 
                                     label="log2FoldChange(max):", 
                                     value="2"),
                           
                           actionButton(inputId="findheat", label="Create Plot"),
                           tags$style("button#findheat{margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 20px; font-family:Andika, Arial, sans-serif; font-size:1.0em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 10px;box-shadow: rgba(0, 0, 0, .55) 0 1px 3px;}"),
                           br(),
                           br(),
                           
                           useShinyalert() 
                           
                    ),
                    column(2,
                           radioButtons(inputId="optionheat",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                           radioButtons(inputId="optionl",label="Select label option",choices=c(ON="1",OFF="2"), selected='2'),
                           uiOutput("checkdan"),
                           numericInput("plot_heightheat", "Plot_Height (# pixels): ",min = 300, max = 800, value = 650),
                           numericInput("plot_widthheat", "Plot_Width (# pixels):", min = 600, max = 1000,value = 800)
                    ),
                    column(2,
                           textInput(inputId="Name", 
                                     label= "Output Filename", 
                                     value="")
                    )
                  ),
                  fluidRow(
                    column(6,
                           actionButton(inputId="FilterAnalysis", label="Submit"),
                           tags$style("button#FilterAnalysis {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                    )  
                  )
              )
            ),
            
            box(title = "Top Significant Genes Heatmap",
                width = 18, 
                height = "1400px",
                solidHeader = T,
                status = "info",
                tabsetPanel(type = "tabs",
                            tabPanel("Heatmap Plot",
                                     column(2,
                                            br(),
                                            textOutput("upregulated"),
                                            tags$head(tags$style("#upregulated{color: red;
                                 font-size: 15px;
                                 font-style: italic;
                                  font-weight: bold;
                                 }"
                                            )),
                                            br(),
                                            textOutput("Downregulated"),
                                            tags$head(tags$style("#Downregulated{color: blue;
                                 font-size: 15px;
                                 font-style: italic;
                                 font-weight: bold;
                                 }"
                                            )),
                                            br(),
                                            p(strong("Plotly Modebar control:")," From left to right"),
                                            p( "1)Download plot as a selected format."),
                                            p( "2) Zoom on selection. Click and hold you mouse button to select a region to zoom on."),
                                            p( "3)Image panning."),
                                            p( "4)Zoom IN."),
                                            p( "5)Zoom OUT."),
                                            p( "6)Reset axis to the original axis ranges."),
                                            p( "7)Toggle spike lines on and off. Shows lines on a graph indicating the exact x-axis and y-axis."),
                                     ),
                                     column(8,
                                            plotlyOutput("heatmapMSGPlot")
                                     )
                            ),
                            
                            tabPanel("Differential Gene Expression Table",
                                     downloadButton("DGEtable", 
                                                    label = "Download txt-file"),
                                     tags$style("#TableDownload {float:left; }"),
                                     
                                     DT::dataTableOutput("DGStable")     

                                     
                            ) 
                )
            )
            
            
    ), #end heatmap
#################################################################################################################
#######################################################################################################################
###Volcano plot
tabItem(tabName="Volcano",
        box(title = "Volcano Analysis",
            width = 16, 
            height = "3000px",
            solidHeader = T,
            status = "info",
            tabsetPanel(type = "tabs",
                        tabPanel("Volcano Plot",
                                 column(3,
                                        br(),
                                        p("The Volcano plot is created using default values as below using Plotly. You can change these values and press Update plot.",style="color:red;"),
                                        p("Afterwards, you can pick any gene names (max=10), choose selected Genes option and then new plot will be created.",style="color:red;"),
                                        p(strong(span("Volcano-Video"),a(strong("Watch it here"),href="https://youtu.be/x42WkcCnG08", target="_blank"),style = "color:green")),
                                        radioButtons(inputId="volc",label="Select Plot type ",choices=c("using-plotly"="1","selectedGenes"="2"), selected='1'),
                                        sliderInput("point_volc", "Size of the points", 0, 10, 2),  
                                        textInput(inputId="volcname", 
                                                  label= "Output Filename", 
                                                  value=""),
                                        
                                        textInput(inputId="fc_cutoff", 
                                                  label="Fold Change threshold(min):", 
                                                  value="-2"),
                                        
                                        textInput(inputId="fc_cutoff1", 
                                                  label="log2FoldChange(min):", 
                                                  value="2"),
                                        
                                        textInput(inputId="p_cutoff", 
                                                  label="Significance threshold:", 
                                                  value="0.05"),
                                        
                                        actionButton(inputId="volcAnalysis", label="Update Plot"),
                                        br(),
                                        br(),
                                        selectizeInput(inputId = 'user_gene_list',
                                                       label = "User selected Genes (max= 10):",
                                                       choices = "",
                                                       selected = "",
                                                       multiple = TRUE, # allow for multiple inputs
                                                       options = list(create = TRUE)), # if TRUE, allows newly created inputs
                                        br(),
                                        downloadButton("volcdownload", 
                                                       label = "Download Selected GenePlot"),
                                        tags$style("#volcdownload {float:left; }"),
                                        br(),
                                        br(),
                                        br(),
                                        radioButtons(inputId="optionvolc",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                                        numericInput("plot_heightvolc", "Plot_Height (# pixels): ",min = 300, max = 800, value = 560),
                                        numericInput("plot_widthvolc", "Plot_Width (# pixels):", min = 600, max = 1000,value = 700),
                                        br(),
                                       
                                        p(strong("Plotly Modebar control:")," From left to right"),
                                        p( "1)Download plot as a selected format."),
                                        p( "2) Zoom on selection. Click and hold you mouse button to select a region to zoom on."),
                                        p( "3)Image panning."),
                                        p( "4)Zoom IN."),
                                        p( "5)Zoom OUT."),
                                        p( "6)Reset axis to the original axis ranges."),
                                        p( "7)Toggle spike lines on and off. Shows lines on a graph indicating the exact x-axis and y-axis."),
                                 ),
                                 column(8,
                                        plotlyOutput("volcanoPlot"),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(), 
                                        br(),
                                        br(), 
                                        br(),
                                        br(), 
                                        br(),
                                        br(), 
                                        plotOutput("volcanoPlot2" )
                                     
                                 )
                        ),
                        
                        tabPanel("Differential Gene Expression Table",
                                 downloadButton("TableDownloadMSG","Download txt-file", 
                                                label = "Download txt-file"),
                                 tags$style("#TableDownload {float:left; }"),
                                 
                                 DT::dataTableOutput("volc")
                                 
                        )
                        
                        
            )
        )
),
    
#########################################################################################################
    #Selected Gene List Heatmap
    tabItem(tabName="heatmap",
            fluidRow(
              box(title = "Selected Gene List Heatmap Plot",
                  width = 12,
                  solidHeader = T,
                  status = "info",
                  fluidRow(
                    column(3,
                           p(strong(span("Heatmap Selected Gene-Video"),a(strong("Watch it here"),href="https://youtu.be/GHQMm5Rzw4Q", target="_blank"),style = "color:green")),
                           p("Please Pick one option for your plots Format:",style="color:red;"),
                           radioButtons(inputId="optiongroupsample",label="Select one Options ",choices=c("Plot by Samples"="1","Plot by z-score and Groups"="2"), selected='1'),
                           radioButtons(inputId="optionsheatplot",label="Select Plot Options ",choices=c(Without_Thresholds="1",Default_Thresholds="2"), selected='1'),
                           p("Please Click on Minimum Two Genes For Clustering:",style="color:red;"),
                           p("Without-Thresholds:This Table is created without applying any Thresholds.",style="color:blue;"),
                           p("Default-Thresholds:This Table is created using  Default values in previous section.",style="color:blue;"),
                           radioButtons(inputId="optionselecheat",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                           radioButtons(inputId="optionlselectheat",label="Select label option",choices=c(ON="1",OFF="2"), selected='1'),
                           uiOutput("checkdanselected"),
                           numericInput("plot_heightselectheat", "Plot_Height (# pixels): ",min = 300, max = 800, value = 650),
                           numericInput("plot_widthselectheat", "Plot_Width (# pixels):", min = 600, max = 1000,value = 800),
                           textInput(inputId="Nameheat", 
                                     label= "Output Filename", 
                                     value="")
                    ),
                    column(8,
                           DT::dataTableOutput("Heatmapselect"))
                  )
              )
            ),
            box(title = "Selected  Genes List Heatmap",
                width = 16, 
                height = "1000px",
                solidHeader = T,
                status = "info",
                column(2,
                       p(strong("Plotly Modebar control:")," From left to right"),
                       p( "1)Download plot as a selected format."),
                       p( "2) Zoom on selection. Click and hold you mouse button to select a region to zoom on."),
                       p( "3)Image panning."),
                       p( "4)Zoom IN."),
                       p( "5)Zoom OUT."),
                       p( "6)Reset axis to the original axis ranges."),
                       p( "7)Toggle spike lines on and off. Shows lines on a graph indicating the exact x-axis and y-axis."),
                ),
                column(8,
                       plotlyOutput("dotplot")
                )
            ) 
    ),
###############################################################################################################
####################################################################################################################
    ## Analysis mMCP-counter
    tabItem(tabName="mMCP",
            fluidRow(    
              box(title = "Estimate The Abundance of Microenvironment Cell Populations In Mouse/Human Tissue",
                  solidHeader = T, status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,
                  radioButtons(inputId="conter",label="Select one Species",choices=c(Mouse="Mouse",Human="Human"), selected='Mouse'),
                  p(strong(span("mMcp-Video"),a(strong("Watch it here"),href="https://youtu.be/sc4y54-FIL8", target="_blank"),style = "color:green")),
                  p(strong("If your Uploaded data is Read count (Integer Numbers) and you didn't run the Differential Gene Expression Analysis yet, first you need to run Differential Gene Expression Analysis (Heatmap section) to get your data normalized then come back to run this analysis."),style="color:red;"),
                  p(strong("mMCP-counter:"),"This function estimates the quantity of several immune and stromal cell populations from heterogeneous transcriptomic data, which has been modified for use specifically murine samples.",a("mMCP-counter",href="https://doi.org/10.1186/s13073-020-00783-w", target="_blank")), 
                  # br(),
                  p(strong("MCP-counter:"),"This function estimates the quantity of several immune and stromal cell populations from absolute abundance cell populations in heterogeneous transcriptomic data, which has been proposed for use specifically with human samples."),
                  a("MCP-counter",href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5", target="_blank"), 
                  br(),
                  actionButton(inputId="bb", label="Submit"),
                  tags$style("button#bb{margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 4px;}"),
                  br(),
                  DT::dataTableOutput("mMCPRes"),
                  br(),
                  downloadButton("TableDownloadmMCP", 
                                 label = "Download"),
                  tags$style("#TableDownload {float:right; }"),
              )   
            )
    ),
    ############################                  
    tabItem(tabName="mMCPplot",  
            fluidRow(  
              box(title = "Analysis Plot",
                  solidHeader = T, status = "info",
                  collapsible = F, collapsed = F,
                  width = 12,
                  height = "1000px",
                  column(3,
                         radioButtons(inputId="optionmcp",label="Select Plot Downloading Format ",choices=c(png="1",svg="2"), selected='1'),
                         br(),
                         numericInput("plot_heightmcp", "Plot height (# pixels): ", min = 300, max = 800, value = 500),
                         numericInput("plot_widthmcp", "Plot width (# pixels):", min = 300, max = 750, value = 600),
                         verbatimTextOutput("console"),
                         downloadButton("PlotDownloadmMCP", 
                                        label = "Download"),
                         br(),
                         tags$style("#PlotDownload {float:left; }")),
                  column(8,  
                         (plotOutput("mMCPlot"))
                         
                  )
              )
              
            )
            
    )
#####################################################################################################
###################################################################################################################        
######################################################   
  ) #tabItems(
) #body <- dashboardBody(

ui <- dashboardPage(header, sidebar, body, skin = "purple")

