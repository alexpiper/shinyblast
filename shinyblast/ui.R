library(shinythemes)
library(DT)
library(shiny)
library(shinyalert)

custom_db <- c("local_db")

ui <- fluidPage(theme = shinytheme("lumen"),
                tagList(
                  tags$head(
                    tags$link(rel="stylesheet", type="text/css",href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js")
                  )
                ),
                useShinyalert(),
                #This block gives us all the inputs:
                mainPanel(
                  headerPanel(
                    list(HTML('<img src="agvic.jpg" width="240" height="70"/>'), "shinyBLAST"),
                    windowTitle="shinyBLAST"),
                  textAreaInput('query', 'Input sequence:', value = "", placeholder = "", width = "600px", height="200px"),
                  radioButtons("db", label = "", choices = c("Local", "NCBI"), selected = "Local", width="120px"),
                  fileInput("custom_db_path", "Select fasta files to blast against", multiple=TRUE,
                            accept = c(".fa", ".fasta", ".fa.gz", ".fasta.gz"), placeholder="Select FASTA file"),
                  #selectInput("db", "Databse:", choices = c("local", "nr"), width="120px"),
                  div(style="display:inline-block",
                      selectInput("program", "Program:", choices=c("blastn","tblastn"), width="100px")),
                  div(style="display:inline-block",
                      selectInput("eval", "e-value:", choices=c(1,0.001,1e-4,1e-5,1e-10), width="120px")),
                  actionButton("blast", "BLAST!")
                ),
                
                #this snippet generates a progress indicator for long BLASTs
                div(class = "busy",  
                    p("Calculation in progress.."), 
                    img(src="https://i.stack.imgur.com/8puiO.gif", height = 100, width = 100,align = "center")
                ),
                
                #Basic results output
                mainPanel(
                  h4("Results"),
                  DT::dataTableOutput("blastResults"),
                  p("Alignment:", tableOutput("clicked") ),
                  verbatimTextOutput("alignment")
                ),
                hr(),
                print("App maintainer: alexander.piper@agriculture.vic.gov.au")
)