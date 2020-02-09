require(XML)
library(plyr)
library(dplyr)
library(DT)
library(stringr)
library(shiny)
library(shinyalert)




# Start server code -------------------------------------------------------

server <- function(input, output, session){
  # Read in BLAST functions
  source("blast_utilities.R", local = TRUE)
  
  
#bserveEvent(input$preview, {
#   if (condition == TRUE) {
#     do something }
#   shinyalert("Oops!", "Something went wrong.", type = "error")
# })
#
#   need(Sys.which("blastn") != "", shinyalert(
#       title = "Hello",
#       text = "This is a modal",
#       closeOnEsc = TRUE,
#       closeOnClickOutside = TRUE,
#       html = TRUE,
#       type = "error",
#       showConfirmButton = TRUE,
#       showCancelButton = FALSE,
#       confirmButtonText = "OK",
#       confirmButtonCol = "#AEDEF4",
#       timer = 0,
#       imageUrl = "",
#       animation = TRUE
# )
  
  blastresults <- eventReactive(input$blast, {
    if(Sys.which("blastn") == ""){
      shinyalert("Error!", "Could not locate blast install.", type = "error")
      #shinyalert(
      #  title = "Error!",
      #  text = "Please manually input blast location",
      #  type = "input",
      #  inputType="text",
      #  callbackR = function(x) {
      #    message("Attempting to load blast from ", normalizePath(x)) 
      #    blast_setpath(normalizePath(x))
      #    message(Sys.which("blastn"))
      #    if (Sys.which("blastn") == ""){
      #      shinyalert("Error!", "Could not locate blast install.", type = "error")
      #    } else if (Sys.which("blastn") != ""){
      #      shinyalert("Sucess!", "Blast sucessfully located", type = "success")
      #    }
      #  })
    } else {
    #gather input and set up temp file
    query <- input$query
    tmp <- tempfile(fileext = ".fa")
    tmpdb <- tempfile(fileext = ".fa")
    
  #if else chooses the right database
  if (input$db == "Local" && length(input$custom_db_path$datapath) == 1){
      makeblastdb(input$custom_db_path$datapath)
      db <-  input$custom_db_path$datapath
      remote <- c("")
  } else if (input$db == "Local" && length(input$custom_db_path$datapath) > 1){
    test <- Biostrings::readDNAStringSet(filepath=input$custom_db_path$datapath, format="fasta")
    Biostrings::writeXStringSet(test, filepath=tmpdb, format="fasta")
    makeblastdb(tmpdb)
    db <-  tmpdb
    remote <- c("")
   } else {
    db <- c("nr")
    remote <- c("-remote")
  }
  
    # ensure query is formatted as fasta
    if (startsWith(query, ">")){
      writeLines(query, tmp)
    } else {
      writeLines(paste0(">Query\n",query), tmp)
    }
    
    # Run BLAST 
    data <- system(paste0(input$program," -query ",tmp," -db ",db," -evalue ",input$eval," -outfmt 5 -max_hsps 1 -max_target_seqs 10 ",remote), intern = T)
    xmlParse(data)
    }
  }, ignoreNULL= T)
  
  # Parse Results
  parsedresults <- reactive({
    if (is.null(blastresults())){}
    else {
      shinyalert("Success!", "BLAST complete.", type = "success", timer=2000, animation = TRUE)
      xmltop = xmlRoot(blastresults())
      
      #the first chunk is for multi-fastas
      results <- xpathApply(blastresults(), '//Iteration',function(row){
        query_ID <- getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
        hit_IDs <- getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
        hit_def <- getNodeSet(row, 'Iteration_hits//Hit//Hit_def') %>% sapply(., xmlValue)
        hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
        bitscore <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score') %>% sapply(., xmlValue)
        eval <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue') %>% sapply(., xmlValue)
        Hsp_identity <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_identity') %>% sapply(., xmlValue)
        Hsp_align <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_align-len') %>% sapply(., xmlValue)
        perc_id <- signif((as.numeric(Hsp_identity) / as.numeric(Hsp_align)) * 100,digits=4 )
        cbind(query_ID, hit_IDs, hit_def, hit_length, bitscore, eval, perc_id)
      })
      #this ensures that NAs get added for no hits
      results <-  rbind.fill(lapply(results,function(y){as.data.frame((y), stringsAsFactors=FALSE)}))
    }
  })
  
  #makes the datatable
 #output$blastResults <- renderDataTable({
 #  if (is.null(blastresults())){
 #  } else {
 #   
 #  }
 #}, selection="single",class = "white-space: nowrap")
 #
  output$blastResults <- renderDataTable({
    if(is.null(blastresults())){return ()}
    DT::datatable( parsedresults() , options = list(lengthMenu = c(5,10),
                                                        columnDefs = list(list(
                                                          targets = "_all",
                                                          render = JS(
                                                            "function(data, type, row, meta) {",
                                                            "return type === 'display' && data != null && data.length > 50 ?",
                                                            "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
                                                            "}")
                                                        ))),
                   selection="single", class = "display")
  })
  
  
  
  #Get the alignemnt information from a clicked row
  output$clicked <- renderTable({
    if(is.null(input$blastResults_rows_selected)){}
    else{
      xmltop <- xmlRoot(blastresults())
      clicked <- input$blastResults_rows_selected
      tableout <- data.frame(parsedresults()[clicked,])

      tableout <- t(tableout)
      names(tableout) <- c("")
      rownames(tableout) <- c("Query ID","Hit ID","Hit Def", "Length", "Bit Score", "e-value", "% ID")
      colnames(tableout) <- NULL
      data.frame(tableout)
    }
  }, rownames =T,colnames =F)
  
  #this chunk makes the alignments for clicked rows
  output$alignment <- renderText({
    if(is.null(input$blastResults_rows_selected)){}
    else{
      xmltop = xmlRoot(blastresults())
      
      clicked = input$blastResults_rows_selected
      
      #loop over the xml to get the alignments
      align <- xpathApply(blastresults(), '//Iteration',function(row){
        top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
        mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
        bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
        rbind(top, mid, bottom)
      })
      
      #split the alignments every 60 carachters to get a "wrapped look"
      alignx <- do.call("cbind", align)
      splits <- strsplit(gsub("(.{60})", "\\1,", alignx[1:3, clicked]),",")
      
      #paste them together with returns '\n' on the breaks
      split_out <- lapply(1:length(splits[[1]]),function(i){
        rbind(paste0("Q-", splits[[1]][i], "\n"),
              paste0("M-", splits[[2]][i], "\n"),
              paste0("H-", splits[[3]][i], "\n"))
      })
      unlist(split_out)
    }
  })
}