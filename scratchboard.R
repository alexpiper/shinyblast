
# testing -----------------------------------------------------------------

source("blast_utilities.R", local = TRUE)

query <- "TTTTGGAAACTGACTAATCCCCCTTATAATTAATTCTACAGATATAGCATTCCCACGAATAAATAATCTTAGTTTTTGATTTCTTCCTCCCTCTTTATTTC"
tmp <- tempfile(fileext = ".fa")
writeLines(paste0(">Query\n",query), tmp)

tmpdb <- tempfile(fileext = ".fa")


makeblastdb("/Users/alexanderpiper/Dropbox/R/shinyblast/Anisembiidae_COI-5P.fa")
db <-  "/Users/alexanderpiper/Dropbox/R/shinyblast/Anisembiidae_COI-5P.fa"
remote <- c("")

# BBLAST
data <- system(paste0("blastn -query ",tmp ," -db ",db," -evalue 1 -outfmt 5 -max_hsps 1 -max_target_seqs 10 ",remote), intern = T)

blastresults <- xmlParse(data)

xmltop = xmlRoot(blastresults)

#the first chunk is for multi-fastas
results <- xpathApply(blastresults, '//Iteration',function(row){
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


# Parse alignments to get hits - write out

#loop over the xml to get the alignments
hits <- xpathApply(blastresults, '//Iteration',function(row){
  #Get query subset
  que <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
  print(que)
  q_from <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_query-from') %>% sapply(., xmlValue)
  q_to <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_query-to') %>% sapply(., xmlValue)

  q_seq <-  Biostrings::DNAStringSet(que)
  for (i in 1:length(que)){
    print(as.numeric(q_from)[i])
    print(as.numeric(q_to)[i])
    q_seq[i] <- subseq(q_seq[i], start=as.numeric(q_from)[i], end=as.numeric(q_to)[i])
  }
  #q_seq <- Biostrings::DNAStringSet(q_seq)

    # Get hits
  hit_def <- getNodeSet(row, 'Iteration_hits//Hit//Hit_def') %>% sapply(., xmlValue)
  hit <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
  h_from <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hit-from') %>% sapply(., xmlValue)
  h_to <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hit-to') %>% sapply(., xmlValue)
  #rbind(query, q_from, q_to, hit_def, hit, h_from, h_to)
  
  h_seq <-  Biostrings::DNAStringSet(hit)
  for (i in 1:length(hit)){
    print(as.numeric(q_from)[i])
    print(as.numeric(q_to)[i])
    h_seq[i] <- subseq(h_seq[i], start=as.numeric(h_from )[i], end=as.numeric(h_to)[i])
  }
  
  out <- c(q_seq, h_seq)
  return(out)
})

# Download full sequences

# Or align just trimmed sequences.

# Should be able to use the q_from and q_to as an index to the initial vector to add relevant padding 
# In order to make the sequences align. ie if q-from = 10, add 10 gaps to start of hit

ape::nj(hits[[1]])


align <- xpathApply(blastresults, '//Iteration',function(row){
  top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
  mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
  bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
  rbind(top, mid, bottom)
})
