
# .findExecutable ---------------------------------------------------------


.findExecutable <- function(exe, interactive=TRUE) {
  path <- Sys.which(exe)
  if(all(path=="")) {
    if(interactive) stop("Executable for ", paste(exe, collapse=" or "), " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.", call.=FALSE)
    return(character(0))
  }
  
  path[which(path!="")[1]]
}


# Install BLAST -----------------------------------------------------------

#' Install BLAST
#'
#' @param url A URL for the specific BLAST release to download. If missing, automatically download latest release.
#' @param destdir Destination to install BLAST to
#'
#' @return
#' @export
#'
#' @examples
blast_install <- function(url, destdir = "bin") {
  # get start time
  time <- Sys.time()
  # get OS
  localos <- Sys.info()["sysname"]
  
  if (missing(url)) {
    # find the latest version of BLAST
    url <- 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/'
    filenames <- RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE) %>%
      stringr::str_split("\r*\n") %>%
      unlist()
    
    if(localos == "Windows"){
    url <- filenames[str_detect(filenames,"win64.tar.gz$")] %>%
      paste0(url,.)
    } else if(localos == "Darwin"){
      url <- filenames[str_detect(filenames,"macosx.tar.gz$")] %>%
        paste0(url,.)
    } else if(localos == "unix"){
      url <- filenames[str_detect(filenames,"linux.tar.gz$")] %>%
        paste0(url,.)
    }
    
  }
  
  if (!dir.exists(destdir)) {
    dir.create(destdir) # Create first directory
  }
  
  version <- basename(url) %>% str_replace("(-x64)(.*?)(?=$)", "")
  if (dir.exists(paste0(destdir, "/",version))) {
    unlink(paste0(destdir, "/",version), recursive = TRUE) # Remove old version
  }
  
  destfile <- file.path(destdir, basename(url))
  if (exists(destfile)) {
    file.remove(destfile) # Remove old zip file
  }
  httr::GET(url, httr::write_disk(destfile, overwrite=TRUE))
  
  #unzip file and remove download
  utils::untar(destfile, exdir = destdir)
  file.remove(destfile)
  
  #Set new $Paths variable for mac & linux
  if(localos == "Darwin" | localos == "unix"){
    old_path <- Sys.getenv("PATH")
    install_path <- list.dirs(destdir, full.names = TRUE)[str_detect(list.dirs(destdir, full.names = TRUE),"/bin$")]
    Sys.setenv(PATH = paste(old_path, normalizePath(install_path), sep = ":"))
  }
  
  time <- Sys.time() - time
  message(paste0("Downloaded ", version, " in ", format(time, digits = 2)))
}


# Set blast PATH ----------------------------------------------------------

#Set new $Paths variable for mac & linux
blast_setpath <- function(path) {
  old_path <- Sys.getenv("PATH")
  install_path <- list.dirs(path, full.names = TRUE)[str_detect(list.dirs(path, full.names = TRUE),"/bin$")]
  Sys.setenv(PATH = paste(old_path, normalizePath(install_path), sep = ":"))
}


# Make Blast DB -----------------------------------------------------------

#' Make blast Database
#'
#' @param file A fasta file to create a database from
#' @param dbtype Type of database, default is nucleotide
#' @param args Extra arguments passed to BLAST
#'
#' @return
#' @export
#'
#' @examples
makeblastdb <- function (file, dbtype = "nucl", args= NULL) {
  if (is.null(args)){args <- ""}
  results <- system2(command = .findExecutable("makeblastdb"),
                     args = c("-in", file, "-dbtype", dbtype, args),
                     wait = TRUE,
                     stdout = TRUE)
}


# BLAST -------------------------------------------------------------------

#' Run BLAST search
#'
#' @param query Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database
#' @param type type of search to conduct, default 'blastn'
#' @param evalue Minimum evalue from search
#' @param args Extra arguments passed to BLAST
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#'
#' @examples
blast <- function (query, db, type="blastn", evalue = 1e-6, args=NULL, quiet=FALSE){
  
  time <- Sys.time() # get time
  # Create temp files
  tmp <- tempdir()
  tmpquery <- paste0(tmp, "/tmpquery.fa")
  tmpdb <- paste0(tmp, "/tmpquery.fa")
  
  
  # Database
  if(inherits(db, "DNAbin")){
    if (!quiet) { message("Database input is DNAbin: Creating temporary blast database") }
    insect::writeFASTA(db, tmpdb)
    makeblastdb(tmpdb)
    db <- tmpdb
  } else if (inherits(db, "DNAString") | inherits(db, "DNAStringSet")){
    if (!quiet) { message("Database input is DNAStringSet: Creating temporary blast database") }
    Biostrings::writeXStringSet(db, tmpdb)
    makeblastdb(tmpdb)
    db <- tmpdb
  } else if (inherits(db, "character") &&  all(str_to_upper(str_split(db,"")[[1]]) %in% Biostrings::DNA_ALPHABET)) { # Handle text input
    if (!quiet) { message("Database input is character string: Creating temporary blast database") }
    if (nchar(db[1]) == 1) {db <- paste0(db, collapse = "")}
    db <- insect::char2dna(db)
    insect::writeFASTA(db, tmpdb)
    makeblastdb(tmpdb)
    db <- tmpdb
  } else if (inherits(db, "character") &&  file.exists(file.path(db))){ # Handle filename
    db <- db
  }
  
  # Query
  if(inherits(query, "DNAbin")){
    if (!quiet) { message("Query is DNAbin: Creating temporary fasta file") }
    insect::writeFASTA(query, tmpquery)
    input <- tmpquery
  } else if (inherits(query, "DNAString") | inherits(query, "DNAStringSet")){
    if (!quiet) { message("Query is DNAString: Creating temporary fasta file") }
    Biostrings::writeXStringSet(query, tmpquery)
    input <- tmpquery
  }else if (inherits(query, "character") &&  all(str_to_upper(str_split(query,"")[[1]]) %in% Biostrings::DNA_ALPHABET)) { # Handle text query
    if (!quiet) { message("Query is character string: Creating temporary fasta file") }
    if (nchar(query[1]) == 1) {query <- paste0(query, collapse = "")}
    query <- insect::char2dna(query)
    insect::writeFASTA(query, tmpquery)
    input <- tmpquery
  } else if (inherits(query, "character") &&  file.exists(file.path(query))){ # Handle filenames
    input <- query
  }
  
  
  #  Conduct BLAST search
  if (!quiet) { message("Running BLAST") }
  results <- system2(command = .findExecutable(type),
                     args = c("-db", db,
                              "-query", input,
                              "-outfmt 6",
                              "-evalue", evalue,
                              "-ungapped", args),
                     wait = TRUE,
                     stdout = TRUE)
  
  # Parse BLAST results
  out <- results %>%
    enframe() %>%
    separate(col = value,
             into = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
             sep = "\t",
             convert = TRUE)
  time <- Sys.time() - time
  if (!quiet) (message(paste0("finished BLAST in ", format(time, digits = 2))))
  
  if(file.exists(tmpdb)){file.remove(tmpdb)}
  if(file.exists(tmpquery)){file.remove(tmpquery)}
  return(out)
}


# BLAST_top_hit -----------------------------------------------------------

#' BLAST Top Hit
#'
#' @description Conduct BLAST search and return top hit
#' @param query Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database
#' @param type type of search to conduct, default 'blastn'
#' @param threshold Minimum identity threshold to accept
#' @param taxranks The taxonomic ranks contained in the fasta headers
#' @param delim The delimiter between taxonomic ranks in fasta headers
#' @param args Extra arguments passed to BLAST
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#'
#' @examples
blast_top_hit <- function(query, db, type="blastn", threshold=90, taxranks=c("Kingdom", "Phylum","Class", "Order", "Family", "Genus", "Species"), delim=";", args="-max_target_seqs 5", quiet=FALSE ){
  #Conduct BLAST
  result <- blast(query=query, db=db, args=args)
  #Subset to top hit
  top_hit <- result %>%
    filter(pident > threshold) %>%
    group_by(qseqid) %>%
    top_n(1, bitscore) %>%
    top_n(1, row_number(name)) %>% # Break ties by position
    separate(sseqid, c("acc",taxranks), delim)
  
}