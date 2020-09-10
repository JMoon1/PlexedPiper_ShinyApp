
filter_masic_data <- function (x, interference_score_threshold = 0.9, s2n_threshold = 4) 
{
  x <- x %>% filter(InterferenceScore >= interference_score_threshold)
  selected <- x %>% dplyr::select(Dataset, ScanNumber, contains("SignalToNoise")) %>% 
    gather(channel, s2n, -c(Dataset, ScanNumber)) %>% mutate(s2n = ifelse(is.na(s2n), 
                                                                          0, s2n)) %>% filter(s2n >= s2n_threshold) %>% dplyr::select(-s2n) %>% 
    mutate(channel = sub("_SignalToNoise", "", channel))
  x <- x %>% dplyr::select(Dataset, ScanNumber, starts_with("Ion"), 
                           -contains("SignalToNoise")) %>% gather(channel, intensity, 
                                                                  -c(Dataset, ScanNumber)) %>% inner_join(selected, by = c("Dataset", 
                                                                                                                           "ScanNumber", "channel")) %>% spread(channel, intensity)
}

get_driver <- function(){
  if(.Platform$OS.type == "unix"){
    return("FreeTDS")
  }else if(.Platform$OS.type == "windows"){
    return("SQL Server")
  }else{
    stop("Unknown OS type.")
  }
}

get_auth <- function(){
  if(.Platform$OS.type == "unix"){
    return("PORT=1433;UID=dmsreader;PWD=dms4fun;")
  }else if(.Platform$OS.type == "windows"){
    return("")
  }else{
    stop("Unknown OS type.")
  }
}

## util tools ==================================================================
tool2suffix <- list("MSGFPlus"="_msgfplus_syn.txt",
                    "MSGFPlus_MzML"="_msgfplus_syn.txt",
                    "MSGFPlus_DTARefinery"="_msgfplus_syn.txt",
                    "MSGFDB_DTARefinery"="_msgfdb_syn.txt",
                    "MASIC_Finnigan"="_ReporterIons.txt",
                    "TopPIC" = "_TopPIC_PrSMs.txt")

progress_shiny <-function (progress, step = 1){
  list(
    init = function(n){},
    step = function() {
      progress$set( progress$getValue() + step )
    }, 
    term = function(){}
  )
}

get_results_for_multiple_jobs.dt.shiny <- function(jobRecords, progress){
  toolName = unique(jobRecords[["Tool"]])
  if (length(toolName) > 1){
    stop("Contains results of more then one tool.")
  }
  results = llply(jobRecords[["Folder"]],
                  get_results_for_single_job.dt.http,
                  fileNamePttrn=tool2suffix[[toolName]],
                  .progress = progress)
  results.dt <- rbindlist(results)
  return( as.data.frame(results.dt) ) # in the future I may keep it as data.table
}

get_results_for_single_job.dt.http <- function(pathToFile, fileNamePttrn){
  
  pathToFile <- as.character(pathToFile)
  
  remote_folder <- gsub("\\\\","/",pathToFile)
  dataset <- sub("^//.*/.*/.*/(.*)/.*", "\\1", remote_folder)
  url_string <- paste0("http:", remote_folder)
  
  pathToFile <- paste0(url_string, "/", dataset, fileNamePttrn)
  
  if(length(pathToFile) == 0){
    stop("can't find the results file")
  }
  if(length(pathToFile) > 1){
    stop("ambiguous results files")
  }
  
  results <- read_tsv(pathToFile, col_types=readr::cols(), progress=FALSE)
  
  out <- data.table(Dataset=dataset, results)
  return(out)
}



read_masic_data_from_DMS_shiny <- function (data_pkg, interference_score = FALSE, progress) 
{
  if (length(data_pkg) > 1) {
    job_rec_ls <- lapply(data_pkg, get_job_records_by_dataset_package)
    jobRecords <- Reduce(rbind, job_rec_ls)
  }
  else {
    jobRecords <- get_job_records_by_dataset_package(data_pkg)
  }
  jobRecords <- jobRecords[grepl("MASIC", jobRecords$Tool), 
                           ]
  masicData <- get_results_for_multiple_jobs.dt.shiny(jobRecords, progress)
  if (interference_score) {
    results = llply(jobRecords[["Folder"]], get_results_for_single_job.dt.http, 
                    fileNamePttrn = "_SICstats.txt", .progress = progress)
    results.dt <- rbindlist(results)
    masicStats <- as.data.frame(results.dt)
    masicStats <- masicStats[, -2]
    masicData <- masicData[, -2]
    x <- select(masicData, Dataset, ScanNumber, starts_with("Ion"), 
                -contains("Resolution"))
    y <- select(masicStats, Dataset, ScanNumber = FragScanNumber, 
                contains("InterferenceScore"))
    z <- inner_join(x, y)
    return(z)
  }
  masicData <- masicData[, -2]
  masicData <- select(masicData, Dataset, ScanNumber, starts_with("Ion"), 
                      -contains("Resolution"))
  return(masicData)
}

read_msms_data_from_DMS_shiny <- function(DataPkgNumber, progress) 
{
  msnid <- MSnID(".")
  if (!is.null(DataPkgNumber)) {
    if (length(DataPkgNumber) > 1) {
      job_rec_ls <- lapply(DataPkgNumber, get_job_records_by_dataset_package)
      jobRecords <- Reduce(rbind, job_rec_ls)
    }
    else {
      jobRecords <- get_job_records_by_dataset_package(DataPkgNumber)
    }
    jobRecords <- jobRecords[grepl("MSGFPlus", jobRecords$Tool), 
                             ]
    x <- get_results_for_multiple_jobs.dt.shiny(jobRecords, progress)
    x <- x %>% mutate(accession = Protein, calculatedMassToCharge = (MH + 
                                                                       (Charge - 1) * MSnID:::.PROTON_MASS)/Charge, chargeState = Charge, 
                      experimentalMassToCharge = PrecursorMZ, isDecoy = grepl("^XXX", 
                                                                              Protein), peptide = Peptide, spectrumFile = Dataset, 
                      spectrumID = Scan)
    x <- mutate(x, pepSeq = MSnID:::.get_clean_peptide_sequence(peptide))
    psms(msnid) <- x
    return(msnid)
  }
}

path_to_FASTA_used_by_DMS.http <- function(data_package_number){
  
  # make sure it was the same fasta used for all msgf jobs
  # at this point this works only with one data package at a time
  jobRecords <- get_job_records_by_dataset_package(data_package_number)
  jobRecords <- jobRecords[grepl("MSGFPlus", jobRecords$Tool),]
  if(length(unique(jobRecords$`Organism DB`)) != 1){
    stop("There should be exactly one FASTA file per data package!")
  }
  
  strSQL <- sprintf("Select [Organism DB],
                             [Organism DB Storage Path]
                     From V_Analysis_Job_Detail_Report_2
                     Where JobNum = %s", jobRecords$Job[1])
  
  con_str <- sprintf("DRIVER={%s};SERVER=gigasax;DATABASE=dms5;%s",
                     get_driver(),
                     get_auth())
  
  con <- dbConnect(odbc(), .connection_string=con_str)
  qry <- dbSendQuery(con, strSQL)
  res <- dbFetch(qry)
  dbClearResult(qry)
  dbDisconnect(con)
  
  remote_folder <- gsub("\\\\","/", res['Organism DB Storage Path'])
  dataset <- sub("^//.*/.*/.*/(.*)/.*", "\\1", remote_folder)
  url_string <- paste0("http:", remote_folder)
  
  path_to_FASTA <- file.path(url_string, res['Organism DB'])
  
  return(path_to_FASTA)
}


get_study_design_by_dataset_package.http <- function(dataPkgNumber) {
  
  con_str <- sprintf("DRIVER={%s};SERVER=gigasax;DATABASE=DMS_Data_Package;%s",
                     get_driver(),
                     get_auth())
  con <- dbConnect(odbc(), .connection_string=con_str)
  
  ## fetch table with path to DataPackage
  strSQL <- sprintf("
                    SELECT *
                    FROM V_Data_Package_Detail_Report
                    WHERE ID = %s",
                    dataPkgNumber)
  qry <- dbSendQuery(con, strSQL)
  dataPkgReport <- dbFetch(qry)
  dbClearResult(qry)

  remote_folder <- gsub("\\\\","/", dataPkgReport$`Share Path`)
  remote_folder <- paste0("http:", remote_folder)

  ## fetch samples.txt
  samples <- readr::read_tsv(paste0(remote_folder, "/samples.txt"), 
                             col_types=readr::cols(), progress=FALSE)
  if (!setequal(colnames(samples), c("PlexID",
                                     "QuantBlock",
                                     "ReporterAlias",
                                     "ReporterName",
                                     "MeasurementName"))) {
    stop("There are incorrect column names or missing columns in the 'samples'
         study design table.")
  }
  
  ## fetch fractions.txt
  fractions <- readr::read_tsv(paste0(remote_folder, "/fractions.txt"), 
                             col_types=readr::cols(), progress=FALSE)
  if (!setequal(colnames(fractions), c("PlexID",
                                       "Dataset"))) {
    stop("There are incorrect column names or missing columns in the 'fractions'
         study design table.")
  }
  
  ## fetch references.txt
  references <- readr::read_tsv(paste0(remote_folder, "/references.txt"), 
                             col_types=readr::cols(), progress=FALSE)
  if (!setequal(colnames(references), c("PlexID",
                                        "QuantBlock",
                                        "Reference"))) {
    stop("There are incorrect column names or missing columns in the 'references'
         study design table.")
  }

  study_des <- list(samples = samples,
                    fractions = fractions,
                    references = references)
  return(study_des)
}

remap_accessions_refseq_to_gene_fasta_shiny <- function(path_to_FASTA,
                                                        organism_name,
                                                        conversion_table){
  
  is_compressed <- FALSE
  if(grepl("[.]gz$", path_to_FASTA)){
    is_compressed <- TRUE
  }else if(grepl("[.](bz2|xz|zip)$", path_to_FASTA)){
    stop("The only supported compression is gzip!")
  }
  
  mySequences <- Biostrings::readAAStringSet(path_to_FASTA)
  names(mySequences) <- sub("^(.P_\\d+)(\\.\\d+)?\\s.*","\\1",names(mySequences))
  prot_lengths <- data.frame(REFSEQ = names(mySequences),
                             Length = width(mySequences),
                             stringsAsFactors = FALSE)
  
  if(missing(conversion_table)){
    ah <- AnnotationHub()
    orgs <- subset(ah, ah$rdataclass == "OrgDb")
    db <- query(orgs, organism_name)
    db <- db[[1]]
    conversion_table <- AnnotationDbi::select(db,
                                              keys=names(mySequences),
                                              columns="SYMBOL",
                                              keytype="REFSEQ")
    conversion_table <- conversion_table[,c("REFSEQ","SYMBOL")]
  }
  # I don't know what column names are, but I require that refseq is first,
  # gene second
  colnames(conversion_table) <- c("REFSEQ","SYMBOL")
  # resolving mapping ambiguity by selecting longest RefSeq per gene
  conversion_table <- conversion_table %>%
    filter(!is.na(SYMBOL)) %>%
    inner_join(prot_lengths, by = "REFSEQ") %>%
    group_by(SYMBOL) %>%
    dplyr::slice(which.max(Length))
  conversion_vec <- conversion_table %>% pull(2)
  names(conversion_vec) <- conversion_table %>% pull(1)
  mySequences <- mySequences[names(conversion_vec)]
  names(mySequences) <- conversion_vec[names(mySequences)]
  
  file_no_ext <- tools::file_path_sans_ext(path_to_FASTA, compression=TRUE)
  ext <- sub(file_no_ext, "", path_to_FASTA, fixed=TRUE)
  
  path_to_FASTA_gene <- paste0(file_no_ext, '_gene', ext)
  
  path_to_FASTA_gene <- sub(".*/(.*)$", "\\1", path_to_FASTA_gene)
  Biostrings::writeXStringSet(mySequences, path_to_FASTA_gene, compress = is_compressed)
  return(path_to_FASTA_gene)
}
