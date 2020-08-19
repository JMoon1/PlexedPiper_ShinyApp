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
                  PlexedPiper:::get_results_for_single_job.dt,
                  fileNamePttrn=tool2suffix[[toolName]],
                  .progress = progress)
  results.dt <- rbindlist(results)
  return( as.data.frame(results.dt) ) # in the future I may keep it as data.table
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
    results = llply(jobRecords[["Folder"]], get_results_for_single_job.dt, 
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