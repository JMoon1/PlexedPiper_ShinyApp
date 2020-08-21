library(AnnotationDbi)
library(PlexedPiper)
library(progress)
library(zeallot)
library(plyr)
library(odbc)
library(tidyverse)
library(data.table)

source("shiny_utils.R")

server <- function(input, output, session) {
  
  progress <- reactiveValues(text = "MASIC not downloaded",
                             col = "red",
                             icon = "times-circle",
                             text2 = "MS-GF+ not downloaded",
                             col2 = "red",
                             text3 = "Data not quantified",
                             col3 = "red")
  
  
  
  ## 1st event
  observeEvent(input$submit, {
    
    jr <- PlexedPiper:::get_job_records_by_dataset_package(input$DataPkgNumber)
    max_steps <- length(unique(jr$Dataset))
    
    progress_bar <- shiny::Progress$new(session, min=0, max=max_steps)
    on.exit(progress_bar$close())
    
    # use the main progress outside of llply
    progress_bar$set( value = 0, "Getting MASIC Data...")
    
    if (input$DataPkgNumber == "") {
      showNotification("Please select Data Package ID", type = "warning")
    }
    
    req(input$DataPkgNumber)
    validate(need(input$DataPkgNumber, "Please select a Data Package"))
    
    progress[['text']] <<- "MASIC data being fetched."
    # withCallingHandlers({
    #   shinyjs::html("text", "")
    # withProgress(message = 'Fetching MASIC data', value = 0, {
    
    
    masic_data <<- read_masic_data_from_DMS_shiny(input$DataPkgNumber, 
                                                  interference_score = T, progress_shiny(progress_bar))
    
    
    if (!is.null(masic_data$InterferenceScore)) {
      progress_bar$set( value = 0.5, "Filtering MASIC data")
      masic_data <<- filter_masic_data(masic_data, interference_score_threshold = input$interference_score,
                                       s2n_threshold = input$s2n)
    }
    
    # masic_progress[[1]] <<- "MASIC filtered"
    progress[['text']] <<- "MASIC Processed" # do i really need <<- here?
    progress[['col']] <<- "green"
    
  })
  
  output$ibox <- renderInfoBox({
    infoBox("MASIC", value = progress[['text']], color = progress[["col"]])
  })
  
  output$masic_out <- renderDataTable({
    masic_data
  })
  
  observeEvent(input$submit2, {
    
    if (input$DataPkgNumber == "") {
      showNotification("Please select Data Package ID", type = "warning")
    }
    
    req(input$DataPkgNumber)
    validate(
      need(input$DataPkgNumber, "Please select a Data Package")
    )
    
    progress[['text2']] <<- "MS-GF+ data being fetched."
    
    jr <- PlexedPiper:::get_job_records_by_dataset_package(input$DataPkgNumber)
    max_steps <- length(unique(jr$Dataset))
    org_name <- gsub("_", " ", jr$Organism)
    ## TODO if org_name isn't in list of supported species, throw error
    
    progress_bar <- shiny::Progress$new(session, min=0, max=max_steps)
    
    # use the main progress outside of llply
    progress_bar$set( value = 0, "Getting MS-GF+ Data...")
    
    msnid <<- read_msms_data_from_DMS_shiny(input$DataPkgNumber, progress = progress_shiny(progress_bar))
    
    ## TODO close progress bar
    ## open new progress bar
    
    # withProgress(message = 'Processing MS-GF+ data', value = 0, {
      
      msnid <<- correct_peak_selection(msnid)
      
      progress_bar$close()
      
      ## start new progress bar
      # progress_bar <- shiny::Progress$new(session, min=0, max=1)
      # on.exit(progress_bar$close())
      # incProgress(amount = 0.25, message = "Filtering unique peptide ID FDR")
      progress_bar <- shiny::Progress$new(session, min=0, max=1)
      on.exit(progress_bar$close())
      
      progress_bar$set(value = 0.25, "Filtering unique peptide ID FDR")
      msnid <<- filter_msgf_data_peptide_level(msnid, input$pep_fdr)
      
      if (input$remap_genes == T) {
        progress_bar$set(value = 0.4, "Remapping to genes")
        
        msnid <<- remap_accessions_refseq_to_gene(msnid,
                                                  organism_name=org_name)
        path_to_FASTA <- path_to_FASTA_used_by_DMS(input$DataPkgNumber)
        path_to_FASTA <- remap_accessions_refseq_to_gene_fasta(
          path_to_FASTA,
          organism_name = org_name)
      }
      else {
        path_to_FASTA <- path_to_FASTA_used_by_DMS(input$DataPkgNumber)
      }
      # incProgress(amount = 0.5, message = "Filtering accession-level FDR")
      progress_bar$set(value = 0.5, "Filtering accession-level FDR")
      msnid <<- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
      msnid <<- filter_msgf_data_protein_level(msnid, input$prot_fdr)
      
      # incProgress(amount = 0.75, message = "Performing parsimonious inference")
      progress_bar$set(value = 0.75, "Performing parsimonious inference")
      
      if (input$protein_inference == T) {
        msnid <<- infer_parsimonious_accessions(msnid, unique_only=input$unique_only)
      }
      
      msnid <<- apply_filter(msnid, "!isDecoy")
      ## TODO remove contaminants?
      # incProgress(amount = 100, message = "MS/MS data complete")
      progress_bar$set(value = 1, "MS/MS data processing complete")
    # })
    
    
    progress[['text2']] <<- "MS-GF+ data processed." 
    progress[['col2']] <<- "green"
    
  })
  
  output$ibox2 <- renderInfoBox({
    infoBox("MS-GF+", value = progress[['text2']], color = progress[['col2']], icon = icon("list"))
  })
  
  output$msms_out <- renderDataTable({
    msnid
  })
  
  observeEvent(input$submit3, {
    
    if (!exists("masic_data") | !exists("msnid")) {
      showNotification("Please download both MASIC and MS-GF+ data before this step.",
                       type = "error")
    }
    
    validate(
      need(exists('masic_data'), "Please download MASIC")
    )
    validate(
      need(exists('msnid'), "Please download MS-GF+")
    )
    
    
    progress[['text3']] <<- "Data not quantified"
    
    withProgress(message = 'Creating crosstab', value = 0, {
      
      samples <- data.frame()
      fractions <- data.frame()
      references <- data.frame()
      
      c(samples, fractions, references) %<-% PlexedPiper:::get_study_design_by_dataset_package(input$DataPkgNumber)
      
      aggregation_level <- c("accession") # hard coded for now
      incProgress(amount = 0.25, message = "Rolling up to accession")
      crosstab <<- PlexedPiper::create_crosstab(msnid, masic_data, aggregation_level,
                                                fractions, samples, references)
      
      incProgress(amount = 1, message = "Crosstab created")
    })
    
    progress[['text3']] <<- "Data quantified and ready for download" # do i really need <<- here?
    progress[['col3']] <<- "green"
    
  })
  
  output$ibox3 <- renderInfoBox({
    infoBox("Quant", value = progress[['text3']], color = progress[['col3']], icon = icon("file-alt"))
  })
  
  output$crosstab_out <- renderDataTable({
    if (input$submit3 == 0)
      return()
    
    validate(
      need(exists('crosstab'), "Crosstab doesn't exist")
    )
    isolate(
      as.data.frame(crosstab)
    )
  })
  
  ## handle downloading data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("crosstab", ".txt")
    },
    content = function(file) {
      write.table(signif(crosstab, 3), file, sep = "\t", col.names = NA, quote = F, 
                  na = "")
    }
  )
  
}
