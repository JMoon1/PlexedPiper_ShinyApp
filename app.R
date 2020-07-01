library(shiny)
library(DT)
# source("PNNL_DMS_utils.R")
library(PlexedPiper)
library(zeallot)
library(odbc)
library(shinycssloaders)

source("PNNL_DMS_utils.R")

# umount_cmd <- sprintf("umount %s", "~/temp_msms_results")
# system(umount_cmd)
# unlink("~/temp_msms_results", recursive = T)


# 3 components
ui <- navbarPage("PlexedPiper",
                   tabPanel("Obtain Data",
                            # Sidebar panel ----
                            sidebarPanel(
                              
                              textInput(inputId = "DataPkgNumber",
                                        label = "Data Package ID:", 
                                        value = "3606 or 3442"),
                              
                              actionButton("submit", "Submit"),
                              downloadButton("downloadData", "Download")
                            ),
                            # Main panel for displaying outputs ----
                            mainPanel(
                              dataTableOutput("out_crosstab") %>% 
                                withSpinner(color="#34b1eb")
                              ))
                   
                   # tabPanel("MS-GF+",
                   #          dataTableOutput(outputId = "msgf_data"))
                   

                   
)


server <- function(input, output, session) {
  crosstab <- eventReactive(input$submit, {
    
    ## MS/MS
    # get_job_records_by_dataset_package(input$DataPkgNumber)
    msnid <- read_msms_data_from_DMS(input$DataPkgNumber)
    msnid <- correct_peak_selection(msnid)
    msnid <- filter_msgf_data_peptide_level(msnid, 0.01)

    path_to_FASTA <- path_to_FASTA_used_by_DMS(input$DataPkgNumber)
    msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
    msnid <- filter_msgf_data_protein_level(msnid, 0.01)
    
    # msnid <- infer_parsimonious_accessions(msnid, unique_only=TRUE)
    
    msnid <- apply_filter(msnid, "!isDecoy")
    psms(msnid)$isContam <- grepl("^Contam", psms(msnid)$accession)
    msnid <- apply_filter(msnid, "!isContam")
    
    
    ## MASIC
    masic_data <- read_masic_data_from_DMS(input$DataPkgNumber, 
                                           interference_score = T)
    masic_data <- filter_masic_data(masic_data, s2n_threshold = 0, 
                                    interference_score_threshold = 0.5)
    
    ## QUANT
    c(samples, fractions, references) %<-% get_study_design_by_dataset_package(input$DataPkgNumber)
    
    aggregation_level <- c("accession")
    create_crosstab(msnid, 
                    masic_data, 
                    aggregation_level, 
                    fractions, samples, references)
    
  }) 

  ## render crosstab
  output$out_crosstab <- renderDataTable({crosstab()})
  
  ## handle downloading data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("crosstab", ".txt")
    },
    content = function(file) {
      write.table(crosstab(), file, sep = "\t", col.names = NA, quote = F)
    }
  )
}


shinyApp(ui, server)
