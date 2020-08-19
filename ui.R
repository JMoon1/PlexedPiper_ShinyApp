library(shinydashboard)

ui <- dashboardPage(
  
  dashboardHeader(title = "PlexedPiper"),
  
  # Sidebar panel ----
  dashboardSidebar(
    
    textInput(inputId = "DataPkgNumber",
              label = "Data Package ID:", 
              value = "3606"),
    fluidRow(
      column(2,
             downloadButton("downloadData", "Download"), style = "padding-left:30px;")),
    tags$style(".skin-blue .sidebar a { color: #444; }")
  ),
  # Main panel for displaying outputs ----
  dashboardBody(
    tags$head(
      tags$style(HTML("hr {border-top: 1px solid #000000;}")) # make hr() visible
    ),
    tags$head(
      tags$style("label{font-size: 15px}")
    ),
    fluidRow(
      
      # tags$div(class = "row",
      #          tags$div(class = "span")),
      column(5, uiOutput("ibox")),
      column(2, textInput(inputId = "s2n",
                          "S/N Ratio:",
                          value = 0),
             textInput(inputId = "interference_score",
                       "Interference Score:",
                       value = 0.5), style = 'padding-left:0px; padding-right:50px; padding-top:0px; padding-bottom:0px'),
      column(2, actionButton("submit", "MASIC"))
    ),
    br(),
    hr(),
    br(),
    fluidRow(
      
      column(5, uiOutput("ibox2")),
      column(2, 
             # selectInput("org_name", "Organism:", choices = list("Rattus norvegicus" = 1,
             #                                                     "Homo sapiens" = 2)),
             textInput(inputId = "pep_fdr",
                       label = "Peptide-level FDR:",
                       value = 0.01),
             textInput("prot_fdr",
                       label = "Accession-level FDR:",
                       value = 0.01),
             checkboxInput("remap_genes",
                           "Remap to Genes"),
             checkboxInput("protein_inference",
                           "Parsimonious Inference"), 
             checkboxInput("unique_only", "Drop Shared Peptides"),
             style = 'padding-left:0px; padding-right:50px; padding-top:0px; padding-bottom:0px'
      ),
      column(2, actionButton("submit2", "MS-GF+"))
      
    ),
    br(),
    hr(),
    br(),
    fluidRow(
      column(5, uiOutput("ibox3")),
      column(2, selectInput(inputId = "aggregate_to",
                            label = "Roll up to:",
                            choices = list("Accession" = 1, 
                                           "Peptide" = 2),
                            selected = 1), style = 'padding-left:0px; padding-right:50px; padding-top:0px; padding-bottom:0px'),
      column(2, actionButton("submit3", "Quant"))
    ),
    fluidRow(
      textOutput('textout')
    )
  )
)
