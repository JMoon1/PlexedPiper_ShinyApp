library(shiny)

source("ui.R", local = T)
source("server.R")

shinyApp(
  ui = ui,
  server = server
)

