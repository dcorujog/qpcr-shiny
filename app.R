library(shiny)

#### Source scripts ####

source("qPCR_analysis_shiny.R")
source("/Users/dcorujo/Documents/Lab/Scripts/RScripts/theme_paper.R")

#### UI ####
ui <- fluidPage(
  
  # App title ----
  titlePanel("qPCR analysis"),
  
  # Layout ----
  sidebarLayout(
    # Sidebar panel ----
    sidebarPanel(
      # Input file ----
      fileInput("input_file", label = "Upload data", buttonLabel = "Choose file"),
      
      # Input design ----
      fileInput("input_design", label = "Upload design", buttonLabel = "Choose file"),
      
      # qPCR function arguments ----
      textInput("exp_name", "Experiment name", "my_experiment"),
      textInput("calibsample", "Reference sample"),
      textInput("hkg", "Normalization genes (comma separated, no spaces)", "GAPDH,RPLP0"),
      checkboxGroupInput("check_args", "Other arguments",
                         choices = c("Fix names" = "fix_names",
                                     "Exclude" = "exclude")),
      
      # Action button to run analysis ----
      actionButton("run_button", "Run analysis"),
      
      # Download button ----
      downloadButton("download_data", "Download results"),
    
      # Download button for QC plots ----
      downloadButton("download_qc", "Download QC plots")),
    
      # Main panel ----
      # Tabset with Results table, QC plots and (to-do) plot of the results
      mainPanel(
        tabsetPanel(
          tabPanel("Results table", dataTableOutput("output_table")),
          tabPanel("QC", 
                   fluidRow(column(6, plotOutput("hkg_scatter")),
                            column(6, plotOutput("ct_sd_qc"))),
                   fluidRow(column(6, plotOutput("ct_avg_qc")),
                            column(6, plotOutput("eff_avg_qc"))))
          )
        )
      )
)
plotOutput("hkg_scatter", width = 300, height = 300)

#### SERVER ####

server <- function(input, output) {
  
  # Save input file as a data frame ----
  input_data <- reactive({
    req(input$input_file)
    input_file <- input$input_file
    return(read.delim(input_file$datapath))
    })
  
  # Save design table as data frame ----
  design_table <- reactive({
    req(input$input_design)
    design_file <- input$input_design
    return(read.delim(design_file$datapath))
    })
  
  # Render input file as table ----
  output$input_file <- renderTable(input_data())
  
  # Call qPCR analysis function with arguments supplied as input by the user when the run button is clicked ----
  # Plots and download handler also depend on the run button
  observeEvent(input$run_button, {
    output_data <- qpcr_analysis(ct_data = input_data(),
                                 design_table = design_table(),
                                 calibsample = input$calibsample,
                                 hkg = unlist(strsplit(input$hkg, ",")),
                                 exp_name = input$exp_name,
                                 fix_names = "fix_names" %in% input$check_args,
                                 exclude = "exclude" %in% input$check_args
                                )
    # Render output data as table ----
    output$output_table <- renderDataTable(output_data$norm_data)
    
    # Render QC plots ----
    output$hkg_scatter <- renderPlot(output_data$hkg_scatter)
    output$ct_sd_qc <- renderPlot(output_data$ct_sd_qc)
    output$ct_avg_qc <- renderPlot(output_data$ct_avg_qc)
    output$eff_avg_qc <- renderPlot(output_data$eff_avg_qc)
    
    # Save QC plots as PDF ----
    ggsave("hkg_scatter.pdf", plot = output_data$hkg_scatter)
    ggsave("ct_sd.pdf", plot = output_data$ct_sd_qc)
    ggsave("ct_avg.pdf", plot = output_data$ct_avg)
    ggsave("eff_avg.pdf", plot = output_data$eff_avg)
    qc_pdfs <- c("hkg_scatter.pdf", "ct_sd.pdf", "ct_avg.pdf", "eff_avg.pdf")
    
    # Download handler for QC plots ----
    output$download_qc <- downloadHandler(
      filename = function() { 
        paste(input$exp_name, '_qc_plots.zip', sep='') 
      },
      content = function(file) {
        zip(zipfile = file,
            files = qc_pdfs)
      },
      contentType = "application/zip"
      )
    
    # Download handler for output data ----
    output$download_data <- downloadHandler(filename = function() {
        paste(input$exp_name, "_norm_data.csv", sep = "")
      },
      content = function(file) {
        write.csv(output_data$norm_data, file, row.names = FALSE)
      }
    )
    
  })
}

shinyApp(ui = ui, server = server)