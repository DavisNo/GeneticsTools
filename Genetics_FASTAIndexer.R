# Load libraries
library(shiny)
library(shinythemes)
library(DT)
library(Rsamtools)
library(Biostrings)
library(dplyr)

# Set the maximum upload file size to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)

# Define UI
ui <- fluidPage(

  tags$head(
    tags$style(HTML("
      body, h1, h2, h3, h4, h5, h6, p, label, input, textarea, select, button, table {
        color: black !important;
      }
    "))
  ),
  
  titlePanel("FASTA Indexing w/ Rsamtools"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_file", "Upload FASTA File (.fasta, .fa):",
                accept = c(".fasta", ".fa", ".txt")),
      
      actionButton("indexBtn", "Index FASTA"),
      
      br(), br(),
      
      h5("Indexing Log"),
      verbatimTextOutput("logOutput"),
      
      hr(),
      
      downloadButton("downloadIndex", "Download Indexed File (.fai)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          "FASTA Sequences",
          h4("Uploaded Sequences"),
          DT::dataTableOutput("fastaTable")
        ),
        tabPanel(
          "FAI Index",
          h4("Parsed .fai Contents"),
          tableOutput("faiTable"),
          br(),
          verbatimTextOutput("summaryInfo")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Store the path of the uploaded FASTA
  fasta_path <- reactiveVal(NULL)
  
  # Store indexing messages or errors
  process_log <- reactiveVal("")
  
  # Store the FASTA sequences
  fasta_data <- reactiveVal(NULL)
  
  # Store the .fai data
  fai_data <- reactiveVal(NULL)
  
  # Observe the file input
  observeEvent(input$fasta_file, {
    req(input$fasta_file)
    
    #Check extension
    valid_ext <- c("fa", "fasta", "txt")
    uploaded_ext <- tools::file_ext(input$fasta_file$name)
    if (!uploaded_ext %in% valid_ext) {
      showNotification(
        "Unsupported file extension. Please upload a .fa, .fasta, or .txt file.",
        type = "error"
      )
      return(NULL)
    }
    
    # Save
    fasta_path(input$fasta_file$datapath)
    
    process_log("")
    fasta_data(NULL)
    fai_data(NULL)
  })
  
  # Read and display the FASTA sequences in a table
  output$fastaTable <- DT::renderDataTable({
    req(fasta_path())
    
    tryCatch({
      # Read the FASTA file
      fasta <- readDNAStringSet(fasta_path())

      df <- data.frame(
        Gene = names(fasta),
        Sequence = as.character(fasta),
        stringsAsFactors = FALSE
      )
      
      # Update reactive value
      fasta_data(df)
      
      # Render as a data table
      DT::datatable(
        df,
        options = list(
          pageLength = 10,
          lengthMenu = c(10, 25, 50, 100),
          searchHighlight = TRUE,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    }, error = function(e) {
      showNotification("Error reading FASTA file. Please ensure it is properly formatted.", type = "error")
      return(NULL)
    })
  })
  
  # Observe the "Index FASTA" button
  observeEvent(input$indexBtn, {
    req(fasta_path())
  
    process_log("Indexing started...\n")
    
    # Index with Rsamtools
    tryCatch({
      indexFa(fasta_path())
      
      # Read the .fai file
      fai_file <- paste0(fasta_path(), ".fai")
      fai_df <- read.table(fai_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      colnames(fai_df) <- c("Sequence", "Length", "Offset", "LineBases", "LineWidth")
      
      # Update reactive value with indexed data
      fai_data(fai_df)
      
      # Update log
      process_log(paste0(process_log(), "Indexing completed successfully!\n"))
    }, error = function(e) {
      # Log errors
      process_log(paste0(process_log(), "Error during indexing: ", e$message, "\n"))
      showNotification("Error during indexing. See log for details.", type = "error")
    })
  })
  
  # Display the .fai file contents
  output$faiTable <- renderTable({
    req(fai_data())
    fai_data()
  })
  
  # Provide summary information about the indexing
  output$summaryInfo <- renderText({
    req(fai_data(), fasta_data())
    
    # Total sequences in the uploaded FASTA
    total_sequences <- nrow(fasta_data())
    
    # Total sequences indexed
    indexed_sequences <- nrow(fai_data())
    
    # Determine if any sequences failed
    failed_sequences <- total_sequences - indexed_sequences
    if (failed_sequences > 0) {
      paste(
        "Summary Information:\n",
        "Total Sequences in FASTA: ", total_sequences, "\n",
        "Sequences Indexed: ", indexed_sequences, "\n",
        "Sequences Failed to Index: ", failed_sequences, "\n"
      )
    } else {
      paste(
        "Summary Information:\n",
        "Total Sequences in FASTA: ", total_sequences, "\n",
        "Sequences Indexed: ", indexed_sequences, "\n",
        "All sequences were successfully indexed."
      )
    }
  })
  
  # Show log
  output$logOutput <- renderText({
    process_log()
  })
  
  # Download handler for the indexed .fai file
  output$downloadIndex <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$fasta_file$name), ".fai")
    },
    content = function(file) {
      req(fai_data())
      
      # Write the .fai file to the download path
      fai_file <- paste0(fasta_path(), ".fai")
      file.copy(fai_file, file)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)
