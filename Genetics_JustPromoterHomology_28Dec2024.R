# Load libraries
library(shiny)
library(shinythemes)
library(DT)           
library(Biostrings)   
library(dplyr)        
library(stringr)      

# Set the maximum upload file size to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  
  tags$head(
    tags$style(HTML("
      body, h1, h2, h3, h4, h5, h6, p, label, input, textarea, select, button, table {
        color: black !important;
      }
    "))
  ),
  
  titlePanel("Homology Target Search on Uploaded FASTA Sequence"),
  
  sidebarLayout(
    sidebarPanel(
      # Input: File upload for FASTA sequences
      fileInput("fasta_file", "Upload FASTA File:",
                accept = c(".txt", ".fasta", ".fa")),
      
      hr(),
      
      # Input: Text area for the target promoter sequence
      textAreaInput("target_seq", "Input Target Promoter Sequence:",
                    placeholder = "Enter raw DNA sequence here...",
                    rows = 3),
      
      # Input: Number of allowed mismatches
      numericInput("max_mismatches", "Allowed Mismatches:",
                   value = 0, min = 0, max = 5, step = 1),
      
      # Button to run homology search
      actionButton("searchBtn", "Search for Target Sequence"),
      
      hr(),
      
      # Button for the homology search results
      downloadButton("downloadHomology", "Download Homology Results as CSV")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "FASTA Sequences",
          h4("Uploaded FASTA Sequences"),
          DT::dataTableOutput("fastaTable")
        ),
        tabPanel(
          "Homology Results",
          h4("Occurrences of Target Sequence in Each FASTA Sequence"),
          DT::dataTableOutput("homologyTable")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression to read the uploaded FASTA file
  fasta_data <- reactive({
    req(input$fasta_file)  # Ensure a file is uploaded
    
    ext <- tools::file_ext(input$fasta_file$name)
    if (!ext %in% c("txt", "fasta", "fa")) {
      showNotification("Unsupported file format. Please upload a TXT or FASTA file.", type = "error")
      return(NULL)
    }
    
    # Read the FASTA file
    tryCatch({
      fasta <- readDNAStringSet(input$fasta_file$datapath)
      if (length(fasta) == 0) {
        showNotification("FASTA file is empty or improperly formatted.", type = "error")
        return(NULL)
      }
      # Create a data frame with sequence names and sequences
      df <- data.frame(
        Header = names(fasta),
        Sequence = as.character(fasta),
        stringsAsFactors = FALSE
      )
      # Extract Gene Name from Header
      df$`Gene Name` <- sapply(df$Header, function(x) {
        # Split by space or pipe and take the first element
        str_split(x, "[|\\s]+")[[1]][1]
      })
      return(df)
    }, error = function(e) {
      showNotification("Error reading FASTA file. Please ensure it is properly formatted.", type = "error")
      return(NULL)
    })
  })
  
  # Display uploaded FASTA sequences in a table
  output$fastaTable <- DT::renderDataTable({
    req(fasta_data())
    df <- fasta_data()
    
    # Reorder columns
    df <- df %>% dplyr::select(`Gene Name`, Header, Sequence)
    
    DT::datatable(
      df,
      options = list(
        pageLength = 10,
        lengthMenu = c(10, 25, 50, 100),
        searchHighlight = TRUE,
        scrollX = TRUE
      ),
      rownames = FALSE,
      colnames = c("Gene Name", "Header", "Sequence")
    )
  })
  
  # Reactive value to store homology search results
  homology_res <- reactiveVal(NULL)
  
  # Reactive value to store match positions
  match_positions <- reactiveVal(NULL)
  
  # "Search for Target Sequence" button
  observeEvent(input$searchBtn, {
    req(fasta_data(), input$target_seq, input$max_mismatches)
    
    # Clean target promoter sequence
    target_clean <- toupper(gsub("[^ATCG]", "", input$target_seq))
    if (nchar(target_clean) == 0) {
      showNotification(
        "Invalid target promoter sequence. Please enter a valid DNA sequence containing only A, T, C, G.",
        type = "error"
      )
      return(NULL)
    }
    
    # Ensure that the number of mismatches is a non-negative integer
    max_mismatches <- input$max_mismatches
    if (!is.numeric(max_mismatches) || max_mismatches < 0 || max_mismatches != floor(max_mismatches)) {
      showNotification("Allowed mismatches must be a non-negative integer.", type = "error")
      return(NULL)
    }
    
    # Prepare the target as a DNAString
    target_dna <- Biostrings::DNAString(target_clean)
    # Use nchar() for width()
    target_length <- nchar(target_clean)
    
    # Get FASTA data
    fasta_df <- fasta_data()
    
    # Results data frame with Matched Snippets column
    homology_results <- data.frame(
      `Gene Name` = character(),
      Header = character(),
      Occurrences = numeric(),
      `Matched Snippets` = character(),
      stringsAsFactors = FALSE
    )
    
    # Initialize a list to store match positions
    matches_list <- list()
    
    # Progress bar
    progress <- Progress$new(session, min = 0, max = nrow(fasta_df))
    on.exit(progress$close())
    progress$set(message = "Performing search...", value = 0)
    
    # Perform (substring) search with allowed mismatches
    for (i in 1:nrow(fasta_df)) {
      gene_name <- fasta_df$`Gene Name`[i]
      header <- fasta_df$Header[i]
      seq_text <- fasta_df$Sequence[i]
      
      seq_clean <- toupper(gsub("[^ATCG]", "", seq_text))
      
      if (nchar(seq_clean) < target_length) {
        count_hits <- 0
        match_positions_seq <- integer(0)
        snippets <- NA_character_
      } else {
        gene_dna <- Biostrings::DNAString(seq_clean)
        matches <- Biostrings::matchPattern(
          pattern      = target_dna,
          subject      = gene_dna,
          max.mismatch = max_mismatches
        )
        count_hits <- length(matches)
        match_positions_seq <- Biostrings::start(matches)
        
        if (count_hits > 0) {
          snippets <- sapply(match_positions_seq, function(pos) {
            # Ensure position does not exceed sequence boundaries
            end_pos <- pos + target_length - 1
            if (end_pos > nchar(seq_clean)) {
              return(NA_character_)
            } else {
              return(as.character(Biostrings::subseq(gene_dna, start = pos, end = end_pos)))
            }
          })
          # Remove any NAs
          snippets <- snippets[!is.na(snippets)]
          #Seperate 
          snippets <- paste(snippets, collapse = "; ")
        } else {
          snippets <- NA_character_
        }
      }
      
      homology_results <- rbind(
        homology_results,
        data.frame(
          `Gene Name` = gene_name,
          Header = header,
          Occurrences = count_hits,
          `Matched Snippets` = snippets,
          stringsAsFactors = FALSE
        )
      )
      
      # Store match positions
      matches_list[[header]] <- match_positions_seq
      
      progress$inc(1, detail = paste("Processed", i, "of", nrow(fasta_df)))
    }
    
    # Filter to only include sequences with at least one occurrence
    homology_results_filtered <- homology_results %>%
      dplyr::filter(Occurrences > 0)
    
    # Store the homology results
    homology_res(homology_results_filtered)
    
    # Store match positions
    match_positions(matches_list)
    
    showNotification("Homology search completed!", type = "message")
  })
  
  # Render the homology results
  output$homologyTable <- DT::renderDataTable({
    req(homology_res())
    res <- homology_res() %>% dplyr::arrange(desc(Occurrences))
    
    DT::datatable(
      res,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      colnames = c("Gene Name", "Header", "Occurrences", "Matched Snippets")
    )
  })
  
  # Download handler for homology search
  output$downloadHomology <- downloadHandler(
    filename = function() {
      paste("homology_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(homology_res())
      write.csv(homology_res(), file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
