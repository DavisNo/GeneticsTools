# Load libraries
library(shiny)
library(DT)
library(dplyr)

# Set the maximum upload file size to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)

# Function to Label Up/Down Degulation
direction_label <- function(lfc) {
  ifelse(lfc > 0, "Up",
         ifelse(lfc <= 0, "Down", NA))
}

ui <- fluidPage(
  titlePanel("DESeq2 Results Filter and Downloader"),
  
  sidebarLayout(
    sidebarPanel(
      # Dataset 1 Input
      h4("Dataset 1"),
      fileInput("deseq1", "Upload DESeq2 Results Set 1", 
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      textInput("title1", "Title for Dataset 1", value = "Dataset_1"),
      hr(),
      
      # Dataset 2 Input
      h4("Dataset 2"),
      fileInput("deseq2", "Upload DESeq2 Results Set 2", 
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      textInput("title2", "Title for Dataset 2", value = "Dataset_2"),
      hr(),
      
      # Filtering Criteria Inputs
      h4("Filtering Criteria"),
      
      # Wald statistic
      checkboxInput("use_stat", "Use Wald Statistic Filtering", value = FALSE),
      conditionalPanel(
        condition = "input.use_stat == true",
        checkboxInput("twoSided_stat", "Two-Sided Range", value = FALSE),
        textInput("stat_min", "Minimum Wald Statistic", placeholder = "Leave blank for no min"),
        textInput("stat_max", "Maximum Wald Statistic", placeholder = "Leave blank for no max")
      ),
      
      # baseMean
      checkboxInput("use_baseMean", "Use baseMean Filtering", value = FALSE),
      conditionalPanel(
        condition = "input.use_baseMean == true",
        checkboxInput("twoSided_baseMean", "Two-Sided Range", value = FALSE),
        textInput("baseMean_min", "Minimum baseMean", placeholder = "Leave blank for no min"),
        textInput("baseMean_max", "Maximum baseMean", placeholder = "Leave blank for no max")
      ),
      
      # log2FoldChange
      checkboxInput("use_log2FC", "Use log2FoldChange Filtering", value = FALSE),
      conditionalPanel(
        condition = "input.use_log2FC == true",
        checkboxInput("twoSided_log2FC", "Two-Sided Range", value = FALSE),
        textInput("log2FC_min", "Minimum log2FoldChange", placeholder = "Leave blank for no min"),
        textInput("log2FC_max", "Maximum log2FoldChange", placeholder = "Leave blank for no max")
      ),
      
      # p-value
      checkboxInput("use_pvalue", "Use p-value Filtering", value = FALSE),
      conditionalPanel(
        condition = "input.use_pvalue == true",
        textInput("pvalue_min", "Minimum p-value", placeholder = "Leave blank for no min"),
        textInput("pvalue_max", "Maximum p-value", placeholder = "Leave blank for no max")
      ),
      
      hr(),
      
      # Button to run filtering
      actionButton("run_filter", "Run Filtering"),
      
      # Panel for download buttons (visible after filtering)
      conditionalPanel(
        condition = "output.filtering_done == true",
        hr(),
        h4("Download Filtered Results"),
        uiOutput("download_buttons")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Filtered Results",
                 h3("Filtered Data"),
                 uiOutput("results_ui")
        ),
        tabPanel("DE Gene List",
                 h3("Intersection of DE Genes (Both Datasets)"),
                 DT::dataTableOutput("de_gene_table"),
                 br(),
                 h4("Contradictory Regulation: Up in One, Down in the Other"),
                 DT::dataTableOutput("contradictory_table"),
                 br(),
                 
                 # Download Options in DE List Tab
                 checkboxInput("exclude_contradictory", 
                               "Exclude contradictory genes from the Full DE List download",
                               value = FALSE),
                 br(),
                 downloadButton("download_upregulated", "Upregulated Genes"),
                 downloadButton("download_downregulated", "Downregulated Genes"),
                 downloadButton("download_full_de_list", "Full DE List")
        ),
        tabPanel("Gene Search",
                 h3("Gene Search and Download"),
                 selectInput("selected_dataset", "Select Dataset:", choices = NULL),
                 br(),
                 downloadButton("download_gene_search", "Download Selected Genes with Data")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive values to store data
  data <- reactiveValues(
    deseq1 = NULL,
    deseq2 = NULL,
    filtered1 = NULL,
    filtered2 = NULL,
    filtering_done = FALSE,
    de_genes = NULL,       
    de_table = NULL,
    diff_table = NULL 
  )
  
  # Function to load and validate DESeq2 data
  load_deseq <- function(input_file, dataset_num) {
    req(input_file)
    tryCatch({
      # Handle scientific notation (e.g., 1.92E-225)
      deseq_data <- read.csv(input_file$datapath, row.names = 1, check.names = FALSE)
      required_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
      if (!all(required_cols %in% colnames(deseq_data))) {
        showModal(modalDialog(
          title = "Error",
          paste0("Dataset ", dataset_num, " is missing required DESeq2 columns."),
          easyClose = TRUE,
          footer = NULL
        ))
        return(NULL)
      } else {
        return(deseq_data)
      }
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste0("Error loading Dataset ", dataset_num, ": ", e$message),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    })
  }
  
  # Observe file uploads
  observeEvent(input$deseq1, {
    data$deseq1 <- load_deseq(input$deseq1, 1)
  })
  observeEvent(input$deseq2, {
    data$deseq2 <- load_deseq(input$deseq2, 2)
  })
  
  # Update dataset choices for gene search
  observe({
    titles <- c()
    if (!is.null(data$deseq1)) titles <- c(titles, input$title1)
    if (!is.null(data$deseq2)) titles <- c(titles, input$title2)
    updateSelectInput(session, "selected_dataset", choices = titles)
  })
  
  # Function to parse numeric text inputs 
  parse_num <- function(txt) {
    val <- suppressWarnings(as.numeric(txt))
    if (is.na(val)) return(NULL) else return(val)
  }
  
  # Function to apply filtering to one dataset
  filter_dataset <- function(deseq_data,
                             use_stat, twoSided_stat, stat_min, stat_max,
                             use_baseMean, twoSided_baseMean, baseMean_min, baseMean_max,
                             use_log2FC, twoSided_log2FC, log2FC_min, log2FC_max,
                             use_pvalue, pvalue_min, pvalue_max) {
    
    if (is.null(deseq_data)) return(NULL)
    df <- deseq_data
    
    # Convert text inputs to numeric
    min_stat <- parse_num(stat_min)
    max_stat <- parse_num(stat_max)
    min_baseMean <- parse_num(baseMean_min)
    max_baseMean <- parse_num(baseMean_max)
    min_log2FC <- parse_num(log2FC_min)
    max_log2FC <- parse_num(log2FC_max)
    min_pval <- parse_num(pvalue_min)
    max_pval <- parse_num(pvalue_max)
    
    # Wald statistic filtering
    if (use_stat) {
      if (!is.null(min_stat)) {
        if (twoSided_stat) {
          df <- df %>% filter(abs(stat) >= min_stat)
        } else {
          df <- df %>% filter(stat >= min_stat)
        }
      }
      if (!is.null(max_stat)) {
        if (twoSided_stat) {
          df <- df %>% filter(abs(stat) <= max_stat)
        } else {
          df <- df %>% filter(stat <= max_stat)
        }
      }
    }
    
    # baseMean filtering
    if (use_baseMean) {
      if (!is.null(min_baseMean)) {
        if (twoSided_baseMean) {
          df <- df %>% filter(abs(baseMean) >= min_baseMean)
        } else {
          df <- df %>% filter(baseMean >= min_baseMean)
        }
      }
      if (!is.null(max_baseMean)) {
        if (twoSided_baseMean) {
          df <- df %>% filter(abs(baseMean) <= max_baseMean)
        } else {
          df <- df %>% filter(baseMean <= max_baseMean)
        }
      }
    }
    
    # log2FoldChange filtering
    if (use_log2FC) {
      if (!is.null(min_log2FC)) {
        if (twoSided_log2FC) {
          df <- df %>% filter(abs(log2FoldChange) >= min_log2FC)
        } else {
          df <- df %>% filter(log2FoldChange >= min_log2FC)
        }
      }
      if (!is.null(max_log2FC)) {
        if (twoSided_log2FC) {
          df <- df %>% filter(abs(log2FoldChange) <= max_log2FC)
        } else {
          df <- df %>% filter(log2FoldChange <= max_log2FC)
        }
      }
    }
    
    # p-value filtering
    if (use_pvalue) {
      if (!is.null(min_pval)) {
        df <- df %>% filter(pvalue >= min_pval)
      }
      if (!is.null(max_pval)) {
        df <- df %>% filter(pvalue <= max_pval)
      }
    }
    
    return(df)
  }
  
  # When user clicks "Run Filtering"
  observeEvent(input$run_filter, {
    # Ensure at least one dataset is uploaded
    if (is.null(data$deseq1) && is.null(data$deseq2)) {
      showModal(modalDialog(
        title = "Error",
        "Please upload at least one DESeq2 results file before running the filter.",
        easyClose = TRUE,
        footer = NULL
      ))
      data$filtering_done <- FALSE
      return()
    }
    
    # Filter each dataset separately
    data$filtered1 <- filter_dataset(
      deseq_data = data$deseq1,
      use_stat = input$use_stat, twoSided_stat = input$twoSided_stat,
      stat_min = input$stat_min, stat_max = input$stat_max,
      use_baseMean = input$use_baseMean, twoSided_baseMean = input$twoSided_baseMean,
      baseMean_min = input$baseMean_min, baseMean_max = input$baseMean_max,
      use_log2FC = input$use_log2FC, twoSided_log2FC = input$twoSided_log2FC,
      log2FC_min = input$log2FC_min, log2FC_max = input$log2FC_max,
      use_pvalue = input$use_pvalue, pvalue_min = input$pvalue_min, pvalue_max = input$pvalue_max
    )
    
    data$filtered2 <- filter_dataset(
      deseq_data = data$deseq2,
      use_stat = input$use_stat, twoSided_stat = input$twoSided_stat,
      stat_min = input$stat_min, stat_max = input$stat_max,
      use_baseMean = input$use_baseMean, twoSided_baseMean = input$twoSided_baseMean,
      baseMean_min = input$baseMean_min, baseMean_max = input$baseMean_max,
      use_log2FC = input$use_log2FC, twoSided_log2FC = input$twoSided_log2FC,
      log2FC_min = input$log2FC_min, log2FC_max = input$log2FC_max,
      use_pvalue = input$use_pvalue, pvalue_min = input$pvalue_min, pvalue_max = input$pvalue_max
    )
    
    if (is.null(data$filtered1) || nrow(data$filtered1) == 0 ||
        is.null(data$filtered2) || nrow(data$filtered2) == 0) {
      showModal(modalDialog(
        title = "No intersection",
        "No genes were retained in both datasets",
        easyClose = TRUE,
        footer = NULL
      ))
      data$filtering_done <- FALSE
      return()
    }
    
    # Identify intersection of rownames
    common_genes <- intersect(rownames(data$filtered1), rownames(data$filtered2))
    if (length(common_genes) == 0) {
      showModal(modalDialog(
        title = "No intersection",
        "No genes passed the filtering criteria in both datasets",
        easyClose = TRUE,
        footer = NULL
      ))
      data$filtering_done <- FALSE
      return()
    }
    
    # Build DE gene table only intersection
    data$de_genes <- common_genes
    df_de <- data.frame(Geneid = common_genes, stringsAsFactors = FALSE)
    
    # Fill regulation columns
    df_de$Direction1 <- direction_label(data$filtered1[common_genes, "log2FoldChange"])
    df_de$Direction2 <- direction_label(data$filtered2[common_genes, "log2FoldChange"])
    
    data$de_table <- df_de
    
    # Build contradictory table
    df_diff <- df_de %>%
      filter(!is.na(Direction1), !is.na(Direction2)) %>%
      filter(
        (Direction1 == "Up" & Direction2 == "Down") |
          (Direction1 == "Down" & Direction2 == "Up")
      )
    data$diff_table <- df_diff
    
    data$filtering_done <- TRUE
  })
  
  output$filtering_done <- reactive({
    data$filtering_done
  })
  outputOptions(output, "filtering_done", suspendWhenHidden = FALSE)
  
  # Render download buttons for each dataset
  output$download_buttons <- renderUI({
    buttons <- list()
    if (!is.null(data$filtered1) && nrow(data$filtered1) > 0) {
      buttons[[length(buttons) + 1]] <- downloadButton("download1", input$title1)
    }
    if (!is.null(data$filtered2) && nrow(data$filtered2) > 0) {
      buttons[[length(buttons) + 1]] <- downloadButton("download2", input$title2)
    }
    do.call(tagList, buttons)
  })
  
  # Download Handlers for Filtered Results (Datasets 1 & 2)
  output$download1 <- downloadHandler(
    filename = function() {
      paste0(gsub(" ", "_", input$title1), "_",
             format(Sys.Date(), "%d%b%Y"), "_filtered_genes.csv")
    },
    content = function(file) {
      req(data$filtered1)
      write.csv(data$filtered1, file, row.names = TRUE)
    }
  )
  
  output$download2 <- downloadHandler(
    filename = function() {
      paste0(gsub(" ", "_", input$title2), "_",
             format(Sys.Date(), "%d%b%Y"), "_filtered_genes.csv")
    },
    content = function(file) {
      req(data$filtered2)
      write.csv(data$filtered2, file, row.names = TRUE)
    }
  )
  
  # Show filtered results in the main panel
  output$results_ui <- renderUI({
    if (!data$filtering_done) {
      h5("Run the filtering to see the results here.")
    } else {
      tabs <- list()
      if (!is.null(data$filtered1) && nrow(data$filtered1) > 0) {
        tabs[[length(tabs) + 1]] <- tabPanel(
          input$title1, DT::dataTableOutput("table1")
        )
      }
      if (!is.null(data$filtered2) && nrow(data$filtered2) > 0) {
        tabs[[length(tabs) + 1]] <- tabPanel(
          input$title2, DT::dataTableOutput("table2")
        )
      }
      do.call(tabsetPanel, tabs)
    }
  })
  
  output$table1 <- DT::renderDataTable({
    req(data$filtered1)
    datatable(data$filtered1, options = list(pageLength = 10))
  })
  
  output$table2 <- DT::renderDataTable({
    req(data$filtered2)
    datatable(data$filtered2, options = list(pageLength = 10))
  })
  
  # DE Gene List Table (intersection), color-coded
  output$de_gene_table <- DT::renderDataTable({
    req(data$de_table)
    
    df_display <- data$de_table
    
    # Rename columns to specify dataset
    colnames(df_display)[2] <- paste("Direction in", input$title1)
    colnames(df_display)[3] <- paste("Direction in", input$title2)
    
    # Apply color coding
    df_display[[2]] <- ifelse(
      df_display[[2]] == "Up", "<span style='color:red'>Up</span>",
      ifelse(df_display[[2]] == "Down", "<span style='color:blue'>Down</span>", "")
    )
    df_display[[3]] <- ifelse(
      df_display[[3]] == "Up", "<span style='color:red'>Up</span>",
      ifelse(df_display[[3]] == "Down", "<span style='color:blue'>Down</span>", "")
    )
    
    datatable(df_display,
              escape = FALSE,
              options = list(pageLength = 10))
  })
  
  # Contradictory table
  output$contradictory_table <- DT::renderDataTable({
    req(data$diff_table)
    
    df_display <- data$diff_table
    
    # Rename columns for clarity
    colnames(df_display)[2] <- paste("Direction in", input$title1)
    colnames(df_display)[3] <- paste("Direction in", input$title2)
    
    # Color them
    df_display[[2]] <- ifelse(
      df_display[[2]] == "Up", "<span style='color:red'>Up</span>",
      ifelse(df_display[[2]] == "Down", "<span style='color:blue'>Down</span>", "")
    )
    df_display[[3]] <- ifelse(
      df_display[[3]] == "Up", "<span style='color:red'>Up</span>",
      ifelse(df_display[[3]] == "Down", "<span style='color:blue'>Down</span>", "")
    )
    
    datatable(df_display,
              escape = FALSE,
              options = list(pageLength = 10))
  })
  
  # (A) Download Upregulated Genes (Up in BOTH datasets)
  output$download_upregulated <- downloadHandler(
    filename = function() {
      "UpregulatedGenes.csv"
    },
    content = function(file) {
      req(data$de_table)
      # Subset to Up in both
      df_up <- data$de_table %>%
        filter(Direction1 == "Up", Direction2 == "Up")
      write.csv(df_up, file, row.names = FALSE)
    }
  )
  
  # (B) Download Downregulated Genes (Down in BOTH datasets)
  output$download_downregulated <- downloadHandler(
    filename = function() {
      "DownregulatedGenes.csv"
    },
    content = function(file) {
      req(data$de_table)
      # Subset to Down in both
      df_down <- data$de_table %>%
        filter(Direction1 == "Down", Direction2 == "Down")
      write.csv(df_down, file, row.names = FALSE)
    }
  )
  
  # (C) Download Full DE List (optionally exclude contradictory genes)
  output$download_full_de_list <- downloadHandler(
    filename = function() {
      "FullDEList.csv"
    },
    content = function(file) {
      req(data$de_table)
      df_full <- data$de_table
      
      if (input$exclude_contradictory) {
        # Exclude genes that appear in diff_table
        if (!is.null(data$diff_table) && nrow(data$diff_table) > 0) {
          contras <- data$diff_table$Geneid
          df_full <- df_full[!(df_full$Geneid %in% contras), ]
        }
      }
      
      write.csv(df_full, file, row.names = FALSE)
    }
  )
  
  # Gene Search Download Handler
  output$download_gene_search <- downloadHandler(
    filename = function() {
      paste0("Gene_search_",
             gsub(" ", "_", input$selected_dataset), "_",
             format(Sys.Date(), "%d%b%Y"), ".csv")
    },
    content = function(file) {
      req(data$de_genes)
      selected_title <- input$selected_dataset
      if (!is.null(data$deseq1) && selected_title == input$title1) {
        selected_data <- data$filtered1
      } else if (!is.null(data$deseq2) && selected_title == input$title2) {
        selected_data <- data$filtered2
      } else {
        selected_data <- NULL
      }
      req(selected_data)
      gene_subset <- selected_data[rownames(selected_data) %in% data$de_genes, ]
      write.csv(gene_subset, file, row.names = TRUE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)

