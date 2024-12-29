# Load libraries
library(shiny)
library(readr)
library(dplyr)
library(biomaRt)
library(DT)
library(Biostrings)
library(httr)
library(jsonlite)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body, h1, h2, h3, h4, h5, h6, p, label, input, textarea, select, button, table {
        color: black !important;
      }
    "))
  ),
  
  titlePanel("FASTA & BED Fetcher w/ Promoter Using Ensembl REST"),
  
  sidebarLayout(
    sidebarPanel(
     
      fileInput("file1", "Upload Gene List (CSV/TSV):",
                accept = c(".csv", ".tsv", ".txt")),
      
      numericInput("extension_length", "Extension (bases) on each side:", 
                   value = 250, min = 0, step = 50),
      
      # Reference Genome 
      tags$h5("Reference Genome: GRCh38"),
      
      # FASTA Section
      tags$h4("FASTA"),
      
      # Exon notation FASTA checkbox
      checkboxInput("exon_notation", "Include Exons in FASTA Headers?", TRUE),
      
      br(),
      
      # BED Section
      tags$h4("BED"),
      # Exon checkbox
      checkboxInput("include_exon", "Include Exon as BED sub-blocks", FALSE),
      
      # Promoter checkbox
      checkboxInput("include_promoter", "Include Promoter as BED sub-blocks", TRUE),
      
      # Checkbox to separate promoter and exon entries
      conditionalPanel(
        condition = "input.include_promoter == true",
        checkboxInput("separate_promoter_exon", "Separate Promoter Entries", FALSE)
      ),
      
      actionButton("process", "Fetch FASTA & BED")

    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "FASTA Results",
          h4("Genes and Extended FASTA Sequences"),
          DT::dataTableOutput("displayResults"),
          br(),
          h4("FASTA Retrieval Summary"),
          verbatimTextOutput("retrievalSummary"),
          br(),
          h4("Download FASTA Files"),
          downloadButton("downloadDataFASTA", "Download All FASTAs as TXT")
        ),
        tabPanel(
          "BED Results",
          h4("Extended BED Data (Exons, Promoter)"),
          DT::dataTableOutput("bedTable"),
          br(),
          h4("Download BED File"),
          downloadButton("downloadBED", "Download BED File")
        ),
        tabPanel(
          "Specific Genes",
          h4("Selected Genes' FASTA and BED Data"),
          textInput("specific_genes", 
                    "Download Specific Genes' FASTA & BED (space-separated):",
                    placeholder = "e.g., GENE1 GENE2 GENE3"),
          actionButton("downloadSpecificBtn", "Prepare FASTA & BED Download for Specific Genes"),
          br(), br(),
          h4("Selected Genes' FASTA Sequences"),
          DT::dataTableOutput("specificFastaTable"),
          br(),
          h4("Selected Genes' BED Data"),
          DT::dataTableOutput("specificBedTable"),
          br(),
          h4("Download Selected FASTA and BED Files"),
          downloadButton("downloadSpecificFASTA", "Download Selected FASTAs as TXT"),
          downloadButton("downloadSpecificBED", "Download Selected BED as BED6 or BED12")
        ),
        tabPanel(
          "Visualization",
          h4("Basic Exon/Promoter Visualization"),
          textInput("visual_genes", "Enter Gene(s) to Visualize (space-separated)", ""),
          actionButton("visualizeBtn", "Visualize"),
          br(),
          plotOutput("exonPlot", height = "600px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(

    gene_col       = "Geneid",
    attribute_gene = "external_gene_name" 
  )
  
  # Read uploaded data
  uploaded_data <- reactive({
    req(input$file1)
    ext <- tools::file_ext(input$file1$name)
    
    # Read the file and column named 'Geneid'
    if (ext == "csv") {
      df <- tryCatch({
        read_csv(input$file1$datapath, col_names = TRUE, show_col_types = FALSE)
      }, error = function(e) {
        showNotification("Error reading CSV file.", type = "error")
        return(NULL)
      })
    } else if (ext %in% c("tsv", "txt")) {
      df <- tryCatch({
        read_delim(input$file1$datapath, delim = "\t", col_names = TRUE, show_col_types = FALSE)
      }, error = function(e) {
        showNotification("Error reading TSV/TXT file.", type = "error")
        return(NULL)
      })
    } else {
      showNotification("Unsupported format. Use CSV.", type = "error")
      return(NULL)
    }
    
    # Check for 'Geneid' column
    if (!"Geneid" %in% colnames(df)) {
      showNotification("Input file must contain a 'Geneid' column.", type = "error")
      return(NULL)
    }
    
    genes <- unique(df[["Geneid"]])
    genes <- genes[!is.na(genes) & genes != ""]
    if (length(genes) == 0) {
      showNotification("No valid gene IDs found in 'Geneid' column.", type = "error")
      return(NULL)
    }
    return(genes)
  })
  
  results <- reactiveVal(NULL)
  missing_genes <- reactiveVal(NULL)
  retrieval_summary_val <- reactiveVal(NULL)
  
  bed_data <- reactiveVal(NULL)
  
  specific_fasta_data <- reactiveVal(NULL)
  specific_bed_data <- reactiveVal(NULL)
  
  get_promoter_positions <- function(start_pos, end_pos, strand, promoter_size = 2000) {
    # Define TSS based on strand
    if (strand == 1) {
      promoter_start <- max(start_pos - promoter_size, 1)
      promoter_end   <- max(start_pos - 1, 1)
    } else if (strand == -1) {
      promoter_start <- end_pos + 1
      promoter_end   <- end_pos + promoter_size
    } else {
      # If no valid strand, return NULL
      return(NULL)
    }
    if (promoter_start > promoter_end) {
      return(NULL)
    }
    return(c(promoter_start, promoter_end))
  }
  
  # Process
  observeEvent(input$process, {
    genes <- uploaded_data()
    req(genes)
    
    ensembl <- tryCatch({
      useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        host = "https://www.ensembl.org"
      )
    }, error = function(e) {
      showNotification("Failed to connect to Ensembl BioMart. Please try again later.", type = "error")
      return(NULL)
    })
    
    if (is.null(ensembl)) {
      return(NULL)
    }
    
    withProgress(message = "Mapping genes to Ensembl", value = 0, {
      
      attributes <- c(
        rv$attribute_gene, 
        "ensembl_gene_id", 
        "ensembl_transcript_id", 
        "transcript_is_canonical",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand"
      )
      
      # Fetch gene mappings
      mapped <- tryCatch({
        biomaRt::getBM(
          attributes = attributes,
          filters    = "external_gene_name",  # Assuming 'Geneid' corresponds to 'external_gene_name'
          values     = genes,
          mart       = ensembl
        )
      }, error = function(e) {
        showNotification("Error fetching data from BioMart.", type = "error")
        return(NULL)
      })
      
      if (is.null(mapped) || nrow(mapped) == 0) {
        showNotification("No data retrieved from BioMart.", type = "error")
        return(NULL)
      }
      
      required_cols <- c(rv$attribute_gene, "ensembl_gene_id", "ensembl_transcript_id",
                         "transcript_is_canonical", "chromosome_name", 
                         "start_position", "end_position", "strand")
      missing_cols <- setdiff(required_cols, colnames(mapped))
      if (length(missing_cols) > 0) {
        showNotification(
          paste("Missing columns in BioMart data:", paste(missing_cols, collapse = ", ")),
          type = "error"
        )
        return(NULL)
      }
      
      incProgress(0.3)
      
      mapped_genes <- unique(mapped[[rv$attribute_gene]])
      unmapped     <- setdiff(genes, mapped_genes)
      if (length(unmapped) > 0) {
        showNotification(
          paste(length(unmapped), "genes failed to map."),
          type = "warning"
        )
      }
      
      mapped <- dplyr::distinct(mapped)
      
      # Filter for canonical transcripts
      mapped_canonical <- tryCatch({
        dplyr::filter(mapped, transcript_is_canonical == TRUE)
      }, error = function(e) {
        showNotification("Error filtering canonical transcripts.", type = "error")
        return(NULL)
      })
      
      if (is.null(mapped_canonical) || nrow(mapped_canonical) == 0) {
        showNotification("No canonical transcripts found.", type = "error")
        return(NULL)
      }
      
      mapped_canonical_unique <- tryCatch({
        dplyr::group_by(mapped_canonical, !!sym(rv$attribute_gene)) %>%
          dplyr::slice_head(n = 1) %>%
          dplyr::ungroup()
      }, error = function(e) {
        showNotification("Error ensuring unique canonical transcripts.", type = "error")
        return(NULL)
      })
      
      if (is.null(mapped_canonical_unique)) {
        return(NULL)
      }
      
      incProgress(0.1, detail = "Fetching sequences from REST")
      
      # Initialize the main data frame
      extended_sequences <- data.frame(
        Gene                  = character(),
        ensembl_gene_id       = character(),
        ensembl_transcript_id = character(),
        Chromosome            = character(),
        Start                 = numeric(),
        End                   = numeric(),
        Strand                = numeric(),
        FASTA                 = character(),
        block_positions_list  = I(list()),
        stringsAsFactors      = FALSE
      )
      
      total_genes <- nrow(mapped_canonical_unique)
      extension   <- input$extension_length
      
      # Loop over each gene
      for (i in seq_len(total_genes)) {
        gene_info <- mapped_canonical_unique[i, , drop = FALSE]
        
        gene_symbol <- if (!is.na(gene_info[[rv$attribute_gene]]) &&
                           nzchar(gene_info[[rv$attribute_gene]])) {
          gene_info[[rv$attribute_gene]]
        } else {
          gene_info[["ensembl_gene_id"]]
        }
        
        chrom    <- as.character(gene_info[["chromosome_name"]])
        startpos <- gene_info[["start_position"]]
        endpos   <- gene_info[["end_position"]]
        strand   <- gene_info[["strand"]]
        
        # Extend coordinates (ensure ext_start >= 1)
        ext_start <- max(startpos - extension, 1)
        ext_end   <- endpos + extension
        
        # 1) Main sequence
        region_str <- paste0(chrom, ":", ext_start, "..", ext_end)
        seq_url    <- paste0(
          "https://rest.ensembl.org/sequence/region/human/",
          region_str,
          "?assembly=GRCh38"
        )
        
        seq_text <- NA_character_
        seq_resp <- tryCatch({
          GET(seq_url, add_headers("Accept" = "text/plain"))
        }, error = function(e) {
          showNotification(paste("Failed to fetch sequence for", gene_symbol), type = "warning")
          return(NULL)
        })
        
        if (!is.null(seq_resp) && status_code(seq_resp) == 200) {
          seq_tmp <- content(seq_resp, as = "text", encoding = "UTF-8")
          seq_tmp <- gsub("\n", "", seq_tmp)
          if (strand == -1) {
            seq_rc <- reverseComplement(DNAString(seq_tmp))
            seq_tmp <- as.character(seq_rc)
          }
          seq_text <- gsub("[^ATCG]", "", toupper(seq_tmp))
        } else {
          showNotification(paste("Sequence not found for", gene_symbol), type = "warning")
        }
        
        # Initialize list to hold all sub-blocks (promoter, exons)
        block_positions <- list()
        
        # 2) Promoter region: using extension_length
        if (input$include_promoter) {
          promoter_c <- get_promoter_positions(startpos, endpos, strand, extension)
          if (!is.null(promoter_c)) {
            # Ensure it’s within [ext_start..ext_end]
            pr_start <- max(promoter_c[1], ext_start)
            pr_end   <- min(promoter_c[2], ext_end)
            if (pr_start <= pr_end) {
              block_positions[[length(block_positions) + 1]] <- list(
                type = "promoter",
                coords = c(pr_start, pr_end)
              )
            }
          }
        }
        
        # 3) Exons
        exon_info <- ""
        exon_blocks <- list()
        if (input$include_exon) {
          exon_url <- paste0(
            "https://rest.ensembl.org/overlap/region/human/",
            region_str,
            "?feature=exon;assembly=GRCh38"
          )
          
          exon_resp <- tryCatch({
            GET(exon_url, add_headers("Accept" = "application/json"))
          }, error = function(e) {
            showNotification(paste("Failed to fetch exon data for", gene_symbol), type = "warning")
            return(NULL)
          })
          
          if (!is.null(exon_resp) && status_code(exon_resp) == 200) {
            exon_data <- tryCatch({
              fromJSON(content(exon_resp, as = "text", encoding = "UTF-8"))
            }, error = function(e) {
              NULL
            })
            
            # Filter exons for canonical transcript
            if (!is.null(exon_data) && length(exon_data) > 0) {
              exon_df <- as.data.frame(exon_data, stringsAsFactors = FALSE)
              if ("Parent" %in% colnames(exon_df)) {
                exon_data_filtered <- dplyr::filter(exon_df, Parent == gene_info[["ensembl_transcript_id"]])
                
                if (nrow(exon_data_filtered) > 0) {
                  exon_coords <- sapply(seq_len(nrow(exon_data_filtered)), function(j) {
                    row_j <- exon_data_filtered[j, ]
                    s_val <- row_j[["start"]]
                    e_val <- row_j[["end"]]
                    if (!is.na(s_val) && !is.na(e_val)) {
                      paste0("exon(", s_val, "-", e_val, ")")
                    } else {
                      "exon(NA-NA)"
                    }
                  })
                  exon_info <- paste(exon_coords, collapse = ";")
                  
                  # Add exons to block_positions
                  for (j in seq_len(nrow(exon_data_filtered))) {
                    row_j <- exon_data_filtered[j, ]
                    s_val <- max(row_j[["start"]], ext_start)
                    e_val <- min(row_j[["end"]], ext_end)
                    if (s_val <= e_val) {
                      block_positions[[length(block_positions) + 1]] <- list(
                        type = "exon",
                        coords = c(s_val, e_val)
                      )
                      exon_blocks <- c(exon_blocks, list(c(s_val, e_val)))
                    }
                  }
                }
              }
            }
          }
        }
        
        # Build FASTA header
        if (!is.na(seq_text)) {
          # Possibly annotate exons
          exon_part <- ""
          if (input$exon_notation && nchar(exon_info) > 0) {
            exon_part <- paste0(" |exons=", exon_info)
          }
          header_info <- paste0(
            ">", gene_symbol, 
            " |chr=", chrom, 
            ":", ext_start, "-", ext_end, 
            " |strand=", strand,
            exon_part
          )
          seq_text <- paste0(header_info, "\n", seq_text)
        }
        
        # Build a single-row data frame
        temp_row <- data.frame(
          Gene                  = gene_symbol,
          ensembl_gene_id       = gene_info[["ensembl_gene_id"]],
          ensembl_transcript_id = gene_info[["ensembl_transcript_id"]],
          Chromosome            = chrom,
          Start                 = ext_start,
          End                   = ext_end,
          Strand                = strand,
          FASTA                 = if (!is.na(seq_text)) seq_text else NA_character_,
          stringsAsFactors      = FALSE
        )
        
        # Attach sub-blocks
        temp_row$block_positions_list <- list(block_positions)
        
        # Now rbind safely using dplyr::bind_rows to ensure consistency
        extended_sequences <- tryCatch({
          dplyr::bind_rows(extended_sequences, temp_row)
        }, error = function(e) {
          showNotification(paste("Failed to append data for", gene_symbol), type = "warning")
          return(extended_sequences)
        })
        
        incProgress(1 / total_genes, detail = paste("Processed", i, "of", total_genes))
        
        # Short pause 
        Sys.sleep(0.1)
      }
      
      # Handle unmapped genes
      if (length(unmapped) > 0) {
        for (um_gene in unmapped) {
          temp_unmapped <- data.frame(
            Gene                  = um_gene,
            ensembl_gene_id       = NA_character_,
            ensembl_transcript_id = NA_character_,
            Chromosome            = NA_character_,
            Start                 = NA_real_,
            End                   = NA_real_,
            Strand                = NA_real_,
            FASTA                 = NA_character_,
            stringsAsFactors      = FALSE
          )
          temp_unmapped$block_positions_list <- list(list())
          
          extended_sequences <- dplyr::bind_rows(extended_sequences, temp_unmapped)
        }
      }
      
      # Mark status
      final_df <- extended_sequences %>%
        mutate(
          FASTA_display = dplyr::case_when(
            is.na(FASTA) & !is.na(ensembl_gene_id) ~ "Sequence Unavailable",
            is.na(FASTA) & is.na(ensembl_gene_id)  ~ "Failed to Map",
            TRUE                                   ~ FASTA
          ),
          Status = dplyr::case_when(
            FASTA_display == "Sequence Unavailable" ~ "Sequence Unavailable",
            FASTA_display == "Failed to Map"        ~ "Failed to Map",
            TRUE                                    ~ "Sequence Retrieved"
          )
        )
      
      missing_seq <- final_df %>% dplyr::filter(Status != "Sequence Retrieved")
      if (nrow(missing_seq) > 0) {
        showNotification(
          paste(nrow(missing_seq), "genes had issues with FASTA retrieval or mapping."),
          type = "warning"
        )
      }
      missing_genes(missing_seq %>% dplyr::select(Gene, Status))
      
      retrieved_genes <- final_df %>% dplyr::filter(Status == "Sequence Retrieved")
      
      results(final_df)
      
      retrieval_summary_val(
        list(
          total     = nrow(final_df),
          retrieved = nrow(retrieved_genes),
          missing   = nrow(missing_seq)
        )
      )
      
      # Prepare the BED data
      df_for_bed <- final_df %>%
        dplyr::filter(Status == "Sequence Retrieved") %>%
        dplyr::distinct(Gene, Chromosome, Start, End, Strand, block_positions_list)
      
      # Initialize BED dataframe
      bed_entries <- data.frame(
        chrom       = character(),
        chromStart  = numeric(),
        chromEnd    = numeric(),
        name        = character(),
        score       = numeric(),
        strand      = character(),
        thickStart  = numeric(),
        thickEnd    = numeric(),
        itemRgb     = character(),
        blockCount  = numeric(),
        blockSizes  = character(),
        blockStarts = character(),
        stringsAsFactors = FALSE
      )
      
      # Define colors
      promoter_color <- "255,0,0"      # Red for promoters
      exon_color     <- "0,0,255"      # Blue for exons
      default_color  <- "0,0,0"        # Black as default
      mixed_color    <- "128,128,128"  # Gray for mixed features
      
      # Define feature colors without UTR
      feature_colors <- c("promoter" = "red", "exon" = "blue")
      
      # Determine if any BED-related checkboxes are selected
      bed12_selected <- input$include_exon || input$include_promoter || input$separate_promoter_exon
      
      if (bed12_selected) {
        # Construct BED entries
        for (i in seq_len(nrow(df_for_bed))) {
          row_i <- df_for_bed[i, ]
          gene_name <- row_i$Gene
          chrom <- row_i$Chromosome
          bed_start <- row_i$Start - 1  # BED is 0-based
          bed_end <- row_i$End
          strand <- row_i$Strand
          
          # Initialize lists to hold promoters and exons
          promoter_blocks <- list()
          exon_blocks <- list()
          
          # Extract block positions
          blocks <- row_i$block_positions_list[[1]]
          
          # Categorize blocks
          for (block in blocks) {
            if (block$type == "promoter") {
              promoter_blocks <- c(promoter_blocks, list(block$coords))
            } else if (block$type == "exon") {
              exon_blocks <- c(exon_blocks, list(block$coords))
            }
          }
          
          # Handle promoter inclusion
          if (input$include_promoter && length(promoter_blocks) > 0) {
            if (input$separate_promoter_exon) {
              # Separate promoter entry
              promoter_coords <- promoter_blocks[[1]]
              promoter_start <- promoter_coords[1] - 1  # 0-based
              promoter_end <- promoter_coords[2]
              
              bed_entries <- dplyr::bind_rows(bed_entries, data.frame(
                chrom       = chrom,
                chromStart  = promoter_start,
                chromEnd    = promoter_end,
                name        = paste0(gene_name, "_promoter"),
                score       = 0,
                strand      = ifelse(strand == 1, "+", ifelse(strand == -1, "-", ".")),
                thickStart  = promoter_start,
                thickEnd    = promoter_end,
                itemRgb     = promoter_color,
                blockCount  = 1,
                blockSizes  = as.character(promoter_end - promoter_start),
                blockStarts = "0",
                stringsAsFactors = FALSE
              ))
            }
            # If not separating, promoter will be included in exons blocks
          }
          
          # Handle exons as sub-blocks
          if (input$include_exon) {
            if (length(exon_blocks) > 0) {
              # Combine promoter and exon blocks if not separating
              combined_blocks <- exon_blocks
              if (input$include_promoter && !input$separate_promoter_exon && length(promoter_blocks) > 0) {
                combined_blocks <- c(combined_blocks, promoter_blocks)
              }
              
              # Sort blocks by start position
              combined_blocks <- combined_blocks[order(sapply(combined_blocks, function(x) x[1]))]
              
              # Calculate block sizes and starts
              block_sizes <- paste(sapply(combined_blocks, function(x) x[2] - x[1] + 1), collapse = ",")
              block_starts <- paste(sapply(combined_blocks, function(x) x[1] - 1 - bed_start), collapse = ",")
              block_count <- length(combined_blocks)
              
              # Determine if multiple feature types are present
              feature_types <- sapply(blocks, function(b) b$type)
              unique_types <- unique(feature_types)
              
              if (input$include_promoter && !input$separate_promoter_exon && length(unique_types) > 1) {
                # Mixed feature types (promoter and exon)
                item_rgb <- mixed_color
              } else {
                # Single feature type
                if (input$include_promoter && length(promoter_blocks) > 0 && !input$separate_promoter_exon) {
                  item_rgb <- promoter_color
                } else {
                  item_rgb <- exon_color
                }
              }
              
              bed_entries <- dplyr::bind_rows(bed_entries, data.frame(
                chrom       = chrom,
                chromStart  = bed_start,
                chromEnd    = bed_end,
                name        = gene_name,
                score       = 0,
                strand      = ifelse(strand == 1, "+", ifelse(strand == -1, "-", ".")),
                thickStart  = bed_start,  # Can be adjusted if needed
                thickEnd    = bed_end,    # Can be adjusted if needed
                itemRgb     = item_rgb,
                blockCount  = block_count,
                blockSizes  = block_sizes,
                blockStarts = block_starts,
                stringsAsFactors = FALSE
              ))
            }
          }
          
          # If only promoter is included and no exons, ensure BED12 format
          if (input$include_promoter && !input$separate_promoter_exon && length(promoter_blocks) > 0 && !input$include_exon) {
            # Add promoter as the only block
            promoter_coords <- promoter_blocks[[1]]
            promoter_start <- promoter_coords[1] - 1  # 0-based
            promoter_end <- promoter_coords[2]
            
            bed_entries <- dplyr::bind_rows(bed_entries, data.frame(
              chrom       = chrom,
              chromStart  = promoter_start,
              chromEnd    = promoter_end,
              name        = gene_name,
              score       = 0,
              strand      = ifelse(strand == 1, "+", ifelse(strand == -1, "-", ".")),
              thickStart  = promoter_start,
              thickEnd    = promoter_end,
              itemRgb     = promoter_color,
              blockCount  = 1,
              blockSizes  = as.character(promoter_end - promoter_start),
              blockStarts = "0",
              stringsAsFactors = FALSE
            ))
          }
        }
      } else {
        # No BED-related checkboxes selected; generate BED6
        for (i in seq_len(nrow(df_for_bed))) {
          row_i <- df_for_bed[i, ]
          gene_name <- row_i$Gene
          chrom <- row_i$Chromosome
          bed_start <- row_i$Start - 1  # BED is 0-based
          bed_end <- row_i$End
          strand <- row_i$Strand
          
          bed_entries <- dplyr::bind_rows(bed_entries, data.frame(
            chrom       = chrom,
            chromStart  = bed_start,
            chromEnd    = bed_end,
            name        = gene_name,
            score       = 0,
            strand      = ifelse(strand == 1, "+", ifelse(strand == -1, "-", ".")),
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Final BED data frame
      final_bed_df <- bed_entries
      
      bed_data(final_bed_df)
    })  
  })
  
  # Download full FASTA
  output$downloadDataFASTA <- downloadHandler(
    filename = function() {
      paste("extended_fasta_sequences_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      req(results())
      fasta_sequences <- results() %>% dplyr::filter(Status == "Sequence Retrieved")
      
      complete_fasta <- c()
      for (i in seq_len(nrow(fasta_sequences))) {
        if (!is.na(fasta_sequences$FASTA[i])) {
          lines_i <- unlist(strsplit(fasta_sequences$FASTA[i], "\n"))
          complete_fasta <- c(complete_fasta, lines_i, "")
        }
      }
      writeLines(complete_fasta, con = file)
    }
  )
  
  # Display table of results
  output$displayResults <- DT::renderDataTable({
    req(results())
    df <- results()
    DT::datatable(
      df %>% dplyr::select(
        Gene, ensembl_gene_id, ensembl_transcript_id,
        Chromosome, Start, End, Strand, FASTA_display, Status
      ),
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      colnames = c(
        "Gene Symbol/ID", "Ensembl Gene ID", "Ensembl TX ID",
        "Chromosome", "Start", "End", "Strand",
        "FASTA Header + Sequence", "Status"
      )
    ) %>%
      DT::formatStyle(
        'Status',
        target = 'row',
        backgroundColor = DT::styleEqual(
          c("Sequence Unavailable", "Failed to Map"),
          c('lightcoral', 'lightyellow')
        )
      )
  })
  
  # Summary
  output$retrievalSummary <- renderText({
    req(retrieval_summary_val())
    s <- retrieval_summary_val()
    paste0(
      "Total Genes Processed: ", s$total, "\n",
      "Successfully Retrieved: ", s$retrieved, "\n",
      "Issues: ", s$missing
    )
  })
  
  # Specific Genes FASTA and BED
  observeEvent(input$downloadSpecificBtn, {
    req(input$specific_genes)
    gene_names <- unlist(strsplit(input$specific_genes, "\\s+"))
    gene_names <- toupper(gene_names)
    
    sel <- tryCatch({
      results() %>%
        dplyr::mutate(Gene_upper = toupper(Gene)) %>%
        dplyr::filter(Gene_upper %in% gene_names) %>%
        dplyr::filter(Status == "Sequence Retrieved")
    }, error = function(e) {
      showNotification("Error filtering specific genes.", type = "error")
      return(NULL)
    })
    
    if (is.null(sel) || nrow(sel) == 0) {
      showNotification("No matching genes found.", type = "error")
      return(NULL)
    }
    
    # Generate FASTA data for specific genes
    fasta_selected <- sel %>%
      dplyr::select(Gene, FASTA)
    
    specific_fasta_data(fasta_selected)
    
    # Generate BED entries for specific genes
    bed_all <- bed_data()
    req(bed_all)
    
    # Filter BED entries for selected genes
    if (input$separate_promoter_exon && input$include_promoter) {
      specific_bed <- bed_all %>%
        dplyr::filter(grepl(paste(sel$Gene, collapse = "|"), name))
    } else {
      specific_bed <- bed_all %>%
        dplyr::filter(name %in% sel$Gene)
    }
    
    specific_bed_data(specific_bed)
    
    showNotification("Selected genes' FASTA & BED prepared.", type = "message")
  })
  
  # Display specific genes FASTA table
  output$specificFastaTable <- DT::renderDataTable({
    req(specific_fasta_data())
    df <- specific_fasta_data()
    DT::datatable(
      df,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      colnames = c("Gene Symbol/ID", "FASTA Sequence")
    )
  })
  
  # Display specific genes BED table
  output$specificBedTable <- DT::renderDataTable({
    req(specific_bed_data())
    df <- specific_bed_data()
    
    # Determine BED format based on BED checkboxes
    bed6_selected <- !(input$include_exon || input$include_promoter || input$separate_promoter_exon)
    
    if (bed6_selected) {
      DT::datatable(
        df[, c("chrom", "chromStart", "chromEnd", "name", "score", "strand")],
        options = list(pageLength = 25, scrollX = TRUE),
        rownames = FALSE,
        colnames = c("Chromosome", "Start", "End", "Name", "Score", "Strand")
      )
    } else {
      # BED12 format
      DT::datatable(
        df[, c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")],
        options = list(pageLength = 25, scrollX = TRUE),
        rownames = FALSE,
        colnames = c("Chromosome", "Start", "End", "Name", "Score", "Strand",
                     "Thick Start", "Thick End", "Item RGB", "Block Count", "Block Sizes", "Block Starts")
      )
    }
  })
  
  # Download Specific Genes FASTA
  output$downloadSpecificFASTA <- downloadHandler(
    filename = function() {
      paste("selected_fasta_sequences_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      req(specific_fasta_data())
      fasta_sequences <- specific_fasta_data()
      
      complete_fasta <- c()
      for (i in seq_len(nrow(fasta_sequences))) {
        if (!is.na(fasta_sequences$FASTA[i])) {
          lines_i <- unlist(strsplit(fasta_sequences$FASTA[i], "\n"))
          complete_fasta <- c(complete_fasta, lines_i, "")
        }
      }
      writeLines(complete_fasta, con = file)
    }
  )
  
 
  
  # BED Results Table Rendering
  output$bedTable <- DT::renderDataTable({
    req(bed_data())
    df <- bed_data()
    
    # Determine BED format based on BED checkboxes
    bed6_selected <- !(input$include_exon || input$include_promoter || input$separate_promoter_exon)
    
    if (bed6_selected) {
      DT::datatable(
        df[, c("chrom", "chromStart", "chromEnd", "name", "score", "strand")],
        options = list(pageLength = 25, scrollX = TRUE),
        rownames = FALSE,
        colnames = c("Chromosome", "Start", "End", "Name", "Score", "Strand")
      )
    } else {
      # BED12 format
      DT::datatable(
        df[, c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")],
        options = list(pageLength = 25, scrollX = TRUE),
        rownames = FALSE,
        colnames = c("Chromosome", "Start", "End", "Name", "Score", "Strand",
                     "Thick Start", "Thick End", "Item RGB", "Block Count", "Block Sizes", "Block Starts")
      )
    }
  })
  
  # Download Genes BED
  output$downloadBED <- downloadHandler(
    filename = function() {
      paste("genes_bed_", Sys.Date(), ".bed", sep = "")
    },
    content = function(file) {
      req(bed_data())
      df <- bed_data()
      
      # Determine BED format based on BED checkboxes
      bed6 <- !(input$include_exon || input$include_promoter || input$separate_promoter_exon)
      
      if (bed6) {
        # BED6 format
        write.table(
          df[, c("chrom", "chromStart", "chromEnd", "name", "score", "strand")],
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE
        )
      } else {
        # BED12 format
        required_bed_cols <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
        if (!all(required_bed_cols %in% colnames(df))) {
          showNotification(" BED data is missing required columns for BED12 format.", type = "error")
          return(NULL)
        }
        
        write.table(
          df[, required_bed_cols],
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE
        )
      }
    }
  )
  
  # Visualization
  observeEvent(input$visualizeBtn, {
    req(input$visual_genes)
    output$exonPlot <- renderPlot({
      genes_to_plot <- unlist(strsplit(input$visual_genes, "\\s+"))
      genes_to_plot <- toupper(gene_names <- genes_to_plot[genes_to_plot != ""])
      
      df <- results()
      req(df)
      
      df_plot <- tryCatch({
        df %>%
          dplyr::mutate(Gene_upper = toupper(Gene)) %>%
          dplyr::filter(Status == "Sequence Retrieved") %>%
          dplyr::filter(Gene_upper %in% genes_to_plot)
      }, error = function(e) {
        showNotification("Error processing visualization data.", type = "error")
        return(NULL)
      })
      
      if (is.null(df_plot) || nrow(df_plot) == 0) {
        plot.new()
        title("No matching retrieved genes found.")
        return()
      }
      
      df_plot <- df_plot %>%
        dplyr::mutate(
          StrandSymbol = dplyr::case_when(
            Strand ==  1 ~ "+",
            Strand == -1 ~ "-",
            TRUE         ~ "."
          )
        )
      
      par(mar = c(5, 5, 4, 2) + 0.1)
      
      n_genes    <- nrow(df_plot)
      max_length <- max(df_plot$End - df_plot$Start + 1, na.rm = TRUE)
      
      plot(
        1, 1,
        type = "n",
        xlim = c(0, max_length),
        ylim = c(0.5, n_genes + 0.5),
        xlab = "Relative Position (Extended Region)",
        ylab = "",
        yaxt = "n"
      )
      title("Exon/Promoter Visualization")
      axis(
        side = 2, 
        at = seq_len(n_genes), 
        labels = paste0(df_plot$Gene, " (", df_plot$StrandSymbol, ")"),
        las = 1
      )
      
      # Define colors for features
      feature_colors <- c("promoter" = "red", "exon" = "blue")
      
      for (i in seq_len(n_genes)) {
        row_i <- df_plot[i, ]
        gene_length <- (row_i$End - row_i$Start + 1)
        y_center <- i
        
        # Draw main line
        segments(
          x0 = 0, 
          x1 = gene_length, 
          y0 = y_center, 
          y1 = y_center, 
          col = "gray50", lwd = 2
        )
        
        # Draw directionality arrows
        if (row_i$StrandSymbol == "+") {
          arrows(x0 = gene_length - 10, y0 = y_center, x1 = gene_length, y1 = y_center,
                 length = 0.1, angle = 30, code = 2, col = "black", lwd = 2)
        } else if (row_i$StrandSymbol == "-") {
          arrows(x0 = 10, y0 = y_center, x1 = 0, y1 = y_center,
                 length = 0.1, angle = 30, code = 2, col = "black", lwd = 2)
        }
        
        sub_list <- row_i$block_positions_list[[1]]
        exons_sorted <- sub_list[sapply(sub_list, function(b) b$type == "exon")]
        exons_sorted <- exons_sorted[order(sapply(exons_sorted, function(b) b$coords[1]))]
        
        # Plot exons
        for (block in exons_sorted) {
          block_type <- block$type
          block_coords <- block$coords
          if (length(block_coords) != 2 || any(is.na(block_coords))) next
          
          # Determine color based on block type
          block_color <- feature_colors[block_type]
          if (is.na(block_color)) {
            block_color <- "black"  #Default
          }
          
          block_start_rel <- block_coords[1] - row_i$Start + 1
          block_end_rel   <- block_coords[2] - row_i$Start + 1
          if (block_start_rel < 0) block_start_rel <- 0
          if (block_end_rel   < 0) block_end_rel   <- 0
          if (block_start_rel > gene_length) next
          if (block_end_rel   > gene_length) block_end_rel <- gene_length
          
          rect(
            xleft   = block_start_rel,
            xright  = block_end_rel,
            ybottom = y_center - 0.25,
            ytop    = y_center + 0.25,
            border  = block_color,
            col     = block_color
          )
        }
        
        # Plot promoters
        promoters_sorted <- sub_list[sapply(sub_list, function(b) b$type == "promoter")]
        promoters_sorted <- promoters_sorted[order(sapply(promoters_sorted, function(b) b$coords[1]))]
        
        for (block in promoters_sorted) {
          block_type <- block$type
          block_coords <- block$coords
          if (length(block_coords) != 2 || any(is.na(block_coords))) next
          
          # Determine color based on block type
          block_color <- feature_colors[block_type]
          if (is.na(block_color)) {
            block_color <- "black" #Default
          }
          
          block_start_rel <- block_coords[1] - row_i$Start + 1
          block_end_rel   <- block_coords[2] - row_i$Start + 1
          if (block_start_rel < 0) block_start_rel <- 0
          if (block_end_rel   < 0) block_end_rel   <- 0
          if (block_start_rel > gene_length) next
          if (block_end_rel   > gene_length) block_end_rel <- gene_length
          
          rect(
            xleft   = block_start_rel,
            xright  = block_end_rel,
            ybottom = y_center - 0.25,
            ytop    = y_center + 0.25,
            border  = block_color,
            col     = block_color
          )
        }
        
        # Plot intron lines between exons
        if (input$include_exon && length(exons_sorted) > 1) {
          exons_sorted <- exons_sorted[order(sapply(exons_sorted, function(b) b$coords[1]))]
          for (j in 1:(length(exons_sorted) - 1)) {
            current_exon <- exons_sorted[[j]]$coords
            next_exon <- exons_sorted[[j + 1]]$coords
            current_end <- current_exon[2] - row_i$Start + 1
            next_start <- next_exon[1] - row_i$Start + 1
            # Draw intron line
            segments(
              x0 = current_end,
              y0 = y_center,
              x1 = next_start,
              y1 = y_center,
              col = "blue",
              lwd = 6,
              lty = 2
            )
          }
        }
      }
      
      # Plot Legend
      legend("topright", 
             legend = c("Promoter", "Exon", "Intron"),
             fill = c(feature_colors["promoter"], feature_colors["exon"], NA),
             border = c(feature_colors["promoter"], feature_colors["exon"], NA),
             lty = c(NA, NA, 2),
             lwd = c(NA, NA, 1),
             col = c(feature_colors["promoter"], feature_colors["exon"], "black"),
             pch = c(15, 15, NA),
             pt.cex = 1.5,
             bty = "n",
             title = "Genomic Features")
    
    })
  })
}

shinyApp(ui = ui, server = server)
