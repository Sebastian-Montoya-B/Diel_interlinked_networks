packs<-c("emln",
         "igraph",
         "magrittr",
         "tidyverse",
         "tnet", 
         "dendextend")

for (i in 1:length(packs)){
  
  if(!require(packs[i],character.only=TRUE)){
    install.packages(packs[i])
  }
  library(packs[i],character.only=TRUE)
}





#------------------------------------------------------------------------------#
#                                                                              #
#  1. edge_list_to_bipartite: transforms an edgelist into a biadjacency matrix #
#                                                                              #
#------------------------------------------------------------------------------#

edge_list_to_bipartite <- function(edge_df) {
  edge_df %>%
    pivot_wider(names_from = to, values_from = weight, values_fill = 0) %>%
    column_to_rownames(var = "from") %>%
    as.matrix()
}


#------------------------------------------------------------------------------#
#                                                                              #
#  2. reorder_vector: reorders vectors to create interlayer links in a cycle   #
#                                                                              #
#------------------------------------------------------------------------------#

reorder_vector <- function(v) {
  if (length(v) > 1) {
    return(c(v[-1], v[1]))
  } else {
    return(v)  # If the vector has only one element, return it unchanged
  }
}



#------------------------------------------------------------------------------#
#                                                                              #
#                          3. as weighted biadjacency                          #
#                                                                              #
#------------------------------------------------------------------------------#


# Function to transform edge list into biadjacency matrix
as_weighted_biadjacency_matrix <- function(graph) {
  require(tidyverse)
  require(igraph)
  edge_list<-data.frame(
    as_edgelist(graph),
    weight= E(graph)$weight)
  
  # Get unique nodes for X1 and X2
  set1 <- unique(edge_list$X1)
  set2 <- unique(edge_list$X2)
  
  # Create a data frame for all possible combinations of nodes
  biadj_matrix_df <- expand.grid(X1 = set1, X2 = set2)
  
  # Merge the edge list with the combination of all pairs, adding weights
  biadj_matrix_df <- biadj_matrix_df %>%
    left_join(edge_list, by = c("X1", "X2")) %>%
    mutate(weight = ifelse(is.na(weight), 0, weight))  # Fill NA with 0
  
  # Pivot the data frame to wide format to create the biadjacency matrix
  biadj_matrix <- biadj_matrix_df %>%
    pivot_wider(names_from = X2, values_from = weight, values_fill = 0) %>%
    column_to_rownames(var = "X1") %>%
    as.matrix()
  
  return(biadj_matrix)
}


#------------------------------------------------------------------------------#
#                                                                              #
#                             4. pref_Shannon                                  #
#                                                                              #
#------------------------------------------------------------------------------#

pref_Shannon<-function(data,q){
  
  if (class(data)[1]!="matrix"){
    data<-matrix(data, ncol=length(data))
  }
  
  if (length(q)!=NCOL(data)){
    stop("The number of resources and their abundance vector do not match")
  }
  
  q<-q/sum(q)
  
  sha_fun<-function(p,q){
      p<-p/sum(p)
      pre<-p/q
      di<-ifelse(pre>0, (pre/sum(pre))*log(pre/sum(pre)), 0)
      d<--sum(di)
      return(exp(d))
  }

  
  psha<-apply(data, 1, sha_fun, q)
  return(psha)
}

#------------------------------------------------------------------------------#
#                                                                              #
#                5. Generate null models based on Sauve et al                  #
#                                                                              #
#------------------------------------------------------------------------------#
generate_null_models <- function(data, columns_to_shuffle = NULL, num_models = 1) {
  # If no columns specified, shuffle all columns
  if (is.null(columns_to_shuffle)) {
    columns_to_shuffle <- colnames(data)
  }
  
  # Convert matrix to data frame if needed
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  
  # Validate columns exist
  if (!all(columns_to_shuffle %in% colnames(data))) {
    stop("Some specified columns don't exist in the data frame")
  }
  
  # Handle single-row case
  if (nrow(data) == 1) {
    warning("Data has only 1 row - shuffling will have no effect")
    return(replicate(num_models, data, simplify = FALSE))
  }
  
  # Generate null models
  null_models <- lapply(1:num_models, function(i) {
    shuffled_data <- data
    for (col in columns_to_shuffle) {
      shuffled_data[[col]] <- sample(shuffled_data[[col]])
    }
    return(shuffled_data)
  })
  
  return(null_models)
}


#------------------------------------------------------------------------------#
#                                                                              #
#                             6. partial correlation                           #
#                                                                              #
#------------------------------------------------------------------------------#
calculate_partial_cor <- function(degree_list, 
                                  ant_col = "antagonistic", 
                                  mut_col = "mutualistic", 
                                  abund_col = "abundance",
                                  method = "spearman",
                                  verbose = TRUE) {
  
  # Check if ppcor package is available
  if (!requireNamespace("ppcor", quietly = TRUE)) {
    stop("Package 'ppcor' is required but not installed. Please install it first.")
  }
  
  # Validate input
  if (!is.list(degree_list)) stop("Input must be a list")
  if (length(degree_list) == 0) stop("Input list is empty")
  
  # Initialize empty list to store results
  result_list <- list()
  
  for (i in seq_along(degree_list)) {
    net_name <- names(degree_list)[i] %||% as.character(i)
    current_df <- degree_list[[i]]
    
    # Create a default result template
    result_template <- data.frame(
      name = net_name,
      correlation = NA_real_,
      n_species = NA_integer_,
      status = NA_character_,
      message = NA_character_,
      stringsAsFactors = FALSE
    )
    
    tryCatch({
      # Check if required columns exist
      if (!all(c(ant_col, mut_col, abund_col) %in% names(current_df))) {
        stop("Missing required columns")
      }
      
      # Extract and prepare data
      data <- data.frame(
        antagonistic = current_df[[ant_col]],
        mutualistic = current_df[[mut_col]],
        log_abundance = log1p(current_df[[abund_col]])  # Using log1p to handle zeros
      )
      
      # Remove rows with NA values
      complete_data <- na.omit(data)
      n <- nrow(complete_data)
      
      # Update species count
      result_template$n_species <- n
      
      # Check for constant variables
      variances <- apply(complete_data, 2, var, na.rm = TRUE)
      if (any(variances == 0, na.rm = TRUE)) {
        result_template$status <- "warning"
        result_template$message <- "Constant variable(s) found"
        if (verbose) message("Warning in ", net_name, ": Constant variable(s) found")
        result_list[[i]] <- result_template
        next
      }
      
      # Check for sufficient sample size
      p <- 1  # Number of variables being controlled for
      if (n <= p + 2) {
        result_template$status <- "warning"
        result_template$message <- paste("Insufficient sample size (n =", n, ")")
        if (verbose) message("Warning in ", net_name, ": Insufficient sample size (n = ", n, ")")
        result_list[[i]] <- result_template
        next
      }
      
      # Compute partial correlation
      pcor_result <- ppcor::pcor(complete_data, method = method)
      cor_value <- pcor_result$estimate[2, 1]  # Correlation between ant and mut
      
      # Prepare successful result
      result_template$correlation <- cor_value
      result_template$status <- "success"
      result_template$message <- "Calculation successful"
      result_list[[i]] <- result_template
      
    }, error = function(e) {
      result_template$status <- "error"
      result_template$message <- gsub("\n", " ", e$message)
      if (verbose) message("Error in ", net_name, ": ", e$message)
      result_list[[i]] <<- result_template
    })
  }
  
  # Combine results safely
  tryCatch({
    final_result <- do.call(rbind, result_list)
    rownames(final_result) <- NULL
    return(final_result)
  }, error = function(e) {
    if (verbose) message("Error combining results: ", e$message)
    return(result_list)  # Return as list if binding fails
  })
}

# Helper function for NULL coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a


calculate_null_partial_cor <- function(null_model_list, 
                                       ant_col = "antagonistic",
                                       mut_col = "mutualistic",
                                       abund_col = "abundance",
                                       method = "spearman",
                                       verbose = TRUE) {
  
  # Process each null model in the list
  results <- lapply(null_model_list, function(single_null_model) {
    tryCatch({
      # Calculate partial correlations for all replicates in this null model
      replicate_results <- lapply(single_null_model, function(df) {
        calculate_partial_cor(list(temp = df),  # Wrap in list to work with our function
                              ant_col, mut_col, abund_col, 
                              method, verbose = FALSE)
      })
      
      # Combine all replicates
      do.call(rbind, replicate_results)
    }, error = function(e) {
      if (verbose) message("Error processing null model: ", e$message)
      NULL
    })
  })
  
  # Name the results
  names(results) <- names(null_model_list)
  
  # Return only non-NULL results
  results[!sapply(results, is.null)]
}
