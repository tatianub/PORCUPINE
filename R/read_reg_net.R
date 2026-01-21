#' Read Regulatory Network File
#' 
#' Reads regulatory network data from either RData files (containing R object)
#' or text files (tab-delimited). The function automatically detects the file
#' type and handles extraction appropriately.
#'
#' @param reg_net_file Character string, path to the regulatory network file.
#'   Can be either .RData/.rda file or a text file (.txt, .tsv, .csv)
#' @param object_name Character string, name of the object to extract from 
#'   RData files. If NULL, will attempt to find a matrix/data.frame object
#'   automatically. Ignored for text files.
#' @param ... Additional arguments passed to fread() for text files
#'
#' @return A numeric matrix where rows are features (genes/edges) and columns
#'   are samples
#'
#' @examples
#' \dontrun{
#' # Read from RData file
#' net1 <- read_networks("network_data.RData", object_name = "net_norm")
#' 
#' # Read from text file  
#' net2 <- read_networks("network_data.txt")
#'
#' }
#' @export
read_networks <- function(reg_net_file,
                          object_name = NULL,
                          ...) {
    # Input validation
    if (!file.exists(reg_net_file)) {
        stop("File does not exist: ", reg_net_file)
    }
    # Determine file type based on extension
    file_ext <- tolower(tools::file_ext(reg_net_file))
    if (file_ext %in% c("rdata", "rda")) {
        # Handle RData files
        net <- read_rdata_network(reg_net_file, object_name)
    } else if (file_ext %in% c("txt", "tsv", "csv")) {
        # Handle text files
        net <- read_text_network(reg_net_file, ...)
    } else {
        stop("Unsupported file format. Use .RData, .rda, .txt, .tsv, or .csv")
    }
    # Ensure it's a numeric matrix
    if (!is.numeric(net)) {
        tryCatch({
            net <- as.matrix(net)
            storage.mode(net) <- "numeric"
        }, error = function(e) {
            stop("Could not convert network data to numeric matrix: ", 
                 e$message)
        })
    }
    # Validate dimensions
    if (nrow(net) == 0 || ncol(net) == 0) {
        stop("Network matrix has zero dimensions")
    }
    message("Loaded network with ", nrow(net), " features and ", 
            ncol(net), " samples")
    return(net)
}

#' Helper function to read RData network files
#' @keywords internal
read_rdata_network <- function(file_path, object_name = NULL) {
    # Load RData file into a new environment
    env <- new.env()
    loaded_objects <- load(file_path, envir = env)
    if (length(loaded_objects) == 0) {
        stop("No objects found in RData file")
    }
    # If object name is specified, extract it directly
    if (!is.null(object_name)) {
        if (!object_name %in% loaded_objects) {
            stop("Object '", object_name, "' not found in RData file. ",
                 "Available objects: ", 
                 paste(loaded_objects, collapse = ", "))
        }
        net <- get(object_name, envir = env)
    } else {
        # Auto-detect: find the first matrix or data.frame
        candidates <- c()
        for (obj_name in loaded_objects) {
            obj <- get(obj_name, envir = env)
            if (is.matrix(obj) || is.data.frame(obj)) {
                candidates <- c(candidates, obj_name)
            }
        }
        if (length(candidates) == 0) {
            stop("No matrix or data.frame objects found in RData file. ",
                 "Available objects: ", 
                 paste(loaded_objects, collapse = ", "))
        }
        if (length(candidates) > 1) {
            warning("Multiple matrix/data.frame objects found: ", 
                    paste(candidates, collapse = ", "), 
                    ". Using: ", candidates[1])
        }
        net <- get(candidates[1], envir = env)
        message("Extracted object: ", candidates[1])
    }
    return(net)
}

#' Helper function to read text network files
#' @keywords internal
read_text_network <- function(file_path, ...) {
    # Try to read with data.table::fread for speed
    if (requireNamespace("data.table", quietly = TRUE)) {
        net <- data.table::fread(file_path, ...)
    } else {
        # Fallback to base R
        net <- read.table(file_path, header = TRUE, sep = "\t", 
                         stringsAsFactors = FALSE, ...)
    }
    # Convert to matrix if it's a data.frame
    if (is.data.frame(net)) {
        # Check if first column contains row names
        if (is.character(net[[1]]) && !any(duplicated(net[[1]]))) {
            rownames(net) <- net[[1]]
            net <- net[, -1, drop = FALSE]
        }
        net <- as.matrix(net)
    }
    return(net)
}