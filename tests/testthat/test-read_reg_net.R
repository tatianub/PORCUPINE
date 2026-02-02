# Helper functions to create test files
setup_test_files <- function() {
    # Use testthat's standard test data directory
    test_dir <- file.path("tests/testdata")
    # Create directory if it doesn't exist
    if (!dir.exists(test_dir)) {
        dir.create(test_dir, recursive = TRUE)
    }
    # Create test data
    test_matrix <- matrix(rnorm(20), nrow = 4, ncol = 5)
    rownames(test_matrix) <- paste0("gene_", 1:4)
    colnames(test_matrix) <- paste0("sample_", 1:5)
    test_df <- data.frame(
        gene_id = paste0("gene_", 1:4),
        sample_1 = rnorm(4),
        sample_2 = rnorm(4),
        sample_3 = rnorm(4)
    )
    test_edges <- data.frame(
        reg = paste0("TF_", 1:6),
        tar = paste0("gene_", rep(1:3, 2)),
        score = runif(6)
    )
    # Create RData files
    net_norm <- test_matrix
    other_obj <- "not a matrix"
    save(net_norm, other_obj, file = file.path(test_dir, "network.RData"))
    edges_data <- test_edges
    save(edges_data, file = file.path(test_dir, "edges.RData"))
    # Create text files
    write.table(test_matrix, file.path(test_dir, "network.txt"), 
                sep = "\t", quote = FALSE, col.names = NA)
    write.table(test_df, file.path(test_dir, "network_with_rownames.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(test_edges, file.path(test_dir, "edges.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    # Create CSV and TSV files for network data
    write.table(test_matrix, file.path(test_dir, "network.tsv"), 
                sep = "\t", quote = FALSE, col.names = NA)
    write.table(test_matrix, file.path(test_dir, "network.csv"), 
                sep = ",", quote = FALSE, col.names = NA)
    # Create files with issues
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    save(empty_matrix, file = file.path(test_dir, "empty_network.RData"))
    bad_edges <- data.frame(regulator = "TF1", target = "gene1")
    save(bad_edges, file = file.path(test_dir, "bad_edges.RData"))
    # Create non-numeric text file
    char_data <- data.frame(
        gene = paste0("edge_", 1:3),
        sample1 = c("high", "low", "medium"),
        sample2 = c("up", "down", "same")
    )
    write.table(char_data, file.path(test_dir, "non_numeric.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    return(test_dir)
}

test_that("read_networks handles file existence validation", {
    expect_error(
        read_networks("nonexistent_file.txt"),
        "File does not exist:"
    )
})

test_that("read_networks handles unsupported file formats", {
    temp_file <- tempfile(fileext = ".xlsx")
    file.create(temp_file)
    on.exit(unlink(temp_file))
    expect_error(
        read_networks(temp_file),
        "Unsupported file format"
    )
})

test_that("read_networks reads RData files correctly", {
    temp_dir <- setup_test_files()
    # Test with specific object name
    result <- read_networks(
        file.path(temp_dir, "network.RData"),
        object_name = "net_norm"
    )
    expect_true(is.matrix(result))
    expect_true(is.numeric(result))
    expect_equal(nrow(result), 4)
    expect_equal(ncol(result), 5)
    # Test auto-detection
    result_auto <- read_networks(file.path(temp_dir, "network.RData"))
    expect_equal(result, result_auto)
})

test_that("read_networks handles RData file errors", {
    temp_dir <- setup_test_files()
    on.exit(unlink(temp_dir, recursive = TRUE))
    # Test non-existent object name
    expect_error(
        read_networks(
            file.path(temp_dir, "network.RData"),
            object_name = "nonexistent"
        ),
        "Object 'nonexistent' not found"
    )
    # Test empty network
    expect_error(
        read_networks(file.path(temp_dir, "empty_network.RData")),
        "Network matrix has zero dimensions"
    )
})

test_that("read_networks reads text files correctly", {
    temp_dir <- setup_test_files()
    on.exit(unlink(temp_dir, recursive = TRUE))

    # Test basic text file reading
    result <- read_networks(file.path(temp_dir, "network.txt"))
    expect_true(is.matrix(result))
    expect_true(is.numeric(result))
    expect_equal(nrow(result), 4)
    expect_equal(ncol(result), 5)
    # Test file with row names in first column
    result_rownames <- read_networks(
        file.path(temp_dir, "network_with_rownames.txt")
    )
    expect_true(is.matrix(result_rownames))
    expect_equal(nrow(result_rownames), 4)
    expect_equal(ncol(result_rownames), 3)
})

test_that("read_networks handles data conversion", {
    temp_dir <- setup_test_files()
    on.exit(unlink(temp_dir, recursive = TRUE))
    # Test non-convertible data throws error
    expect_error(
        read_networks(file.path(temp_dir, "non_numeric.txt")),
        "Data contains non-numeric values that 
            cannot be converted to numeric matrix"
    )
})

test_that("read_edges handles file existence validation", {
    expect_error(
        read_edges("nonexistent_file.txt"),
        "File does not exist:"
    )
})

test_that("read_edges handles unsupported file formats", {
    temp_file <- tempfile(fileext = ".json")
    file.create(temp_file)
    on.exit(unlink(temp_file))
    expect_error(
        read_edges(temp_file),
        "Unsupported file format"
    )
})

test_that("read_edges reads RData files correctly", {
    temp_dir <- setup_test_files()
    on.exit(unlink(temp_dir, recursive = TRUE))
    # Test with specific object name
    result <- read_edges(
        file.path(temp_dir, "edges.RData"),
        object_name = "edges_data"
    )
    expect_s3_class(result, "data.frame")
    expect_true(all(c("reg", "tar") %in% colnames(result)))
    expect_equal(nrow(result), 6)
    # Test auto-detection
    result_auto <- read_edges(file.path(temp_dir, "edges.RData"))
    expect_equal(result, result_auto)
})

test_that("read_edges validates required columns", {
    temp_dir <- setup_test_files()
    on.exit(unlink(temp_dir, recursive = TRUE))
    expect_error(
        read_edges(file.path(temp_dir, "bad_edges.RData")),
        "Edges file must contain columns 'reg' and 'tar'"
    )
})

test_that("read_edges reads text files correctly", {
    temp_dir <- setup_test_files()
    on.exit(unlink(temp_dir, recursive = TRUE))
    result <- read_edges(file.path(temp_dir, "edges.txt"))
    expect_s3_class(result, "data.frame")
    expect_true(all(c("reg", "tar") %in% colnames(result)))
    expect_equal(nrow(result), 6)
    expect_equal(ncol(result), 3)
})

test_that("read_edges handles empty files", {
    temp_file <- tempfile(fileext = ".txt")
    write.table(data.frame(reg = character(0), 
        tar = character(0)), temp_file, row.names = FALSE, )
    on.exit(unlink(temp_file))
    expect_error(
        read_edges(temp_file),
        "Edges file is empty"
    )
})

test_that("read_networks reads other extensions correctly", {
     temp_dir <- setup_test_files()
     # Test .tsv extension
    tsv_result <- read_networks(
        file.path(temp_dir, "network.tsv")
    )
    expect_true(is.matrix(tsv_result))
    expect_true(is.numeric(tsv_result))
    expect_equal(nrow(tsv_result), 4)
    expect_equal(ncol(tsv_result), 5)
    # Test .csv extension
    csv_result <- read_networks(
        file.path(temp_dir, "network.csv")
    )
    expect_true(is.matrix(csv_result))
    expect_true(is.numeric(csv_result))
    expect_equal(nrow(csv_result), 4)
    expect_equal(ncol(csv_result), 5)
})
