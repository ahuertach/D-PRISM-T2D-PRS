#' DPRISM Type 2 Diabetes PRS Calculator
#'
#' @description
#' Combines PRS scores from five ancestries (AFR, AMR, EAS, EUR, SAS) into a single meta-score
#' using ancestry-specific weights optimized for Type 2 Diabetes.
#'
#' @details
#' **IMPORTANT:** You must explicitly name all file path arguments to ensure the correct weights
#' are applied (e.g., \code{AFR_profile_path = "file.txt"}).
#'
#' @param targetcohort_ancestry Target ancestry: "AFR", "AMR", "EAS", "EUR", or "SAS".
#' @param output_file_path Optional. Path to save output file. If NULL, returns data frame only.
#' @param ... Forces subsequent arguments to be named.
#' @param AFR_profile_path Path to African ancestry score file.
#' @param AMR_profile_path Path to Admixed American ancestry score file.
#' @param EAS_profile_path Path to East Asian ancestry score file.
#' @param EUR_profile_path Path to European ancestry score file.
#' @param SAS_profile_path Path to South Asian ancestry score file.
#'
#' @return Data frame with columns: FID, IID, METASCORE, ZMETASCORE.
#' @export
dprism_t2dprs <- function(targetcohort_ancestry,
                          output_file_path = NULL,
                          ...,
                          AFR_profile_path,
                          AMR_profile_path,
                          EAS_profile_path,
                          EUR_profile_path,
                          SAS_profile_path) {

  # --- 1. Safety Check: if the user tries to put file paths in the `...` (by just listing them without names),stop ---
  if (length(list(...)) > 0) {
    stop("Error: You must provide the file paths as NAMED arguments. \nExample: dprism_t2dprs(..., AFR_profile_path = 'path/to/file', ...)")
  }

  # --- Definition of Ancestry-specific trained coefficients ---
  all_weights <- list(
    AFR = c(AFR=0.26258, AMR=0.05058, EAS=0.01167, EUR=0.25675, SAS=0.09701),
    AMR = c(AFR=0.06117, AMR=0.26515, EAS=0.23973, EUR=0.62115, SAS=0.0563),
    EAS = c(AFR=0.05231, AMR=0.0835,  EAS=0.55889, EUR=0.30485, SAS=0.05935),
    EUR = c(AFR=0.04457, AMR=0.05245, EAS=0.113,   EUR=0.62345, SAS=0.00089),
    SAS = c(AFR=0.10939, AMR=0.03381, EAS=0.26606, EUR=0.53909, SAS=0.32257)
  )

  # --- Setup ---
  if (!targetcohort_ancestry %in% names(all_weights)) {
    stop(paste("Invalid target ancestry. Must be one of:", paste(names(all_weights), collapse = ", ")))
  }

  # Get the weights for the target
  weights <- all_weights[[targetcohort_ancestry]]

  # create a NAMED list.
  # Makes sure 'AFR_profile_path' holds the file the user INTENDED to be AFR.
  file_paths <- list(
    AFR = AFR_profile_path,
    AMR = AMR_profile_path,
    EAS = EAS_profile_path,
    EUR = EUR_profile_path,
    SAS = SAS_profile_path
  )

  # --- 2. Load and Scale ---
  # We iterate through the NAMES (AFR, AMR, etc.) to keep track of identity
  data_list <- lapply(names(file_paths), function(ancestry_name) {
    path <- file_paths[[ancestry_name]]
    df <- read.table(path, header = TRUE, stringsAsFactors = FALSE)

    if (!all(c("FID", "IID", "SCORE") %in% names(df))) {
      stop(paste("File", path, "is missing required columns (FID, IID, SCORE)."))
    }

    df$SCALED_SCORE <- as.numeric(scale(df$SCORE, center = TRUE, scale = TRUE))
    return(df[, c("FID", "IID", "SCALED_SCORE")])
  })

  # Assign names to the list so data_list$AFR is definitely the AFR data
  names(data_list) <- names(file_paths)

  # --- 3. Consistency Check ---
  for (i in 2:length(data_list)) {
    if (!all(data_list[[1]]$FID == data_list[[i]]$FID) || !all(data_list[[1]]$IID == data_list[[i]]$IID)) {
      stop("FID/IID columns do not match or are not in the same order across all input files.")
    }
  }

  # --- 4. Linear Combination ---

  # We ensure the weights align with the data list.
  # "Give me the weights in the same order as my data list names"
  ordered_weights <- weights[names(data_list)]

  # Final Double Check:
  # If data_list is (AFR, AMR, EAS...), ordered_weights must be (0.26, 0.05, 0.01...)

  scaled_score_matrix <- do.call(cbind, lapply(data_list, function(df) df$SCALED_SCORE))
  combined_score_vector <- scaled_score_matrix %*% ordered_weights

  # --- 5. Post-combination Scaling ---
  scaled_combined_score_vector <- as.numeric(scale(combined_score_vector[, 1], center = TRUE, scale = TRUE))

  # --- 6. Output ---
  output_df <- data.frame(
    FID = data_list[[1]]$FID,
    IID = data_list[[1]]$IID,
    METASCORE = combined_score_vector[, 1],
    ZMETASCORE = scaled_combined_score_vector
  )

  if (!is.null(output_file_path)) {
    write.table(output_df, output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message(paste("Successfully wrote results to:", output_file_path))
  }

  return(output_df)
}
