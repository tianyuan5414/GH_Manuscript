#' Run end-to-end PLS pipeline for spectra and p-coumaric acid response
 #'
 #' Builds a treatment-level dataset by combining spectra and response data,
 #' fits a PLS model, selects the optimal number of components via the 1-SE
 #' rule, and returns cross-validation and diagnostic outputs.
 #'
 #' @param spectra_matrix Numeric matrix of spectra with samples in rows and
 #'   wavenumbers in columns. Row names must be `unique_id` values.
#' @param pca_df Data frame containing response values and grouping variable.
#'   Must include the grouping column (default `sample_level`) and `pCA_ng_grain`.
#' @param sample_data Data frame containing sample metadata with at least
#'   `unique_id` and the grouping column (default `sample_level`).
#' @param group_var Character. Name of the grouping column used to aggregate
#'   spectra and response values (default: `"sample_level"`).
 #' @param crossval Character cross-validation strategy passed to
 #'   `run_pls_mod()`. Default is `"LOO"`.
 #'
 #' @return A named list with:
 #' \describe{
 #'   \item{pls_model}{PLS model object returned by `run_pls_mod()`.}
 #'   \item{selected_components}{Output from `pick_1se_comp()`, including
 #'   selected component count.}
 #'   \item{cv_results}{Cross-validated RMSEP object from [pls::RMSEP()].}
 #'   \item{diagnostics}{Model diagnostic table from `get_diagnostic_table()`.}
 #' }
 #'
 #' @details
 #' The function averages both spectra and response values to treatment/sample
 #' level before model fitting, then joins by `sample_level` (response) and
 #' `sample_treatment` (spectra metadata).
 #'
 #' @examples
 #' # res <- run_pls_pipeline(
 #' #   spectra_matrix = spectra_matrix,
 #' #   pca_df = pca_summary,
 #' #   sample_data = sample_metadata,
 #' #   crossval = "LOO"
 #' # )
run_pls_pipeline <- function(
  spectra_matrix,
  pca_df,
  sample_data,
  crossval = "LOO",
  group_var = "sample_level"
) {
  # Ensure rownames of spectra_matrix match unique_id in sample_data
  if (!all(rownames(spectra_matrix) %in% sample_data$unique_id)) {
    stop("Rownames of spectra_matrix must match unique_id in sample_data")
  }

  # validate grouping column and response column
  if (!is.character(group_var) || length(group_var) != 1) {
    stop("'group_var' must be a single character string naming the grouping column.")
  }
  if (!group_var %in% names(sample_data)) {
    stop(paste0("Grouping column '", group_var, "' not found in sample_data."))
  }
  if (!group_var %in% names(pca_df)) {
    stop(paste0("Grouping column '", group_var, "' not found in pca_df."))
  }
  if (!"pCA_ng_grain" %in% names(pca_df)) {
    stop("'pca_df' must contain a 'pCA_ng_grain' column.")
  }

  # Combine spectra with pCA data
  spectra_df <- as.data.frame(spectra_matrix) |>
    tibble::rownames_to_column("unique_id") |>
    dplyr::left_join(sample_data, by = "unique_id")

  # calculate mean spectra per grouping variable to match pCA data
  spectra_df <- spectra_df |>
    dplyr::group_by(.data[[group_var]]) |>
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::matches("^\\d+"),
        .fns = ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  # summarise pCA data to grouping level to match spectra data
  pca_df <- pca_df |>
    dplyr::group_by(.data[[group_var]]) |>
    dplyr::summarise(
      pCA_ng_grain = mean(pCA_ng_grain, na.rm = TRUE),
      .groups = "drop"
    )

  # Combine data for PLS regression, matching on sample level (e.g., treatment)
  combined_df <- pca_df |>
    dplyr::left_join(
      spectra_df,
      by = group_var
    )

  # Run PLS regression using the combined data
  mod_pls <- run_pls_mod(
    combined_df,
    pred_col = "pCA_ng_grain",
    crossval = crossval
  )

  # Pick components based on cross-validation
  mod_comp <- pick_1se_comp(mod_pls)

  # run cross-validated PLS regression with selected components
  mod_cv <- pls::RMSEP(
    mod_pls,
    estimate = "CV"
  )

  mod_diag <- get_diagnostic_table(
    mod_pls,
    nopt = mod_comp$nopt[1]
  )

  list(
    pls_model = mod_pls,
    selected_components = mod_comp,
    cv_results = mod_cv,
    diagnostics = mod_diag
  )

}
