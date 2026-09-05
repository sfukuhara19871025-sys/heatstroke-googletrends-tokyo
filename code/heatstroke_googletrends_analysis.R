# =============================================================================
# Step 2 (FINAL INTEGRATED): Year-specific HAC + random-effects pooling + four pre-submission revisions
# =============================================================================
# Study:
#   Tokyo, May-September 2021-2025
#   Primary Google Trends query: "熱中症"
#
# IMPORTANT:
#   - This script DOES NOT download Google Trends data.
#   - The final raw Google Trends source is raw_google_trends.zip.
#   - The ZIP is validated, extracted automatically, and all valid session_* folders are read.
#   - If multiple sessions were acquired on the same calendar date, they are
#     first combined within that acquisition date so that one date does not
#     receive extra statistical weight.
#   - This FINAL script requires exactly 10 unique acquisition dates.
#   - No parametric AR/ARMA error model is selected in the primary analysis.
#
# Main workflow:
#   A. Reconstruct each session-year independently from overlapping windows.
#   B. Standardize the reconstructed series WITHIN each session and study year.
#   C. If >1 session exists on the same acquisition date, take the daily median
#      within that acquisition date and re-standardize within study year.
#   D. Across acquisition dates, take the daily median and re-standardize within
#      study year. This becomes the primary multi-session outcome.
#   E. Fit each study year as a separate OLS time-series regression with the same
#      prespecified fixed-effect structure.
#   F. Use Newey-West HAC covariance (primary lag 14; sensitivity lags 7 and 21)
#      for year-specific inference without requiring residual whitening.
#   G. Pool the five year-specific Alert coefficients using REML random-effects
#      with Knapp-Hartung inference as the PRIMARY pooled analysis.
#   H. Re-run acquisition-date robustness and covariate-specification sensitivity
#      analyses using the same HAC + random-effects framework.
#
#
# MANUSCRIPT-TABLE FIX (v3):
#   - Table 1 now exactly follows the manuscript structure (2021-2025 + Overall).
#   - Table 2 now combines year-specific and pooled estimates in one output.
#   - The obsolete separate manuscript "Table 3" output has been removed.
#   - Numeric audit copies are saved under diagnostics/.
#
# INTEGRATED PRE-SUBMISSION ADDITIONS:
#   I.   Explicit temperature-adjustment flexibility trend (linear, df=2, df=3, df=4)
#   II.  Leave-one-study-year-out influence analysis
#   III. Post hoc exploratory +/-1 calendar-day Alert timing model
#   IV.  GitHub/Zenodo-safe public-release package using an explicit whitelist
#
# SUPPLEMENTARY-OUTPUT FIX (v4):
#   - Former Supplementary Table S7 (temperature-curve source data) is no longer
#     treated as a manuscript table; its numeric source data are kept in diagnostics/.
#   - The redundant temperature-flexibility trend table (formerly treated as an
#     extra/S11-type output) is not a Supplementary Table; source data are kept in diagnostics/.
#   - Supplementary Tables are renumbered consecutively S1-S9.
#   - All Supplementary Tables/Figures are written directly to tables/ and figures/.
#   - All formal supplementary outputs are stored directly under tables/ and figures/.
#
# These additions do NOT redefine the prespecified primary analysis.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
required_packages <- c(
  "dplyr", "purrr", "readr", "tibble", "lubridate", "stringr",
  "tidyr", "ggplot2", "splines", "readxl", "MASS",
  "sandwich", "lmtest", "metafor"
)

installed <- rownames(installed.packages())
need_install <- setdiff(required_packages, installed)
if (length(need_install) > 0L) {
  install.packages(need_install)
}

# Do NOT attach MASS because MASS::select masks dplyr::select.
packages_to_attach <- setdiff(required_packages, "MASS")
invisible(lapply(packages_to_attach, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------
# Final Google Trends raw source.
# Recommended: put raw_google_trends.zip in the same project folder as this script,
# data_f(1).txt, and the heat-stroke alert Excel file.
# If needed, specify an absolute ZIP path before source():
#   Sys.setenv(HEATSTROKE_GT_ZIP = "/absolute/path/raw_google_trends.zip")
RAW_GOOGLE_TRENDS_ZIP_EXPLICIT <- Sys.getenv("HEATSTROKE_GT_ZIP", unset = "")

# Optional fallback: use an already-extracted raw_google_trends directory instead.
#   Sys.setenv(HEATSTROKE_GT_RAW_DIR = "/absolute/path/raw_google_trends")
RAW_GOOGLE_TRENDS_DIR_EXPLICIT <- Sys.getenv("HEATSTROKE_GT_RAW_DIR", unset = "")

WEATHER_FILE_EXPLICIT <- Sys.getenv("HEATSTROKE_WEATHER_FILE", unset = "")
ALERT_FILE_EXPLICIT <- Sys.getenv("HEATSTROKE_ALERT_FILE", unset = "")

GT_ZIP_CANDIDATE_NAMES <- c(
  "raw_google_trends.zip",
  "raw_google_trends(1).zip"
)

WEATHER_CANDIDATE_NAMES <- c(
  "data_f.txt",
  "data_f(1).txt",
  "data_f(2).txt",
  "data_f_revise.txt"
)
ALERT_CANDIDATE_NAMES <- c(
  "check_heatstroke_alert_2021_2025_wide_by_region.xlsx",
  "check_heatstroke_alert_2021_2025_wide_by_region(1).xlsx",
  "check_heatstroke_alert_2021_2025_wide_by_region(2).xlsx",
  "check_heatstroke_alert_2021_2025_wide_by_region(3).xlsx"
)

# Paths are resolved later, before the analysis starts.
RAW_SOURCE_TYPE <- NA_character_
RAW_SOURCE_ZIP <- NA_character_
RAW_DIR <- NA_character_
PROJECT_DIR <- NA_character_
WEATHER_FILE <- NA_character_
ALERT_FILE <- NA_character_
ANALYSIS_OUT <- NA_character_
DIR_TABLES <- NA_character_
DIR_FIGURES <- NA_character_
DIR_DIAGNOSTICS <- NA_character_

# Optional: analyze only selected sessions.
# NULL = use all session_* folders found in the ZIP/extracted RAW_DIR.
# Example:
# TARGET_SESSION_IDS <- c("20260802_022825", "20260803_090000")
TARGET_SESSION_IDS <- NULL

STUDY_YEARS <- 2021:2025
STUDY_START_MMDD <- "05-01"
STUDY_END_MMDD <- "09-30"
EXPECTED_DAYS_PER_YEAR <- 153L
EXPECTED_TOTAL_DAYS <- 765L
EXPECTED_WINDOWS_PER_YEAR <- 3L

# Primary model specification.
TEMP_DF_PRIMARY <- 3L
SEASON_DF_PRIMARY <- 5L

# Final analysis is locked to exactly 10 unique acquisition dates.
TARGET_ACQUISITION_DATES <- 10L

# The main temperature curve is restricted to observed central 95% of temperature
# values to avoid over-interpreting sparse extreme tails.
TEMP_CURVE_QUANTILES <- c(0.025, 0.975)

# Primary HAC specification.
# NeweyWest() uses the Bartlett kernel. We explicitly disable prewhitening so that
# "lag = 14" has a simple interpretation as a prespecified 14-day HAC truncation.
# A finite-sample n/(n-k) adjustment is applied.
HAC_PRIMARY_LAG <- 14L
HAC_SENSITIVITY_LAGS <- c(7L, 21L)
HAC_PREWHITE <- FALSE
HAC_FINITE_SAMPLE_ADJUST <- TRUE

# Primary second-stage pooling.
# Five year-specific estimates are combined with REML random-effects and
# Knapp-Hartung inference. A common/equal-effect inverse-variance result is
# reported only as a sensitivity analysis.
POOL_METHOD_PRIMARY <- "REML"
POOL_TEST_PRIMARY <- "knha"

# -----------------------------------------------------------------------------
# 2. General helpers
# -----------------------------------------------------------------------------
log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

# Locate the final raw Google Trends ZIP (preferred) or an extracted raw directory.
find_candidate_file <- function(explicit, candidate_names, label) {
  if (nzchar(trimws(explicit))) {
    p <- path.expand(explicit)
    if (!file.exists(p)) {
      stop(label, " not found: ", p, call. = FALSE)
    }
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  ancestors <- unique(c(wd, dirname(wd), dirname(dirname(wd))))
  dirs <- unique(c(
    ancestors,
    file.path(ancestors, "data"),
    file.path(ancestors, "input"),
    file.path(ancestors, "inputs"),
    path.expand("~/Downloads"),
    path.expand("~/Desktop")
  ))
  dirs <- dirs[dir.exists(dirs)]

  hits <- unique(unlist(lapply(dirs, function(d) file.path(d, candidate_names)), use.names = FALSE))
  hits <- hits[file.exists(hits)]

  # Prefer a file in the current working directory.
  if (length(hits) > 1L) {
    in_wd <- hits[normalizePath(dirname(hits), winslash = "/", mustWork = TRUE) == wd]
    if (length(in_wd) == 1L) hits <- in_wd
  }

  if (length(hits) == 1L) {
    return(normalizePath(hits, winslash = "/", mustWork = TRUE))
  }
  if (length(hits) > 1L) {
    stop(
      "Multiple candidate ", label, " files were found. Set an explicit path with HEATSTROKE_GT_ZIP.",
      call. = FALSE
    )
  }
  NA_character_
}

inspect_gt_zip <- function(zip_path) {
  listing <- utils::unzip(zip_path, list = TRUE)
  members <- listing$Name
  pattern <- paste0(
    "^raw_google_trends/session_[0-9]{8}_[0-9]{6}/",
    "gt_heatstroke_[0-9]{4}_window_[0-9]{2}\\.csv$"
  )
  data_members <- members[grepl(pattern, members)]

  if (length(data_members) == 0L) {
    stop(
      "The ZIP does not contain the expected raw_google_trends/session_*/gt_heatstroke_*.csv structure.",
      call. = FALSE
    )
  }

  session_id <- sub(
    "^raw_google_trends/(session_[0-9]{8}_[0-9]{6})/.*$",
    "\\1",
    data_members
  )
  file_name <- basename(data_members)
  year <- as.integer(sub("^gt_heatstroke_([0-9]{4})_window_.*$", "\\1", file_name))
  window_id <- as.integer(sub("^.*_window_([0-9]{2})\\.csv$", "\\1", file_name))

  detail <- tibble::tibble(
    member = data_members,
    session_id = session_id,
    acquisition_date = as.Date(substr(sub("^session_", "", session_id), 1, 8), format = "%Y%m%d"),
    year = year,
    window_id = window_id
  )

  expected_pairs <- tidyr::expand_grid(
    year = STUDY_YEARS,
    window_id = seq_len(EXPECTED_WINDOWS_PER_YEAR)
  )

  session_qc <- detail %>%
    dplyr::group_by(.data$session_id, .data$acquisition_date) %>%
    dplyr::summarise(
      n_files = dplyr::n(),
      n_years = dplyr::n_distinct(.data$year),
      n_year_window_pairs = dplyr::n_distinct(paste(.data$year, .data$window_id)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      complete = .data$n_files == length(STUDY_YEARS) * EXPECTED_WINDOWS_PER_YEAR &
        .data$n_years == length(STUDY_YEARS) &
        .data$n_year_window_pairs == length(STUDY_YEARS) * EXPECTED_WINDOWS_PER_YEAR
    )

  # Check exact expected year x window set in every session.
  set_ok <- vapply(unique(detail$session_id), function(sid) {
    observed <- detail %>%
      dplyr::filter(.data$session_id == sid) %>%
      dplyr::distinct(.data$year, .data$window_id) %>%
      dplyr::arrange(.data$year, .data$window_id)
    identical(observed$year, expected_pairs$year) &&
      identical(observed$window_id, expected_pairs$window_id)
  }, logical(1))
  session_qc$exact_expected_set <- unname(set_ok[session_qc$session_id])
  session_qc$complete <- session_qc$complete & session_qc$exact_expected_set

  if (any(!session_qc$complete)) {
    bad <- session_qc$session_id[!session_qc$complete]
    stop(
      "One or more sessions in raw_google_trends.zip are incomplete: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  n_dates <- dplyr::n_distinct(session_qc$acquisition_date)
  if (n_dates != TARGET_ACQUISITION_DATES) {
    stop(
      "FINAL ZIP must contain exactly ", TARGET_ACQUISITION_DATES,
      " unique acquisition dates, but found ", n_dates, ".",
      call. = FALSE
    )
  }

  list(detail = detail, session_qc = session_qc)
}

extract_gt_zip_once <- function(zip_path) {
  md5 <- unname(tools::md5sum(zip_path))
  extract_root <- file.path(
    dirname(zip_path),
    paste0(".gt_final_input_", substr(md5, 1, 12))
  )
  raw_dir <- file.path(extract_root, "raw_google_trends")

  if (!dir.exists(raw_dir)) {
    if (dir.exists(extract_root)) {
      unlink(extract_root, recursive = TRUE, force = TRUE)
    }
    dir.create(extract_root, recursive = TRUE, showWarnings = FALSE)
    log_message("Extracting final raw Google Trends ZIP: ", zip_path)
    utils::unzip(zip_path, exdir = extract_root)
  } else {
    log_message("Using previously extracted ZIP contents: ", raw_dir)
  }

  if (!dir.exists(raw_dir)) {
    stop(
      "ZIP extraction completed, but raw_google_trends directory was not found.",
      call. = FALSE
    )
  }

  list(
    raw_dir = normalizePath(raw_dir, winslash = "/", mustWork = TRUE),
    extract_root = normalizePath(extract_root, winslash = "/", mustWork = TRUE),
    md5 = md5
  )
}

resolve_raw_google_trends_source <- function() {
  # 1) Explicit extracted directory, if requested.
  if (nzchar(trimws(RAW_GOOGLE_TRENDS_DIR_EXPLICIT))) {
    p <- path.expand(RAW_GOOGLE_TRENDS_DIR_EXPLICIT)
    if (!dir.exists(p)) {
      stop("HEATSTROKE_GT_RAW_DIR not found: ", p, call. = FALSE)
    }
    return(list(
      source_type = "extracted_directory",
      zip_path = NA_character_,
      zip_md5 = NA_character_,
      raw_dir = normalizePath(p, winslash = "/", mustWork = TRUE),
      project_dir = normalizePath(dirname(p), winslash = "/", mustWork = TRUE)
    ))
  }

  # 2) ZIP is the preferred final source.
  zip_path <- find_candidate_file(
    RAW_GOOGLE_TRENDS_ZIP_EXPLICIT,
    GT_ZIP_CANDIDATE_NAMES,
    "raw Google Trends ZIP"
  )

  if (!is.na(zip_path)) {
    # Validate before extracting.
    zip_qc <- inspect_gt_zip(zip_path)
    extracted <- extract_gt_zip_once(zip_path)
    return(list(
      source_type = "zip",
      zip_path = zip_path,
      zip_md5 = extracted$md5,
      raw_dir = extracted$raw_dir,
      project_dir = normalizePath(dirname(zip_path), winslash = "/", mustWork = TRUE),
      zip_qc = zip_qc
    ))
  }

  # 3) Simple fallback: raw_google_trends directory in working directory.
  fallback <- file.path(getwd(), "raw_google_trends")
  if (dir.exists(fallback)) {
    return(list(
      source_type = "extracted_directory",
      zip_path = NA_character_,
      zip_md5 = NA_character_,
      raw_dir = normalizePath(fallback, winslash = "/", mustWork = TRUE),
      project_dir = normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    ))
  }

  stop(
    paste0(
      "Could not find raw_google_trends.zip. Put the ZIP in the working directory, or run:
",
      'Sys.setenv(HEATSTROKE_GT_ZIP = "/absolute/path/raw_google_trends.zip")'
    ),
    call. = FALSE
  )
}

resolve_named_input <- function(explicit, candidate_names, base_dir, label) {
  if (nzchar(trimws(explicit))) {
    p <- path.expand(explicit)
    if (!file.exists(p)) {
      stop(label, " file not found: ", p, call. = FALSE)
    }
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  project_dir <- dirname(base_dir)
  dirs <- unique(c(
    wd,
    base_dir,
    project_dir,
    file.path(wd, "data"),
    file.path(base_dir, "data"),
    file.path(project_dir, "data"),
    path.expand("~/Downloads"),
    path.expand("~/Desktop")
  ))
  dirs <- dirs[dir.exists(dirs)]

  hits <- unique(unlist(lapply(
    dirs,
    function(d) file.path(d, candidate_names)
  ), use.names = FALSE))
  hits <- hits[file.exists(hits)]

  if (length(hits) == 1L) {
    return(normalizePath(hits, winslash = "/", mustWork = TRUE))
  }

  if (length(hits) > 1L) {
    normalized_hits <- normalizePath(hits, winslash = "/", mustWork = TRUE)
    normalized_base <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
    base_hits <- hits[dirname(normalized_hits) == normalized_base]

    if (length(base_hits) == 1L) {
      return(normalizePath(base_hits, winslash = "/", mustWork = TRUE))
    }

    stop(
      paste0(
        "Multiple candidate ", label, " files were found. ",
        "Specify the desired file with an absolute path."
      ),
      call. = FALSE
    )
  }

  stop(
    paste0(
      label, " file was not found. Expected one of: ",
      paste(candidate_names, collapse = ", ")
    ),
    call. = FALSE
  )
}

median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}

safe_scale <- function(x, label = "series") {
  if (anyNA(x)) stop(label, " contains missing values before standardization.", call. = FALSE)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) {
    stop(label, " has zero/non-finite SD and cannot be standardized.", call. = FALSE)
  }
  as.numeric((x - mean(x)) / s)
}

expected_study_calendar <- purrr::map_dfr(STUDY_YEARS, function(y) {
  tibble::tibble(
    year = y,
    date = seq.Date(
      as.Date(paste0(y, "-", STUDY_START_MMDD)),
      as.Date(paste0(y, "-", STUDY_END_MMDD)),
      by = "day"
    ),
    day_of_season = seq_len(EXPECTED_DAYS_PER_YEAR)
  )
})

# -----------------------------------------------------------------------------
# 3. Session discovery and acquisition-date handling
# -----------------------------------------------------------------------------
parse_session_date <- function(session_id, session_dir) {
  first8 <- stringr::str_extract(session_id, "^[0-9]{8}")
  parsed <- suppressWarnings(as.Date(first8, format = "%Y%m%d"))
  if (!is.na(parsed)) return(parsed)

  # Fallback for manually named sessions.
  as.Date(file.info(session_dir)$mtime)
}

list_sessions <- function(raw_dir, target_ids = NULL) {
  stop_if_missing(raw_dir, "Step 1 raw Google Trends directory")

  dirs <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^session_", basename(dirs))]

  if (length(dirs) == 0L) {
    stop("No session_* directories were found under: ", raw_dir, call. = FALSE)
  }

  out <- tibble::tibble(
    session_dir = dirs,
    session_id = sub("^session_", "", basename(dirs))
  ) %>%
    dplyr::mutate(
      acquisition_date = as.Date(
        purrr::map2_chr(
          .data$session_id,
          .data$session_dir,
          ~ as.character(parse_session_date(.x, .y))
        )
      )
    ) %>%
    dplyr::arrange(.data$acquisition_date, .data$session_id)

  if (!is.null(target_ids)) {
    missing_ids <- setdiff(target_ids, out$session_id)
    if (length(missing_ids) > 0L) {
      stop(
        "TARGET_SESSION_IDS contains session IDs that were not found: ",
        paste(missing_ids, collapse = ", "),
        call. = FALSE
      )
    }
    out <- out %>% dplyr::filter(.data$session_id %in% target_ids)
  }

  out
}

# -----------------------------------------------------------------------------
# 4. Read one Step 1 session
# -----------------------------------------------------------------------------
read_session_windows <- function(session_id, session_dir) {
  files <- list.files(
    session_dir,
    pattern = "^gt_heatstroke_.*\\.csv$",
    full.names = TRUE
  )

  if (length(files) == 0L) {
    stop("No Step 1 Google Trends CSV files found in ", session_dir, call. = FALSE)
  }

  dat <- purrr::map_dfr(files, function(f) {
    readr::read_csv(f, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::mutate(source_file = basename(f))
  })

  needed <- c(
    "date", "hits", "year", "window_id", "window_start", "window_end"
  )
  missing_cols <- setdiff(needed, names(dat))
  if (length(missing_cols) > 0L) {
    stop(
      "Session ", session_id, " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  dat %>%
    dplyr::transmute(
      session_id = session_id,
      date = as.Date(.data$date),
      hits = as.numeric(.data$hits),
      year = as.integer(.data$year),
      window_id = as.integer(.data$window_id),
      window_start = as.Date(.data$window_start),
      window_end = as.Date(.data$window_end),
      source_file = .data$source_file
    ) %>%
    dplyr::arrange(.data$year, .data$window_id, .data$date)
}

# -----------------------------------------------------------------------------
# 5. Robust simultaneous overlap alignment within one session-year
# -----------------------------------------------------------------------------
reconstruct_session_year <- function(dat_year) {
  y <- unique(dat_year$year)
  sid <- unique(dat_year$session_id)

  if (length(y) != 1L || length(sid) != 1L) {
    stop("reconstruct_session_year() requires one session and one year.", call. = FALSE)
  }

  window_ids <- sort(unique(dat_year$window_id))
  if (length(window_ids) != EXPECTED_WINDOWS_PER_YEAR) {
    stop(
      "Session ", sid, ", year ", y, ": expected ",
      EXPECTED_WINDOWS_PER_YEAR, " windows but found ", length(window_ids), ".",
      call. = FALSE
    )
  }

  # One row per day in each window.
  window_data <- purrr::map(window_ids, function(w) {
    dat_year %>%
      dplyr::filter(.data$window_id == w) %>%
      dplyr::select(.data$date, .data$hits) %>%
      dplyr::distinct(.data$date, .keep_all = TRUE) %>%
      dplyr::arrange(.data$date)
  })
  names(window_data) <- as.character(window_ids)

  # Build all available overlap equations:
  # log(h_i) + a_i ~= log(h_j) + a_j
  # => log(h_i)-log(h_j) ~= a_j-a_i
  pair_list <- utils::combn(window_ids, 2, simplify = FALSE)

  overlap_equations <- purrr::map_dfr(pair_list, function(pair) {
    wi <- pair[1]
    wj <- pair[2]

    di <- window_data[[as.character(wi)]] %>%
      dplyr::rename(hits_i = .data$hits)
    dj <- window_data[[as.character(wj)]] %>%
      dplyr::rename(hits_j = .data$hits)

    dplyr::inner_join(di, dj, by = "date") %>%
      dplyr::filter(
        is.finite(.data$hits_i), is.finite(.data$hits_j),
        .data$hits_i > 0, .data$hits_j > 0
      ) %>%
      dplyr::transmute(
        date = .data$date,
        window_i = wi,
        window_j = wj,
        y_log_ratio = log(.data$hits_i) - log(.data$hits_j)
      )
  })

  if (nrow(overlap_equations) < 10L) {
    stop(
      "Session ", sid, ", year ", y,
      ": insufficient positive overlap observations for robust alignment (n=",
      nrow(overlap_equations), ").",
      call. = FALSE
    )
  }

  unknown_windows <- window_ids[-1]
  X <- matrix(
    0,
    nrow = nrow(overlap_equations),
    ncol = length(unknown_windows),
    dimnames = list(NULL, paste0("a_w", unknown_windows))
  )

  for (k in seq_along(unknown_windows)) {
    w <- unknown_windows[k]
    X[, k] <-
      as.numeric(overlap_equations$window_j == w) -
      as.numeric(overlap_equations$window_i == w)
  }

  if (qr(X)$rank < ncol(X)) {
    stop(
      "Session ", sid, ", year ", y,
      ": overlap design matrix is not connected/full rank.",
      call. = FALSE
    )
  }

  fit <- MASS::rlm(
    x = X,
    y = overlap_equations$y_log_ratio,
    psi = MASS::psi.huber,
    maxit = 100
  )

  coefs <- as.numeric(fit$coefficients)
  if (length(coefs) != length(unknown_windows) || any(!is.finite(coefs))) {
    stop(
      "Session ", sid, ", year ", y,
      ": robust overlap alignment failed to estimate finite scale factors.",
      call. = FALSE
    )
  }

  log_scale <- c(0, coefs)
  names(log_scale) <- as.character(window_ids)

  aligned_long <- dat_year %>%
    dplyr::mutate(
      log_scale = unname(log_scale[as.character(.data$window_id)]),
      multiplier = exp(.data$log_scale),
      aligned_hits = .data$hits * .data$multiplier
    )

  expected_dates <- expected_study_calendar %>%
    dplyr::filter(.data$year == y) %>%
    dplyr::select(.data$date)

  daily <- aligned_long %>%
    dplyr::group_by(.data$date) %>%
    dplyr::summarise(
      gt_reconstructed = median_or_na(.data$aligned_hits),
      n_windows_contributing = sum(is.finite(.data$aligned_hits)),
      .groups = "drop"
    ) %>%
    dplyr::right_join(expected_dates, by = "date") %>%
    dplyr::arrange(.data$date)

  if (nrow(daily) != EXPECTED_DAYS_PER_YEAR) {
    stop("Session ", sid, ", year ", y, ": reconstructed day count != 153.", call. = FALSE)
  }
  if (anyNA(daily$gt_reconstructed)) {
    stop("Session ", sid, ", year ", y, ": reconstructed series contains missing days.", call. = FALSE)
  }

  daily <- daily %>%
    dplyr::mutate(
      session_id = sid,
      year = y,
      gt_z_session = safe_scale(
        .data$gt_reconstructed,
        paste0("Session ", sid, ", year ", y)
      )
    )

  fitted_log_ratio <- as.numeric(X %*% coefs)
  alignment_diag <- tibble::tibble(
    session_id = sid,
    year = y,
    n_positive_overlap_equations = nrow(overlap_equations),
    median_abs_log_alignment_residual = stats::median(
      abs(overlap_equations$y_log_ratio - fitted_log_ratio)
    ),
    max_abs_log_alignment_residual = max(
      abs(overlap_equations$y_log_ratio - fitted_log_ratio)
    )
  )

  scale_table <- tibble::tibble(
    session_id = sid,
    year = y,
    window_id = window_ids,
    log_scale = unname(log_scale),
    multiplier = exp(unname(log_scale))
  )

  list(daily = daily, alignment_diag = alignment_diag, scale_table = scale_table)
}

# -----------------------------------------------------------------------------
# 6. Reconstruct a complete session safely
# -----------------------------------------------------------------------------
reconstruct_one_session <- function(session_row) {
  sid <- session_row$session_id[[1]]
  sdir <- session_row$session_dir[[1]]
  acq_date <- session_row$acquisition_date[[1]]

  log_message("Reading/reconstructing session: ", sid)

  raw <- tryCatch(
    read_session_windows(sid, sdir),
    error = function(e) e
  )

  if (inherits(raw, "error")) {
    return(list(
      success = FALSE,
      session_id = sid,
      acquisition_date = acq_date,
      error = conditionMessage(raw)
    ))
  }

  # Basic session completeness check.
  session_summary <- raw %>%
    dplyr::distinct(.data$year, .data$window_id) %>%
    dplyr::count(.data$year, name = "n_windows")

  if (!all(STUDY_YEARS %in% session_summary$year) ||
      any(session_summary$n_windows != EXPECTED_WINDOWS_PER_YEAR)) {
    return(list(
      success = FALSE,
      session_id = sid,
      acquisition_date = acq_date,
      error = "Session does not contain exactly 3 windows for every study year."
    ))
  }

  reconstructed <- tryCatch({
    pieces <- purrr::map(STUDY_YEARS, function(y) {
      reconstruct_session_year(raw %>% dplyr::filter(.data$year == y))
    })

    list(
      daily = purrr::map_dfr(pieces, "daily"),
      alignment_diag = purrr::map_dfr(pieces, "alignment_diag"),
      scale_table = purrr::map_dfr(pieces, "scale_table")
    )
  }, error = function(e) e)

  if (inherits(reconstructed, "error")) {
    return(list(
      success = FALSE,
      session_id = sid,
      acquisition_date = acq_date,
      error = conditionMessage(reconstructed)
    ))
  }

  daily <- reconstructed$daily %>%
    dplyr::mutate(acquisition_date = acq_date)

  if (nrow(daily) != EXPECTED_TOTAL_DAYS) {
    return(list(
      success = FALSE,
      session_id = sid,
      acquisition_date = acq_date,
      error = paste0("Reconstructed session has ", nrow(daily), " days instead of 765.")
    ))
  }

  list(
    success = TRUE,
    session_id = sid,
    acquisition_date = acq_date,
    daily = daily,
    alignment_diag = reconstructed$alignment_diag,
    scale_table = reconstructed$scale_table,
    error = NA_character_
  )
}

# -----------------------------------------------------------------------------
# 7. Meteorological and alert data
# -----------------------------------------------------------------------------
read_weather_data <- function(path) {
  stop_if_missing(path, "Weather file")

  weather_names <- c(
    "prefecture", "date",
    "b_mean_temp", "b_max_temp", "b_min_temp", "b_humidity", "b_rain", "b_wind",
    "p_mean_temp", "p_max_temp", "p_min_temp", "p_humidity", "p_rain", "p_wind"
  )

  dat <- readr::read_table(
    path,
    skip = 2,
    col_names = weather_names,
    col_types = readr::cols(
      prefecture = readr::col_character(),
      date = readr::col_date(format = "%Y-%m-%d"),
      .default = readr::col_double()
    ),
    locale = readr::locale(encoding = "UTF-8"),
    show_col_types = FALSE,
    progress = FALSE
  )

  tokyo <- dat %>%
    dplyr::filter(.data$prefecture %in% c("東京", "東京都")) %>%
    dplyr::transmute(
      date = .data$date,
      pop_max_temp = .data$p_max_temp,
      pop_humidity = .data$p_humidity,
      pop_rain = .data$p_rain,
      pop_wind = .data$p_wind,
      built_max_temp = .data$b_max_temp,
      built_humidity = .data$b_humidity,
      built_rain = .data$b_rain,
      built_wind = .data$b_wind
    ) %>%
    dplyr::arrange(.data$date)

  if (nrow(tokyo) != EXPECTED_TOTAL_DAYS) {
    stop("Expected 765 Tokyo weather rows, but found ", nrow(tokyo), ".", call. = FALSE)
  }
  if (anyDuplicated(tokyo$date)) stop("Duplicate dates in Tokyo weather data.", call. = FALSE)
  if (anyNA(tokyo)) stop("Missing values in Tokyo weather data.", call. = FALSE)

  tokyo
}

read_alert_data <- function(path) {
  stop_if_missing(path, "Alert file")

  sheets <- readxl::excel_sheets(path)
  sheet_use <- if ("Alert_Wide" %in% sheets) "Alert_Wide" else sheets[1]
  raw <- readxl::read_excel(path, sheet = sheet_use)

  tokyo_col <- intersect(c("東京", "東京都"), names(raw))
  if (length(tokyo_col) != 1L) {
    stop("Could not uniquely identify Tokyo alert column.", call. = FALSE)
  }
  if (!("date" %in% names(raw))) {
    stop("Column 'date' not found in alert file.", call. = FALSE)
  }

  out <- raw %>%
    dplyr::transmute(
      date = as.Date(.data$date),
      alert = as.integer(.data[[tokyo_col]])
    ) %>%
    dplyr::filter(.data$date %in% expected_study_calendar$date) %>%
    dplyr::arrange(.data$date)

  if (nrow(out) != EXPECTED_TOTAL_DAYS) {
    stop("Expected 765 Tokyo alert rows, but found ", nrow(out), ".", call. = FALSE)
  }
  if (anyDuplicated(out$date)) stop("Duplicate dates in Tokyo alert data.", call. = FALSE)
  if (anyNA(out)) stop("Missing values in Tokyo alert data.", call. = FALSE)
  if (!all(out$alert %in% c(0L, 1L))) stop("Alert variable must be 0/1.", call. = FALSE)

  out
}

# -----------------------------------------------------------------------------
# 8. Combine repeated sessions without pseudo-replicating same-day acquisitions
# -----------------------------------------------------------------------------
combine_sessions <- function(session_daily) {
  # A) Within acquisition date: median across same-day sessions.
  by_acquisition_date <- session_daily %>%
    dplyr::group_by(.data$acquisition_date, .data$year, .data$date) %>%
    dplyr::summarise(
      gt_z_date_raw = stats::median(.data$gt_z_session),
      n_sessions_same_acquisition_date = dplyr::n_distinct(.data$session_id),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$acquisition_date, .data$year) %>%
    dplyr::mutate(
      gt_z_date = safe_scale(
        .data$gt_z_date_raw,
        paste0("Acquisition date ", unique(.data$acquisition_date),
               ", study year ", unique(.data$year))
      )
    ) %>%
    dplyr::ungroup()

  # B) Across acquisition dates: median daily standardized value.
  aggregate_daily <- by_acquisition_date %>%
    dplyr::group_by(.data$year, .data$date) %>%
    dplyr::summarise(
      gt_z_multisession_raw = stats::median(.data$gt_z_date),
      gt_z_q25_raw = stats::quantile(.data$gt_z_date, 0.25, names = FALSE),
      gt_z_q75_raw = stats::quantile(.data$gt_z_date, 0.75, names = FALSE),
      gt_z_min_raw = min(.data$gt_z_date),
      gt_z_max_raw = max(.data$gt_z_date),
      n_acquisition_dates = dplyr::n_distinct(.data$acquisition_date),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$year) %>%
    dplyr::mutate(
      aggregate_year_mean = mean(.data$gt_z_multisession_raw),
      aggregate_year_sd = stats::sd(.data$gt_z_multisession_raw),
      gt_z = (.data$gt_z_multisession_raw - .data$aggregate_year_mean) / .data$aggregate_year_sd,
      gt_z_q25 = (.data$gt_z_q25_raw - .data$aggregate_year_mean) / .data$aggregate_year_sd,
      gt_z_q75 = (.data$gt_z_q75_raw - .data$aggregate_year_mean) / .data$aggregate_year_sd,
      gt_z_min = (.data$gt_z_min_raw - .data$aggregate_year_mean) / .data$aggregate_year_sd,
      gt_z_max = (.data$gt_z_max_raw - .data$aggregate_year_mean) / .data$aggregate_year_sd
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$date)

  if (any(!is.finite(aggregate_daily$gt_z))) {
    stop("Aggregated multi-session outcome could not be standardized.", call. = FALSE)
  }

  list(
    by_acquisition_date = by_acquisition_date,
    aggregate_daily = aggregate_daily
  )
}

# -----------------------------------------------------------------------------
# 9. Reproducibility: pairwise Spearman correlations across acquisition dates
# -----------------------------------------------------------------------------
calculate_pairwise_agreement <- function(by_acquisition_date) {
  acq_dates <- sort(unique(by_acquisition_date$acquisition_date))

  if (length(acq_dates) < 2L) {
    return(tibble::tibble(
      study_year = integer(),
      acquisition_date_1 = as.Date(character()),
      acquisition_date_2 = as.Date(character()),
      n_days = integer(),
      spearman_rho = numeric()
    ))
  }

  # Use index pairs so the Date class is never lost by combn().
  index_pairs <- utils::combn(seq_along(acq_dates), 2, simplify = FALSE)

  purrr::map_dfr(c(STUDY_YEARS, NA_integer_), function(y) {
    dat_y <- if (is.na(y)) {
      by_acquisition_date
    } else {
      by_acquisition_date %>% dplyr::filter(.data$year == y)
    }

    purrr::map_dfr(index_pairs, function(pair_index) {
      pair <- acq_dates[pair_index]
      d1 <- dat_y %>%
        dplyr::filter(.data$acquisition_date == pair[1]) %>%
        dplyr::select(.data$date, z1 = .data$gt_z_date)
      d2 <- dat_y %>%
        dplyr::filter(.data$acquisition_date == pair[2]) %>%
        dplyr::select(.data$date, z2 = .data$gt_z_date)

      joined <- dplyr::inner_join(d1, d2, by = "date") %>%
        dplyr::filter(is.finite(.data$z1), is.finite(.data$z2))

      rho <- if (nrow(joined) >= 10L &&
                 stats::sd(joined$z1) > 0 && stats::sd(joined$z2) > 0) {
        suppressWarnings(stats::cor(joined$z1, joined$z2, method = "spearman"))
      } else {
        NA_real_
      }

      tibble::tibble(
        study_year = ifelse(is.na(y), 0L, y),
        acquisition_date_1 = pair[1],
        acquisition_date_2 = pair[2],
        n_days = nrow(joined),
        spearman_rho = rho
      )
    })
  })
}

# -----------------------------------------------------------------------------
# 10. Prepare covariates for statistical model
# -----------------------------------------------------------------------------
prepare_analysis_data <- function(gt_daily, weather, alert) {
  dat <- gt_daily %>%
    dplyr::select(.data$date, .data$year, .data$gt_z) %>%
    dplyr::left_join(weather, by = "date") %>%
    dplyr::left_join(alert, by = "date") %>%
    dplyr::left_join(
      expected_study_calendar %>% dplyr::select(.data$date, .data$day_of_season),
      by = "date"
    ) %>%
    dplyr::arrange(.data$date) %>%
    dplyr::mutate(
      month = lubridate::month(.data$date),
      dow = factor(
        lubridate::wday(.data$date, week_start = 1, label = TRUE),
        ordered = FALSE
      ),
      year_f = factor(.data$year),
      pop_humidity_c = .data$pop_humidity - mean(.data$pop_humidity),
      pop_rain_log_c = log1p(.data$pop_rain) - mean(log1p(.data$pop_rain)),
      pop_wind_c = .data$pop_wind - mean(.data$pop_wind),
      built_humidity_c = .data$built_humidity - mean(.data$built_humidity),
      built_rain_log_c = log1p(.data$built_rain) - mean(log1p(.data$built_rain)),
      built_wind_c = .data$built_wind - mean(.data$built_wind)
    )

  if (nrow(dat) != EXPECTED_TOTAL_DAYS) {
    stop("Expected 765 analysis days, but found ", nrow(dat), ".", call. = FALSE)
  }

  required <- c(
    "gt_z", "alert", "pop_max_temp", "pop_humidity_c", "pop_rain_log_c",
    "pop_wind_c", "built_max_temp", "built_humidity_c", "built_rain_log_c",
    "built_wind_c", "dow", "year_f", "day_of_season"
  )
  if (anyNA(dat[required])) {
    stop("One or more required model variables contain missing values.", call. = FALSE)
  }

  dat
}

# -----------------------------------------------------------------------------
# 11. Build an aggregate series from any subset of acquisition dates
# -----------------------------------------------------------------------------
aggregate_selected_dates <- function(by_acquisition_date, selected_dates) {
  out <- by_acquisition_date %>%
    dplyr::filter(.data$acquisition_date %in% selected_dates) %>%
    dplyr::group_by(.data$year, .data$date) %>%
    dplyr::summarise(
      gt_z_raw = stats::median(.data$gt_z_date),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$year) %>%
    dplyr::mutate(
      gt_z = safe_scale(
        .data$gt_z_raw,
        paste0("Subset aggregate, study year ", unique(.data$year))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$date)

  out
}

# -----------------------------------------------------------------------------
# 12. YEAR-SPECIFIC OLS + NEWEY-WEST HAC HELPERS
# -----------------------------------------------------------------------------
# Each summer is a separate 153-day time series.
# The regression coefficients are ordinary least-squares estimates; uncertainty
# is calculated with a Bartlett-kernel Newey-West HAC covariance matrix.
# HAC inference does NOT require the regression residuals themselves to become
# white noise, so residual ACF/Ljung-Box statistics below are descriptive only.

resolve_weather_terms <- function(weather_source = c("population", "building")) {
  weather_source <- match.arg(weather_source)
  if (weather_source == "population") {
    list(
      temp = "pop_max_temp",
      humidity = "pop_humidity_c",
      rain = "pop_rain_log_c",
      wind = "pop_wind_c"
    )
  } else {
    list(
      temp = "built_max_temp",
      humidity = "built_humidity_c",
      rain = "built_rain_log_c",
      wind = "built_wind_c"
    )
  }
}

fit_year_hac_model <- function(dat_year,
                               hac_lag = HAC_PRIMARY_LAG,
                               temp_df = TEMP_DF_PRIMARY,
                               season_df = SEASON_DF_PRIMARY,
                               weather_source = c("population", "building"),
                               linear_temperature = FALSE,
                               model_label = "year-specific HAC model") {
  weather_source <- match.arg(weather_source)

  if (dplyr::n_distinct(dat_year$year) != 1L) {
    stop("fit_year_hac_model() requires exactly one study year.", call. = FALSE)
  }
  if (nrow(dat_year) != EXPECTED_DAYS_PER_YEAR) {
    stop(
      "Each year-specific HAC model must contain exactly ",
      EXPECTED_DAYS_PER_YEAR, " study days.", call. = FALSE
    )
  }
  if (!is.numeric(hac_lag) || length(hac_lag) != 1L ||
      !is.finite(hac_lag) || hac_lag < 0 || hac_lag >= nrow(dat_year)) {
    stop("Invalid HAC lag: ", hac_lag, call. = FALSE)
  }

  dat_model <- dat_year %>%
    dplyr::arrange(.data$day_of_season) %>%
    droplevels()

  wt <- resolve_weather_terms(weather_source)

  temp_term <- if (isTRUE(linear_temperature)) {
    wt$temp
  } else {
    paste0("splines::ns(", wt$temp, ", df = ", as.integer(temp_df), ")")
  }
  season_term <- paste0(
    "splines::ns(day_of_season, df = ", as.integer(season_df), ")"
  )

  rhs <- c(
    "alert",
    temp_term,
    season_term,
    wt$humidity,
    wt$rain,
    wt$wind,
    "dow"
  )

  formula_use <- stats::as.formula(
    paste("gt_z ~", paste(rhs, collapse = " + "))
  )

  fit <- stats::lm(
    formula = formula_use,
    data = dat_model,
    na.action = stats::na.fail,
    x = TRUE,
    y = TRUE,
    model = TRUE
  )

  # Explicit day ordering protects against accidental row-order changes.
  V_hac <- sandwich::NeweyWest(
    fit,
    lag = as.integer(hac_lag),
    order.by = dat_model$day_of_season,
    prewhite = HAC_PREWHITE,
    adjust = HAC_FINITE_SAMPLE_ADJUST
  )

  if (any(!is.finite(V_hac))) {
    stop(
      "Newey-West covariance contains non-finite values for year ",
      unique(dat_model$year), ".", call. = FALSE
    )
  }

  attr(fit, "hac_vcov") <- V_hac
  attr(fit, "hac_lag") <- as.integer(hac_lag)
  attr(fit, "study_year") <- as.integer(unique(dat_model$year))
  attr(fit, "model_label") <- model_label
  attr(fit, "weather_source") <- weather_source
  attr(fit, "temp_df") <- as.integer(temp_df)
  attr(fit, "season_df") <- as.integer(season_df)
  attr(fit, "linear_temperature") <- isTRUE(linear_temperature)
  attr(fit, "temp_var") <- wt$temp
  attr(fit, "analysis_data") <- dat_model

  fit
}

extract_year_hac_alert_result <- function(model) {
  if (!inherits(model, "lm")) {
    stop("extract_year_hac_alert_result() requires an lm object.", call. = FALSE)
  }

  V <- attr(model, "hac_vcov")
  if (is.null(V)) stop("HAC covariance was not stored in the model.", call. = FALSE)

  beta <- stats::coef(model)
  if (!("alert" %in% names(beta))) {
    stop("Alert coefficient not found in year-specific HAC model.", call. = FALSE)
  }

  est <- as.numeric(beta[["alert"]])
  se <- sqrt(as.numeric(V["alert", "alert"]))
  df_resid <- stats::df.residual(model)
  stat <- est / se

  # Use the residual df for a modest finite-sample t reference distribution.
  pval <- 2 * stats::pt(abs(stat), df = df_resid, lower.tail = FALSE)
  crit <- stats::qt(0.975, df = df_resid)

  tibble::tibble(
    year = as.integer(attr(model, "study_year")),
    method = "OLS with Newey-West HAC",
    hac_lag = as.integer(attr(model, "hac_lag")),
    hac_prewhite = HAC_PREWHITE,
    hac_finite_sample_adjust = HAC_FINITE_SAMPLE_ADJUST,
    estimate = est,
    std_error = se,
    lower_95 = est - crit * se,
    upper_95 = est + crit * se,
    t_value = stat,
    df_residual = as.integer(df_resid),
    p_value = pval,
    n = stats::nobs(model)
  )
}

extract_year_hac_all_coefficients <- function(model) {
  V <- attr(model, "hac_vcov")
  beta <- stats::coef(model)
  se <- sqrt(diag(V))
  df_resid <- stats::df.residual(model)
  stat <- beta / se
  crit <- stats::qt(0.975, df = df_resid)

  tibble::tibble(
    year = as.integer(attr(model, "study_year")),
    hac_lag = as.integer(attr(model, "hac_lag")),
    term = names(beta),
    estimate = as.numeric(beta),
    std_error = as.numeric(se),
    lower_95 = as.numeric(beta - crit * se),
    upper_95 = as.numeric(beta + crit * se),
    t_value = as.numeric(stat),
    df_residual = as.integer(df_resid),
    p_value = as.numeric(2 * stats::pt(abs(stat), df = df_resid, lower.tail = FALSE))
  )
}

diagnose_year_ols_residuals <- function(model, lag_max = 21L) {
  y <- as.integer(attr(model, "study_year"))
  r <- stats::residuals(model)
  n <- length(r)
  lag_max <- min(as.integer(lag_max), n - 2L)

  acf_obj <- stats::acf(r, lag.max = lag_max, plot = FALSE, na.action = na.pass)
  acf_tbl <- tibble::tibble(
    year = y,
    lag = seq_len(lag_max),
    acf = as.numeric(acf_obj$acf)[seq_len(lag_max) + 1L]
  )

  lb_lags <- c(7L, 14L, 21L)
  lb_lags <- lb_lags[lb_lags < n]
  lb_tbl <- purrr::map_dfr(lb_lags, function(L) {
    tibble::tibble(
      year = y,
      lag = L,
      p_value = stats::Box.test(r, lag = L, type = "Ljung-Box")$p.value
    )
  })

  lag1 <- if (n >= 3L) {
    suppressWarnings(stats::cor(r[-n], r[-1], use = "complete.obs"))
  } else {
    NA_real_
  }

  list(
    summary = tibble::tibble(
      year = y,
      lag1_raw_residual_correlation = lag1
    ),
    acf = acf_tbl,
    ljung_box = lb_tbl,
    residuals = tibble::tibble(
      year = y,
      day_of_season = seq_along(r),
      residual = as.numeric(r)
    )
  )
}

fit_all_years_hac <- function(dat,
                              hac_lag = HAC_PRIMARY_LAG,
                              temp_df = TEMP_DF_PRIMARY,
                              season_df = SEASON_DF_PRIMARY,
                              weather_source = "population",
                              linear_temperature = FALSE) {
  fits <- purrr::map(STUDY_YEARS, function(y) {
    fit_year_hac_model(
      dat %>% dplyr::filter(.data$year == y),
      hac_lag = hac_lag,
      temp_df = temp_df,
      season_df = season_df,
      weather_source = weather_source,
      linear_temperature = linear_temperature,
      model_label = paste0(y, " HAC(", hac_lag, ")")
    )
  })
  names(fits) <- as.character(STUDY_YEARS)

  year_results <- purrr::map_dfr(fits, extract_year_hac_alert_result)
  all_coefficients <- purrr::map_dfr(fits, extract_year_hac_all_coefficients)

  list(
    fits = fits,
    year_results = year_results,
    all_coefficients = all_coefficients
  )
}

# -----------------------------------------------------------------------------
# 13. SECOND-STAGE POOLING
# -----------------------------------------------------------------------------
pool_year_hac_estimates <- function(year_results) {
  dat <- year_results %>%
    dplyr::filter(
      is.finite(.data$estimate),
      is.finite(.data$std_error),
      .data$std_error > 0
    ) %>%
    dplyr::arrange(.data$year)

  if (nrow(dat) != length(STUDY_YEARS)) {
    stop(
      "Primary pooling requires all ", length(STUDY_YEARS),
      " valid year-specific estimates; found ", nrow(dat), ".", call. = FALSE
    )
  }

  re <- metafor::rma.uni(
    yi = dat$estimate,
    sei = dat$std_error,
    method = POOL_METHOD_PRIMARY,
    test = POOL_TEST_PRIMARY,
    slab = dat$year
  )

  ee <- metafor::rma.uni(
    yi = dat$estimate,
    sei = dat$std_error,
    method = "EE",
    test = "z",
    slab = dat$year
  )

  pred_re <- tryCatch(
    stats::predict(re),
    error = function(e) NULL
  )

  re_pi_lb <- if (!is.null(pred_re) && !is.null(pred_re$pi.lb)) {
    as.numeric(pred_re$pi.lb[1])
  } else NA_real_
  re_pi_ub <- if (!is.null(pred_re) && !is.null(pred_re$pi.ub)) {
    as.numeric(pred_re$pi.ub[1])
  } else NA_real_

  tibble::tibble(
    pooling_model = c(
      "Random-effects REML + Knapp-Hartung (PRIMARY)",
      "Common-effect inverse-variance (sensitivity)"
    ),
    estimate = c(as.numeric(re$b[1]), as.numeric(ee$b[1])),
    std_error = c(as.numeric(re$se[1]), as.numeric(ee$se[1])),
    lower_95 = c(as.numeric(re$ci.lb[1]), as.numeric(ee$ci.lb[1])),
    upper_95 = c(as.numeric(re$ci.ub[1]), as.numeric(ee$ci.ub[1])),
    p_value = c(as.numeric(re$pval[1]), as.numeric(ee$pval[1])),
    tau2 = c(as.numeric(re$tau2), 0),
    Q = c(as.numeric(re$QE), as.numeric(ee$QE)),
    Q_p_value = c(as.numeric(re$QEp), as.numeric(ee$QEp)),
    I2_percent = c(as.numeric(re$I2), as.numeric(ee$I2)),
    prediction_interval_lower = c(re_pi_lb, NA_real_),
    prediction_interval_upper = c(re_pi_ub, NA_real_),
    k_years = nrow(dat)
  )
}

fit_all_years_hac_and_pool <- function(dat,
                                       hac_lag = HAC_PRIMARY_LAG,
                                       temp_df = TEMP_DF_PRIMARY,
                                       season_df = SEASON_DF_PRIMARY,
                                       weather_source = "population",
                                       linear_temperature = FALSE) {
  yr <- fit_all_years_hac(
    dat,
    hac_lag = hac_lag,
    temp_df = temp_df,
    season_df = season_df,
    weather_source = weather_source,
    linear_temperature = linear_temperature
  )

  list(
    fits = yr$fits,
    year_results = yr$year_results,
    all_coefficients = yr$all_coefficients,
    pooled = pool_year_hac_estimates(yr$year_results)
  )
}

# -----------------------------------------------------------------------------
# 14. HAC-BASED TEMPERATURE CONTRASTS (DESCRIPTIVE/ADJUSTED)
# -----------------------------------------------------------------------------
# These curves use the same OLS coefficients and HAC covariance as the primary
# year-specific models. Each year is displayed against its own median temperature
# and only across that year's central 95% observed temperature range.

model_matrix_for_newdata <- function(model, newdata) {
  tt <- stats::delete.response(stats::terms(model))
  stats::model.matrix(
    tt,
    data = newdata,
    contrasts.arg = model$contrasts,
    xlev = model$xlevels
  )
}

make_year_temperature_curve_hac <- function(model, points = 100L) {
  dat <- attr(model, "analysis_data")
  temp_var <- attr(model, "temp_var")
  V <- attr(model, "hac_vcov")
  beta <- stats::coef(model)
  y <- as.integer(attr(model, "study_year"))

  if (isTRUE(attr(model, "linear_temperature"))) {
    stop("Temperature curve requires the nonlinear spline model.", call. = FALSE)
  }

  qs <- stats::quantile(
    dat[[temp_var]],
    probs = TEMP_CURVE_QUANTILES,
    names = FALSE,
    na.rm = TRUE
  )
  grid <- seq(qs[1], qs[2], length.out = as.integer(points))
  ref_temp <- stats::median(dat[[temp_var]], na.rm = TRUE)

  ref_dat <- dat
  ref_dat[[temp_var]] <- ref_temp
  ref_dat$alert <- 0L
  X_ref <- model_matrix_for_newdata(model, ref_dat)
  x_ref <- colMeans(X_ref[, names(beta), drop = FALSE])

  df_resid <- stats::df.residual(model)
  crit <- stats::qt(0.975, df = df_resid)

  purrr::map_dfr(grid, function(tv) {
    nd <- dat
    nd[[temp_var]] <- tv
    nd$alert <- 0L
    X <- model_matrix_for_newdata(model, nd)
    x_tv <- colMeans(X[, names(beta), drop = FALSE])
    delta <- x_tv - x_ref

    est <- as.numeric(sum(delta * beta))
    se <- sqrt(as.numeric(t(delta) %*% V %*% delta))

    tibble::tibble(
      year = y,
      temperature = tv,
      reference_temperature = ref_temp,
      adjusted_difference_sd = est,
      std_error = se,
      lower_95 = est - crit * se,
      upper_95 = est + crit * se
    )
  })
}



# 15. START ANALYSIS
# -----------------------------------------------------------------------------
log_message("Starting FINAL year-specific HAC + random-effects Step 2 analysis.")

# Resolve the final raw Google Trends source first.
RAW_SOURCE <- resolve_raw_google_trends_source()
RAW_SOURCE_TYPE <- RAW_SOURCE$source_type
RAW_SOURCE_ZIP <- RAW_SOURCE$zip_path
RAW_DIR <- RAW_SOURCE$raw_dir
PROJECT_DIR <- RAW_SOURCE$project_dir

WEATHER_FILE <- resolve_named_input(
  WEATHER_FILE_EXPLICIT, WEATHER_CANDIDATE_NAMES, PROJECT_DIR, "Weather"
)
ALERT_FILE <- resolve_named_input(
  ALERT_FILE_EXPLICIT, ALERT_CANDIDATE_NAMES, PROJECT_DIR, "Alert"
)

# Put final outputs next to raw_google_trends.zip (or the extracted directory).
ANALYSIS_OUT <- file.path(PROJECT_DIR, "outputs_gt_step2_hac_twostage_final_PRIVATE")
DIR_TABLES <- file.path(ANALYSIS_OUT, "tables")
DIR_FIGURES <- file.path(ANALYSIS_OUT, "figures")
DIR_DIAGNOSTICS <- file.path(ANALYSIS_OUT, "diagnostics")
dir.create(ANALYSIS_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_DIAGNOSTICS, showWarnings = FALSE, recursive = TRUE)

log_message("Raw Google Trends source type: ", RAW_SOURCE_TYPE)
if (!is.na(RAW_SOURCE_ZIP)) log_message("Resolved raw ZIP: ", RAW_SOURCE_ZIP)
log_message("Resolved raw session directory: ", RAW_DIR)
log_message("Resolved weather file: ", WEATHER_FILE)
log_message("Resolved alert file: ", ALERT_FILE)
log_message("Step 2 output directory: ", ANALYSIS_OUT)

input_source_manifest <- tibble::tibble(
  item = c("google_trends_source_type", "google_trends_zip", "google_trends_zip_md5",
           "raw_google_trends_directory", "weather_file", "alert_file"),
  value = c(
    RAW_SOURCE_TYPE,
    ifelse(is.na(RAW_SOURCE_ZIP), "", RAW_SOURCE_ZIP),
    ifelse(is.null(RAW_SOURCE$zip_md5) || is.na(RAW_SOURCE$zip_md5), "", RAW_SOURCE$zip_md5),
    RAW_DIR, WEATHER_FILE, ALERT_FILE
  )
)
readr::write_csv(
  input_source_manifest,
  file.path(DIR_DIAGNOSTICS, "final_input_source_manifest.csv")
)

if (identical(RAW_SOURCE_TYPE, "zip")) {
  readr::write_csv(
    RAW_SOURCE$zip_qc$session_qc,
    file.path(DIR_DIAGNOSTICS, "raw_zip_session_qc.csv")
  )
  readr::write_csv(
    RAW_SOURCE$zip_qc$detail,
    file.path(DIR_DIAGNOSTICS, "raw_zip_file_inventory.csv")
  )
}

session_manifest <- list_sessions(RAW_DIR, TARGET_SESSION_IDS)
readr::write_csv(session_manifest, file.path(DIR_DIAGNOSTICS, "session_manifest_detected.csv"))

log_message("Detected session folders: ", nrow(session_manifest))
log_message(
  "Detected unique acquisition dates: ",
  dplyr::n_distinct(session_manifest$acquisition_date)
)

# Reconstruct each session independently.
reconstructed_objects <- purrr::map(
  seq_len(nrow(session_manifest)),
  function(i) reconstruct_one_session(session_manifest[i, ])
)

session_qc <- purrr::map_dfr(reconstructed_objects, function(x) {
  tibble::tibble(
    session_id = x$session_id,
    acquisition_date = x$acquisition_date,
    reconstruction_success = x$success,
    error = x$error
  )
})
readr::write_csv(session_qc, file.path(DIR_DIAGNOSTICS, "session_reconstruction_qc.csv"))

valid_objects <- reconstructed_objects[purrr::map_lgl(reconstructed_objects, "success")]
if (length(valid_objects) == 0L) {
  stop("No complete valid sessions were available for multi-session analysis.", call. = FALSE)
}

session_daily <- purrr::map_dfr(valid_objects, "daily")
alignment_diag <- purrr::map_dfr(valid_objects, "alignment_diag")
scale_table <- purrr::map_dfr(valid_objects, "scale_table")

readr::write_csv(
  session_daily,
  file.path(DIR_TABLES, "google_trends_reconstructed_all_valid_sessions.csv")
)
readr::write_csv(
  alignment_diag,
  file.path(DIR_DIAGNOSTICS, "overlap_alignment_diagnostics.csv")
)
readr::write_csv(
  scale_table,
  file.path(DIR_DIAGNOSTICS, "overlap_scale_factors.csv")
)

valid_session_ids <- sort(unique(session_daily$session_id))
valid_acquisition_dates <- sort(unique(session_daily$acquisition_date))

log_message("Valid complete sessions: ", length(valid_session_ids))
log_message("Valid unique acquisition dates: ", length(valid_acquisition_dates))

if (length(valid_acquisition_dates) != TARGET_ACQUISITION_DATES) {
  stop(
    "FINAL analysis requires exactly ", TARGET_ACQUISITION_DATES,
    " unique acquisition dates, but found ", length(valid_acquisition_dates), ". ",
    "If more than 10 dates now exist, lock the intended 10 sessions with TARGET_SESSION_IDS before final reporting.",
    call. = FALSE
  )
}
message("FINAL target reached and locked: exactly 10 unique acquisition dates.")

# Combine same-day sessions first, then combine acquisition dates.
combined <- combine_sessions(session_daily)
by_acquisition_date <- combined$by_acquisition_date
multi_gt_daily <- combined$aggregate_daily

readr::write_csv(
  by_acquisition_date,
  file.path(DIR_TABLES, "google_trends_by_acquisition_date.csv")
)
readr::write_csv(
  multi_gt_daily,
  file.path(DIR_TABLES, "google_trends_multisession_daily_median.csv")
)

# Pairwise reproducibility.
pairwise_agreement <- calculate_pairwise_agreement(by_acquisition_date)
readr::write_csv(
  pairwise_agreement,
  file.path(DIR_DIAGNOSTICS, "pairwise_spearman_acquisition_dates.csv")
)

agreement_summary <- pairwise_agreement %>%
  dplyr::group_by(.data$study_year) %>%
  dplyr::summarise(
    n_pairs = sum(is.finite(.data$spearman_rho)),
    median_spearman = median_or_na(.data$spearman_rho),
    min_spearman = ifelse(any(is.finite(.data$spearman_rho)),
                          min(.data$spearman_rho, na.rm = TRUE), NA_real_),
    max_spearman = ifelse(any(is.finite(.data$spearman_rho)),
                          max(.data$spearman_rho, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    period = ifelse(.data$study_year == 0L, "All 2021-2025 summer days", as.character(.data$study_year))
  )
readr::write_csv(
  agreement_summary,
  file.path(DIR_DIAGNOSTICS, "pairwise_spearman_summary.csv")
)

# Read non-Google data.
weather <- read_weather_data(WEATHER_FILE)
alert <- read_alert_data(ALERT_FILE)

# Primary multi-session analysis dataset.
analysis_df <- prepare_analysis_data(multi_gt_daily, weather, alert)
readr::write_csv(
  analysis_df,
  file.path(DIR_TABLES, "analysis_dataset_multisession.csv")
)

# -----------------------------------------------------------------------------
# 16. MANUSCRIPT TABLE 1: Study characteristics by summer
# -----------------------------------------------------------------------------
# IMPORTANT:
#   This object is intentionally formatted to match the manuscript Table 1:
#   2021-2025 plus an Overall row. Google Trends acquisition-reproducibility
#   metrics are diagnostics and are NOT mixed into manuscript Table 1.

Table1_numeric_year <- analysis_df %>%
  dplyr::group_by(.data$year) %>%
  dplyr::summarise(
    days = dplyr::n(),
    n_alert_days = sum(.data$alert),
    pct_alert_days = 100 * mean(.data$alert),
    mean_pop_max_temp = mean(.data$pop_max_temp),
    sd_pop_max_temp = stats::sd(.data$pop_max_temp),
    mean_pop_humidity = mean(.data$pop_humidity),
    mean_pop_rain = mean(.data$pop_rain),
    mean_pop_wind = mean(.data$pop_wind),
    .groups = "drop"
  )

Table1_numeric_overall <- analysis_df %>%
  dplyr::summarise(
    year = NA_integer_,
    days = dplyr::n(),
    n_alert_days = sum(.data$alert),
    pct_alert_days = 100 * mean(.data$alert),
    mean_pop_max_temp = mean(.data$pop_max_temp),
    sd_pop_max_temp = stats::sd(.data$pop_max_temp),
    mean_pop_humidity = mean(.data$pop_humidity),
    mean_pop_rain = mean(.data$pop_rain),
    mean_pop_wind = mean(.data$pop_wind)
  )

Table1_numeric <- dplyr::bind_rows(Table1_numeric_year, Table1_numeric_overall)

# Structural QC: exactly five summer rows plus Overall, and 153 days per summer.
if (nrow(Table1_numeric) != length(STUDY_YEARS) + 1L) {
  stop("Manuscript Table 1 must contain 5 study-year rows plus Overall.", call. = FALSE)
}
if (!all(Table1_numeric_year$days == EXPECTED_DAYS_PER_YEAR)) {
  stop("One or more Table 1 study-year rows do not contain 153 days.", call. = FALSE)
}
if (Table1_numeric_overall$days[[1]] != EXPECTED_TOTAL_DAYS) {
  stop("Table 1 Overall row does not contain 765 days.", call. = FALSE)
}

Table1 <- Table1_numeric %>%
  dplyr::transmute(
    Year = dplyr::if_else(is.na(.data$year), "Overall", as.character(.data$year)),
    Days = .data$days,
    `Alert days, n (%)` = sprintf("%d (%.1f)", .data$n_alert_days, .data$pct_alert_days),
    `Max temp, °C, mean (SD)` = sprintf("%.1f (%.1f)", .data$mean_pop_max_temp, .data$sd_pop_max_temp),
    `Humidity, %, mean` = sprintf("%.1f", .data$mean_pop_humidity),
    `Precipitation, mm/day mean` = sprintf("%.1f", .data$mean_pop_rain),
    `Wind speed, m/s mean` = sprintf("%.1f", .data$mean_pop_wind)
  )

readr::write_csv(
  Table1,
  file.path(DIR_TABLES, "Table_1_study_characteristics.csv")
)
# Keep unrounded values for audit without creating another manuscript-numbered table.
readr::write_csv(
  Table1_numeric,
  file.path(DIR_DIAGNOSTICS, "table1_study_characteristics_numeric_audit.csv")
)

# -----------------------------------------------------------------------------
# 17. PRIMARY YEAR-SPECIFIC HAC(14) MODELS + MANUSCRIPT TABLE 2
# -----------------------------------------------------------------------------
primary_hac <- fit_all_years_hac_and_pool(
  analysis_df,
  hac_lag = HAC_PRIMARY_LAG,
  temp_df = TEMP_DF_PRIMARY,
  season_df = SEASON_DF_PRIMARY,
  weather_source = "population",
  linear_temperature = FALSE
)

primary_year_results <- primary_hac$year_results %>%
  dplyr::left_join(
    analysis_df %>%
      dplyr::group_by(.data$year) %>%
      dplyr::summarise(
        n_alert_days = sum(.data$alert),
        pct_alert_days = 100 * mean(.data$alert),
        .groups = "drop"
      ),
    by = "year"
  ) %>%
  dplyr::mutate(
    analysis_role = "Primary year-specific inference: HAC lag 14"
  )

# Preserve machine-readable year-specific results as diagnostics. The manuscript
# Table 2 is generated below after pooling so that year-specific and pooled rows
# are written as one table.
readr::write_csv(
  primary_year_results,
  file.path(DIR_DIAGNOSTICS, "primary_year_specific_HAC14_alert_results_numeric.csv")
)
readr::write_csv(
  primary_hac$all_coefficients,
  file.path(DIR_DIAGNOSTICS, "primary_HAC14_all_coefficients_numeric.csv")
)

pooled_results <- primary_hac$pooled %>%
  dplyr::mutate(
    year_specific_method = "OLS + Newey-West HAC lag 14",
    hac_prewhite = HAC_PREWHITE,
    hac_finite_sample_adjust = HAC_FINITE_SAMPLE_ADJUST,
    analysis_status = "FINAL"
  )

# Pooled numeric results remain available for figures/sensitivity analyses, but
# they are NOT written as a separate manuscript Table 3.
readr::write_csv(
  pooled_results,
  file.path(DIR_DIAGNOSTICS, "primary_two_stage_pooled_HAC14_results_numeric.csv")
)

# Formatting helpers for manuscript Table 2.
format_trim_2 <- function(x) {
  out <- sprintf("%.2f", x)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  out
}

format_p_manuscript <- function(p) {
  dplyr::if_else(
    is.finite(p) & p < 0.001,
    "<0.001",
    ifelse(is.finite(p), sprintf("%.3f", p), NA_character_)
  )
}

Table2_year <- primary_year_results %>%
  dplyr::arrange(.data$year) %>%
  dplyr::transmute(
    `Year/model` = as.character(.data$year),
    `Alert days` = as.character(.data$n_alert_days),
    `Estimate (SD)` = format_trim_2(.data$estimate),
    SE = format_trim_2(.data$std_error),
    `95% CI` = paste0(sprintf("%.2f", .data$lower_95), " to ", sprintf("%.2f", .data$upper_95)),
    `p value` = format_p_manuscript(.data$p_value)
  )

Table2_pool <- pooled_results %>%
  dplyr::mutate(
    manuscript_label = dplyr::case_when(
      grepl("PRIMARY", .data$pooling_model) ~ "Random-effects pooled",
      TRUE ~ "Common-effect sensitivity"
    ),
    manuscript_order = dplyr::if_else(grepl("PRIMARY", .data$pooling_model), 1L, 2L)
  ) %>%
  dplyr::arrange(.data$manuscript_order) %>%
  dplyr::transmute(
    `Year/model` = .data$manuscript_label,
    `Alert days` = "—",
    `Estimate (SD)` = format_trim_2(.data$estimate),
    SE = format_trim_2(.data$std_error),
    `95% CI` = paste0(sprintf("%.2f", .data$lower_95), " to ", sprintf("%.2f", .data$upper_95)),
    `p value` = format_p_manuscript(.data$p_value)
  )

Table2 <- dplyr::bind_rows(Table2_year, Table2_pool)

# Structural QC: five year-specific rows + two pooled rows, in manuscript order.
if (nrow(Table2) != length(STUDY_YEARS) + 2L) {
  stop("Manuscript Table 2 must contain 5 year-specific rows plus 2 pooled rows.", call. = FALSE)
}
expected_table2_labels <- c(
  as.character(STUDY_YEARS),
  "Random-effects pooled",
  "Common-effect sensitivity"
)
if (!identical(Table2[["Year/model"]], expected_table2_labels)) {
  stop("Manuscript Table 2 row order/labels are not as expected.", call. = FALSE)
}

readr::write_csv(
  Table2,
  file.path(DIR_TABLES, "Table_2_alert_associations.csv")
)

primary_pooled <- pooled_results %>%
  dplyr::filter(grepl("PRIMARY", .data$pooling_model))

log_message(
  "PRIMARY random-effects pooled Alert coefficient = ",
  round(primary_pooled$estimate, 4),
  " (95% CI ", round(primary_pooled$lower_95, 4),
  " to ", round(primary_pooled$upper_95, 4), ")"
)

# -----------------------------------------------------------------------------
# 18. DESCRIPTIVE RESIDUAL DIAGNOSTICS
# -----------------------------------------------------------------------------
# These diagnostics are intentionally NOT used as pass/fail criteria.
# Newey-West HAC corrects the covariance estimate in the presence of residual
# heteroskedasticity/autocorrelation; it is not a residual-whitening model.
primary_diag <- purrr::map(primary_hac$fits, diagnose_year_ols_residuals)
primary_diag_summary <- purrr::map_dfr(primary_diag, "summary")
primary_diag_acf <- purrr::map_dfr(primary_diag, "acf")
primary_diag_lb <- purrr::map_dfr(primary_diag, "ljung_box") %>%
  dplyr::mutate(
    p_holm_descriptive = stats::p.adjust(.data$p_value, method = "holm")
  )
primary_diag_resid <- purrr::map_dfr(primary_diag, "residuals")

readr::write_csv(
  primary_diag_summary,
  file.path(DIR_DIAGNOSTICS, "primary_HAC_model_raw_residual_summary.csv")
)
readr::write_csv(
  primary_diag_acf,
  file.path(DIR_DIAGNOSTICS, "primary_HAC_model_raw_residual_acf.csv")
)
readr::write_csv(
  primary_diag_lb,
  file.path(DIR_DIAGNOSTICS, "primary_HAC_model_raw_residual_ljung_box_DESCRIPTIVE.csv")
)
readr::write_csv(
  primary_diag_resid,
  file.path(DIR_DIAGNOSTICS, "primary_HAC_model_raw_residuals.csv")
)

# -----------------------------------------------------------------------------
# 19. HAC-LAG SENSITIVITY: 7, 14, 21 DAYS
# -----------------------------------------------------------------------------
hac_lags_all <- sort(unique(c(HAC_SENSITIVITY_LAGS, HAC_PRIMARY_LAG)))

hac_lag_sensitivity <- purrr::map_dfr(hac_lags_all, function(L) {
  out <- tryCatch(
    fit_all_years_hac_and_pool(
      analysis_df,
      hac_lag = L,
      temp_df = TEMP_DF_PRIMARY,
      season_df = SEASON_DF_PRIMARY,
      weather_source = "population",
      linear_temperature = FALSE
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    return(tibble::tibble(
      hac_lag = L,
      pooling_model = "Random-effects REML + Knapp-Hartung (PRIMARY framework)",
      estimate = NA_real_, std_error = NA_real_,
      lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_,
      tau2 = NA_real_, Q = NA_real_, Q_p_value = NA_real_,
      I2_percent = NA_real_,
      prediction_interval_lower = NA_real_,
      prediction_interval_upper = NA_real_,
      fit_status = "failed",
      error = conditionMessage(out)
    ))
  }

  row <- out$pooled %>%
    dplyr::filter(grepl("PRIMARY", .data$pooling_model)) %>%
    dplyr::mutate(
      hac_lag = L,
      fit_status = "success",
      error = NA_character_
    )

  row
})

readr::write_csv(
  hac_lag_sensitivity,
  file.path(DIR_TABLES, "Table_S1_HAC_lag_7_14_21_random_effects_sensitivity.csv")
)

# Also save year-specific estimates under all three HAC lags.
hac_lag_year_specific <- purrr::map_dfr(hac_lags_all, function(L) {
  out <- fit_all_years_hac(
    analysis_df,
    hac_lag = L,
    temp_df = TEMP_DF_PRIMARY,
    season_df = SEASON_DF_PRIMARY,
    weather_source = "population",
    linear_temperature = FALSE
  )
  out$year_results
})
readr::write_csv(
  hac_lag_year_specific,
  file.path(DIR_TABLES, "Table_S2_year_specific_alert_results_by_HAC_lag.csv")
)

# -----------------------------------------------------------------------------
# 20. ACQUISITION-DATE ROBUSTNESS USING THE PRIMARY HAC(14) FRAMEWORK
# -----------------------------------------------------------------------------
acquisition_date_pooled <- purrr::map_dfr(valid_acquisition_dates, function(acq) {
  gt_one <- by_acquisition_date %>%
    dplyr::filter(.data$acquisition_date == acq) %>%
    dplyr::transmute(
      date = .data$date,
      year = .data$year,
      gt_z = .data$gt_z_date
    )
  dat_one <- prepare_analysis_data(gt_one, weather, alert)

  out <- tryCatch(
    fit_all_years_hac_and_pool(
      dat_one,
      hac_lag = HAC_PRIMARY_LAG,
      temp_df = TEMP_DF_PRIMARY,
      season_df = SEASON_DF_PRIMARY,
      weather_source = "population",
      linear_temperature = FALSE
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    return(tibble::tibble(
      acquisition_date = acq,
      estimate = NA_real_, std_error = NA_real_,
      lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_,
      tau2 = NA_real_, I2_percent = NA_real_,
      fit_status = "failed", error = conditionMessage(out)
    ))
  }

  re <- out$pooled %>%
    dplyr::filter(grepl("PRIMARY", .data$pooling_model))

  tibble::tibble(
    acquisition_date = acq,
    estimate = re$estimate,
    std_error = re$std_error,
    lower_95 = re$lower_95,
    upper_95 = re$upper_95,
    p_value = re$p_value,
    tau2 = re$tau2,
    I2_percent = re$I2_percent,
    fit_status = "success",
    error = NA_character_
  )
})

readr::write_csv(
  acquisition_date_pooled,
  file.path(DIR_TABLES, "Table_S3_acquisition_date_specific_random_effects_HAC14.csv")
)

# Sequential convergence: first 1, 2, ..., 10 acquisition dates.
sequential_pooled <- purrr::map_dfr(seq_along(valid_acquisition_dates), function(k) {
  dates_use <- valid_acquisition_dates[seq_len(k)]
  gt_sub <- aggregate_selected_dates(by_acquisition_date, dates_use)
  dat_sub <- prepare_analysis_data(gt_sub, weather, alert)

  out <- tryCatch(
    fit_all_years_hac_and_pool(
      dat_sub,
      hac_lag = HAC_PRIMARY_LAG,
      temp_df = TEMP_DF_PRIMARY,
      season_df = SEASON_DF_PRIMARY,
      weather_source = "population",
      linear_temperature = FALSE
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    return(tibble::tibble(
      n_acquisition_dates = k,
      first_date = min(dates_use),
      last_date = max(dates_use),
      estimate = NA_real_, std_error = NA_real_,
      lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_,
      fit_status = "failed", error = conditionMessage(out)
    ))
  }

  re <- out$pooled %>%
    dplyr::filter(grepl("PRIMARY", .data$pooling_model))

  tibble::tibble(
    n_acquisition_dates = k,
    first_date = min(dates_use),
    last_date = max(dates_use),
    estimate = re$estimate,
    std_error = re$std_error,
    lower_95 = re$lower_95,
    upper_95 = re$upper_95,
    p_value = re$p_value,
    fit_status = "success",
    error = NA_character_
  )
})

readr::write_csv(
  sequential_pooled,
  file.path(DIR_TABLES, "Table_S4_sequential_acquisition_convergence_random_effects_HAC14.csv")
)

# Leave one acquisition date out.
loo_pooled <- purrr::map_dfr(valid_acquisition_dates, function(omit_date) {
  dates_use <- valid_acquisition_dates[valid_acquisition_dates != omit_date]
  gt_sub <- aggregate_selected_dates(by_acquisition_date, dates_use)
  dat_sub <- prepare_analysis_data(gt_sub, weather, alert)

  out <- tryCatch(
    fit_all_years_hac_and_pool(
      dat_sub,
      hac_lag = HAC_PRIMARY_LAG,
      temp_df = TEMP_DF_PRIMARY,
      season_df = SEASON_DF_PRIMARY,
      weather_source = "population",
      linear_temperature = FALSE
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    return(tibble::tibble(
      omitted_acquisition_date = omit_date,
      n_dates_retained = length(dates_use),
      estimate = NA_real_, std_error = NA_real_,
      lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_,
      fit_status = "failed", error = conditionMessage(out)
    ))
  }

  re <- out$pooled %>%
    dplyr::filter(grepl("PRIMARY", .data$pooling_model))

  tibble::tibble(
    omitted_acquisition_date = omit_date,
    n_dates_retained = length(dates_use),
    estimate = re$estimate,
    std_error = re$std_error,
    lower_95 = re$lower_95,
    upper_95 = re$upper_95,
    p_value = re$p_value,
    fit_status = "success",
    error = NA_character_
  )
})

readr::write_csv(
  loo_pooled,
  file.path(DIR_TABLES, "Table_S5_leave_one_acquisition_date_out_random_effects_HAC14.csv")
)

# -----------------------------------------------------------------------------
# 21. PRESPECIFIED MODEL-SPECIFICATION SENSITIVITY ANALYSES
# -----------------------------------------------------------------------------
sensitivity_specs <- tibble::tribble(
  ~model, ~hac_lag, ~temp_df, ~season_df, ~weather_source, ~linear_temperature,
  "Primary: HAC14, temp spline df=3, season spline df=5, population weather",
      14L, 3L, 5L, "population", FALSE,
  "HAC lag 7",  7L, 3L, 5L, "population", FALSE,
  "HAC lag 21", 21L, 3L, 5L, "population", FALSE,
  "Temperature spline df=2", 14L, 2L, 5L, "population", FALSE,
  "Temperature spline df=4", 14L, 4L, 5L, "population", FALSE,
  "Season spline df=3", 14L, 3L, 3L, "population", FALSE,
  "Season spline df=7", 14L, 3L, 7L, "population", FALSE,
  "Building-land-weighted weather", 14L, 3L, 5L, "building", FALSE,
  "Linear maximum temperature", 14L, 3L, 5L, "population", TRUE
)

sensitivity_pooled <- purrr::pmap_dfr(
  sensitivity_specs,
  function(model, hac_lag, temp_df, season_df, weather_source, linear_temperature) {
    out <- tryCatch(
      fit_all_years_hac_and_pool(
        analysis_df,
        hac_lag = hac_lag,
        temp_df = temp_df,
        season_df = season_df,
        weather_source = weather_source,
        linear_temperature = linear_temperature
      ),
      error = function(e) e
    )

    if (inherits(out, "error")) {
      return(tibble::tibble(
        model = model,
        hac_lag = hac_lag,
        estimate = NA_real_, std_error = NA_real_,
        lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_,
        tau2 = NA_real_, I2_percent = NA_real_,
        fit_status = "failed", error = conditionMessage(out)
      ))
    }

    re <- out$pooled %>%
      dplyr::filter(grepl("PRIMARY", .data$pooling_model))

    tibble::tibble(
      model = model,
      hac_lag = hac_lag,
      estimate = re$estimate,
      std_error = re$std_error,
      lower_95 = re$lower_95,
      upper_95 = re$upper_95,
      p_value = re$p_value,
      tau2 = re$tau2,
      I2_percent = re$I2_percent,
      fit_status = "success",
      error = NA_character_
    )
  }
)

readr::write_csv(
  sensitivity_pooled,
  file.path(DIR_TABLES, "Table_S6_model_specification_random_effects_sensitivity.csv")
)

# -----------------------------------------------------------------------------
# 22. YEAR-SPECIFIC TEMPERATURE CONTRAST CURVES
# -----------------------------------------------------------------------------
temperature_curves <- purrr::map_dfr(
  primary_hac$fits,
  make_year_temperature_curve_hac
)
# Numeric source data for Supplementary Fig. S5 are retained for reproducibility,
# but are not designated as a Supplementary Table.
readr::write_csv(
  temperature_curves,
  file.path(DIR_DIAGNOSTICS, "source_data_Supplementary_Figure_S5_temperature_curves_HAC14.csv")
)

# -----------------------------------------------------------------------------
# 23. FIGURES
# -----------------------------------------------------------------------------
# Figure 1: multi-session daily series with acquisition-date IQR and Alert days.
fig1_df <- analysis_df %>%
  dplyr::left_join(
    multi_gt_daily %>%
      dplyr::select(.data$date, .data$gt_z_q25, .data$gt_z_q75),
    by = "date"
  )

fig1 <- ggplot2::ggplot(
  fig1_df,
  ggplot2::aes(x = .data$date, y = .data$gt_z)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data$gt_z_q25, ymax = .data$gt_z_q75),
    fill = "grey75", alpha = 0.35
  ) +
  ggplot2::geom_line(linewidth = 0.45, colour = "black") +
  ggplot2::geom_point(
    data = fig1_df %>% dplyr::filter(.data$alert == 1L),
    size = 1.3, colour = "#D55E00"
  ) +
  ggplot2::facet_wrap(~ year, scales = "free_x", ncol = 1) +
  ggplot2::labs(
    x = NULL,
    y = "Standardized heatstroke search activity"
  ) +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Figure_1_multisession_daily_series.png"),
  fig1, width = 8, height = 10, dpi = 300
)
ggplot2::ggsave(
  file.path(DIR_FIGURES, "Figure_1_multisession_daily_series.tiff"),
  fig1, width = 8, height = 10, dpi = 600, compression = "lzw"
)

# Figure 2: year-specific HAC(14) estimates plus random/common pooled estimates.
forest_year <- primary_year_results %>%
  dplyr::transmute(
    label = as.character(.data$year),
    estimate = .data$estimate,
    lower_95 = .data$lower_95,
    upper_95 = .data$upper_95,
    type = "Year-specific HAC14"
  )

forest_pool <- pooled_results %>%
  dplyr::transmute(
    label = ifelse(
      grepl("PRIMARY", .data$pooling_model),
      "Pooled random-effects",
      "Pooled common-effect"
    ),
    estimate = .data$estimate,
    lower_95 = .data$lower_95,
    upper_95 = .data$upper_95,
    type = ifelse(grepl("PRIMARY", .data$pooling_model), "Primary pooled", "Sensitivity pooled")
  )

forest_df <- dplyr::bind_rows(forest_year, forest_pool) %>%
  dplyr::mutate(
    label = factor(
      .data$label,
      levels = rev(c(
        as.character(STUDY_YEARS),
        "Pooled random-effects",
        "Pooled common-effect"
      ))
    )
  )

fig2 <- ggplot2::ggplot(
  forest_df,
  ggplot2::aes(x = .data$estimate, y = .data$label)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = .data$lower_95, xmax = .data$upper_95),
    height = 0.18
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = .data$type),
    size = 2.4
  ) +
  ggplot2::labs(
    x = "Adjusted Alert coefficient (SD units)",
    y = NULL,
    shape = NULL
  ) +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Figure_2_year_specific_HAC14_random_effects_forest.png"),
  fig2, width = 7, height = 5.5, dpi = 300
)
ggplot2::ggsave(
  file.path(DIR_FIGURES, "Figure_2_year_specific_HAC14_random_effects_forest.tiff"),
  fig2, width = 7, height = 5.5, dpi = 600, compression = "lzw"
)

# Supplementary Figure S1: HAC lag sensitivity.
plot_lag <- hac_lag_sensitivity %>%
  dplyr::filter(.data$fit_status == "success")

fig_s1 <- ggplot2::ggplot(
  plot_lag,
  ggplot2::aes(x = factor(.data$hac_lag), y = .data$estimate)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
    width = 0.12
  ) +
  ggplot2::geom_point(size = 2.3) +
  ggplot2::labs(
    x = "Newey-West maximum lag (days)",
    y = "Random-effects pooled Alert coefficient (SD units)"
  ) +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S1_HAC_lag_sensitivity.png"),
  fig_s1, width = 6.5, height = 4.8, dpi = 300
)

# Supplementary Figure S2: acquisition-date-specific pooled coefficients.
plot_acq <- acquisition_date_pooled %>%
  dplyr::filter(.data$fit_status == "success") %>%
  dplyr::mutate(
    acquisition_date_f = factor(
      .data$acquisition_date,
      levels = sort(unique(.data$acquisition_date))
    )
  )

fig_s2 <- ggplot2::ggplot(
  plot_acq,
  ggplot2::aes(x = .data$acquisition_date_f, y = .data$estimate)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  ggplot2::geom_hline(
    yintercept = primary_pooled$estimate,
    linetype = 3,
    colour = "#0072B2"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
    width = 0.15
  ) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::labs(
    x = "Google Trends acquisition date",
    y = "Random-effects pooled Alert coefficient (SD units)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S2_acquisition_date_coefficients.png"),
  fig_s2, width = 8, height = 5, dpi = 300
)

# Supplementary Figure S3: coefficient convergence as acquisition dates accumulate.
plot_seq <- sequential_pooled %>%
  dplyr::filter(.data$fit_status == "success")

fig_s3 <- ggplot2::ggplot(
  plot_seq,
  ggplot2::aes(x = .data$n_acquisition_dates, y = .data$estimate)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
    alpha = 0.18
  ) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(
    x = "Number of acquisition dates included",
    y = "Random-effects pooled Alert coefficient (SD units)"
  ) +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S3_acquisition_convergence.png"),
  fig_s3, width = 7, height = 5, dpi = 300
)

# Supplementary Figure S4: raw OLS residual ACF.
# This is diagnostic/descriptive only; it is not a pass/fail requirement for HAC.
fig_s4 <- ggplot2::ggplot(
  primary_diag_acf,
  ggplot2::aes(x = .data$lag, y = .data$acf)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
  ggplot2::geom_segment(
    ggplot2::aes(xend = .data$lag, yend = 0),
    linewidth = 0.55
  ) +
  ggplot2::facet_wrap(~ year, ncol = 1) +
  ggplot2::scale_x_continuous(breaks = seq(1, 21, by = 2)) +
  ggplot2::labs(
    x = "Lag (days)",
    y = "ACF of OLS residuals"
  ) +
  ggplot2::theme_bw(base_size = 10)

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S4_raw_residual_ACF.png"),
  fig_s4, width = 7, height = 9, dpi = 300
)

# Supplementary Figure S5: year-specific adjusted temperature curves.
fig_s5 <- ggplot2::ggplot(
  temperature_curves,
  ggplot2::aes(
    x = .data$temperature,
    y = .data$adjusted_difference_sd
  )
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
    alpha = 0.18
  ) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::facet_wrap(~ year, scales = "free_x") +
  ggplot2::labs(
    x = "Daily maximum temperature (°C)",
    y = "Adjusted difference in search activity (SD)"
  ) +
  ggplot2::theme_bw(base_size = 10)

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S5_year_specific_temperature_curves.png"),
  fig_s5, width = 9, height = 7, dpi = 300
)

# -----------------------------------------------------------------------------
# 24. RUN SUMMARY
# -----------------------------------------------------------------------------
lag7_row <- hac_lag_sensitivity %>%
  dplyr::filter(.data$hac_lag == 7L, .data$fit_status == "success")
lag21_row <- hac_lag_sensitivity %>%
  dplyr::filter(.data$hac_lag == 21L, .data$fit_status == "success")

# Overall median pairwise Spearman correlation across all 2021-2025 summer days.
# study_year == 0 is the all-years row created in calculate_pairwise_agreement().
overall_rho <- agreement_summary %>%
  dplyr::filter(.data$study_year == 0L) %>%
  dplyr::pull(.data$median_spearman)

if (length(overall_rho) != 1L || !is.finite(overall_rho)) {
  overall_rho <- NA_real_
}

run_summary <- tibble::tibble(
  item = c(
    "analysis_status",
    "valid_complete_sessions",
    "unique_valid_acquisition_dates",
    "overall_median_pairwise_spearman",
    "year_specific_model",
    "primary_HAC_lag",
    "HAC_prewhite",
    "HAC_finite_sample_adjust",
    "primary_pooling_model",
    "primary_pooled_alert_estimate",
    "primary_pooled_alert_lower95",
    "primary_pooled_alert_upper95",
    "primary_pooled_alert_p_value",
    "primary_tau2",
    "primary_I2_percent",
    "primary_prediction_interval_lower",
    "primary_prediction_interval_upper",
    "HAC7_pooled_estimate",
    "HAC21_pooled_estimate"
  ),
  value = c(
    "FINAL",
    as.character(length(valid_session_ids)),
    as.character(length(valid_acquisition_dates)),
    as.character(overall_rho[1]),
    "OLS with natural splines + Newey-West HAC",
    as.character(HAC_PRIMARY_LAG),
    as.character(HAC_PREWHITE),
    as.character(HAC_FINITE_SAMPLE_ADJUST),
    "Random-effects REML + Knapp-Hartung",
    as.character(primary_pooled$estimate),
    as.character(primary_pooled$lower_95),
    as.character(primary_pooled$upper_95),
    as.character(primary_pooled$p_value),
    as.character(primary_pooled$tau2),
    as.character(primary_pooled$I2_percent),
    as.character(primary_pooled$prediction_interval_lower),
    as.character(primary_pooled$prediction_interval_upper),
    ifelse(nrow(lag7_row) == 1L, as.character(lag7_row$estimate), ""),
    ifelse(nrow(lag21_row) == 1L, as.character(lag21_row$estimate), "")
  )
)

readr::write_csv(
  run_summary,
  file.path(DIR_DIAGNOSTICS, "run_summary.csv")
)
capture.output(
  sessionInfo(),
  file = file.path(DIR_DIAGNOSTICS, "sessionInfo.txt")
)

log_message("Core FINAL HAC + random-effects analysis completed; starting integrated pre-submission additions.")

# =============================================================================
# 25. INTEGRATED PRE-SUBMISSION ANALYSES
# =============================================================================
# These analyses are part of the final pre-submission workflow and are executed
# automatically after the primary analysis. They are therefore stored together
# with the other manuscript tables and figures in the standard output directories.
#
# IMPORTANT:
#   - The primary analysis is NOT changed.
#   - Existing in-memory objects are reused; no second source() call is required.
#   - Supplementary Tables are numbered S1-S9 and Figures S1-S8.
# =============================================================================

FINAL_OUT <- ANALYSIS_OUT

# Reuse the definitive objects created above.
year_results <- primary_year_results
spec_sens <- sensitivity_pooled

HAC_LAG <- HAC_PRIMARY_LAG
TEMP_DF <- TEMP_DF_PRIMARY
SEASON_DF <- SEASON_DF_PRIMARY

log_message("Starting integrated pre-submission revisions.")

# -----------------------------------------------------------------------------
# Shared second-stage pooling helper
# -----------------------------------------------------------------------------
pool_random_effects_hk <- function(dat) {
  dat <- dat %>%
    dplyr::filter(is.finite(.data$estimate), is.finite(.data$std_error), .data$std_error > 0)

  if (nrow(dat) < 3L) stop("At least 3 estimates are required for pooling.", call. = FALSE)

  fit <- metafor::rma.uni(
    yi = dat$estimate,
    sei = dat$std_error,
    method = "REML",
    test = "knha"
  )

  pred <- tryCatch(stats::predict(fit), error = function(e) NULL)

  tibble::tibble(
    k = nrow(dat),
    estimate = as.numeric(fit$b[1]),
    std_error = as.numeric(fit$se[1]),
    lower_95 = as.numeric(fit$ci.lb[1]),
    upper_95 = as.numeric(fit$ci.ub[1]),
    p_value = as.numeric(fit$pval[1]),
    tau2 = as.numeric(fit$tau2),
    Q = as.numeric(fit$QE),
    Q_p_value = as.numeric(fit$QEp),
    I2_percent = as.numeric(fit$I2),
    prediction_interval_lower = if (!is.null(pred) && !is.null(pred$pi.lb)) as.numeric(pred$pi.lb[1]) else NA_real_,
    prediction_interval_upper = if (!is.null(pred) && !is.null(pred$pi.ub)) as.numeric(pred$pi.ub[1]) else NA_real_
  )
}

# =============================================================================
# REVISION 1: Temperature-adjustment flexibility
# =============================================================================
# This does NOT add new model searching. It reorganizes the already-prespecified
# sensitivity results so the systematic attenuation is transparent.

temp_flex <- spec_sens %>%
  dplyr::filter(
    .data$model %in% c(
      "Linear maximum temperature",
      "Temperature spline df=2",
      "Primary: HAC14, temp spline df=3, season spline df=5, population weather",
      "Temperature spline df=4"
    )
  ) %>%
  dplyr::mutate(
    temperature_adjustment = dplyr::case_when(
      .data$model == "Linear maximum temperature" ~ "Linear",
      .data$model == "Temperature spline df=2" ~ "Spline df=2",
      grepl("temp spline df=3", .data$model, fixed = TRUE) ~ "Spline df=3 (primary)",
      .data$model == "Temperature spline df=4" ~ "Spline df=4",
      TRUE ~ .data$model
    ),
    flexibility_order = dplyr::case_when(
      .data$temperature_adjustment == "Linear" ~ 1L,
      .data$temperature_adjustment == "Spline df=2" ~ 2L,
      .data$temperature_adjustment == "Spline df=3 (primary)" ~ 3L,
      .data$temperature_adjustment == "Spline df=4" ~ 4L,
      TRUE ~ 99L
    )
  ) %>%
  dplyr::arrange(.data$flexibility_order) %>%
  dplyr::select(
    .data$temperature_adjustment, .data$flexibility_order,
    .data$estimate, .data$std_error, .data$lower_95, .data$upper_95,
    .data$p_value, .data$tau2, .data$I2_percent
  )

# Numeric source data for Supplementary Fig. S6 are retained for reproducibility,
# but are not designated as a Supplementary Table.
readr::write_csv(
  temp_flex,
  file.path(DIR_DIAGNOSTICS, "source_data_Supplementary_Figure_S6_temperature_flexibility.csv")
)

fig_temp <- ggplot2::ggplot(
  temp_flex,
  ggplot2::aes(
    x = factor(.data$temperature_adjustment, levels = .data$temperature_adjustment),
    y = .data$estimate
  )
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
    width = 0.12
  ) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0("I²=", sprintf("%.0f", .data$I2_percent), "%")),
    vjust = -1.1,
    hjust = -0.2,
    size = 3.2
  ) +
  ggplot2::labs(
    x = "Adjustment for daily maximum temperature",
    y = "Random-effects pooled Alert coefficient (SD units)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S6_temperature_flexibility_sensitivity.png"),
  fig_temp, width = 7.2, height = 5.2, dpi = 300
)
ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S6_temperature_flexibility_sensitivity.tiff"),
  fig_temp, width = 7.2, height = 5.2, dpi = 600, compression = "lzw"
)

# =============================================================================
# REVISION 2: Leave-one-YEAR-out influence analysis
# =============================================================================
loo_year <- purrr::map_dfr(STUDY_YEARS, function(omit_year) {
  dat_sub <- year_results %>% dplyr::filter(.data$year != omit_year)
  pooled <- pool_random_effects_hk(dat_sub)
  pooled %>%
    dplyr::mutate(omitted_year = omit_year, .before = 1)
})

readr::write_csv(
  loo_year,
  file.path(DIR_TABLES, "Table_S7_leave_one_year_out_influence_HAC14.csv")
)

fig_loo <- ggplot2::ggplot(
  loo_year,
  ggplot2::aes(x = factor(.data$omitted_year), y = .data$estimate)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
    width = 0.12
  ) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(
    x = "Omitted study year",
    y = "Random-effects pooled Alert coefficient (SD units)"
  ) +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S7_leave_one_year_out_influence.png"),
  fig_loo, width = 6.6, height = 4.9, dpi = 300
)

# =============================================================================
# REVISION 3: Exploratory Alert timing (+/- 1 calendar day)
# =============================================================================
# IMPORTANT INTERPRETATION:
#   outcome on calendar day t is modeled jointly against:
#     alert_prev_day = Alert status for t-1 (possible carry-over)
#     alert          = Alert status for t   (main calendar-day indicator)
#     alert_next_day = Alert status for t+1 (a lead indicator)
#
# Because the Japanese Alert for t+1 may be announced around 17:00 on day t,
# a positive alert_next_day coefficient can be consistent with anticipatory
# information-seeking. However, daily data cannot separate that mechanism from
# persistent heat conditions, media coverage, or serial clustering of Alerts.
# This is POST HOC / EXPLORATORY and does not replace the primary same-day model.

analysis_timing <- analysis_df %>%
  dplyr::arrange(.data$year, .data$day_of_season) %>%
  dplyr::group_by(.data$year) %>%
  dplyr::mutate(
    alert_prev_day = dplyr::lag(.data$alert, 1L),
    alert_next_day = dplyr::lead(.data$alert, 1L)
  ) %>%
  dplyr::ungroup()

fit_timing_one_year <- function(dat_year) {
  dat_model <- dat_year %>%
    dplyr::filter(!is.na(.data$alert_prev_day), !is.na(.data$alert_next_day)) %>%
    dplyr::arrange(.data$day_of_season) %>%
    droplevels()

  y <- unique(dat_model$year)
  if (length(y) != 1L) stop("Timing model requires one year.", call. = FALSE)

  fit <- stats::lm(
    gt_z ~ alert_prev_day + alert + alert_next_day +
      splines::ns(pop_max_temp, df = TEMP_DF) +
      splines::ns(day_of_season, df = SEASON_DF) +
      pop_humidity_c + pop_rain_log_c + pop_wind_c + dow,
    data = dat_model,
    na.action = stats::na.fail,
    x = TRUE,
    y = TRUE,
    model = TRUE
  )

  V <- sandwich::NeweyWest(
    fit,
    lag = HAC_LAG,
    order.by = dat_model$day_of_season,
    prewhite = FALSE,
    adjust = TRUE
  )

  beta <- stats::coef(fit)
  df_resid <- stats::df.residual(fit)
  crit <- stats::qt(0.975, df = df_resid)

  purrr::map_dfr(c("alert_prev_day", "alert", "alert_next_day"), function(term) {
    est <- as.numeric(beta[[term]])
    se <- sqrt(as.numeric(V[term, term]))
    tibble::tibble(
      year = as.integer(y),
      timing_term = term,
      timing_label = dplyr::case_when(
        term == "alert_prev_day" ~ "Previous-day Alert status (t-1)",
        term == "alert" ~ "Same-day Alert status (t)",
        term == "alert_next_day" ~ "Following-day Alert status (t+1)",
        TRUE ~ term
      ),
      estimate = est,
      std_error = se,
      lower_95 = est - crit * se,
      upper_95 = est + crit * se,
      t_value = est / se,
      df_residual = as.integer(df_resid),
      p_value = 2 * stats::pt(abs(est / se), df = df_resid, lower.tail = FALSE),
      n = stats::nobs(fit)
    )
  })
}

timing_year_specific <- purrr::map_dfr(
  STUDY_YEARS,
  ~ fit_timing_one_year(analysis_timing %>% dplyr::filter(.data$year == .x))
)

readr::write_csv(
  timing_year_specific,
  file.path(DIR_TABLES, "Table_S8_alert_timing_year_specific_HAC14.csv")
)

timing_pooled <- timing_year_specific %>%
  dplyr::group_by(.data$timing_term, .data$timing_label) %>%
  dplyr::group_split(.keep = TRUE) %>%
  purrr::map_dfr(function(dat_term) {
    pooled <- pool_random_effects_hk(dat_term)
    pooled %>%
      dplyr::mutate(
        timing_term = unique(dat_term$timing_term),
        timing_label = unique(dat_term$timing_label),
        analysis_role = "Post hoc exploratory mutually adjusted +/-1-day timing model",
        .before = 1
      )
  }) %>%
  dplyr::mutate(
    timing_order = dplyr::case_when(
      .data$timing_term == "alert_prev_day" ~ 1L,
      .data$timing_term == "alert" ~ 2L,
      .data$timing_term == "alert_next_day" ~ 3L,
      TRUE ~ 99L
    )
  ) %>%
  dplyr::arrange(.data$timing_order)

readr::write_csv(
  timing_pooled,
  file.path(DIR_TABLES, "Table_S9_alert_timing_pooled_HAC14_random_effects.csv")
)

fig_timing <- ggplot2::ggplot(
  timing_pooled,
  ggplot2::aes(
    x = factor(.data$timing_label, levels = .data$timing_label),
    y = .data$estimate
  )
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
    width = 0.12
  ) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(
    x = NULL,
    y = "Random-effects pooled Alert coefficient (SD units)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 18, hjust = 1))

ggplot2::ggsave(
  file.path(DIR_FIGURES, "Supplementary_Figure_S8_alert_timing_exploratory.png"),
  fig_timing, width = 8, height = 5, dpi = 300
)

# -----------------------------------------------------------------------------
# Supplementary output numbering QC
# -----------------------------------------------------------------------------
# Formal Supplementary Tables must be exactly S1-S9. Source-data CSVs used to
# draw Figures S5/S6 are intentionally stored in diagnostics/ and are not tables.
expected_supp_tables <- sprintf("Table_S%d_", 1:9)
actual_supp_table_files <- list.files(DIR_TABLES, pattern = "^Table_S[0-9]+_.*\\.csv$", full.names = FALSE)
actual_supp_numbers <- suppressWarnings(as.integer(sub("^Table_S([0-9]+)_.*$", "\\1", actual_supp_table_files)))
if (length(actual_supp_table_files) != 9L ||
    !identical(sort(unique(actual_supp_numbers)), 1:9)) {
  stop(
    "Supplementary Table numbering QC failed. Expected exactly one file for each of S1-S9; found files: ",
    paste(actual_supp_table_files, collapse = ", "),
    call. = FALSE
  )
}
if (any(grepl("Table_S(10|11)_", actual_supp_table_files))) {
  stop("Obsolete Supplementary Table S10/S11 output detected.", call. = FALSE)
}

# =============================================================================
# 26. FINAL RUN SUMMARY / COMPLETION MARKERS
# =============================================================================
# These files are written to the full local output directory first. The public
# release package created below copies only an explicit whitelist of files.

presubmission_summary <- tibble::tibble(
  item = c(
    "primary_analysis_changed",
    "temperature_flexibility_figure_created",
    "leave_one_year_out_created",
    "exploratory_adjacent_day_timing_created",
    "repository_safe_release_created"
  ),
  value = c("FALSE", "TRUE", "TRUE", "TRUE", "TRUE")
)
readr::write_csv(
  presubmission_summary,
  file.path(DIR_DIAGNOSTICS, "presubmission_analysis_summary.csv")
)

# Refresh session information after all integrated analyses.
capture.output(
  sessionInfo(),
  file = file.path(DIR_DIAGNOSTICS, "sessionInfo.txt")
)

integrated_completion <- tibble::tibble(
  item = c(
    "integrated_analysis_status",
    "primary_analysis",
    "temperature_flexibility_analysis",
    "leave_one_year_out_analysis",
    "adjacent_day_alert_timing_analysis",
    "repository_safe_release"
  ),
  value = c(
    "FINAL_COMPLETE",
    "Year-specific OLS + Newey-West HAC lag 14; REML random-effects + Knapp-Hartung",
    "completed",
    "completed",
    "completed_post_hoc_exploratory",
    "created_from_explicit_whitelist"
  )
)
readr::write_csv(
  integrated_completion,
  file.path(DIR_DIAGNOSTICS, "integrated_analysis_completion.csv")
)

# =============================================================================
# 27. GITHUB / ZENODO REPOSITORY-SAFE PUBLIC RELEASE
# =============================================================================
# IMPORTANT:
#   - The full local analysis directory contains restricted NARO-derived daily
#     meteorological variables and MUST NOT be uploaded as a whole.
#   - The public-release directory below is built from an EXPLICIT WHITELIST.
#   - Raw Google Trends exports are NOT copied.
#   - NARO source data and daily NARO-derived meteorological variables are NOT copied.
#   - The integrated analysis dataset (analysis_dataset_multisession.csv) is NOT copied.
#   - Source data for Supplementary Fig. S5 (daily temperature values) are NOT copied.
#   - The release package contains processed standardized Google Trends series,
#     aggregate manuscript/SI outputs, selected reproducibility diagnostics,
#     figures, metadata, and R session information.
#
# The R source code itself should be committed separately to GitHub as:
#   code/heatstroke_googletrends_analysis.R
# =============================================================================

PUBLIC_RELEASE_VERSION <- "1.0.0"
PUBLIC_RELEASE_DIR <- file.path(
  PROJECT_DIR,
  paste0("github_release_v", PUBLIC_RELEASE_VERSION)
)
PUBLIC_RELEASE_ZIP <- paste0(PUBLIC_RELEASE_DIR, ".zip")

if (dir.exists(PUBLIC_RELEASE_DIR)) {
  unlink(PUBLIC_RELEASE_DIR, recursive = TRUE, force = TRUE)
}
if (file.exists(PUBLIC_RELEASE_ZIP)) {
  unlink(PUBLIC_RELEASE_ZIP, force = TRUE)
}

PUBLIC_DIR_PROCESSED <- file.path(PUBLIC_RELEASE_DIR, "processed_data")
PUBLIC_DIR_METADATA <- file.path(PUBLIC_RELEASE_DIR, "metadata")
PUBLIC_DIR_TABLES <- file.path(PUBLIC_RELEASE_DIR, "results", "tables")
PUBLIC_DIR_FIGURES <- file.path(PUBLIC_RELEASE_DIR, "results", "figures")
PUBLIC_DIR_DIAGNOSTICS <- file.path(PUBLIC_RELEASE_DIR, "diagnostics")
PUBLIC_DIR_ENVIRONMENT <- file.path(PUBLIC_RELEASE_DIR, "environment")

for (d in c(
  PUBLIC_DIR_PROCESSED,
  PUBLIC_DIR_METADATA,
  PUBLIC_DIR_TABLES,
  PUBLIC_DIR_FIGURES,
  PUBLIC_DIR_DIAGNOSTICS,
  PUBLIC_DIR_ENVIRONMENT
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

copy_public_file <- function(src, dest_dir, required = TRUE) {
  if (!file.exists(src)) {
    if (isTRUE(required)) {
      stop("Required public-release file was not found: ", src, call. = FALSE)
    }
    warning("Optional public-release file was not found: ", src, call. = FALSE)
    return(invisible(FALSE))
  }
  dest <- file.path(dest_dir, basename(src))
  ok <- file.copy(src, dest, overwrite = TRUE, copy.mode = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy public-release file: ", src, call. = FALSE)
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 27A. Processed Google Trends data (NO raw exports)
# -----------------------------------------------------------------------------
copy_public_file(
  file.path(DIR_TABLES, "google_trends_by_acquisition_date.csv"),
  PUBLIC_DIR_PROCESSED
)
copy_public_file(
  file.path(DIR_TABLES, "google_trends_multisession_daily_median.csv"),
  PUBLIC_DIR_PROCESSED
)

# -----------------------------------------------------------------------------
# 27B. Clean metadata for independent retrieval / replication
# -----------------------------------------------------------------------------
acquisition_metadata <- session_manifest %>%
  dplyr::group_by(.data$acquisition_date) %>%
  dplyr::summarise(
    n_complete_sessions = dplyr::n_distinct(.data$session_id),
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data$acquisition_date)

readr::write_csv(
  acquisition_metadata,
  file.path(PUBLIC_DIR_METADATA, "acquisition_dates.csv")
)

retrieval_specifications <- tibble::tibble(
  item = c(
    "google_trends_query",
    "google_trends_geography",
    "study_period",
    "study_years",
    "acquisition_date_count",
    "windows_per_study_year",
    "first_two_window_length_days",
    "adjacent_window_overlap_days",
    "final_window_definition",
    "within_session_alignment",
    "zero_value_handling",
    "within_session_standardization",
    "across_acquisition_combination",
    "final_outcome_standardization"
  ),
  value = c(
    "熱中症",
    "Tokyo, Japan (JP-13)",
    "May 1 through September 30",
    "2021-2025",
    as.character(dplyr::n_distinct(session_manifest$acquisition_date)),
    as.character(EXPECTED_WINDOWS_PER_YEAR),
    "75",
    "21",
    "Remaining end-of-summer period",
    "Simultaneous Huber robust regression on positive-valued overlap days",
    "Reported zeros retained in reconstructed series; positive overlap values only used to estimate scale factors",
    "Mean 0 and SD 1 within each session-year",
    "Daily median across acquisition dates after within-date combination",
    "Re-standardized to mean 0 and SD 1 within each study year"
  )
)
readr::write_csv(
  retrieval_specifications,
  file.path(PUBLIC_DIR_METADATA, "google_trends_retrieval_and_reconstruction_specifications.csv")
)

required_external_inputs <- tibble::tibble(
  input = c(
    "Original Google Trends exports",
    "NARO-derived meteorological data",
    "Heat Stroke Alert data"
  ),
  redistribution_status = c(
    "Not redistributed in this repository",
    "Not redistributed in this repository",
    "Not redistributed in this repository"
  ),
  expected_local_input = c(
    "raw_google_trends.zip or raw_google_trends/",
    "data_f.txt (or compatible local filename)",
    "check_heatstroke_alert_2021_2025_wide_by_region.xlsx (or compatible local filename)"
  ),
  note = c(
    "Retrieve independently using the documented query, geography, dates, and window settings.",
    "Obtain access directly from NARO; the analysis requires population- and building-land-weighted daily variables.",
    "Obtain from the official Japanese Heat Illness Prevention Information source."
  )
)
readr::write_csv(
  required_external_inputs,
  file.path(PUBLIC_DIR_METADATA, "required_external_inputs.csv")
)

public_source_manifest <- tibble::tibble(
  item = c(
    "release_version",
    "google_trends_source_type",
    "google_trends_raw_archive_md5",
    "naro_daily_data_included",
    "raw_google_trends_included",
    "integrated_weather_analysis_dataset_included"
  ),
  value = c(
    PUBLIC_RELEASE_VERSION,
    RAW_SOURCE_TYPE,
    ifelse(is.null(RAW_SOURCE$zip_md5) || is.na(RAW_SOURCE$zip_md5), "", RAW_SOURCE$zip_md5),
    "FALSE",
    "FALSE",
    "FALSE"
  )
)
readr::write_csv(
  public_source_manifest,
  file.path(PUBLIC_DIR_METADATA, "public_source_manifest.csv")
)

# -----------------------------------------------------------------------------
# 27C. Aggregate manuscript and Supplementary Tables
# -----------------------------------------------------------------------------
public_table_files <- c(
  "Table_1_study_characteristics.csv",
  "Table_2_alert_associations.csv",
  "Table_S1_HAC_lag_7_14_21_random_effects_sensitivity.csv",
  "Table_S2_year_specific_alert_results_by_HAC_lag.csv",
  "Table_S3_acquisition_date_specific_random_effects_HAC14.csv",
  "Table_S4_sequential_acquisition_convergence_random_effects_HAC14.csv",
  "Table_S5_leave_one_acquisition_date_out_random_effects_HAC14.csv",
  "Table_S6_model_specification_random_effects_sensitivity.csv",
  "Table_S7_leave_one_year_out_influence_HAC14.csv",
  "Table_S8_alert_timing_year_specific_HAC14.csv",
  "Table_S9_alert_timing_pooled_HAC14_random_effects.csv"
)

for (f in public_table_files) {
  copy_public_file(file.path(DIR_TABLES, f), PUBLIC_DIR_TABLES)
}

# -----------------------------------------------------------------------------
# 27D. Formal manuscript / Supplementary figures
# -----------------------------------------------------------------------------
public_figure_files <- list.files(
  DIR_FIGURES,
  pattern = "^(Figure_|Supplementary_Figure_).*[.](png|tiff)$",
  full.names = TRUE,
  ignore.case = TRUE
)
if (length(public_figure_files) == 0L) {
  stop("No formal figure files were found for the public release.", call. = FALSE)
}
for (f in public_figure_files) {
  copy_public_file(f, PUBLIC_DIR_FIGURES)
}

# -----------------------------------------------------------------------------
# 27E. Selected reproducibility / statistical diagnostics
# -----------------------------------------------------------------------------
# Deliberately excluded:
#   - final_input_source_manifest.csv (contains local paths in the private output)
#   - raw_zip_file_inventory.csv
#   - google_trends_reconstructed_all_valid_sessions.csv
#   - analysis_dataset_multisession.csv
#   - source_data_Supplementary_Figure_S5_temperature_curves_HAC14.csv
#   - primary_HAC_model_raw_residuals.csv
# The excluded files are not needed for public verification of the reported results.

public_diagnostic_files <- c(
  "session_reconstruction_qc.csv",
  "overlap_alignment_diagnostics.csv",
  "overlap_scale_factors.csv",
  "pairwise_spearman_acquisition_dates.csv",
  "pairwise_spearman_summary.csv",
  "primary_year_specific_HAC14_alert_results_numeric.csv",
  "primary_two_stage_pooled_HAC14_results_numeric.csv",
  "primary_HAC_model_raw_residual_summary.csv",
  "primary_HAC_model_raw_residual_acf.csv",
  "primary_HAC_model_raw_residual_ljung_box_DESCRIPTIVE.csv",
  "source_data_Supplementary_Figure_S6_temperature_flexibility.csv",
  "run_summary.csv",
  "presubmission_analysis_summary.csv",
  "integrated_analysis_completion.csv"
)

for (f in public_diagnostic_files) {
  copy_public_file(file.path(DIR_DIAGNOSTICS, f), PUBLIC_DIR_DIAGNOSTICS)
}

copy_public_file(
  file.path(DIR_DIAGNOSTICS, "sessionInfo.txt"),
  PUBLIC_DIR_ENVIRONMENT
)

# -----------------------------------------------------------------------------
# 27F. Safety checks: filenames, restricted columns, and local paths
# -----------------------------------------------------------------------------
all_public_files <- list.files(
  PUBLIC_RELEASE_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)

# Deny files that should never enter the public package.
forbidden_name_patterns <- c(
  "analysis_dataset_multisession",
  "data_f",
  "raw_google_trends",
  "raw_zip_file_inventory",
  "google_trends_reconstructed_all_valid_sessions",
  "source_data_Supplementary_Figure_S5_temperature_curves_HAC14",
  "primary_HAC_model_raw_residuals[.]csv"
)

forbidden_name_hits <- purrr::map_dfr(forbidden_name_patterns, function(pat) {
  hits <- all_public_files[grepl(pat, basename(all_public_files), ignore.case = TRUE)]
  if (length(hits) == 0L) return(tibble::tibble())
  tibble::tibble(
    check_type = "forbidden_filename",
    pattern = pat,
    file = substring(hits, nchar(PUBLIC_RELEASE_DIR) + 2L)
  )
})

# Text-level checks. These exact internal column names identify daily NARO-derived
# meteorological variables used by the local analysis dataset and should not appear
# in public CSV/TXT/MD/JSON content.
text_files <- all_public_files[grepl("[.](csv|txt|md|json|R)$", all_public_files, ignore.case = TRUE)]
restricted_content_patterns <- c(
  "pop_max_temp",
  "pop_humidity",
  "pop_rain",
  "pop_wind",
  "built_max_temp",
  "built_humidity",
  "built_rain",
  "built_wind"
)

restricted_content_hits <- purrr::map_dfr(text_files, function(f) {
  lines <- tryCatch(
    readLines(f, warn = FALSE, encoding = "UTF-8"),
    error = function(e) character()
  )
  purrr::map_dfr(restricted_content_patterns, function(pat) {
    hit_lines <- grep(pat, lines, fixed = TRUE)
    if (length(hit_lines) == 0L) return(tibble::tibble())
    tibble::tibble(
      check_type = "restricted_weather_column",
      pattern = pat,
      file = substring(f, nchar(PUBLIC_RELEASE_DIR) + 2L),
      line_number = hit_lines,
      matched_line = lines[hit_lines]
    )
  })
})

# Scan for common absolute local paths / user-profile strings.
current_user <- as.character(Sys.info()[["user"]])
privacy_patterns <- unique(c(
  "/Users/",
  "\\\\Users\\\\",
  "CloudStorage",
  "Dropbox",
  if (nzchar(current_user)) current_user else character()
))

privacy_hits <- purrr::map_dfr(text_files, function(f) {
  lines <- tryCatch(
    readLines(f, warn = FALSE, encoding = "UTF-8"),
    error = function(e) character()
  )
  purrr::map_dfr(privacy_patterns, function(pat) {
    hit_lines <- grep(pat, lines, fixed = TRUE)
    if (length(hit_lines) == 0L) return(tibble::tibble())
    tibble::tibble(
      check_type = "local_path_or_profile",
      pattern = pat,
      file = substring(f, nchar(PUBLIC_RELEASE_DIR) + 2L),
      line_number = hit_lines,
      matched_line = lines[hit_lines]
    )
  })
})

release_safety_scan <- dplyr::bind_rows(
  forbidden_name_hits,
  restricted_content_hits,
  privacy_hits
)

readr::write_csv(
  release_safety_scan,
  file.path(PUBLIC_RELEASE_DIR, "release_safety_scan.csv")
)

if (nrow(release_safety_scan) > 0L) {
  stop(
    "Repository-safety scan detected one or more files that should not be released. ",
    "Review release_safety_scan.csv. The private/full analysis outputs remain unchanged.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 27G. Public file inventory with checksums
# -----------------------------------------------------------------------------
all_public_files <- list.files(
  PUBLIC_RELEASE_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)
all_public_files <- all_public_files[file.info(all_public_files)$isdir %in% FALSE]
all_public_files <- all_public_files[basename(all_public_files) != "public_file_inventory.csv"]

public_file_inventory <- tibble::tibble(
  relative_path = substring(all_public_files, nchar(PUBLIC_RELEASE_DIR) + 2L),
  size_bytes = as.numeric(file.info(all_public_files)$size),
  md5 = unname(tools::md5sum(all_public_files))
) %>%
  dplyr::arrange(.data$relative_path)

readr::write_csv(
  public_file_inventory,
  file.path(PUBLIC_RELEASE_DIR, "public_file_inventory.csv")
)

public_release_manifest <- tibble::tibble(
  item = c(
    "release_version",
    "release_status",
    "release_strategy",
    "raw_google_trends_included",
    "naro_daily_meteorological_data_included",
    "integrated_analysis_dataset_included",
    "processed_google_trends_included",
    "aggregate_statistical_outputs_included",
    "formal_figures_included",
    "selected_reproducibility_diagnostics_included",
    "local_path_scan_passed",
    "restricted_weather_column_scan_passed",
    "note"
  ),
  value = c(
    PUBLIC_RELEASE_VERSION,
    "repository_safety_scan_passed",
    "explicit_whitelist",
    "FALSE",
    "FALSE",
    "FALSE",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    paste0(
      "Upload only this release package plus code/heatstroke_googletrends_analysis.R to GitHub. ",
      "Do not upload the full local *_PRIVATE analysis directory."
    )
  )
)
readr::write_csv(
  public_release_manifest,
  file.path(PUBLIC_RELEASE_DIR, "public_release_manifest.csv")
)

# Rebuild inventory once more so it includes the manifest itself.
all_public_files <- list.files(
  PUBLIC_RELEASE_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)
all_public_files <- all_public_files[file.info(all_public_files)$isdir %in% FALSE]
all_public_files <- all_public_files[basename(all_public_files) != "public_file_inventory.csv"]
public_file_inventory <- tibble::tibble(
  relative_path = substring(all_public_files, nchar(PUBLIC_RELEASE_DIR) + 2L),
  size_bytes = as.numeric(file.info(all_public_files)$size),
  md5 = unname(tools::md5sum(all_public_files))
) %>%
  dplyr::arrange(.data$relative_path)
readr::write_csv(
  public_file_inventory,
  file.path(PUBLIC_RELEASE_DIR, "public_file_inventory.csv")
)

# Make a distributable ZIP using relative paths.
old_wd <- getwd()
tryCatch(
  {
    setwd(dirname(PUBLIC_RELEASE_DIR))
    utils::zip(
      zipfile = basename(PUBLIC_RELEASE_ZIP),
      files = basename(PUBLIC_RELEASE_DIR),
      flags = "-r9X"
    )
  },
  finally = {
    setwd(old_wd)
  }
)

log_message("FINAL INTEGRATED analysis completed successfully.")
log_message("PRIVATE full analysis outputs (do not upload): ", ANALYSIS_OUT)
log_message("GitHub/Zenodo-safe release directory: ", PUBLIC_RELEASE_DIR)
log_message("GitHub/Zenodo-safe release ZIP: ", PUBLIC_RELEASE_ZIP)
log_message(
  "Primary inference remains HAC lag ", HAC_PRIMARY_LAG,
  " with REML random-effects + Knapp-Hartung pooling."
)
