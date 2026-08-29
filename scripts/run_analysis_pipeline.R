library(rstudioapi)
library(here)
library(purrr)
library(later)

# Settings ----------------------------------------------------------------

# Use 1 for the lowest memory usage. Increase to 2 only if enough RAM is
# available to hold two independent R sessions and their data at the same time.
max_parallel_jobs <- 1L
job_check_interval <- 5

# Analysis plan ------------------------------------------------------------

# Stages are run in the order shown below. Scripts within the same stage may
# run in parallel, subject to `max_parallel_jobs`.
pipeline_plan <- list(
  preprocessing = c(
    raids = "scripts/raids/1_Bio-RAIDs_data_preprocessing.R",
    tcga  = "scripts/tcga/1_TCGA-CESC_data_preprocessing.R"
  ),

  cox_lasso = c(
    raids = "scripts/raids/2_Bio-RAIDs_cox_lasso.R",
    tcga  = "scripts/tcga/2_TCGA-CESC_cox_lasso.R"
  ),

  cox_boosting = c(
    raids = "scripts/raids/3_Bio-RAIDs_cox_boosting.R",
    tcga  = "scripts/tcga/3_TCGA-CESC_cox_boosting.R"
  ),

  vimp_survML = c(
    raids    = "scripts/raids/4_Bio-RAIDs_vimp_survML.R",
    tcga     = "scripts/tcga/4_TCGA-CESC_vimp_survML.R",
    combined = "scripts/combined/2_combined_vimp_survML.R"
  ),

  vimp_ranger = c(
    raids = "scripts/raids/5_Bio-RAIDs_vimp_ranger.R",
    tcga  = "scripts/tcga/5_TCGA-CESC_vimp_ranger.R"
  )
)

# Build memory-safe batches ------------------------------------------------

create_execution_plan <- function(plan, max_parallel = 1L) {
  if (
    length(max_parallel) != 1L ||
      is.na(max_parallel) ||
      max_parallel < 1L ||
      max_parallel != as.integer(max_parallel)
  ) {
    stop(
      "`max_parallel` must be a single positive integer.",
      call. = FALSE
    )
  }

  execution_plan <- list()

  for (stage in names(plan)) {
    scripts <- plan[[stage]]
    batch_number <- ceiling(seq_along(scripts) / max_parallel)
    batches <- split(scripts, batch_number)

    for (i in seq_along(batches)) {
      batch_name <- paste0(stage, "_batch_", i)
      execution_plan[[batch_name]] <- batches[[i]]
    }
  }

  execution_plan
}

execution_plan <- create_execution_plan(
  plan = pipeline_plan,
  max_parallel = max_parallel_jobs
)

# Preliminary checks -------------------------------------------------------

rstudioapi::verifyAvailable()

all_scripts <- unlist(execution_plan, use.names = FALSE)

missing_scripts <- purrr::keep(
  all_scripts,
  \(script) !file.exists(here::here(script))
)

if (length(missing_scripts) > 0L) {
  stop(
    paste0(
      "The following scripts could not be found:\n- ",
      paste(missing_scripts, collapse = "\n- ")
    ),
    call. = FALSE
  )
}

# Launch one batch ---------------------------------------------------------

launch_jobs <- function(scripts, batch_name) {
  purrr::imap(
    scripts,
    \(script, dataset) {
      job_name <- paste(batch_name, dataset, sep = " | ")

      message("Starting job: ", job_name)

      rstudioapi::jobRunScript(
        path = here::here(script),
        name = job_name,
        workingDir = here::here(),
        importEnv = FALSE,
        exportEnv = ""
      )
    }
  )
}

# Monitor one batch without blocking the RStudio console ------------------

wait_for_jobs <- function(
    job_ids,
    batch_name,
    on_success,
    interval = 5
) {
  check_jobs <- function() {
    states <- purrr::map_chr(
      job_ids,
      \(job_id) {
        tryCatch(
          rstudioapi::jobGetState(job_id),
          error = \(e) "unknown"
        )
      }
    )

    names(states) <- names(job_ids)
    pipeline_status$last_job_states <- states

    unsuccessful <- states %in% c("failed", "cancelled", "unknown")

    if (any(unsuccessful)) {
      failed_jobs <- paste0(
        names(states)[unsuccessful],
        " [",
        states[unsuccessful],
        "]"
      )

      pipeline_status$status <- "failed"
      pipeline_status$failed_batch <- batch_name
      pipeline_status$current_batch <- NULL

      warning(
        paste0(
          "The pipeline stopped during batch `",
          batch_name,
          "`.\nAffected jobs:\n- ",
          paste(failed_jobs, collapse = "\n- ")
        ),
        call. = FALSE
      )

      return(invisible(NULL))
    }

    if (all(states == "succeeded")) {
      message("Completed batch: ", batch_name)
      on_success()
      return(invisible(NULL))
    }

    later::later(check_jobs, delay = interval)
  }

  later::later(check_jobs, delay = interval)
  invisible(job_ids)
}

# Run all batches in sequence ---------------------------------------------

# Job identifiers and pipeline state remain available for inspection after
# this script has been sourced.
pipeline_jobs <- new.env(parent = emptyenv())

pipeline_status <- new.env(parent = emptyenv())
pipeline_status$status <- "running"
pipeline_status$current_batch <- NULL
pipeline_status$failed_batch <- NULL
pipeline_status$last_job_states <- NULL

run_next_batch <- function(batch_number = 1L) {
  if (batch_number > length(execution_plan)) {
    pipeline_status$status <- "succeeded"
    pipeline_status$current_batch <- NULL
    message("All analysis batches completed successfully.")
    return(invisible(NULL))
  }

  batch_name <- names(execution_plan)[batch_number]
  scripts <- execution_plan[[batch_number]]

  pipeline_status$current_batch <- batch_name

  message(
    "\nStarting batch ",
    batch_number,
    "/",
    length(execution_plan),
    ": ",
    batch_name
  )

  job_ids <- launch_jobs(
    scripts = scripts,
    batch_name = batch_name
  )

  pipeline_jobs[[batch_name]] <- job_ids

  wait_for_jobs(
    job_ids = job_ids,
    batch_name = batch_name,
    interval = job_check_interval,
    on_success = function() {
      run_next_batch(batch_number + 1L)
    }
  )

  invisible(job_ids)
}

# Start the pipeline -------------------------------------------------------

message(
  "Starting the analysis pipeline with a maximum of ",
  max_parallel_jobs,
  " simultaneous job(s)."
)

run_next_batch()

# Monitoring commands -----------------------------------------------------

# pipeline_status$status
# pipeline_status$current_batch
# pipeline_status$failed_batch
# pipeline_status$last_job_states
# as.list(pipeline_jobs)
