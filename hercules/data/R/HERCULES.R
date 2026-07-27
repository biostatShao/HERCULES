#' Run HERCULES from a validated YAML configuration
#'
#' Python owns configuration, paths, stage orchestration, manifests, and
#' subprocess handling. R remains the boundary for the validated SuperLearner
#' procedure used by the installed workflow.
#'
#' @param config Path to hercules.yaml.
#' @param hercules_executable Executable name or path.
#' @param dry_run Validate and print the planned stages without inference.
#' @param extra_args Additional command-line arguments.
#' @return Invisibly returns the captured command output.
#' @export
HERCULES <- function(config,
                     hercules_executable = "hercules",
                     dry_run = FALSE,
                     extra_args = character()) {
  if (length(config) != 1L || !nzchar(config)) {
    stop("config must be one non-empty path", call. = FALSE)
  }
  if (!file.exists(config)) {
    stop(sprintf("HERCULES configuration does not exist: %s", config), call. = FALSE)
  }

  config_path <- normalizePath(config, winslash = "/", mustWork = TRUE)
  args <- c("run", "--config", shQuote(config_path))
  if (isTRUE(dry_run)) {
    args <- c(args, "--dry-run")
  }
  args <- c(args, vapply(extra_args, shQuote, character(1)))

  command <- hercules_executable
  if (file.exists(command)) {
    command <- shQuote(normalizePath(command, winslash = "/", mustWork = TRUE))
  }

  output <- system2(
    command = command,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      paste(c(sprintf("HERCULES exited with status %s", status), output), collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(output)
}
