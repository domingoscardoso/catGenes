#' Run MrBayes from RStudio using an existing NEXUS file
#'
#' @author Domingos Cardoso
#'
#' @description
#' Runs a Bayesian phylogenetic analysis in
#' [MrBayes](https://nbisweden.github.io/MrBayes/) directly from R/RStudio using a
#' pre-configured NEXUS file that already contains a complete \code{begin mrbayes;}
#' command block. The function copies (or moves) the NEXUS file to a dedicated run
#' directory, launches MrBayes, and captures the program output.
#'
#' When \pkg{processx} is installed and \code{use_processx = TRUE}, output can be
#' streamed live to the R console and the returned object includes a \code{stop()}
#' handle that can terminate the run programmatically. Programmatic stopping is
#' only useful when \code{background = TRUE}, because in foreground mode the call
#' blocks until MrBayes finishes.
#'
#' Without \pkg{processx}, the analysis can still be run and live output can be
#' printed, but stopping is interactive (e.g., Esc in RStudio).
#'
#' @usage
#' mrbayesRun(nexus_file,
#'            mrbayes_dir,
#'            run_dir = NULL,
#'            copy_mode = "copy",
#'            quiet = FALSE,
#'            live = TRUE,
#'            use_processx = TRUE,
#'            background = FALSE)
#'
#' @param nexus_file
#' Path to a NEXUS file (\code{.nex} or \code{.nexus}) that already contains a valid
#' MrBayes block (i.e., \code{begin mrbayes;} ... \code{end;}).
#'
#' @param mrbayes_dir
#' Path to the folder containing the MrBayes executable. On macOS/Linux this is
#' typically a folder containing \code{mb}; on Windows, \code{mb.exe}. If the
#' executable is not found directly in \code{mrbayes_dir}, the function also checks
#' \code{file.path(mrbayes_dir, "bin")}.
#'
#' Typical values for \code{mrbayes_dir}:
#' \itemize{
#'   \item macOS (Apple Silicon, Homebrew): \code{"/opt/homebrew/bin"}
#'   \item macOS (Intel, Homebrew): \code{"/usr/local/bin"}
#'   \item macOS/Linux (local install): \code{"~/.local/bin"} (use \code{path.expand("~/.local/bin")})
#'   \item Windows (typical): \code{"C:/Program Files/MrBayes"} or a folder containing \code{mb.exe}
#' }
#'
#' @param run_dir
#' Directory where the analysis will run. If \code{NULL}, a new folder named
#' \code{mrbayes_run_YYYYmmdd_HHMMSS} is created in \code{getwd()}.
#'
#' @param copy_mode
#' Either \code{"copy"} (default) or \code{"move"}. Determines whether \code{nexus_file}
#' is copied to \code{run_dir} or moved there before running MrBayes.
#'
#' @param quiet
#' Logical. If \code{TRUE}, suppresses most informational messages (default: \code{FALSE}).
#'
#' @param live
#' Logical. If \code{TRUE}, prints MrBayes output to the R console (default: \code{TRUE}).
#' Live streaming is most reliable with \pkg{processx}.
#'
#' @param use_processx
#' Logical. If \code{TRUE} (default), uses \pkg{processx} when installed to provide
#' robust output streaming and a programmatic stop handle.
#'
#' @param background
#' Logical. If \code{TRUE}, starts MrBayes and returns immediately (non-blocking).
#' The returned object includes \code{poll()} to print new output on demand and
#' \code{stop()} to terminate the run. Requires \pkg{processx}. Default: \code{FALSE}.
#'
#' @return
#' A list with (at minimum) \code{exe_path}, \code{run_dir}, \code{nexus_path},
#' \code{stdout}, \code{stderr}, \code{exit_status}, and \code{error}. When
#' \pkg{processx} is used, the list additionally includes \code{proc} (a
#' \pkg{processx} process object) and \code{stop}. When \code{background = TRUE},
#' the list also includes \code{poll}.
#'
#' @details
#' Output capture:
#' \itemize{
#'   \item \code{stdout}: file path where standard output is recorded.
#'   \item \code{stderr}: file path where errors/warnings are recorded.
#' }
#' The function runs MrBayes in \code{run_dir}, so all MrBayes output files
#' (e.g., \code{.p}, \code{.t}, \code{.con.tre}) will be written there.
#'
#' Stopping the run:
#' \itemize{
#'   \item With \pkg{processx} + \code{background = TRUE}: call \code{mrbayesStop(res)} or \code{res$stop()}.
#'   \item With \pkg{processx} + \code{background = FALSE}: stop interactively (Esc) because the call blocks.
#'   \item Without \pkg{processx}: stop interactively (Esc in RStudio).
#' }
#'
#' @examples
#' \dontrun{
#' # Background mode (recommended for stop-anytime control)
#' res <- mrbayesRun("analysis.nex",
#'                   mrbayes_dir = "/opt/homebrew/bin",
#'                   background = TRUE)
#' res$poll()          # print new output
#' mrbayesStop(res)    # graceful stop
#' mrbayesStop(res, TRUE) # force kill
#'
#' # Foreground mode (blocking): stop with Esc in RStudio
#' mrbayesRun("analysis.nex",
#'            mrbayes_dir = "/opt/homebrew/bin")
#' }
#'
#' @seealso \code{\link{evomodelTest}}
#'
#' @importFrom processx process
#'
#' @export
#'
mrbayesRun <- function(nexus_file,
                       mrbayes_dir,
                       run_dir = NULL,
                       copy_mode = "copy",
                       quiet = FALSE,
                       live = TRUE,
                       use_processx = TRUE,
                       background = FALSE) {

  if (!copy_mode %in% c("copy", "move")) stop("copy_mode must be 'copy' or 'move'.")
  if (!is.logical(background) || length(background) != 1) stop("background must be TRUE/FALSE.")
  if (!is.logical(live) || length(live) != 1) stop("live must be TRUE/FALSE.")
  if (!is.logical(use_processx) || length(use_processx) != 1) stop("use_processx must be TRUE/FALSE.")

  if (is.null(nexus_file) || !nzchar(nexus_file)) stop("nexus_file is required.")
  if (!file.exists(nexus_file)) stop("nexus_file does not exist: ", nexus_file)
  nexus_file <- normalizePath(nexus_file, winslash = "/", mustWork = TRUE)

  if (is.null(mrbayes_dir) || !nzchar(mrbayes_dir)) stop("mrbayes_dir is required.")
  if (!dir.exists(mrbayes_dir)) stop("mrbayes_dir does not exist: ", mrbayes_dir)
  mrbayes_dir <- normalizePath(mrbayes_dir, winslash = "/", mustWork = TRUE)

  exe_name <- if (identical(.Platform$OS.type, "windows")) "mb.exe" else "mb"

  exe_path <- file.path(mrbayes_dir, exe_name)
  if (!file.exists(exe_path)) {
    exe_alt <- file.path(mrbayes_dir, "bin", exe_name)
    if (!file.exists(exe_alt)) {
      stop("MrBayes executable not found.\nTried: ", exe_path, "\nAlso tried: ", exe_alt)
    }
    exe_path <- exe_alt
  }
  exe_path <- normalizePath(exe_path, winslash = "/", mustWork = TRUE)

  if (is.null(run_dir)) {
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    run_dir <- file.path(getwd(), paste0("mrbayes_run_", stamp))
  }
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

  dest_nexus <- file.path(run_dir, basename(nexus_file))
  if (copy_mode == "copy") {
    ok <- file.copy(nexus_file, dest_nexus, overwrite = TRUE)
    if (!ok) stop("Failed to copy nexus_file into run_dir.")
  } else {
    ok <- file.rename(nexus_file, dest_nexus)
    if (!ok) stop("Failed to move nexus_file into run_dir.")
  }

  stdout_file <- file.path(run_dir, "mrbayes_stdout.txt")
  stderr_file <- file.path(run_dir, "mrbayes_stderr.txt")

  if (!quiet) {
    message("MrBayes executable: ", exe_path)
    message("Run directory: ", run_dir)
    message("Running NEXUS: ", dest_nexus)
  }

  want_processx <- isTRUE(use_processx) && requireNamespace("processx", quietly = TRUE)
  if (isTRUE(background) && !want_processx) stop("background=TRUE requires the processx package.")

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(run_dir)

  # ---- processx path (supports background + programmatic stop) ----
  if (want_processx) {
    out_con <- file(stdout_file, open = "ab")
    err_con <- file(stderr_file, open = "ab")
    on.exit({
      try(close(out_con), silent = TRUE)
      try(close(err_con), silent = TRUE)
    }, add = TRUE)

    buf_write <- function(con, txt) {
      if (length(txt) == 0) return(invisible(NULL))
      if (is.character(txt) && length(txt) > 1) txt <- paste(txt, collapse = "\n")
      if (!is.character(txt)) txt <- as.character(txt)
      if (!nzchar(txt)) return(invisible(NULL))
      writeChar(paste0(txt, "\n"), con, eos = NULL, useBytes = TRUE)
      flush(con)
      invisible(NULL)
    }

    proc <- processx::process$new(
      command = exe_path,
      args = c(basename(dest_nexus)),
      wd = run_dir,
      stdout = "|",
      stderr = "|",
      cleanup = TRUE
    )

    stop_fun <- function(grace = TRUE) {
      if (!proc$is_alive()) return(invisible(FALSE))
      if (grace) {
        if (identical(.Platform$OS.type, "windows")) proc$signal(15) else proc$signal(2)
      } else {
        proc$kill()
      }
      invisible(TRUE)
    }

    poll_fun <- function() {
      out <- proc$read_output_lines()
      if (length(out)) {
        if (isTRUE(live)) cat(paste0(out, "\n"))
        buf_write(out_con, out)
      }

      errl <- proc$read_error_lines()
      if (length(errl)) {
        if (isTRUE(live)) cat(paste0(errl, "\n"), file = stderr())
        buf_write(err_con, errl)
      }
      invisible(NULL)
    }

    if (!quiet) {
      msg <- if (background) {
        "Process started (background). Use res$poll() to follow; res$stop() to stop."
      } else {
        "Process started (foreground). Stop with Esc, or use background=TRUE for stop-anytime control."
      }
      message(msg, " (PID ", proc$get_pid(), ")")
    }

    if (isTRUE(background)) {
      return(list(
        exe_path = exe_path,
        run_dir = run_dir,
        nexus_path = dest_nexus,
        stdout = stdout_file,
        stderr = stderr_file,
        exit_status = NA_integer_,
        error = NULL,
        proc = proc,
        stop = stop_fun,
        poll = poll_fun
      ))
    }

    # Foreground (blocking) streaming
    exit_status <- NA_integer_
    err <- NULL
    tryCatch({
      while (proc$is_alive()) {
        poll_fun()
        Sys.sleep(0.05)
      }
      poll_fun()
      exit_status <- proc$get_exit_status()
    }, interrupt = function(e) {
      err <<- e
      try(proc$kill(), silent = TRUE)
      exit_status <<- NA_integer_
      if (!quiet) message("\nInterrupted: process killed.")
    }, error = function(e) {
      err <<- e
      try(proc$kill(), silent = TRUE)
      exit_status <<- NA_integer_
    })

    return(list(
      exe_path = exe_path,
      run_dir = run_dir,
      nexus_path = dest_nexus,
      stdout = stdout_file,
      stderr = stderr_file,
      exit_status = exit_status,
      error = err,
      proc = proc,
      stop = stop_fun
    ))
  }

  # ---- base R fallback (no process handle) ----
  if (!quiet && isTRUE(live)) {
    message("Tip: install.packages('processx') for background mode and programmatic stop().")
  }

  res <- tryCatch(
    system2(
      exe_path,
      args = c(basename(dest_nexus)),
      stdout = if (live) "" else stdout_file,
      stderr = if (live) "" else stderr_file
    ),
    error = function(e) structure(NA_integer_, error = e)
  )

  list(
    exe_path = exe_path,
    run_dir = run_dir,
    nexus_path = dest_nexus,
    stdout = stdout_file,
    stderr = stderr_file,
    exit_status = if (is.numeric(res)) res else NA_integer_,
    error = attr(res, "error")
  )
}

#' Stop a MrBayes run started by \code{mrbayesRun}
#'
#' @description
#' Convenience wrapper to stop a MrBayes run returned by \code{\link{mrbayesRun}}.
#' For programmatic stopping, the run must have been started with \pkg{processx}
#' (i.e., \code{use_processx = TRUE}) and, for practical "stop anytime" control,
#' \code{background = TRUE}.
#'
#' @param res A result object returned by \code{\link{mrbayesRun}}.
#' @param force Logical. If \code{FALSE} (default), requests a graceful stop.
#' If \code{TRUE}, force-kills the process.
#'
#' @return Invisibly returns \code{TRUE} if a stop signal was sent, otherwise \code{FALSE}.
#'
#' @export
mrbayesStop <- function(res, force = FALSE) {
  if (is.null(res) || !is.list(res)) stop("res must be an object returned by mrbayesRun().")
  if (!is.logical(force) || length(force) != 1) stop("force must be TRUE/FALSE.")

  if (!"stop" %in% names(res) || !is.function(res$stop)) {
    stop(
      "No programmatic stop handle.\n",
      "Use processx and run: mrbayesRun(..., background = TRUE)\n",
      "Otherwise, stop interactively (Esc in RStudio)."
    )
  }

  invisible(res$stop(!isTRUE(force)))
}
