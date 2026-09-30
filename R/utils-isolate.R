# Isolated fan-out ####
#
# Runs a kernel function (R/kernel-*.R) over a list of plain-data tasks on a
# pool of mirai daemons that never load GiottoDisk. A daemon that loads the
# package pays ~7 s and holds its namespace until it exits; one that receives
# only kernel code pays the load time of arrow. The pool is started for the
# call and shut down when it returns, so it does not share workers with the
# user's future plan, and Arrow memory it retains goes with it. adr/0019 has
# the measurement and the alternatives (fork, future, a persistent pool).


# Workers for an isolated fan-out: `.par_workers()` when mirai and codetools
# are installed, otherwise 1, which tells the caller to run serially.
.isolated_workers <- function() {
    n <- .par_workers()
    if (n <= 1L) return(1L)
    if (!requireNamespace("mirai", quietly = TRUE) ||
        !requireNamespace("codetools", quietly = TRUE)) {
        return(1L)
    }
    n
}

# Copy `entry` and every GiottoDisk function it reaches into one environment
# parented on globalenv(), so nothing in the result refers to the namespace and
# unserializing it loads no package. Functions GiottoDisk imports are recorded
# by package and resolved on the worker from their own namespace.
#
# Reaching an S4 generic is an error rather than something to copy: it is
# outside the kernel, and dispatching it needs GiottoDisk loaded. That includes
# base names GiottoDisk turns into generics (`unique`), which a worker would
# silently resolve to base instead.
.kernel_bundle <- function(entry) {
    ns <- environment(.kernel_bundle)
    imp_env <- parent.env(ns)
    env <- new.env(parent = globalenv())
    imports <- list()
    todo <- codetools::findGlobals(entry)
    seen <- character()
    while (length(todo)) {
        f <- todo[[1L]]
        todo <- todo[-1L]
        if (f %in% seen) next
        seen <- c(seen, f)
        if (exists(f, envir = ns, inherits = FALSE)) {
            obj <- get(f, envir = ns)
            if (methods::is(obj, "genericFunction")) {
                stop("[.kernel_bundle] `", f, "` is an S4 generic, so the ",
                     "code reaching it cannot run without GiottoDisk. Kernel ",
                     "code must stay within R/kernel-*.R.", call. = FALSE)
            }
            if (is.function(obj)) {
                environment(obj) <- env
                todo <- c(todo, codetools::findGlobals(obj))
            }
            assign(f, obj, envir = env)
        } else if (exists(f, envir = imp_env, inherits = FALSE)) {
            fenv <- environment(get(f, envir = imp_env))
            if (!is.null(fenv)) imports[[f]] <- getNamespaceName(fenv)
        }
    }
    environment(entry) <- env
    list(entry = entry, env = env, imports = imports)
}

# One task on a worker. Shipped with its environment reset to globalenv(), so
# it must use nothing but base R.
.isolated_run <- function(x, bundle, opts) {
    old <- options(opts)
    on.exit(options(old), add = TRUE)
    for (f in names(bundle$imports)) {
        assign(f, getExportedValue(bundle$imports[[f]], f),
            envir = bundle$env)
    }
    bundle$entry(x)
}

# Map `entry` (a kernel function) over `X` on a pool of `n_workers` daemons
# started for this call, and return the results in order.
#
# Each task carries the bundle and a snapshot of the `giottodisk.*` options,
# so no worker depends on state from an earlier task. The pool runs under its
# own compute profile, which leaves the default one -- where future.mirai puts
# a user's plan -- untouched, and `on.exit` shuts it down on error or
# interrupt. Daemons launch with R_LIBS set to this session's `.libPaths()`, so
# a library prepended at runtime is where they find arrow, and mirai itself.
.isolated_map <- function(X, entry, n_workers, site) {
    bundle <- .kernel_bundle(entry)
    run <- .isolated_run
    environment(run) <- globalenv()
    opts <- options()[grep("^giottodisk[.]", names(options()))]

    profile <- paste0("giottodisk_", .make_uid())
    old_libs <- Sys.getenv("R_LIBS", unset = NA)
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    mirai::daemons(min(n_workers, length(X)), .compute = profile)
    if (is.na(old_libs)) Sys.unsetenv("R_LIBS")
    else Sys.setenv(R_LIBS = old_libs)
    on.exit(mirai::daemons(0L, .compute = profile), add = TRUE)

    # A daemon that cannot start never connects, and mirai_map() would wait
    # on it indefinitely; one timed round trip turns that into an error.
    ping <- mirai::mirai(TRUE, .compute = profile, .timeout = 60000)[]
    if (mirai::is_error_value(ping)) {
        stop("[", site, "] could not start fan-out workers: ",
             as.character(ping), call. = FALSE)
    }

    res <- mirai::mirai_map(X, run,
        .args = list(bundle = bundle, opts = opts),
        .compute = profile)[]
    failed <- vapply(res, mirai::is_error_value, logical(1L))
    if (any(failed)) {
        k <- which(failed)[1L]
        stop("[", site, "] fan-out task ", k, " of ", length(X), " failed: ",
             as.character(res[[k]]), call. = FALSE)
    }
    res
}
