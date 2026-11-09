# The registry (adr/0023): versioned files committed only if the version is
# still free, edits that create a record or add tags and nothing else, and
# writers that lose a race retrying instead of overwriting each other.

# callers remove `td` with on.exit()
.reg_src <- function() {
    td <- tempfile("registry_")
    list(td = td, src = gDirSource(td))
}

# a test hook that runs `f` the first time the code under test reaches `point`
.interleave <- function(point, f) {
    done <- FALSE
    function(p) {
        if (identical(p, point) && !done) {
            done <<- TRUE
            f()
        }
    }
}

describe("edits", {
    create <- list(uid = "u1", add = list(giottosave = NA),
        create = list(time = "t0", store = "parquetStore", hash = "h"))
    tag_a <- list(uid = "u1", add = list(giottosave = "a.aaaaaaaa"))
    tag_b <- list(uid = c("u1", "u2"), add = list(giottosave = "b.bbbbbbbb"))

    it("give the same records in any order, applied once or twice", {
        one <- Reduce(.gdsrc_apply_edit, list(create, tag_a, tag_b), list())
        two <- Reduce(.gdsrc_apply_edit,
            list(tag_b, tag_a, create, tag_a, create, tag_b), list())
        expect_setequal(one$u1$giottosave, c("a.aaaaaaaa", "b.bbbbbbbb"))
        expect_setequal(two$u1$giottosave, one$u1$giottosave)
        expect_identical(two$u1[c("time", "store", "hash")],
            one$u1[c("time", "store", "hash")])
        expect_identical(one$u2$giottosave, "b.bbbbbbbb")
    })

    it("never change a field of a record that exists", {
        again <- create
        again$create$store <- "other"
        rec <- Reduce(.gdsrc_apply_edit, list(create, again), list())$u1
        expect_identical(rec$store, "parquetStore")
        expect_true(is.na(rec$giottosave))
    })

    it("only add snapshot tags to a record once it exists", {
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        u <- sourceWrite(x$src, data.frame(a = 1:3))@uid
        expect_error(x$src[u, "store"] <- "x", "fixed")
        expect_error(x$src[u] <- list(store = "x"), "fixed")
        x$src[u, "giottosave"] <- "qc.aaaaaaaa"
        x$src[u, "giottosave"] <- "qc.bbbbbbbb"
        expect_setequal(x$src[u, "giottosave"], c("qc.aaaaaaaa", "qc.bbbbbbbb"))
        expect_identical(x$src[u, "store"], "parquet")
    })
})

describe("commits", {
    it("leave giottodir.json as a stub stamped with the schema", {
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        stub <- .json_read(.gdsrc_json_path(x$td))
        expect_identical(stub$header$versions$manifest, .GDSRC_SCHEMA_VERSION)
        expect_length(stub$content, 0L)
        expect_s4_class(resolveSource(x$td), "gDirSource")
    })

    it("keep two versions, and remove edit files one commit after applying", {
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        for (i in 1:3) {
            sourceWrite(x$src, data.frame(a = i))
            .gdsrc_json_consolidate(x$td)
        }
        expect_identical(.gdsrc_registry_versions(x$td), 3:4) # 1 is the empty start
        pending <- list.files(.gdsrc_json_pending_dir(x$td))
        expect_length(pending, 1L)
        expect_identical(pending, .gdsrc_registry_load(x$td)$header$applied)
        expect_length(as.list(x$src), 3L)
    })

    it("by a writer that loses the race retry, keeping both writers' edits", {
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        a <- sourceWrite(x$src, data.frame(a = 1))@uid
        b <- NULL
        GiottoUtils::gwith_options(list(giottodisk.test_hook =
            .interleave("commit", function() {
                b <<- sourceWrite(x$src, data.frame(b = 1))@uid
                .gdsrc_json_consolidate(x$td) # commits the version first
            })
        ), {
            reg <- .gdsrc_json_consolidate(x$td)
        })
        expect_setequal(names(reg$content), c(a, b))
        expect_identical(reg$n, 2L) # the loser's content was never committed
    })

    it("fall back to the version below when the newest is cut short", {
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        u <- sourceWrite(x$src, data.frame(a = 1))@uid
        .gdsrc_json_consolidate(x$td)
        writeLines("{\"header\":", .gdsrc_registry_file(x$td, 2L))
        expect_identical(.gdsrc_registry_load(x$td)$n, 1L)
        # the edit file is still there, so the next read restores the record
        reg <- .gdsrc_json_consolidate(x$td)
        expect_identical(reg$n, 3L)
        expect_true(u %in% names(reg$content))
    })

    it("are skipped for a project that cannot be written", {
        skip_on_os("windows")
        skip_if(identical(Sys.info()[["user"]], "root"), "root writes anyway")
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        u <- sourceWrite(x$src, data.frame(a = 1))@uid
        Sys.chmod(x$td, "555")
        on.exit(Sys.chmod(x$td, "755"), add = TRUE, after = FALSE)
        reg <- .gdsrc_json_read(x$td, consolidate = TRUE)
        expect_true(u %in% names(reg$content))
        expect_identical(.gdsrc_registry_versions(x$td), 1L)
    })

    it("explain a filesystem without hard links, or rename when told to", {
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        sourceWrite(x$src, data.frame(a = 1))
        local_mocked_bindings(.gdsrc_link = function(from, to) FALSE)
        expect_error(.gdsrc_json_consolidate(x$td), "hard links")
        GiottoUtils::gwith_options(list(giottodisk.commit_by_rename = TRUE), {
            reg <- .gdsrc_json_consolidate(x$td)
        })
        expect_length(reg$content, 1L)
    })

    it("by two processes tagging one artifact lose no tags", {
        skip_on_cran()
        skip_if_not_installed("callr")
        skip_if_not_installed("pkgload")
        pkg <- normalizePath(test_path("..", ".."))
        skip_if_not(file.exists(file.path(pkg, "DESCRIPTION")), "needs the source tree")
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        u <- sourceWrite(x$src, data.frame(a = 1))@uid
        .gdsrc_json_consolidate(x$td)

        worker <- function(pkg, td, uid, who, n) {
            suppressMessages(pkgload::load_all(pkg, quiet = TRUE))
            for (i in seq_len(n)) {
                .gdsrc_json_add_tags(td, uid, sprintf("%s.%08d", who, i))
                .gdsrc_json_consolidate(td)
            }
            TRUE
        }
        ps <- lapply(c("a", "b"), function(w) {
            callr::r_bg(worker, list(pkg, x$td, u, w, 20L))
        })
        for (p in ps) p$wait(timeout = 120000)
        expect_true(all(vapply(ps, function(p) isTRUE(p$get_result()), NA)))
        tags <- x$src[u, "giottosave"]
        expect_length(tags, 40L)
    })
})

describe("sourcePrune", {
    it("drops tags whose snapshot is gone", {
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        tag <- basename(dirname(.write_snapshot_fixture(x$td, "qc", list(a = 1))))
        u <- sourceWrite(x$src, data.frame(a = 1))@uid
        x$src[u, "giottosave"] <- c(tag, "gone.abcdefgh")
        suppressMessages(sourcePrune(x$src))
        expect_identical(x$src[u, "giottosave"], tag)
    })

    it("that loses the race to a save keeps what the save protects", {
        x <- .reg_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        u <- sourceWrite(x$src, data.frame(a = 1))@uid
        .gdsrc_json_consolidate(x$td) # recorded, untagged: prunable
        GiottoUtils::gwith_options(list(giottodisk.test_hook =
            .interleave("prune", function() {
                f <- .write_snapshot_fixture(x$td, "qc", list(a = 1))
                .gdsrc_json_add_tags(x$td, u, basename(dirname(f)))
                .gdsrc_json_consolidate(x$td)
            })
        ), {
            pruned <- suppressMessages(sourcePrune(x$src))
        })
        expect_length(pruned, 0L)
        expect_true(dir.exists(.gdsrc_artifact_dir(x$td, u)))
    })
})
