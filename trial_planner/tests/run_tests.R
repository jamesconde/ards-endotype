# Rscript tests/run_tests.R   (from the repo root, env trialplanner-r)
# Exits non-zero on any failure so a CI job or a shell && chain cannot mistake a failure for a pass.
res <- testthat::test_dir(file.path("tests", "testthat"), reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n%d tests, %d expectations, %d failed, %d skipped, %.1f s\n",
            nrow(df), sum(df$nb), sum(df$failed), sum(df$skipped), sum(df$real)))
quit(status = if (sum(df$failed) > 0 || any(df$error)) 1 else 0)
