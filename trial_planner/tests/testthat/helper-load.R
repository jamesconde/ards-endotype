ROOT <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in list.files(file.path(ROOT, "app", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
REF <- jsonlite::fromJSON(file.path(ROOT, "reference", "python_reference.json"), simplifyVector = FALSE)
DEFAULTS <- jsonlite::fromJSON(file.path(ROOT, "app", "defaults.json"), simplifyVector = FALSE)
rows_to_df <- function(rows) {
  df <- do.call(rbind, lapply(rows, as.data.frame))
  df[] <- lapply(df, function(x) if (is.numeric(x)) as.numeric(x) else x)
  df
}
