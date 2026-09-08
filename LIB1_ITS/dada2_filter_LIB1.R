library(dada2)

path <- "primersCut_out"

fnFs <- sort(list.files(
  path,
  pattern = "_R1\\.fq\\.gz$",
  full.names = TRUE
))

fnRs <- sort(list.files(
  path,
  pattern = "_R2\\.fq\\.gz$",
  full.names = TRUE
))

filt_path <- "qualFiltered_out"
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(
  filt_path,
  basename(fnFs)
)

filtRs <- file.path(
  filt_path,
  basename(fnRs)
)

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  truncLen = c(0, 0),
  minLen = 100,
  maxLen = 600,
  matchIDs = TRUE,
  compress = TRUE,
  multithread = TRUE
)

summary <- as.data.frame(out)

summary$percent_retained <- round(
  summary$reads.out / summary$reads.in * 100,
  2
)

write.csv(
  summary,
  "qualFiltered_out/seq_count_summary.csv"
)

print(summary)
