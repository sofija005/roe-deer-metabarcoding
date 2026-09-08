library(dada2)

# =========================
# Plants trnLc-trnLh - LIB2
# DADA2 quality filtering
# =========================

path <- "primersCut_out"

# Input files
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

# Check number of paired files
print(paste("Number of R1 files:", length(fnFs)))
print(paste("Number of R2 files:", length(fnRs)))

# Check that R1 and R2 sample names match
samplesF <- sub("_R1\\.fq\\.gz$", "", basename(fnFs))
samplesR <- sub("_R2\\.fq\\.gz$", "", basename(fnRs))

stopifnot(identical(samplesF, samplesR))

# Output directory
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

# Quality filtering
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

# Create summary table
summary <- as.data.frame(out)

summary$sample <- samplesF

summary$percent_retained <- round(
  summary$reads.out / summary$reads.in * 100,
  2
)

summary <- summary[, c(
  "sample",
  "reads.in",
  "reads.out",
  "percent_retained"
)]

write.csv(
  summary,
  "qualFiltered_out/seq_count_summary.csv",
  row.names = FALSE
)

print(summary)

