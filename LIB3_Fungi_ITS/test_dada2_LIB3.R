library(dada2)

# Paths
path <- "/home/sofija/INTERNSHIP/samples/Capreolus_capreolus/Fungi_ITS"

filt_path <- file.path(path, "qualFiltered_out")
out_path  <- file.path(path, "test_dada2_out")

dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# Find filtered reads
fnFs <- sort(list.files(
  filt_path,
  pattern = "_R1\\.fq\\.gz$",
  full.names = TRUE
))

fnRs <- sort(list.files(
  filt_path,
  pattern = "_R2\\.fq\\.gz$",
  full.names = TRUE
))

# Check pairing
sampleF <- sub("_R1\\.fq\\.gz$", "", basename(fnFs))
sampleR <- sub("_R2\\.fq\\.gz$", "", basename(fnRs))

stopifnot(identical(sampleF, sampleR))

# TEST ONLY: first two samples
fnFs <- fnFs[1:2]
fnRs <- fnRs[1:2]

samples <- sampleF[1:2]

cat("TEST samples:\n")
print(samples)

# 1. Learn error rates
errF <- learnErrors(
  fnFs,
  multithread = TRUE
)

errR <- learnErrors(
  fnRs,
  multithread = TRUE
)

saveRDS(
  errF,
  file.path(out_path, "test_errF.rds")
)

saveRDS(
  errR,
  file.path(out_path, "test_errR.rds")
)

# 2. Denoise
dadaFs <- dada(
  fnFs,
  err = errF,
  multithread = TRUE
)

dadaRs <- dada(
  fnRs,
  err = errR,
  multithread = TRUE
)

names(dadaFs) <- samples
names(dadaRs) <- samples

# 3. Merge paired reads
mergers <- mergePairs(
  dadaFs,
  fnFs,
  dadaRs,
  fnRs,
  verbose = TRUE
)

# 4. Construct ASV table
seqtab <- makeSequenceTable(mergers)

saveRDS(
  seqtab,
  file.path(out_path, "test_seqtab.rds")
)

# 5. Remove chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

saveRDS(
  seqtab.nochim,
  file.path(out_path, "test_seqtab_nochim.rds")
)

cat("\nTEST completed successfully.\n")
cat("Sequence table dimensions:\n")
print(dim(seqtab))

cat("Non-chimeric sequence table dimensions:\n")
print(dim(seqtab.nochim))
