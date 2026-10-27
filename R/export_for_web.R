# =====================================================================
# export_for_web.R
# Export coverage profiles for the Stalling Profile Explorer.
#
# Run this in the R session where you built `cov_long`, i.e. after the
# left_join with STR.info.light.df. cov_long needs the columns
#   seq.id, pos, depth, sample
# plus whatever metadata you want to filter on in the website.
#
# It writes three files into OUT_DIR:
#   sequences.tsv.gz  one row per sequence: seq.id + metadata
#   profiles.tsv.gz   one row per sequence x sample: seq.id, sample, p1 ... p150
#                     (the read count at each position, NOT normalised;
#                      the website normalises on the fly)
#   manifest.json     tells the website which files to load
#
# data/local/ is in .gitignore, so these files are never pushed to GitHub.
# =====================================================================

library(data.table)

OUT_DIR   <- "~/Projects/stalling-profile-explorer/data/local"
DATASET   <- "260805-PhenDC3TP"   # name shown in the website header
STR_START <- 100L                 # STR start position; dist = STR_START - pos
POSITIONS <- 1:150                # reference positions to export

# Metadata columns to carry over (any that are missing are skipped,
# and you can add more: every extra column becomes a filter).
META_COLS <- c("seq", "motif.type", "motif.length", "periodicity",
               "STR.motif.ssDNA", "STR.motif.dsDNA", "pred.struct")

OUT_DIR <- path.expand(OUT_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dt <- as.data.table(cov_long)

# ---- 1. Sequence metadata: one row per seq.id -----------------------
meta.cols <- intersect(META_COLS, names(dt))
sequences <- unique(dt[, c("seq.id", meta.cols), with = FALSE], by = "seq.id")
fwrite(sequences, file.path(OUT_DIR, "sequences.tsv.gz"), sep = "\t", na = "NA")

# ---- 2. Profiles: one row per seq.id x sample, one column per position
prof <- dt[pos %in% POSITIONS, .(seq.id, sample, pos, depth)]
prof[is.na(depth), depth := 0]
prof[, pos := factor(pos, levels = POSITIONS)]

# drop = c(TRUE, FALSE) keeps a column for every position, even ones with
# no reads in any sequence; missing positions are filled with 0.
wide <- dcast(
  prof, seq.id + sample ~ pos,
  value.var = "depth", fun.aggregate = sum, fill = 0,
  drop = c(TRUE, FALSE)
)
setnames(wide, as.character(POSITIONS), paste0("p", POSITIONS))
fwrite(wide, file.path(OUT_DIR, "profiles.tsv.gz"), sep = "\t")

# ---- 3. Manifest -----------------------------------------------------
writeLines(sprintf(
'{
  "name": "%s",
  "files": ["sequences.tsv.gz", "profiles.tsv.gz"],
  "origin": %d,
  "created": "%s"
}', DATASET, as.integer(STR_START), format(Sys.time(), "%Y-%m-%d %H:%M")),
  file.path(OUT_DIR, "manifest.json"))

message(
  "Wrote ", nrow(sequences), " sequences and ", nrow(wide),
  " sequence x sample profiles (", uniqueN(wide$sample), " samples) to ", OUT_DIR
)
