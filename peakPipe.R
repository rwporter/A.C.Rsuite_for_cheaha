#!/usr/bin/env Rscript
# peakPipe_noIDR.R — IDR removed; auto-detects/creates group peak sets for peakQuan.R

suppressWarnings(suppressMessages({
  if (!require(optparse)) install.packages("optparse", repos="https://cloud.r-project.org/")
  library(optparse)
}))

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]),
              help="Output directory path [default= %default]"),
  make_option(c("-g", "--genome"), type="character", default="mm10",
              help="Genome: mm10, hg38 [default= %default]"),
  make_option(c("-a", "--assay"), type="character", default="atac",
              help="Assay: atac, chip [default= %default]"),
  make_option(c("-c", "--homerCor"), type="character", default=" -minDist 200 -size 200",
              help="Extra args for homer peak calling (correlation) [default= '%default']"),
  make_option(c("-q", "--homerQuan"), type="character", default=" -size given -pc 3",
              help="Extra args for annotatePeaks.pl (quant) [default= '%default']"),
  make_option(c("-d", "--distal"), type="numeric", default=3000,
              help="Distance from TSS to call enhancer [default= %default]"),
  make_option(c("-f", "--minNormTag"), type="numeric", default=1,
              help="Min normalized tag count per peak [default= %default]"),
  make_option(c("-n", "--minNormSample"), type="numeric", default=1,
              help="Min number of samples passing minNormTag [default= %default]"),
  make_option(c("-m", "--maxPeak"), type="numeric", default=10000,
              help="Max peaks for heatmap [default= %default]"),
  make_option(c("-l", "--logFC"), type="numeric", default=1,
              help="Log2FC threshold [default= %default]"),
  make_option(c("-t", "--track"), type="character",
              default=paste(Sys.info()[["user"]], gsub("\\.txt$","", basename(args[1])), sep="_"),
              help="Track name (optional) [default= '%default']"),
  make_option(c("-p", "--padj"), type="numeric", default=0.05,
              help="Adjusted P-value threshold [default= %default]")
)

opt_parser = OptionParser("\n\t%prog path/to/sample_definition.txt [options]",
                          option_list=option_list, prog="peakPipe_noIDR.R")

if (length(args) < 1) { print_help(opt_parser); stop("path/to/sample_definition.txt is required.\n", call.=FALSE) }

strSample <- args[1]
opt <- parse_args(opt_parser, args[-1])

strOutput <- file.path(opt$out, "")
if (!dir.exists(strOutput)) {
  if (!dir.create(strOutput, recursive=TRUE, showWarnings=FALSE)) {
    stop(sprintf("Cannot create output folder: %s (parent: %s)", strOutput, dirname(strOutput)))
  }
}

## 1) Alignment stats
cat("\nGet the alignment status: --------------------\n")
strCMD <- sprintf("alignStats.R %s > %s", shQuote(strSample), shQuote(file.path(strOutput, "alignStats.txt")))
cat(strCMD, "\n")
if (system(strCMD) != 0) stop("Error in alignment status!")

## 2) Pairwise correlation (also produces replicate/group peaks in tmp)
cat("\nCalculate the pair-wised correlation: --------------------\n")
corOut <- file.path(strOutput, "peakCor")
strCMD <- paste(
  "peakCor.R", shQuote(strSample),
  "-o", shQuote(paste0(corOut, "/")),
  "-g", shQuote(opt$genome), "-a", shQuote(opt$assay),
  paste0("-c '", opt$homerCor, "'")
)
cat(strCMD, "\n")
if (system(strCMD) != 0) stop("Error in Calculating the pair-wised correlation!")

## 3) Quantification — auto-detect or create group peak sets (no IDR)
cat("\nQuantify peaks for each sample (no IDR): --------------------\n")
tmpdir <- file.path(corOut, "peakCor_tmp")
if (!dir.exists(tmpdir)) stop(sprintf("Expected directory not found: %s", tmpdir))

# Helper: create group merged peaks if missing, using replicate peaks for that group
merge_group_if_missing <- function(group, rep_files) {
  target <- file.path(tmpdir, paste0(group, ".peaks"))
  if (!file.exists(target)) {
    cmd <- paste("mergePeaks", paste(shQuote(rep_files), collapse=" "), ">", shQuote(target))
    cat("Creating missing group peaks via mergePeaks for", group, ":\n ", cmd, "\n")
    status <- system(cmd)
    if (status != 0 || !file.exists(target)) stop(sprintf("Failed to create merged group peaks for '%s'", group))
  }
  return(target)
}

all_peaks <- list.files(tmpdir, pattern="\\.peaks$", full.names=TRUE)
if (length(all_peaks) == 0) stop(sprintf("No *.peaks found in %s", tmpdir))

# Infer groups by stripping the first underscore and following suffix (e.g., ctrl_JM1_ATAC.peaks -> group 'ctrl')
bn <- basename(all_peaks)
groups_inferred <- sub("_.+$", "", sub("\\.peaks$", "", bn), perl=TRUE)

# Files already present as group-level merges (no underscore in basename)
is_group_level <- !grepl("_", sub("\\.peaks$", "", bn))
group_level_files <- all_peaks[is_group_level]

# For groups without a group-level file, build one from their replicate peaks
needed_groups <- unique(groups_inferred)
present_groups <- unique(sub("\\.peaks$", "", basename(group_level_files)))
missing_groups <- setdiff(needed_groups, present_groups)

if (length(missing_groups) > 0) {
  for (g in missing_groups) {
    rep_files <- all_peaks[groups_inferred == g & grepl("_", bn)]
    if (length(rep_files) == 0) stop(sprintf("No replicate peak files found for group '%s' to build merged peaks.", g))
    group_level_files <- c(group_level_files, merge_group_if_missing(g, rep_files))
  }
}

# Deduplicate, sort for stability
group_level_files <- sort(unique(group_level_files))
if (length(group_level_files) == 0) stop("No group-level peak files available for quantification.")

cat("Using group peak sets (-p):\n  ", paste(group_level_files, collapse=",\n  "), "\n", sep="")

strCMD <- paste(
  "peakQuan.R", shQuote(strSample),
  "-o", shQuote(file.path(strOutput, "peakQuan/")),
  "-g", shQuote(opt$genome), "-a", shQuote(opt$assay),
  "-p", shQuote(paste(group_level_files, collapse=",")),
  "-d", shQuote(opt$distal),
  "-t", shQuote(opt$track),
  paste0("-c '", opt$homerQuan, "'")
)
cat(strCMD, "\n")
if (system(strCMD) != 0) stop("Error in Quantifying peaks for each sample!")

## 4) Differential peaks
cat("\nFind differential peaks for pair-wised groups: --------------------\n")
strAnno <- file.path(strOutput, "peakQuan", "allRawTags.txt")
if (!file.exists(strAnno)) stop(sprintf("Expected quant matrix not found: %s", strAnno))

strCMD <- paste(
  "peakDiff.R", shQuote(strSample),
  "-o", shQuote(file.path(strOutput, "peakDiff/")),
  "-g", shQuote(opt$genome), "-a", shQuote(opt$assay),
  "-q", shQuote(strAnno), "-d", shQuote(opt$distal),
  "-c", shQuote(opt$minNormTag), "-n", shQuote(opt$minNormSample),
  "-m", shQuote(opt$maxPeak),
  "-l", shQuote(opt$logFC), "-p", shQuote(opt$padj)
)
cat(strCMD, "\n")
if (system(strCMD) != 0) stop("Error in Finding differential peaks for pair-wised groups!")

warnings()
cat("\n", opt$assay, "peak analysis pipeline finished successfully (no IDR).\n", sep="")
