#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library("annotatr")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

REs_file <- file(args[1]) # Update this with Kevin's annotation
RE_regions = read_regions(con = REs_file, genome = 'hg19', format = 'bed')

print(RE_regions)

# Select annotations for intersection with regions
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)

# Intersect the regions we read in with the annotations
RE_annotated = annotate_regions(
  regions = RE_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(RE_annotated)

# Coerce to a data.frame
df_RE_annotated = data.frame(RE_annotated)

# See the GRanges column of dm_annotaed expanded
print(head(df_RE_annotated))

# Output filename
output_filename <- sub("\\.([a-zA-Z]+)$", "_annotated.\\1", args[1])

# Write to file
write.table(df_RE_annotated, output_filename)

# Find the number of regions per annotation type
RE_annsum = summarize_annotations(
  annotated_regions = RE_annotated,
  quiet = TRUE)
print(RE_annsum)
