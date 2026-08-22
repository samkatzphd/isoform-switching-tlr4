#!/usr/bin/env Rscript

########################################################
# About:
# This script identifies and quantifies isoform switches
# using the R package IsoformSwitchAnalyzeR. It requires
# the per-sample isoform quantification files generated
# by RSEM as input.
########################################################


# Load the IsoformSwitchAnalyzeR library
suppressMessages(library(IsoformSwitchAnalyzeR))  # Bioconductor package 
suppressMessages(library(argparse))               # CRAN package


###################################
# Helper functions
###################################
# Add a timestamp to a log message
timestamp <- function(...) { cat("[", format(Sys.time()), "]", ..., "\n") ; }

# Create an switchAnalyzeRlist object
# by import RSEM counts, design, and
# annotation information:
# https://rdrr.io/bioc/IsoformSwitchAnalyzeR/man/importRdata.html
isa_import <- function(
        rsem_input_dir,      # Directory with per-sample RSEM isoform directories
        output_dir,          # Output directory to write the results
        sample_sheet,        # Sample sheet (TSV), needs sampleID and condition cols
        gtf_file,            # GTF file used for alignment and isoform quantification,
        transcripts_fa,      # Path to the transcriptomic fasta file
        condition_1,         # WT group in the 'KO vs. WT' contrast
        condition_2,         # KO group in the 'KO vs. WT' contrast
        alpha_filter = 0.05, # Adjusted p-value filter
        dif_filter = 0.1,    # Differential isoform fraction filter
        iu_filter = 0.01     # Isform usage pre-filter
    ) {

    # Import the sample sheet with
    # sample to group mappings,
    # which is used to create the
    # design matrix, must contain
    # the following columns:
    #     - 'sampleID'
    #     - 'condition'
    design_matrix <- read.table(
        sample_sheet,
        header = TRUE,
        sep = "\t",
        quote = ""
    )

    # Filter the design matrix to only
    # include the samples of interest
    design_matrix <- design_matrix[
        design_matrix$condition == condition_1 | design_matrix$condition == condition_2,
    ]

    # Create a contrast/comparison
    # data frame object, must contain
    # the following columns:
    #     - 'condition_1'
    #     - 'condition_2'
    # dIF = condition_2 - condition_1
    comparison <- data.frame(
        condition_1 = condition_1,   # Baseline/Control group
        condition_2 = condition_2    # Treatment/Case/Exerimental group
    )

    # Import the per-sample RSEM isoform 
    # quantification results
    timestamp('Started running importIsoformExpression step...')
    rsem <- importIsoformExpression(
        parentDir = rsem_input_dir,
        addIsofomIdAsColumn = TRUE,
        normalizationMethod = 'TMM',
        calculateCountsFromAbundance = TRUE,
        showProgress = TRUE
    )

    # Create SwitchAnalyzeRlist object
    # from provided inputs
    timestamp('Started running importRdata step...')
    isa_list <- importRdata(
        isoformCountMatrix = rsem$counts,
        isoformRepExpression = rsem$abundance,
        isoformExonAnnoation = gtf_file,
        designMatrix = design_matrix,
        comparisonsToMake = comparison,
        showProgress = TRUE,
        isoformNtFasta = transcripts_fa
    )


    # Filter the SwitchAnalyzeRlist,
    # extra options are current set
    # to the packages defaults
    timestamp('Started running preFilter step...')
    isa_list <- preFilter(
        switchAnalyzeRlist = isa_list,
        geneExpressionCutoff = 1,
        isoformExpressionCutoff = 0,
        IFcutoff = iu_filter,
        alpha = alpha_filter,
        dIFcutoff = dif_filter,
        quiet = FALSE
    )

    return(isa_list)
}


###################################
# Parse Arguements
###################################

# Pass command line args to main
# create parser object
timestamp('Started running parsing args...')
parser <- ArgumentParser()

# Parent directory containing the 
# per-sample RSEM output files
parser$add_argument(
    "-i", "--input_rsem_isoform_directory",
    help = "Directory containing the per-sample RSEM isoform directories, required",
    type = "character",
    required = TRUE
)

# Output directory to write the
# IsoformSwitchAnalyzeR results
parser$add_argument(
    "-o", "--output_directory",
    help = "Path to output directory, default: IsoformSwitchAnalyzeR_results",
    type = "character",
    required = FALSE,
    default = "IsoformSwitchAnalyzeR_results"
)

# Path to the GRCh38 plus ISO-seq
# combined GTF file
parser$add_argument(
    "-g", "--gtf_file",
    help = "Path to the GTF file used for alignment and isoform quantification, required",
    type = "character",
    required = TRUE
)

# Path to transcriptomic FASTA file
parser$add_argument(
    "-t", "--transcriptome_fa",
    help = "Path to the transcriptomic fasta file, required",
    type = "character",
    required = TRUE
)

# Sample sheet in TSV format
# containing group information
# for each samples
parser$add_argument(
    "-s", "--sample_sheet",
    help = "Path to the sample sheet, must contain sampleID & condition columns, required",
    type = "character",
    required = TRUE
)

# Group1 in a contrast to analyze,
# the resulting contrast will be
# group1 vs. group2
parser$add_argument(
    "-c1", "--condition_1",
    help = "Baseline/control group in the contrast, it will be WT group in the contrast 'KO vs. WT', required",
    type = "character",
    required = TRUE
)

# Group2 in a contrast to analyze,
# the resulting contrast will be
# group1 vs. group2
parser$add_argument(
    "-c2", "--condition_2",
    help = "Experimental/case group in the contrast, it will be KO group in the contrast 'KO vs. WT', required",
    type = "character",
    required = TRUE
)

# Alpha or p-value threshold for
# filtering the results, default
# has no filtering
parser$add_argument(
    "-p", "--pvalue_filter",
    help = "Alpha or adjusted p-value threshold for filtering the results, default: 0.05",
    type = "double",
    required = FALSE,
    default = 0.05
)

# Differential isoform fraction
# threshold for filtering the
# results, analogous to a fold-
# change, default has no filter
parser$add_argument(
    "-f", "--dif_filter",
    help = "Differential isoform fraction threshold, similar to fold-change, for filtering the results, default: 0.1",
    type = "double",
    required = FALSE,
    default = 0.1
)

# Isoform fraction threshold for
# filtering the results, analogous to a fold-
# change, default has no filter
parser$add_argument(
    "-u", "--iu_filter",
    help = "Isoform usage threshold for pre-filtering, default: 0.01",
    type = "double",
    required = FALSE,
    default = 0.01
)

# Runs extra downstream analysis if
# provided, like analyzeAlternativeSplicing
parser$add_argument(
    "--run_extra_analysis",
    help = "Run extra analysis like analyzeAlternativeSplicing, default: FALSE",
    action = "store_true",
    default = FALSE
)

# Parse the command line arguments
args <- parser$parse_args()

# Create output directory
# if it does not exist
outdir <- args$output_directory
# NOTE: dIF = condition_2 - condition_1
comparison <- paste(args$condition_2, args$condition_1, sep = "-")
prefix <- file.path(outdir, comparison)
dir.create(
    file.path(outdir),
    showWarnings = FALSE,
    recursive = TRUE
)

# Main Entry Point of Program
isa_list <- isa_import(
    rsem_input_dir = args$input_rsem_isoform_directory,
    output_dir = args$output_directory,
    sample_sheet = args$sample_sheet,
    gtf_file = args$gtf_file,
    condition_1 = args$condition_1,
    condition_2 = args$condition_2,
    transcripts_fa = args$transcriptome_fa,
    alpha_filter = args$pvalue_filter,
    dif_filter = args$dif_filter,
    iu_filter = args$iu_filter
)

# Print summary
summary(isa_list)

# Run the main analysis, test for 
# isoform switches using DEX-seq, 
# setting reduceToSwitchingGenes
# to TRUE will cause the function
# to subset the switchAnalyzeRlist
# to the genes which each contain 
# at least one differential used 
# isoform, as indicated by the 
# alpha and dIFcutoff cutoffs
timestamp('Started running isoformSwitchTestDEXSeq step...')
isa_list <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = isa_list,
    reduceToSwitchingGenes = TRUE,
    alpha = args$pvalue_filter,
    dIFcutoff = args$dif_filter,
    showProgress = TRUE,
    quiet = FALSE
)

# Predict ORF/CDS regions using
# analyzeORF function, this is
# needed because this project
# contains known and novel 
# isoforms
timestamp('Started running analyzeORF step...')
isa_list <- analyzeORF(isa_list, showProgress = TRUE)

# Run alternative splicing analysis
# to quantifty different alternative
# splicing events, such as exon skipping,
# alternative 5' and 3' splice sites,
# intron retention, etc.
if (args$run_extra_analysis){
    timestamp('Started running analyzeAlternativeSplicing step...')
    isa_list <- analyzeAlternativeSplicing(
        switchAnalyzeRlist = isa_list,
        onlySwitchingGenes = TRUE,
        alpha = args$pvalue_filter,
        dIFcutoff = args$dif_filter,
        showProgress = TRUE,
        quiet = FALSE
    ) 

    # Visualize the results
    # Create splicing summary plot
    pdf(paste(prefix, "_splicing_summary.pdf", sep = ""))
    timestamp('Started running extractSplicingSummary step...')
    extractSplicingSummary(
        isa_list,
        splicingToAnalyze = 'all',
        asFractionTotal = FALSE,
        onlySigIsoforms = FALSE,
        plotGenes = FALSE,
        localTheme = theme_bw(),
    )
    dev.off()

    # Create splicing enrichment plot
    pdf(paste(prefix, "_splicing_enrichment.pdf", sep = ""))
    timestamp('Started running extractSplicingEnrichment step...')
    splicing_enrichment <- extractSplicingEnrichment(
        isa_list,
        splicingToAnalyze = 'all',
        alpha = args$pvalue_filter,
        dIFcutoff = args$dif_filter,
        plot = TRUE,
        localTheme = theme_bw(base_size = 14),
        minEventsForPlotting = 10,
        returnResult = TRUE,
        returnSummary = TRUE
    )
    dev.off()

    # Create genome wide splicing plot
    pdf(paste(prefix, "_genome_wide_splicing.pdf", sep = ""))
    timestamp('Started running extractSplicingGenomeWide step...')
    genome_wide_splicing <- extractSplicingGenomeWide(
        isa_list,
        featureToExtract = 'all',
        splicingToAnalyze = 'all',
        alpha = args$pvalue_filter,
        dIFcutoff = args$dif_filter,
        log2FCcutoff = 1,
        violinPlot = TRUE,
        alphas = c(0.05, 0.001),
        localTheme = theme_bw(),
        plot = TRUE,
        returnResult = TRUE
    )
    dev.off()
}

# Extract the top gene and isoform switches
timestamp('Started running extractTopSwitches step for genes...')
top_gene_switches <- extractTopSwitches(
    switchAnalyzeRlist = isa_list,
    filterForConsequences = FALSE,
    extractGenes = TRUE,    # extract genes
    alpha = args$pvalue_filter,
    dIFcutoff = args$dif_filter,
    sortByQvals = TRUE,
    n = Inf
)
timestamp('Started running extractTopSwitches step for isoforms...')
top_isoform_switches <- extractTopSwitches(
    switchAnalyzeRlist = isa_list,
    filterForConsequences = FALSE,
    extractGenes = FALSE,   # extract isoforms
    alpha = args$pvalue_filter,
    dIFcutoff = args$dif_filter,
    sortByQvals = TRUE,
    n = Inf
)

# Write the results to output files
# Gene switches
timestamp('Started writing output TSV files...')
write.table(
    top_gene_switches,
    file = paste(prefix, "_top_gene_switches.tsv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
# Isoform switches
write.table(
    top_isoform_switches,
    file = paste(prefix, "_top_isoform_switches.tsv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
if (args$run_extra_analysis){
    # Splicing enrichment
    write.table(
        splicing_enrichment,
        file = paste(prefix, "_splicing_enrichment.tsv", sep = ""),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    # Genome wide splicing
    write.table(
        genome_wide_splicing,
        file = paste(prefix, "_genome_wide_splicing.tsv", sep = ""),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}
# Save all R objects to an Rdata
# file for future figures or analysis
timestamp('Started writing output Rdata file...')
save.image(file = paste(prefix, "_IsoformSwitchAnalyzeR.Rdata", sep = ""))
