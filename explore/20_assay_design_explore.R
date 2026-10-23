#!/usr/bin/env Rscript
# ============================================================================
# 20_assay_design_explore.R
#
# Standalone explorer for the assay-design stage. Reads ONLY CSVs and writes
# NOTHING. Run it line by line, change the cuts, re-plot. No bootstrap, no
# pipeline objects, no .rds.
#
#   Rscript explore/20_assay_design_explore.R          # prints summaries
#   or open in RStudio and step through the blocks
#
# Needs: dplyr, tidyr, ggplot2  (patchwork only for the composed funnel)
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

# Point this at the repo root if you are running from elsewhere.
ROOT <- Sys.getenv("ISO_ROOT", unset = ".")
tp <- function(f) file.path(ROOT, "results/tables", f)

assign   <- read.csv(tp("assay_assignment.csv"),       stringsAsFactors = FALSE)
order_s  <- read.csv(tp("assay_order_sheet.csv"),      stringsAsFactors = FALSE)
spec     <- read.csv(tp("assay_specificity.csv"),      stringsAsFactors = FALSE)
pairspec <- read.csv(tp("assay_pair_specificity.csv"), stringsAsFactors = FALSE)
prep     <- read.csv(tp("pair_selection_report.csv"),  stringsAsFactors = FALSE)
regions  <- read.csv(tp("discriminating_regions.csv"), stringsAsFactors = FALSE)
seqver   <- read.csv(tp("sequence_verification.csv"),  stringsAsFactors = FALSE)
detect   <- read.csv(tp("candidate_detectability_by_switch.csv"), stringsAsFactors = FALSE)

cat("candidates:", nrow(assign), "  oligos ordered:", nrow(order_s),
    "  oligos designed:", nrow(spec), "\n\n")

# ---------------------------------------------------------------------------
# 1. What kind of switch are these, and what assay did each one get?
# ---------------------------------------------------------------------------
assign %>% count(primary_event, sort = TRUE) %>% print()
assign %>% count(primary_assay, sort = TRUE) %>% print()

# Why so few are qPCR-addressable: an isoform that differs only at its 5' end
# is CONTAINED in its longer partner, so no short amplicon is specific to it.
assign %>%
  count(primary_event, qpcr_addressable) %>%
  pivot_wider(names_from = qpcr_addressable, values_from = n, values_fill = 0) %>%
  print()

# ---------------------------------------------------------------------------
# 2. The platform question: which switches needed PacBio?
#    This is the project's standing requirement -- novel-vs-known must stay
#    separable at every stage.
# ---------------------------------------------------------------------------
plat <- assign %>%
  left_join(detect %>% select(isoform_id, detect_class, detect_reason,
                              class_code_candidate, class_code_partner,
                              same_reference_transcript, illumina_resolvable),
            by = "isoform_id")
plat %>% count(detect_class) %>% print()

# The non-obvious case: both isoforms map to the SAME reference transcript, so
# no isoform-level short-read analysis can see the switch at any depth.
plat %>% filter(same_reference_transcript == "TRUE") %>%
  select(gene_name, class_code_candidate, class_code_partner, own_dIF) %>%
  arrange(desc(abs(own_dIF))) %>% print()

# ---------------------------------------------------------------------------
# 3. Specificity, measured rather than predicted
#    Binding several isoforms of the same gene is expected. Binding a DIFFERENT
#    gene is the defect. Plot the second, not the first.
# ---------------------------------------------------------------------------
sp <- spec %>%
  mutate(off = suppressWarnings(as.integer(offtarget_loci)),
         n   = suppressWarnings(as.integer(transcripts_full_match)))

cat("\noligos designed with zero off-target loci:",
    sum(sp$off == 0, na.rm = TRUE), "of", nrow(sp), "\n")

sp %>% group_by(oligo_role) %>%
  summarise(designed = n(), clean = sum(off == 0, na.rm = TRUE),
            median_transcripts = median(n, na.rm = TRUE), .groups = "drop") %>%
  print()

p_spec <- ggplot(sp %>% filter(n > 0), aes(n, oligo_role, colour = off > 0)) +
  geom_vline(xintercept = 1, linetype = "22") +
  geom_jitter(height = 0.2, size = 1, alpha = 0.7) +
  scale_x_log10() +
  scale_colour_manual(values = c(`FALSE` = "#0b6fa4", `TRUE` = "#c2531a"),
                      labels = c("own locus only", "binds another locus"), name = NULL) +
  labs(x = "transcripts bound (log)", y = NULL) + theme_minimal(base_size = 9)
print(p_spec)

# ALTERNATIVE CUT: drop the 3'-anchored partial sites and keep only full matches
#   sp %>% filter(sites_full_match > 0)
# ALTERNATIVE CUT: how much did offering 5 alternatives actually buy?
#   sp %>% group_by(isoform_id, oligo_role) %>%
#     summarise(best = min(off), worst = max(off), .groups="drop") %>%
#     filter(worst > best) %>% arrange(desc(worst - best))

# ---------------------------------------------------------------------------
# 4. Which pairs are genuinely discriminating?
#    A qPCR pair must amplify the candidate, must NOT amplify the partner, and
#    must stay inside the locus. Those are three separate failures.
# ---------------------------------------------------------------------------
prep %>%
  mutate(class = case_when(passes == "True" ~ "usable",
                           amplifies_partner == "True" ~ "amplifies the partner too",
                           TRUE ~ "binds outside the locus")) %>%
  count(pair, class) %>% print()

prep %>% filter(grepl("^qPCR", pair)) %>%
  select(gene_name, chosen_alt, passes, offtarget_loci, amplifies_partner,
         ref_record_unspliced) %>% print()

# ---------------------------------------------------------------------------
# 5. Two QC flags worth looking at before ordering anything
# ---------------------------------------------------------------------------
# (a) pairs that share almost no exonic sequence -- easy to design against,
#     hard to believe as two isoforms of one gene
regions %>% filter(low_overlap_flag == "TRUE") %>%
  select(gene_name, candidate_len, partner_len, shared_bp,
         shared_frac_of_candidate, primary_event) %>%
  arrange(shared_bp) %>% print()

# (b) candidate models under 300 nt -- likely assembly fragments
regions %>% filter(short_model_flag == "TRUE") %>%
  select(gene_name, candidate_len, partner_len, shared_bp) %>% print()

# (c) the reference FASTA defect: entries that were unspliced genomic intervals
seqver %>% filter(pipeline_fasta_was_unspliced == "True") %>%
  select(isoform_id, n_exons, gff_exon_sum, pipeline_fasta_len, final_len) %>%
  print()

# ---------------------------------------------------------------------------
# 6. Effect size against assay class -- the figure in the paper
# ---------------------------------------------------------------------------
p_assign <- assign %>%
  mutate(gene_name = reorder(gene_name, abs(own_dIF)),
         assay = sub(" \\(\\+ 5' RACE to fix the TSS\\)", "", primary_assay)) %>%
  ggplot(aes(abs(own_dIF), gene_name, colour = assay,
             shape = is_nested == "TRUE")) +
  geom_segment(aes(x = 0, xend = abs(own_dIF), yend = gene_name),
               colour = "grey80", linewidth = 0.3) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16),
                     labels = c("resolvable alone", "nested in partner"), name = NULL) +
  labs(x = "|dIF| in the discovering dataset", y = NULL) +
  theme_minimal(base_size = 9)
print(p_assign)

# ALTERNATIVE CUT: only the candidates that need PacBio
#   plat %>% filter(illumina_resolvable == "FALSE") %>% arrange(desc(abs(own_dIF)))
# ALTERNATIVE CUT: only the shared (H and T) candidates
#   assign %>% filter(set == "shared") %>% arrange(desc(abs(own_dIF)))

# ---------------------------------------------------------------------------
# 7. Pull the order sheet for one gene
# ---------------------------------------------------------------------------
show_gene <- function(g) {
  order_s %>% filter(gene_name == g) %>%
    select(oligo_role, sequence, length, tm_C, gc_pct, offtarget_loci,
           expected_product_candidate_bp, expected_product_partner_bp, selection)
}
print(show_gene("NCOA7"))
# show_gene("KDM6A"); show_gene("MERTK"); show_gene("IRAK3")

cat("\nDone. Nothing was written.\n")
