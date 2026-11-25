# The UBL5 knockdown line, measured from the RNA-seq itself

*Attached: `ubl5_knockdown_efficacy.png`. Underlying numbers in
`ubl5_sart1_knockdown_stats.csv`, `ubl5_sart1_cpm_libraries.csv`,
`ubl5_isoform_site_map.csv`.*

Before building candidate lists off the UBL5 arm, I went back and measured UBL5 in the
sequencing data directly, rather than relying on the original qPCR validation. Two things
came out of it that I think you'll want to know about, and one question only you can answer.

## The knockdown is about 30%

Using library-size-normalised counts across all twelve THP-1 libraries (n = 3 per arm):

| condition | wildtype | knockdown | ratio | knockdown | p |
|---|---|---|---|---|---|
| unstimulated | 49.8 CPM | 33.6 CPM | 0.671 (CI 0.48–0.93) | 33% | 0.029 |
| LPS | 43.0 CPM | 29.9 CPM | 0.696 (CI 0.65–0.75) | 30% | 0.0019 |

Ratios, intervals and p-values throughout are from a Welch t-test on log2 CPM, so that the
point estimate and its interval are on the same scale.

So the direction and significance of the original qPCR reproduce cleanly — but the magnitude
is around 30%, leaving roughly 70% of UBL5 mRNA in place. That is well below what "shRNA
knockdown" usually implies, and it matters for how we read everything downstream: with that
much residual message, a negative result in the knockdown is weak evidence that UBL5 isn't
involved, rather than good evidence that it isn't.

I want to be careful about one thing: **this is mRNA, and mRNA only.** miR-30-based shRNAs
loaded into RISC can repress translation as well as degrade transcript, so protein-level
knockdown could be considerably deeper than 30%. I have no way to settle that from
sequencing. If a western was ever run on this line, that single piece of data would change
the interpretation more than anything else I can compute.

## The shRNA is on-target, and no isoform escapes it

I checked whether the modest knockdown could be explained by some UBL5 isoform lacking the
target site. It can't.

Searching the sequence recorded in the thesis — `CGGGATGAACCTGGAGCTTTAT` — against all
315,349 transcripts in our long-read annotation returned **zero** matches. Dropping the
leading `CGG` and searching the remaining 19-mer `GATGAACCTGGAGCTTTAT` returns exactly 14.
So those three bases look like a cloning overhang rather than target-complementary sequence;
worth correcting in the records, because searching the 22-mer as written looks like a
specificity failure when it isn't.

Of those 14 hits: 13 are UBL5 transcript models, and the fourteenth is the pseudogene UBL5P2.
Nothing else in the transcriptome contains the site, so the reagent is clean. UBL5P2 sits
below 1 TPM and the other three pseudogenes are flat zero, so that cross-match has no
practical consequence.

One UBL5 model does lack the site — but it carries about 0.2% of the gene's output
(0.2–0.6 TPM against the dominant model's ~200 TPM), so it cannot be hiding the residual
expression. The knockdown is uniform across every isoform that is actually expressed. The
ceiling is the reagent's efficacy, not an escape route.

My best guess, and it is a guess: LTR-driven shRNA cassettes commonly silence over
continued passage. If the qPCR validation was done shortly after selection and the RNA for
sequencing came from a materially later passage, partial silencing would produce exactly
this picture.

## Why SART1 is the gene to watch

SART1 is the human homolog of yeast Snu66. It carries the HIND motifs that Hub1 — UBL5's
yeast counterpart — physically binds, and it sits in the U4/U6.U5 tri-snRNP. It is, as far
as anyone knows, the adaptor through which UBL5 reaches the spliceosome at all. That makes
its level relevant in three separate ways:

1. **Buffering.** More adaptor can partly offset less UBL5. If SART1 rises in the knockdown,
   the functional perturbation is smaller than the 30% figure already suggests.
2. **Specificity.** If SART1 itself moved substantially, any splicing phenotype we see could
   not be attributed to UBL5 cleanly — we'd have perturbed the partner too.
3. **Complex stability.** Direction matters. A fall would raise the worry that the
   tri-snRNP is destabilised; a rise is the more benign reading.

Measured: SART1 is up about 10–14% in the knockdown, in the same direction in both
conditions — but the baseline comparison is not significant (p = 0.105) and the LPS one only
marginally so (p = 0.049, CI lower bound 1.001, sitting right on the boundary). With n = 3
and no correction for multiple testing, **I would not claim SART1 is changed.** The
consistent direction is suggestive enough to be worth a targeted measurement, and it is the
one partner whose level would change how we read a UBL5 splicing effect, which is why I'm
flagging a non-significant result at all rather than dropping it.

## What would help most, in order

1. Was UBL5 protein ever blotted in this line? This is the single highest-value missing
   number.
2. Is the sh-Luciferase control line still available? The sequenced comparison was against
   untransduced wildtype, so for validation work we need the transduced control to separate
   UBL5 loss from transduction and selection.
3. Roughly what passage was the line at when RNA was harvested, relative to the original
   qPCR validation? This tests the cassette-silencing idea directly.
4. Was knockdown re-checked at the time of harvest, or only after selection?

If the protein knockdown turns out to be deep, the arm is a solid hypomorph and everything
downstream stands. If protein tracks mRNA at ~30%, we should describe this as a UBL5
hypomorph rather than a loss-of-function throughout, and treat negative results as
uninformative rather than as evidence of no requirement.
