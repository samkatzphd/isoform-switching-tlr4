# Deep dive: the novel *ATF5* isoform, macrophages, and what it might be

*Attached: `atf5_biology.png`. Tables: `atf5_isoform_expression_all_libraries.csv`,
`atf5_lps_induction.csv`, `atf5_oligo_order_table.csv`. Literature retrieved from PubMed;
DOIs given inline.*

---

## 1. Yes — it is expressed in wildtype, and that changes the case entirely

The short answer to "do we see this isoform in wildtype at all" is **yes, in every single
library**, and the numbers are better than anything the interaction analysis produced.

Novel isoform **TCONS_00159376**, absolute abundance (TPM):

| system | unstimulated | + LPS | fold | p |
|---|---|---|---|---|
| THP-1 wildtype | 26.8 | 44.5 | **1.66×** | **0.004** |
| primary MDM (3 donors) | 12.1 | 34.6 | **2.87×** | 0.06 |
| THP-1 UBL5 knockdown | 23.8 | 23.9 | **1.00×** | 0.96 |

Detected in all 18 libraries, all 3 primary donors, both conditions. It is not a knockdown
artefact, not a THP-1 artefact, and not a low-abundance curiosity — at 27–45 TPM it is a
well-expressed transcript, and in primary macrophages it carries **20–25% of all *ATF5*
output**, a higher share than in THP-1 (11–14%).

### The cleanest statement of the UBL5 effect — and it is not the one we started with

The interaction analysis framed this as a change in isoform *fraction*. The absolute numbers
tell a sharper story:

- In wildtype, LPS induces the novel isoform **1.66-fold** (p = 0.004).
- In the knockdown, LPS induces it **not at all** (1.00-fold, p = 0.96).
- At baseline the two genotypes are indistinguishable (23.8 vs 26.8 TPM, p = 0.26).
- After LPS they differ two-fold (23.9 vs 44.5 TPM, p = 0.002).

And the specificity is striking: in the knockdown, LPS induces the ***ATF5* gene as a whole
harder** than in wildtype (1.74× vs 1.38×) and the dominant isoform harder still (1.88× vs
1.33×) — while this one isoform stays flat. **UBL5 is not required for *ATF5* expression; it
is required for the LPS-induced production of this particular isoform.**

That is a specific, splicing-level claim, and it rests on marginal contrasts with
p = 0.002–0.004 rather than on the underpowered interaction (q = 0.020, no genome-wide
excess over null). Everything below should be read as following from the absolute-abundance
result, not from the interaction.

---

## 2. What the transcript actually is

Two structural features, both verified from sequence rather than inferred.

### It ends at a genuine intronic poly(A) site

The novel isoform terminates at chr19:49,931,121 — inside what is intron 2 of the dominant
isoform — rather than splicing on to the large terminal exon. Searching its 3′ end finds a
canonical **AATAAA at 15 nt from the transcript end**, which is textbook poly(A)-signal
placement. The dominant isoform carries its own AATAAA 23 nt from its end, as a positive
control for the method.

This makes the transcript an **intronic polyadenylation (IPA) isoform** — a cryptic poly(A)
site inside an intron, used to terminate transcription early and produce a C-terminally
truncated protein. IPA is a recognised, regulated mechanism rather than a splicing error, and
it is prominent in activated immune cells. The AATAAA is strong evidence that this is a real
3′ end and not an assembly artefact or a degradation fragment.

### It would encode an 88-aa protein lacking the entire bZIP

Translating from the canonical ATG gives **88 aa terminating in a real stop codon**, against
282 aa for full-length ATF5 (my translation reproduces the Ensembl length exactly, which
validates the frame). The first 62 residues are shared.

Critically, the stop sits in the final exon, 383 nt downstream of the last exon–exon
junction — so this is **not** predicted to be an NMD substrate. If made, it persists.

### It uses the uORF-rich 5′ leader

*ATF5* is a translationally controlled gene of the ATF4 type. <cite>The ATF5 gene transcribes
mRNAs with two alternative 5′-UTRs, 5′-UTRα and 5′-UTRβ, derived from exon 1α and exon 1β, and
5′-UTRα contains the upstream open reading frames uORF1 and uORF2</cite> (PMID 23876217,
doi:10.1111/febs.12440). <cite>ATF5 translational control involves two uORFs in the 5′-leader:
uORF1 is positive-acting, allowing scanning ribosomes to reinitiate, while uORF2 is inhibitory;
phosphorylation of eIF2 during stress delays reinitiation and allows ribosomes to bypass
uORF2</cite> (PMID 18195013, doi:10.1074/jbc.M708530200).

Our two isoforms differ exactly here. The novel isoform has a **457-nt leader with 4 upstream
AUGs** (three terminating within the leader) against the dominant isoform's **315-nt leader
with 2**. Both retain an overlapping uORF that terminates downstream of the main start codon —
the configuration that mediates eIF2α-dependent control — and the novel isoform adds three
more upstream.

So the novel isoform differs from the dominant one at *both* ends, and both differences land
on documented regulatory machinery: a different uORF architecture at the 5′ end, and an
intronic poly(A) site at the 3′ end.

---

## 3. ATF5 in macrophages and inflammatory signalling

This is the part that makes the candidate worth bench time, and it is better established than
I expected.

**ATF5 restrains inflammatory cytokine output in macrophage-lineage cells.** In BV-2 microglia,
<cite>LPS treatment increased pro-inflammatory cytokine secretion alongside decreased oxygen
consumption and mitochondrial membrane potential and up-regulation of the UPRmt; knockdown of
ATF5 by siRNA increased production of IL-6, IL-1β and TNF-α and decreased membrane potential,
indicating that ATF5-dependent induction of the UPRmt acts as a protective mechanism during
neuroinflammation</cite> (PMID 36834738, doi:10.3390/ijms24043322).

Read that carefully, because it sets the functional prediction: **less ATF5 function → more
IL-6, IL-1β and TNF-α.** ATF5 is a brake on inflammatory output, operating through the
mitochondrial unfolded protein response.

**ATF5 is also a node in myeloid differentiation.** <cite>mTORC2 promotes monocyte development
by downregulating the AKT–Foxo4–ATF5–mtUPR pathway, and overexpression of ATF5 rescued the
monocyte/macrophage differentiation defects of Rictor-deficient myeloid progenitors</cite>
(PMID 37649699, doi:10.1016/j.isci.2023.107540). Relevant here because THP-1 are monocytic and
the protocol uses PMA differentiation.

**And ATF5 is the mammalian UPRmt regulator.** <cite>ATF5 is regulated similarly to C. elegans
ATFS-1 and induces a similar transcriptional response; ATF5 expression rescues UPRmt signalling
in atfs-1-deficient worms, and mammalian cells require ATF5 to maintain mitochondrial activity
during stress</cite> (PMID 27426517, doi:10.1016/j.cub.2016.06.002).

The convergence with UBL5 is the reason this candidate was taken seriously in the first place:
ATFS-1, DVE-1 and UBL-5 are the three nuclear components of the same pathway in *C. elegans*,
and ATF5 and UBL5 are the mammalian counterparts of two of them.

---

## 4. What UBL5 would be doing, mechanistically

The yeast literature makes a specific prediction that fits.

<cite>Spliceosomes that lack Hub1, or are defective in Hub1–HIND interaction, cannot use certain
non-canonical 5′ splice sites and are defective in alternative SRC1 splicing</cite>
(PMID 21614000, doi:10.1038/nature10143). <cite>Hub1 binds the DEAD-box helicase Prp5 and
stimulates its ATPase activity, enhancing splicing and relaxing fidelity; high Hub1 levels
enhance splicing efficiency but also cause missplicing by tolerating suboptimal splice sites
and branchpoint sequences</cite> (PMID 28712727, doi:10.1016/j.molcel.2017.06.021).

So UBL5/Hub1 is a **permissiveness factor** — it licenses the use of suboptimal sites. An
isoform that depends on terminating at a cryptic intronic poly(A) site instead of splicing to
the canonical downstream exon is precisely the kind of marginal, non-canonical outcome that
would require permissive conditions. Lose UBL5, lose the permissiveness, and the cell defaults
to the canonical splice — which is what we see: the novel isoform fails to be induced while the
canonical isoform is induced harder than normal.

**One honest tension.** <cite>UBL5 depletion decreases pre-mRNA splicing efficiency, leading to
globally enhanced intron retention</cite> (PMID 25092792, doi:10.15252/embr.201438679). Reduced
splicing efficiency would, naively, *favour* intronic polyadenylation — IPA competes
kinetically with splicing of the downstream intron. We observe the opposite. Either the effect
here is not kinetic competition, or this locus behaves against the global trend, or the ~30%
knockdown is too mild to produce the global phenotype. This is the clearest thing the bench
work could discriminate, and it should not be smoothed over.

---

## 5. What the transcript might be about — ranked hypotheses

### H1. A translational-control device rather than a protein (most likely)

ATF5's output is governed by uORFs and eIF2α phosphorylation, and LPS/TLR4 engages the
integrated stress response. A transcript that is LPS-induced, carries extra upstream uORFs,
and terminates early would divert *ATF5* transcriptional output into a form that is poorly
translated under basal conditions and differently translated under stress. On this reading the
isoform is a **rheostat on how much functional ATF5 protein a given amount of *ATF5*
transcription yields** — and the "function" is the diversion itself, not the 88-aa peptide.

This fits best because it uses the regulatory machinery the gene is already known to use, it
explains why the transcript is abundant and stable (no NMD), and it requires no novel protein
activity.

### H2. A truncated protein with a residual non-transcriptional activity

The 88-aa product retains the N-terminal region. In ATFS-1/ATF5 the mitochondrial targeting
sequence is N-terminal while the bZIP is C-terminal, so a short N-terminal fragment could still
be targeted to mitochondria while being transcriptionally inert — a decoy for import, or a
fragment with some other activity.

### H3. A dominant negative — but **not** by the published mechanism

This needs stating precisely, because "dominant-negative ATF5" is an established and
well-characterised reagent and our isoform is **not** one. <cite>Dominant-negative ATF5 forms
were created lacking DNA-binding ability but retaining the ATF5 leucine zipper, so that they
associate with and sequester ATF5's requisite leucine-zipper-binding partners; DN-ATF5
associates with CEBPB, CEBPD and CCDC6 and interferes with CEBPB and CEBPD transcriptional
activity</cite> (PMID 31676720, doi:10.1158/1541-7786.MCR-19-0631; see also PMID 27126996,
doi:10.1158/1078-0432.CCR-15-2827).

The published DN retains the leucine zipper and loses only the basic region. **Our 88-aa form
loses the basic region *and* the leucine zipper**, because it terminates before the exon that
encodes them. It therefore cannot dimerise with CEBPB/CEBPD and cannot work by sequestration.
Anyone reading "truncated ATF5" as "dominant negative" will reach for the wrong mechanism and
the wrong experiment. Rank this third, and if tested, test it as competition for something
other than bZIP partners.

### What each hypothesis predicts at the bench

| | predicts |
|---|---|
| H1 translational device | 88-aa peptide absent or barely detectable; isoform abundance responds to eIF2α manipulation (e.g. ISRIB, thapsigargin); cytokine phenotype only via total ATF5 protein |
| H2 truncated protein | ~10 kDa band on an N-terminal antibody blot; mitochondrial localisation by fractionation |
| H3 dominant negative | ~10 kDa band **and** a cytokine phenotype on overexpression, but no CEBPB/CEBPD co-IP |

**The western blot separates H1 from H2/H3 in one experiment**, which is why it is first in the
protocol.

---

## 6. The functional prediction worth testing

Putting the macrophage literature together with the structure: if the novel isoform reduces
the amount of functional ATF5 protein made per unit of *ATF5* transcription, and ATF5 restrains
IL-6/IL-1β/TNF-α through the UPRmt (PMID 36834738), then **LPS-induced production of this
isoform is a feed-forward step that limits ATF5's own brake on inflammation.**

That would make it a self-limiting switch on the anti-inflammatory arm of the macrophage LPS
response — induced by LPS, reducing functional ATF5, de-restraining cytokine output. It also
predicts that the UBL5 knockdown, which cannot make the isoform on stimulation, should have
*more* functional ATF5 under LPS and therefore a *blunted* cytokine response — which is
consistent with the globally blunted LPS response already documented in that arm, though that
blunting has many possible causes and this is not evidence for the mechanism.

This is the hypothesis the assays are designed to test, and it is a hypothesis, not a finding.

---

## 7. What is solid and what is not

**Solid, from these data:** the isoform exists; it is expressed in every library including all
three primary donors; it is LPS-induced 1.66× in wildtype THP-1 (p = 0.004) and 2.87× in
primary MDM; it is completely flat in the UBL5 knockdown (p = 0.96) while the gene is induced
harder; it carries a canonical poly(A) signal 15 nt from its 3′ end; it would encode 88 aa and
is not an NMD target; its 5′ leader carries 4 uAUGs against the dominant isoform's 2.

**Not established:** that the 88-aa protein is made; that the transcript is capped at the
annotated 5′ end; that the poly(A) site is used in vivo rather than inferred from assembly;
that UBL5 acts directly on this event; any functional consequence whatsoever. The domain
assignments (N-terminal MTS, C-terminal bZIP) are from the literature, not measured here.

**Also unresolved:** the sequencing comparison used untransduced wildtype rather than the
sh-Luciferase control, so the genotype contrast carries transduction and selection. The
LPS-induction contrast within wildtype does not — it is internally controlled — which is
another reason to lead with the induction result.
