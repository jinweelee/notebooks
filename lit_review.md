1) clip mentions of esi in papers

2) learn about those hpRNA loci:
  
Ok it seems that in order of apparent significance
CR46342 > CR18854 > CR32205 > CR32207 aka 
esi-2/CG4068-B > CG18854/esi-1 > esi-4.1 > CG32207

Ranking in terms of avg norm_counts for PD rescue:
2k > 300 > 9~ > 7~  

So how are these covered in the papers so far: 

1) Czech
Analysis of the most abundant siRNA from esi-2 in flies mutant for Dcr-2, AGO2, r2d2 or
loqs extended our findings from cell culture (Supplementary Fig. 10). To examine the
unexpected requirement for loqs more broadly, we sequenced small RNAs from loqs-mutant
ovaries and observed a near complete loss of endo-siRNAs from structured loci

2) GSE17171 
we depleted all Loqs isoforms
in S2 cells using dsRNA targeting shared 59-UTR sequences
The re-expression of Loqs-PD restored
levels of esi-2.1 in Loqs-depleted cells

3) GSE26230
The only thing they examine is CG4068B/esi-2


4) Marques 
The only thing they examine is CG4068B/esi-2



Ok so good news: we at least have CR18854.
Bad news: Its all we have. 




Now now, not too fast. First, lets think about whats in the Kawamura index and why it soaks all the TE reads

##Kawamura 2008
To test whether AGO2 exists in a complex
with endogenous siRNAs, we immunopurified AGO2 with a specific
monoclonal antibody from a cultured Drosophila somatic S2 cell line
(Supplementary Fig. 1a, b) and examined its associated RNAs. AGO2
in S2 cells was predominantly associated with small RNAs of about 21
nucleotides in length (Supplementary Fig. 1c).

Among the 77,327 reads, 52,348 sequences
matched the Drosophila genome with 100% identity over their entire
length (Fig. 1). In addition, a large number of the AGO2-associated
small RNAs showed single mismatches


Database searching revealed that a large number of the AGO2-
associated small RNAs corresponded to transposons and other repetitive
elements in the genome **around 60~% of reads**


##GSE26230

General literature:
exo-siRNA precursors areprocessed by Dicer-2 (Dcr-2), then loaded by a complex of
Dcr-2 and R2D2 into Ago2 (17–19). Endo-siRNA biogenesis
depends on Dcr-2 paired with a different isoform
of Loqs, Loqs-PD (20–22).
 Loqs-PD and

R2D2 appear as functional antagonists during both
endo- and exo-siRNA mediated silencing in S2-cells,
arguing that they compete for a common factor. Both
Loqs-PD and R2D2 contribute to silencing by invertedrepeat
transgenes, but a quantitative analysis suggests
that neither protein is absolutely required; instead, their
effects appear to be additive suggesting parallel pathways
rather than sequential action

In these cells, the
efficiency of endo-siRNA silencing can be measured via GFP-fluorescence. 
This is because they transfect cells with GFP construct which is itself a target of endo-siRNA silencing
So a higher flouresence level indicates lower siRNA activity, while lower = stronger 

Right so they found that a loqsPD-KO weakened siRNA silencing but when coupled with a r2d2-KO actually brought silencing closer to normal strength i.e an additional depletion of R2D2
partially restored silencing of the reporter only when
the Loqs-PD levels had been concomitantly lowered.

Interesting.. **a clear impairment of RNAi following depletion of
Dcr-2, Ago2 and R2D2, whereas depletion of Loqs-PD
increased RNAi efficiency**

**Simultaneous knock down of R2D2
did not reduce the levels of mature endo-siRNA further
in all combinations; if anything, endo-siRNA biogenesis
improved upon knock down of R2D2**

hmnn They found that there was not much variation in reads mapping to Transposon-targeting endo-siRNAs 

**In
summary, we could validate the previous finding that the
nucleolytic processing of both hairpin-derived and
transposon-targeting endo-siRNAs does not depend on
r2d2.**

**In addition, it appears that Loqs-PD may not be
essential for the dicing of transposon-derived
endo-siRNAs either**

Hmnn seems like my clear cut-distinction hypothesis might not be too feasible **Okamura and colleagues
present a thorough description about how the r2d2
mutation affects processing of the CG4068 hairpin**

**r2d2 mutant flies do not show a gross reduction of
transposon-matching endo-siRNAs relative to heterozygous
animals, consistent with the results obtained from our own libraries**


Now heres the big question, did they do their homework and map to the kawamura index or did they just see which mapped to transposons in general? And hmnn despite them claiming roughly comparable amounts, is that really a solid claim? I could say that the loqs-KO library actually had a higher percentage of reads mapping to the transposons. 

Based on this result, we
favor the hypothesis that 1) Loqs-PD and R2D2 interact
with the same site on Dcr-2 and that 2) ratio of Loqs-PD to other dsRBD
proteins is critical for optimal endo-siRNA biogenesis, 3)Loqs-PD and R2D2 function independently during RNAi.


Hmnn ok initial impressions of renyi's external library mapping.

Firstly its weird that in 17, a r2d2KO slighltly decreases hpRNA levels (ok reasonable), but in 23, r2d2-KO HAD almost a complete silencing? Mis-handling of data? 

Next, it seems that for GSE17171, we see a mirrored trend in almost all the features, that a dcr2/loqs-ORF KO leads to a drop, while there is almost no discernable change in the r2d2 KO. This suggests that r2d2 dosent really have a role to play, kinda un-validating my hypothesis? 

Or could it be that with no r2d2, loqs-PD is free to do its thing and restore things to normal levels. BUT in OUR libraries, even when loqs-PD is restored, other than in hpRNA, its rescue effect is generally muted because of the presence of abundant r2d2. 


Thirdly, GSE37443 supports the competetion hypothesis / aka too many cooks spoil the soup. And also simultanously disproves my earlier hypothesis because it seems that a PD rescue does affect the counts of other loci. 

What is strange, and very strange though is how the `wt` library maps almost entirely to the miRNA. And in fact, i dare say that the `Fly Ovary` library does seem promising? No not really? this one might be the most complex one to handle?

its funny because the FlyHead libraries show a good `competetion` trend, while the FlyOvary libraries show a good hpRNA trend. 



##GSE17171

###Background info: 

The canonical Dcr-2 partner R2D2 seems not to be
required for the production of siRNAs. Instead, it was
found to impact the loading of siRNA duplexes into the
RNA-induced silencing complex (RISC) and proper guide
strand selection (Liu et al. 2003; Tomari et al. 2004).



###Summary of findings: 

We show that coordinated depletion of all Loqs isoforms in cultured cells affects the biogenesis of both miRNAs and endo-siRNAs,
whereas cells **singly depleted of Loqs-PB or Loqs-PD show
an impact only on the miRNA or on the endo-siRNA
pathway, respectively**. While the **re-expression of Loqs-PD
restored endo-siRNA levels** in cultured cells that had been
depleted of all Loqs isoforms, Loqs-PD was incapable of
rescuing miRNA processing defects. Moreover, we show
that **Loqs-PD preferentially interacts with Dcr-2**, the
enzyme responsible for the processing of all endo-siRNA
species.

###Isoform specific effects on biogenesis:

Investigated the impact of isoform-specific knockdowns
on the biogenesis of **a prevalent endo-siRNA (Fig.
1D, esi-2.1) and a miRNA (Fig. 1D, miR-bantam)**.


We observed miRNA processing defects
in all cells treated with dsRNAs that cotargeted Loqs-
PB, whereas the biogenesis of
endo-siRNA was affected only if cells were treated with
dsRNAs that target Loqs-PD, either singly or together with
other isoforms.

Upon depleting all loqs isoforms, The
expression of Loqs-PA, Loqs-PC, or R2D2 failed to rescue
any observed biogenesis defect (this i assume with reference to esi2 and bantam). **However, The re-expression of Loqs-PD restored
levels of esi-2.1 in Loqs-depleted cells but
was incapable of rescuing miRNA-processing defects**.

SUPER IMPT ATTEMPT AT R2D2 RESCUE:
**Our observation that the
expression of R2D2 was incapable of rescuing any small
RNA-processing defect caused by depletion of Loqs
strongly suggests that R2D2 and all Loqs isoforms cannot
function in a redundant manner.**

**Considered together,
these data indicate that only the Loqs-PB isoform is
required for the biogenesis of miRNAs and suggest that
only Loqs-PD is essential for endo-siRNA production.**

Loqs-PA? The expression of Loqs-PA in loqs f00791
homozygous mutant flies was incapable of restoring either
miRNA or endo-siRNA processing. This is in line with cell-based studies showing that neither Loqs-PA nor
PB is required for the endo-siRNA pathway.The
function of Loqs-PA, remains elusive.


###Affinity between loqs-isoforms and Dcr through IP experiments (note dimer formation)
Immunoprecipitation of Flag-tagged Dcr-2 followed
by immunoblotting with an antibody specific to the
N-terminus of endogenous Loqs, and thus recognizing all
known Loqs isoforms, revealed a strong signal corresponding
to the molecular weight of Loqs-PC/PD

Next, we examined
the interactions between various Loqs isoforms and Dcr-2. **All four tagged Loqs isoforms and the positive control T7-
R2D2 were able to interact with Dcr-2**

WAIT WHAT, HOMO/HETERODIMERS?
**Interestingly, we also found that Loqs isoforms
are capable of forming homo- and heterodimers (or even
oligomers) with each other, as Flag-tagged Loqs isoforms
PA, PB, and PC are capable of pulling down T7-tagged
Loqs-PA in coimmunoprecipitation assays**

Loqs-PB and to a lesser
extent Loqs-PA, are capable of interacting with Dcr-1, a
result consistent with findings reported by Fo¨rstemann
et al. (2005). In contrast, **Dcr-2 predominantly interacts
with endogenous Loqs-PD, while other isoforms interact
with Dcr-2 to a significant extent only when overexpressed.**


###Sequencing experiments involving loqs-depletion

small RNA
libraries from S2 cells treated with loqs-ORF dsRNAs and
from control knockdown cells, including lacZ, dcr-1, dcr-2,
r2d2, and untreated cells (‘‘mock’’).

The indicated small RNA categories were isolated
from the total library bioinformatically (Czech et al. 2008). **Wait, he uses Czech's method of profiling?**

Endo-siRNAs mapping to
overlapping transcripts (exonic antisense) were strongly
reduced in dcr-2 and loqs-ORF knockdowns

Finally,
knockdowns of dcr-2 or loqs caused a substantial (z50%)
reduction in endo-siRNAs corresponding to repeat and
transposon sequences, while the levels of these small RNAs
remained unchanged in all other knockdowns tested



###Validating sequencing results with northern blotting
We saw decreased levels of three
independent endo-siRNAs derived from structured loci
(esi-1.2, esi-2.1, esi-4.1) in dcr-2 and loqs-ORF knockdowns.

We also noted a reduction of klarsicht
siRNAs and sequences derived from the transposon
DM1731 upon loqs-ORF knockdown, but not upon Loqs-
BC depletion

###Effect of loqs-ORF depletion on transposon silencing
Endo-siRNAs that
correspond to repeats and that were 21 nt in length were
extracted bioinformatically

Knockdown of dcr-2 caused a reduction of endosiRNA
sequences for the majority of transposable elements.
**r2d2 knockdowns, showed only minor, if any, effects. And of course, depletion of all
Loqs isoforms reduced levels of repeat endo-siRNAs**.


OK THIS IS INTERESTING **Interestingly,
we noted that knockdown of r2d2 showed a weak
impact on those repetitive elements that did not depend on
Loqs depletion (2 Loci)**


###Final thoughts
Ok crucially they did not do a isoform specific rescue on a loqs-ORF to pinpoint whether or not Loqs-PD is responsible for transposable element silencing. 



##GSE37443

###Background info:
R2D2 and free phosphate
suppress the inherent ability of Dcr-2 to process premiRNAs
into 21 nt duplexes, restricting it to the longer dsRNA
triggers associated with RNAi.

Production of esiRNAs (hpRNAs) by Dcr-2 in cultured Drosophila S2 cells
requires the alternative partner protein, Loqs-PD, rather than
R2D2. **However, it is not known whether the production of exo-siRNAs,
cis-NAT-endo-siRNAs, and mobile element-derived endosiRNAs
also requires Loqs-PD. Moreover, an in vivo role for
Loqs-PD in small RNA production has not been established.**

Interesting stuff on the different structure/ domains of loqs-isoforms:  
The largest Loqs protein isoform, **Loqs-PB, comprises three
double-stranded RNA-binding domains (dsRBDs); dsRBD3 is
required for Loqs-PB to bind to Dcr-1** (Fo¨ rstemann et al., 2005;
Ye et al., 2007). **Loqs-PA, which also binds to Dcr-1, lacks
a part of the third dsRBD of Loqs-PB**. Loqs-PD has a unique carboxy
terminal sequence in place of dsRBD3 and binds to Dcr-2 rather than Dcr-1.

Interesting stuff on the general theme of sRNA processing:  
Partnerships between dicer-like proteins anddsRBD
proteins are a general theme in RNA silencing pathways. For example, plants produce four distinct
Dicer enzymes, each with its own specialized
dsRBD partner



###Summary of findings
We find
that Loqs-PA and Loqs-PB, but not
Loqs-PD, can rescue the lethality of
loqs null mutant flies.

Loqs-PB is crucial
for female fertility and for producing
specific subsets of miRNAs. **Loqs-PD
enhances the production of both endo and
exo-siRNAs.**


In the absence of
Loqs-PB, Dcr-1 produces aberrant products
from pre-miR-307a, pre-miR-87,
and pre-miR-316.


##Set up and background info

loqsKO mutant flies are embryonic lethal (Park et al.,
2007). Loqs-PA or Loqs-PB, but not Loqs-PD, restored embryonic
viability, suggesting that a defect in miRNA biogenesis
underlies the lethality observed in loqsKO embryos.

Oh wow its an entirely in vivo study 


### Normal Exo-siRNA Accumulation Requires Loqs-PD In Vivo

In cultured Drosophila S2 cells, Loqs-PD is required for efficient
esiRNA production, but does does Loqs-PD play a role in RNAi in adult flies? 

**Okay interesting, the GMR-wIR transgene is an exo-siRNA reporter**

The GMR-wIR transgene produces during eye development
an inverted repeat hairpin RNA corresponding to white exon 3.
Dcr-2 processes the wIR hairpin into siRNAs (Lee and Carthew,
2003; Vagin et al., 2006), which in turn silence **white expression,
causing the eye to be white or orange instead of red**. As such, **reduction in the absorbance at 480 nm of pigment extracted from the eye provides a quantitative
measure of white silencing.**


They found that in WT and HT-KO for dcr-2, r2d2, or
ago2, the wIR transgene produced white eyes. In contrast, the
eyes from wIR transgenic female flies **homozygous for dcr-2,
r2d2, or ago2 mutations (A480 = 0.9–1.0) were similar to those
of wild-type** (mean A480 = 1.2).

AH ok, i looked at the diagram wrongly, so it turns out it all makes sense. HT show silencing levels similar to WT, while HM show almost no silencing. 

Relation to Loqs? HT loqs-ORF KO no diff, but HM loqs-ORF KO shows slightly lower levels of silencing (but silencing not entirely attenuated). **Then, the rescue of PA,PB/PA+PB seems to decrease exo-siRNA processing? But adding PD to this combination restores siRNA processing**

Loss of Loqs-PD decreased
the abundance of white siRNA reads across the entire wIR
hairpin sequence for both sense and antisense siRNAs, suggesting
a sequence-independent role for Loqs-PD in facilitating exosiRNA
production by Dcr-2 (


### Oh shit Normal Exo-siRNA Accumulation Requires Loqs-PD In Vivo

Cis-NAT endo-siRNA accumulation also required Loqs-PD in
heads and ovaries (Figures 2D and S2D). Our data were more
equivocal for transposon-derived endo-siRNAs: these
decreased when Loqs-PB was the only isoform present, but
were unaltered compared to loqsKO/CyO heterozygotes when
Loqs-PA was expressed (Figure 2D). In contrast, overexpression
of Loqs-PD increased both cis-NAT and transposon-derived endo-
siRNAs, suggesting that Loqs-PD acts in the production of
transposon-derived endo-siRNAs as it does for exo-siRNA,
esiRNA, and cis-NAT endo-siRNA biogenesis.

In contrast, overexpression
of Loqs-PD increased both cis-NAT and transposon-derived endo-
siRNAs, suggesting that Loqs-PD acts in the production of
transposon-derived endo-siRNAs as it does for exo-siRNA,
esiRNA, and cis-NAT endo-siRNA biogenesis.

Key points / statements they make that we can challenge: 
**suggesting that Loqs-PD acts in the production of
transposon-derived endo-siRNAs as it does for exo-siRNA,
esiRNA, and cis-NAT endo-siRNA bioge
nesis.** 

**Endo-siRNA Accumulation *Requires* Loqs-PD** 


We conclude
that Loqs-PD functions in both endo-siRNA and exo-siRNA
production by decreasing the concentration of substrate
required for Dcr-2 to efficiently produce siRNA.


###Loqs-PB Tunes Where Dcr-1 Cleaves Some Pre-miRNAs
In the **absence of Loqs-PB, miR-307a was shorter and miR-316 and
miR-9b were longer** (Figure 3A–3D). miR-307a was also shorter
and miR-316 was longer in loqsf00791 hypomorphic mutants.

flies lacking Loqs-PB accumulate a 21 nt long miR-307a isomir whose seed sequence
begins at position 4 of the canonical 23 nt long miR-307a isomir.
The different lengths of these miRNAs in loqsKO;Loqs-PA/TM3
and loqsKO;Loqs-PB/TM3 flies likely reflects a direct influence of
Loqs-PB on the choice of cleavage site by Dcr-1.

Our data suggest that Loqs-PB repositions
pre-miR-307a on Dcr-1 so as to favor production of
miR-307a23-mer over miR-307a21-mer

Also note that in vitro: 
We conclude that Dcr-1 and Loqs-PB
form a 1:1 complex.


##GSE45290


To differentiate between Ago2-loaded and other small RNAs,
we made use of the fact that the 3’-terminal nucleotide of Piwi-/
Aub-/Ago3- as well as Ago2-loaded small RNAs is 2’-O-methyl
modified [16,17]. This modification renders the small RNAs
resistant to oxidation of vicinal diols with sodium periodate and
subsequent b-elimination that will shorten the un-modified RNAs
by one nucleotide and prevent them from participating in the 3’-
end ligation reaction required for deep sequencing library
generation. The technique is highly efficient and specific since
the b-elimination resistant small RNAs essentially disappear in
libraries prepared from ago2 null mutant flies filtered out.


A striking observation was that a large proportion ofreads (0.6% to 14.5% of genome matching reads, 5.1% to 80.6%of transposons matching reads) could be attributed to only fourtransposable elements (roo,297,TNFBandblood). 

Loss of Loqs-PD resulted in a 1.8-fold reduction of transposonmatching
endo-siRNAs in libraries without b-elimination, consistent
with the notion that its role is predominantly in siRNA
production (Figure 1 C, left panel). While this was true for the
analysis of all transposons in bulk, some individual exceptions to
this rule exist. For example, the transposons F-element, 412 and Doc
were only slightly affected by loss of Loqs-PD (Figure S8 in File
S1). The overall reduction of transposon-targeting endo-siRNAs
was also observed after b-elimination, in agreement with the
notion that small RNAs must be produced before they can be

We observed no overall reduction of transposon
matching endo-siRNAs between heterozygous and homozygous
mutants before b-elimination (Figure 1 C, left panel) with only one
major exception: The endo-siRNAs directed against F-element were
strongly reduced in the absence of r2d2 but only mildly affected by
the absence of loqs-D (Figure S8 in File S1).



Is there any common principle that could explain why certain transposons differ from the bulk in their requirements for Loqs-PDand R2D2? This distinction is not based on their abundance sincepreference of Loqs-PD for biogenesis or R2D2 for Ago2-loadingdoes not correlate with the absolute amount of small RNAs (Figure3). Furthermore, when transposons were classified into longterminal repeats (LTRs), long interspersed elements (LINEs) andinverted  repeats  (IRs),  we  did  not  observe  any  consistent correlation that could explain R2D2 versus Loqs-PD preference.