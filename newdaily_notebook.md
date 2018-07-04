## 290517

Revised API for os, str and list 

thought about how to implement the other external_library snakefile.

thought about architecture of project folder, decided to create a Init.Snakefile folder to get the necessary architecture out, can possibly add/remove in the future. 

can we rewrite sequential mapping? I think we can. 


## 300517

Reached office late, need to fix that. 

Doublechecked that repbase account was registered, revised some pandas, got kenneth's RNA-seq pipeline.

modified the init.snakefile abit thanks to kenneth.

Going to start re-writing seqeuntial mapping, nope.

Think we defintely need to start from index generation and work from there. 

tmr, need to find out how to ssh to atlas.


## 310517

Started using markdown pad, mobaxterm.

Got Prerequisite files for repeatmasker set up, yet to run them. 

Somewhat finalized architecture for index_gen,

decided to compile all runscripts to `default_config.yaml` extension 


Hmnn realized that the flybase source has alot of the base features but its the full length features.would we want to use those instead of the merge exon we have been using?
Ask greg tomorrow. Hmnn realized there was this all along.


Hmnn esp for miRNA, i think we cant use the fasta file becuase it contains pre miRNA details. 

wa shit, seems like 6 miRNA loci found in rDNA, when rDNA index is mapped for us first. 



##010618

Seems like only jeremy, ahjung and mingchen have rights to install programs on atlas. Need to ask about Repeatmasker ASAP. 

k holy shit, it seems that we can just edit from mobaxterm, we are moving everything server side from now on. 

both `testing grounds` and `layout` effectively moved server side.  

decided to go back to my old method of doing cis nat, seems legit.  

When we get back to work on sunday, start with gffutils



##040618
Note, we have to `git push origin master` after committing and we gotta add files to be tracked before commiting. 


for "Name=hpRNA", note that while `gene` and `ncRNA` features span the same region, `ncRNA` features have repeats because they represent the various hpRNA-isoforms.  

For example, CR33940 has 3 isoforms, D,E and F. 


 
##080618

Ok we seem to having some issues with `bedtools merge -s`, i suspect its a result of having more than 6 fields, look into this. 

Yes thats the issue, so, after gff2bed we have to:
1) Remove additional fields 
2) Add parenttype information to "$4" 


Also, after gff2bed, `eID` represents exon


Note, we have to refine bedtools -intersect parameters for more information.


##110618

well learnt how to deal with multiple whitespace delimited bullshit in pandas. 


`os.chdir()` is life changing 


##120618

Think we should allow the use of flags to link up seemingly separate rules. 

Note, we incorporated the flags already but will only observe them when we re-run the entire thing.

ok, 05_CR14033.bed is buggy, i cant seem to `get fasta` it. Double check coordinates, seems ok, very strange. 

so i think we'll just add the CR14033 fasta to downloads 

Also we've found out that if you remove one file from the output, the entire rule runs to reproduce everything else as well.  
 


##130618

Things we have yet to settle,

got to get the index list sorted out. 



##140618

Hmnnm seems like index0 should have Uxtra excluded? 

Yup seems like i might have made a minor mistake by not excluding Uxtra before creating all the indexes. oh well. 

Hmnn okay, right recall, sequential mapping vs all_index mapping

might want to do 1 additional rule at the end, like a libstats rule:



Ahhh ok, so unmapped files are always in the fasta format, while mapped files are always in the bowtie.txt format. 



what details do we want to include? 

1) Mapped/ Unmapped to Uextra reference genome 
2) Unique/ Multiple mappers    

to do do list for tmr:

1) rerun RM on genome -> DMr6.21excludeUextra.fasta
2) re-configure naming in configfile
3) re-run index generation in layout.
4) re-run sequential mapping/ specifically mapping to index0


<<<<<<< HEAD
Ok, we are also going to scrap the entire dictionary index idea and just position based on a list. we just ensure that the genome index is first, spikein is last. 



##210618
Writing here again after about 1 week hiatus. A combination of forgetulness, laziness and lack of not getting into the habit of commiting daily.

So spent the last few days reconfiguring all Snakefiles to how they actually should be written using wildcards.

I think over the weekend i will reconfigure the directories in a more compact Input/Ouput structure. 

Now moving on to the R phase of my life. 

Also, finding clojure.


##220618

Now that we generate each mapped/ unmapped combo, i think we can rewrite the maplog rule.


ok, over the weekend we are goign to reconfigure architure, 3 impt points: 

1) higher order input, output structure 

2) make it such that its `Okamura_libs`, `GSE17171`.. essentially remove the `External_dataset` folder.

3) make structure more friendly for alternative runs. 
ie, go 1 level deeper, and specify a `run_ID`, take reference from R script if needed, it is most recent structure.

4) change `genome_vers` to just `genome`


Added an index dictionary for downstream R data gathering

Note, once data vizualization set up, i have several runs in mind
1) thickveins only
2) othersiRNA only
3) Each type of transposon only
 
=======
Ok, we are also going to scrap the entire dictionary index idea and just position based on a list. we just ensure that the genome index is first, spikein is last.





##030718
So we are starting a new daily notebook. Because im an idiot and i messed up my git commits. 

K i think we make the DE.snakefile as general as possible? Then we do abit more processing in the R script.

Yup, text processing of loci names will occur in the DE R script.

Note to self, change the name of 'NewCisnat' to 'NewTransnat', maybe when we rerun everything in one shot tomorrow.

nah nah, a better name "ExonicAntisense"

And also check out that Watanabe paper, well shit turns out its a mouse paper

and alright. tmr, we give the R script a full run through with the snakefile and if thats a good to go, we can try a few mapping RUN-IDs tmr. 



##040718

we have to sync up with the index generation snakefile if we want to introduce new stuff
or else, today, we can just use what we have 

K i think we can try a transposon mapping run and a citation mapping run today. 

RM indexes we still keep as a filter: 
             
SatelliteRM, Low_complexityRM, RCRM, Simple_repeatRM, OtherRM, UnknownRM, ARTEFACTRM

We shall call it, `GffTE`

note check extent of ovelap between 4 citation indexes
Also extent of overlap between LTRRM, LINERM and DNARM

`GffTE`: Ok so what have we learned. There is some serious significant overlap between GffTE and RMTransposon, that is quite clearly reflected in the count difference vs `Liling`. Other than that, nothing much significant changed. 

`GffTEmoveRM`: ok, literally nothing changed and not much seems really to be mapped to `otherRM`

`RMTE`: k so this one involves me changing some variables in the index_dict and understanding the overlap between the 3 transposonRM indexes.

About 900+ features between LTR and LINE intersect, which is unsruprising considering how LINEs are a type of LTR. 
And about 100+ and 200+ features between LINE and LTR intersect with DNA respectively. 

So i suppose, the 'lazy' way would be to just go hierchically/ based on the number of features, so 
LTR>LINE>DNARM, which is exactly what Liling did.

Note, removed citation index and exon antisense from this run, because the results for those wouldve been the same anyway. 

Also, changed the `index_dict` pairings. 

Well so... maybe only a slight... very slight hint for DNA transposon? Maybe if we do a gene-DE, we can get something out of it? 



`Citation`: So for this one, we just move all the GffTE/ RM TE features after our citation index. Here all the citation indexes are grouped together 

Now what we notice is that the highest RMtransposon drops by around 100k, presumably this all goes into the citation indexes. (This is when we group all citation indexes together) 


`CitationByType`: this one we just modify our `index_dict` to split our Citation indexes. Nope, nothing 


`KlarsVeins`: aka, `FBgn0001316` and `FBgn0046776` respectively 

k seems like klarsicht and thickveins are both under genes 

Ok so i added CG17046 to the index_dict.

Hmnn, ok, lets see if transposon mapping affects it at all

K lets put the transposons behind first, to see if theres even anything significant to begin with

Ok, so both thickveins and klars dosent amount to anything


`OnlyTE`: well it seems like removing background didnt do anything. 	


`Antisense`: Ok, seems like it gets even more amorphous when we try to push the exonic antisense up, by removing TE




ok final things do do: 

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

Kawamura 2008
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



K now the serious stuff GSE26230

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
