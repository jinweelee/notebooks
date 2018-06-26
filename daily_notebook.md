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



#wait what 
