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

  