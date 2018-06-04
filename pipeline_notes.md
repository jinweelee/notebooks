#Dependencies on atlas
Installed: 
 
- gffutils 

# Overall changes 

I think we should have a Init.snakefile that is gonna initialize the general architecture of the project.
Then for each snakefile, we store the general path as a variable right at the start, so we minimalize clutter. 

Also we need to think of all the paths we need to initialize for each library now and, 
i think for each run we first check if the directory for the library has been created, if not we then move the file to that folder.   

We also really need to think about how to run indvidual libraries/adding libraries without going through the whole damn pipeline again.

We could just change the runlist, then delete the flags

Also we need to find out what happens when we delete a "middle flag", aka trying to re rerun 1 part. 




K i think for sequential mapping we can do `{lib}_{bt_index}.ebwt, expand`



#Index Generation

Ok, i think we defintely need to reconfigure this first, before going into sequential mapping 

K, so all accessory files will be directed to `Index_output`, with `Bowtie.indexes` and `Final_index_fasta` storing both the bowtie indexes and the fasta files required to make the bowtie indexes respectively. 

Hmnn should we do an index run_list, think about it. we can store both 1)indexes/ 2)corresponding type/name in a file,  
read them both in separately as lists, then make a dictionary from them to use later.  


OK so we are just gonna lump all our runlists into 1 config file.

general_feature_types=['exon','rRNA','RNA','tRNA','snRNA','snoRNA','pseudogene','pre_miRNA','mRNA', 'ncRNA','transposable_element']
might also want to include FT= gene 

exon_feature_types =
FBEXPT1=${FBDMELBASE}.txt.exon.mRNA.merge
FBEXPT2=${FBDMELBASE}.txt.exon.ncRNA.merge
FBEXPT3=${FBDMELBASE}.txt.exon.rRNA.merge
FBEXPT4=${FBDMELBASE}.txt.exon.tRNA.merge
FBEXPT5=${FBDMELBASE}.txt.exon.snRNA.merge
FBEXPT6=${FBDMELBASE}.txt.exon.snoRNA.merge
FBEXPT7=${FBDMELBASE}.txt.exon.pseudogene.merge


bedtools merge -s : same strand, -S opposite strands

Ok, so we need to: 

1) split.gff to parent types  
1b) Label exon parent types  
2) split exon.gff based on parent types  
3)Extract miRNA, hpRNA information    
4) bedtools merge -s, forcing same strandedness/ also perform withouy strandedness? 

Then for cis-nat  

4) bedtools intersect -s exon.MRNA vs all other types   
5) bedtools intersect -S gene.bed  
6) bedtools map exon overlap onto gene overlap.  
7) And i think we will leave out the cis/nat side for now, can always incorporate next time.  


k so what did renyi do:  
`db.featuretypes()`
`db.features_of_type`

Note exon  names, CG.. 

k fuck it, we are doing things my way from now onwards but we learn from what he did

So i think we will:  
  
1) For each parenttype (mRNA, etc) we call out all the features using `db.features_of_type`.  

2) Generate doc/list of all ParentIDs eg,  for each parentype `Parent=FBgn0031208`

3) Then use db.childern(type=exon) to get all the exons for each parentID, append them to a txt file. 


##Gffutils
K several things with the GFF file. I think if we use GFFutils we dont need a cleanup

`gffutils.FeatureDB.features_of_type(featuretype)`  
Returns an iterator of gffutils.Feature objects.  
Could be useful. 





# External Libs 

Ok, so i think we are gonna split: 

no theres not much point, we should just run custom scripts, and we make sure all converge at "Lib_final_col_filter.fasta"


## GSE26230 

1. Extract .tar file 
2. Extract *.txt.gz files, store everything else as misc
3. .txt files are your fastq files
4. use cut adapt to convert from fastq -> standardized fasta 
5. enter sequential mapping 
6. count mapped/unmapped by library 
7. r scripts 


## GSE17171

Seems like the libraries are in .csv format and to be confirmed with renyi what he did with the mappers/ non mappers, whether he used both or not  
