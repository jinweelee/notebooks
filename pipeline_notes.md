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


for hpRNA we just `grep "Name=hpRNA" gene.gff`, then settle the isoform related stuff downstream. 


K i am also going to omit the PTID file, i think its unecessary. So we will only extract the information pertaining to `Parent=FBgn0031208` and then subsequently drop duplicates. 

Also note, *`transposable_element` have no exons*, so we just merge them later. 

wow actually the `exon.gff` file is pretty useless. 
Note, when printing we also add the parent type to the end 

Okay, we are also going to try and stop hard coding the lib-vers in there 


##GFF2bed

Seems that after this, geneID is set at column "4"
 
```
2L      901490  901654  eID=FBgn0003916 .       -       FlyBase exon    .       Parent=FBtr0078028;geneID=FBgn0003916;parenttype=snRNA
```

I suspect this is a property of me adding the geneID information to that field at the end. 

ok, it seems that we have an issue with the gene.bed file that makes it unable for us to use gff2bed. I think we scrap it. 

##bedtools

K seems like we have to remove the additional fields from the `gff2bed` output and format it abit, because apparently `bedtools merge -s` cant work properly when theres more than 6 fields present. 

For some really weird reason, `bedtools intersect -S` gives a different output but same number of features with and without `wa & wb` 

Ahhhhhh i see, its because bedtools intersect only reports the feature in -a file, it does not actually merge the details. this is confirmed ? There should be a 1 -> 1 correspondence 

*AHHH THE KEY IS -wb** 

##Gffutils
K several things with the GFF file. I think if we use GFFutils we dont need a cleanup

`gffutils.FeatureDB.features_of_type(featuretype)`  
Returns an iterator of gffutils.Feature objects.  
Could be useful. 


##Repeatmasker

ok... its not tab delimited and to prevent tearing hair out delete all " *" 


##Citation Indexes
* `04_Czech` features taken from Supp info, S2 cells.

* `05_CR14033` discontinued as citation index, now directly downloaded as fasta from flybase 

##mirBASE

* Ok i think its the right call to use mature *miRNA* features. 

#Mapping 

##Mapping algorithmn

Ok, so both our external libraries and OKamura libraries are going to use 1 central `sequential_mapping` snakefile, with the necessary specification of **2 additional parameters** at the command line. 

* `lib_type` = `okamura` or `external`
* `dataset` = `GSE...`

**Note, second config parameter only necessary if `lib_type == external`

the purpose of this custom configS will be to:
1) which external dataset we map from

2) specify whether we map to spikein or not, using a conditional wrapped around sequential mapping rules 


Maybe we will apply this variable passing thing when we do RPM / SPIKein normalzation later.

`runlibstatmapping.sh` is the liling script 

##Logfiles
Multiple vs Single mappers(read identifers): filter 1 occurence vs multiple occurence.
Multiple vs Single mappers(counts) : filter 1 occurence vs multiple occurence. -> gather counts for each

Mapped vs Unmapped (read identifiers): .drop_duplicates vs extract .startswith(">") in unmapped. 
Mapped vs Unmapped (counts): -> gather counts for each. 



##RPM and spikein normalization
RPM and spikein correspond to `mapcolfasta.sh` and `countfeaturesbysensereads.sh` 

## External Libs 

Ok, so i think we are gonna split: 

no theres not much point, we should just run custom scripts, and we make sure all converge at "Lib_final_col_filter.fasta".

And `External_libraries` folder will have the same architecture with 3 internal folders,
`final_preprocessed`, `original_libs`, `temp_libs`

Additionally, within each of these 3 folders theres going to be an additional folder for each dataset, i.e
`final_preprocessed/GSE17171/` -> containing all the final preprocessed GSMs/libs of that dataset. 


Same for output `Mapping_output/External_libs/Mapped_fasta/GSE17171`


We are going to have to set up an additional library parameter ot handle each dataset and its libaries. 


## GSE26230 
Perl scripts were used to sort according to barcodes, remove the 3'-linker sequences and select for 19-25 nt long inserts (available upon request). Bowtie was used for mapping against the target sequences specified in the manuscript.

Ok so in addition theres also this hairpin/transposon mapped file. 

1. Extract .tar file 
2. Extract *.txt.gz files, store everything else as misc
3. .txt files are your fastq files
4. use cut adapt to convert from fastq -> standardized fasta 
5. enter sequential mapping 
6. count mapped/unmapped by library 
7. r scripts 


## GSE17171

K seems like this is already mapped and normalized. 

Renyi then took these sequences and re-mapped them. 

The 36-nt long Illumina reads were stripped of the 3’ linker, collapsed, and the resulting small RNA sequences were mapped to the Release 5 genome excluding the Uextra portion.For annotation, we used Flybase for protein coding genes, UCSC for non-coding RNAs and transposons/repeats (repeat masker tracks) and the most recent miRNA catalog. For endo-siRNAs from structured loci we used custom annotations. To control for sequencing depth, small RNA counts were normalized to the same total number after subtracting those matching abundant cellular RNAs such as rRNAs and tRNAs. SEQ = the sequence of the small RNA. MAPNUM = the number of times the small RNA matches the respective genome sequence. TYPE = annotation category of the small RNA. READS = read count in respective library. All small RNAs were converted into DNA for sequencing, therefore U's are represented as T's.



## GSE37443

Supplementary files format and content: Sequences + number of reads. Barcode and adapter sequences are removed, and the same reads are pooled together.

Seems like this is raw, adapter trimmed filtered reads. 



## GSE11019
Seems like this is raw, adapter trimmed filtered reads. 

## GSE11086
Also seems like raw reads 