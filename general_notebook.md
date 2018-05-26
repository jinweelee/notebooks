|<General Notebook> 

# Pandoc/ markdown

* 2 spaces indicate a line break  
* Asterisk for bullet point with indentation  
* For code lines and blocks use apostrophe   

````
def fn(){

}
````

Inserting a link <http://google.com>  

# Ubuntu/ Unix Changelog 

Reconfigured sources.list and .bashrc

Installed Packages in loqs:  

1. Pandoc  
3. Miniconda  
5. Snakemake  

# Ubuntu/ Unix Help

* To view installed packages: `dpkg -l` or `conda list`  



# Snakemake

## Wildcards 

- If we need to use the file prefix in our shell command `{wildcards.sample}` that is not in the already defined output. 

````
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

````

## Targets 
````
SAMPLES = ["A", "B"]

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

````
bam can be interpreted as ["sorted_reads/A.bam", "sorted_reads/B.bam"]   

## External Scripts 

External script path is always relative to the Snakefile


##Misc 
- Try this way of using `input` as a variable to iterate through 

````
rule NAME:
    input: "path/to/inputfile", "path/to/other/inputfile"
    output: "path/to/outputfile", somename = "path/to/another/outputfile"
    run:
        for f in input:
            ...
            with open(output[0], "w") as out:
                out.write(...)
````


# Git 

When we `git add` a file, it means we are tracking it and git knows about its existence and will be included in the next commit  
But changes to it will not be automatically staged for the next commit.  


# Python

* `int()` and `str()` 

# R  