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

- When we `git add` a file, it means we are tracking it and git knows about its existence and will be included in the next commit  
But changes to it will not be automatically staged for the next commit.    

- To do so, we will have to either `git add` again or when we commit, `git commit -a -m "title"` to automatically stage and commit all changes.  

- `git log` seems super useful for tracking past commits and probably retrieving them

- `git mv` and `git rm` are used so git can keep track of renaming /deletion   

- `.gitignore` file should be set up in each main repo to specify which files shouldnt be tracked/ staged.  

- `git reset file` to unstage a modified file which was staged/ added  

- `git remote add origin git@github.com:jinweelee/repo` -> linking local repo to remote repo  

- `git remote -v`, note how `origin` is the convention for the main repo 

# Python

* `int()` and `str()` 

# R  