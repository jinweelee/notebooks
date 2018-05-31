# Pandoc/ markdown

- 2 spaces indicate a line break  
- Dash for bullet point with indentation  
- ** double asterisk to wrap bold and _ to wrap italics. **_we can combine both wrappers_**
- Direct link <http://google.com>
- [Inserting an inline link](http://youtube.com)
- [Inserting a reference link][someplace]

[someplace]: http://youtube.com

- Making a nested list is as simple as adding a 1 space indentation 

- Inserting an image is kinda just like inserting a link, just that we add an ! infront of the link

![And this text will appear below the image like a captio](https://i2.wp.com/outlookaub.com/wp-content/uploads/2015/04/kafka-on-the-shore.jpg)



# Git 

- 

- When we `git add` a file, it means we are tracking it and git knows about its existence and will be included in the next commit  
But changes to it will not be automatically staged for the next commit.    

- To do so, we will have to either `git add` again or when we commit, `git commit -a -m "title"` to automatically stage and commit all changes.  

- `git log` seems super useful for tracking past commits and probably retrieving them

- `git mv` and `git rm` are used so git can keep track of renaming /deletion   

- `.gitignore` file should be set up in each main repo to specify which files shouldnt be tracked/ staged.  

- `git reset file` to unstage a modified file which was staged/ added  

- `git remote add origin git@github.com:jinweelee/repo` -> linking local repo to remote repo -> `git push origin master`  

- `git remote -v`, note how `origin` is the convention for the main repo   

- `git checkout -- file` reverts the said file to the most recent commit or if a specific version in mind `git checkout SHA-1 ID -- file`  

- `git revert` and ` git reset` note . 

- `git checkout branch` to swap branches.  

- `git push origin branch` just like pushing master and `git fetch origin` to get any new branches that might be added  




# Snakemake

## Wildcards 

- If we need to use the file prefix in our shell command `{wildcards.sample}` that is not in the already defined output. 

```

rule `samtools_sort`:  
    input:   
      "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

```

## Targets 
```  

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

```  

bam can be interpreted as ["sorted_reads/A.bam", "sorted_reads/B.bam"]   

## External Scripts 

External script path is always relative to the Snakefile


## Flags

* Hmmnn, okay it seems that if we an intermediate `.done`, the pipeline will not rerun as long as the final `.done` is still there. 




##Misc 
* Try this way of using `input` as a variable to iterate through 

```

`rule NAME`:  
  
    input: "path/to/inputfile", "path/to/other/inputfile"
    output: "path/to/outputfile", somename = "path/to/another/outputfile"
    run:
        for f in input:
            ...
            with open(output[0], "w") as out:
                out.write(...)
```

* We now store all "runlist" type information, in `default_config.yaml` file. It also appears that by convention we have to read in our file with the variable `configfile` and later on its named as `config`.

# Bioinfo/Tech packages

1. Markdown notepad
2. MobaXterm
3.  

# Ubuntu/ Unix Changelog 

Reconfigured sources.list and .bashrc

Installed Packages in loqs:  

1. Pandoc  
3. Miniconda  
5. Snakemake
1. Git
1. GffUtils   
2. Tandem Repeats Finder (executable binary moved to /usr/local/bin)
3. RMBlast 


Updated Packages in loqs:

1. date-util 




# Ubuntu/ Unix Help

* To view installed packages: `dpkg -l` or `conda list`  

* `tar xopf foo.tar` to unzip tar files 

# Python

* `int()` and `str()` 

* Ah for string slicing, its positive 0-based and negative 1-based. For both, end index not inclusive 

* `type()`

## os

* `os.listdir(path)`: Returns list of strings of all files in directory specified by path    
* `os.mkdir(path)` : Creates directory basd on path  
* `os.getcwd()`: Returns current working directory   
 
### os.path

* `os.path.isdir(path)`: Returns boolean of whether directory specified by path exists.   
* `os.path.join(path, *paths)`: Returns the concatenation of 2 paths, both are strings.   

## str

* `str.endswith/ startswith(suffix)`  
* `str.split(sep)`: Note, this returns a list of strings.  

## list

* `len(lst)`, `lst.count(x)`, `lst.clear()`, `lst.remove(i)`, `lst.pop(i)`, `lst.extend(lst2)`

## dictionary
*  Note we can already iterate through `dict.keys()` but it will return us keys /values in some weird uncomprehensible order.
 
* To get **in the order they were added,** essentialy their position in the dictionary, convert to list.`list(dict.keys())` / `list(dict.values())`
* `dict[key]`
* 



# Pandas 


## Series

This is essentially just a horizontal array, initialized with `s = pd.Series(data,index)` , data can be a dictionary, index is a list

* `s[:3]`, slicing goes by rows and `s[row]` we can retrieve rows by index 


## Dataframes

columns (column labels) , index (row labels)  

### Creating a Dataframe 

1. Passing a dictionary as input. Keys are treated as columns.

`df2 = pd.DataFrame({ 'A' : 1., 'B' : pd.Timestamp('20130102'),'C' : pd.Series(1,index=list(range(4)),dtype='float32'),})` 

2. Specifications in the `columns` field as a list   



### Column and row manipulation

* `df['one']`: indexing is done via columns  
* `df['new']`: adding a new column by calling a new header  
* `df['new'] = df['one'] + df['two']`: this works for text no problem. 


**`.loc[]` works on labels of your index, `.iloc[]` works on the positions in your index.**  
**Both methods return a series and try to use .loc[]**   
* Essentially its the same as working with columns, just throw a `.loc[]` at the back.   
* When adding/setting rows with `.loc[]`, added entry must fit current dimensions


* `df.drop(df.index[], axis = 1)`: For dropping entire rows/columns Optionally axis=1 ensures that we refer to a column and not row.  



### Data entry manipulation 

* `df.loc[file[0] == 'goodbye']`: Removes rows where entry in target column is not 'goodbye'
* `df.drop_duplicates(0 /[0,1], keep = '' )`: Removes rows with duplicate entries in target column(s).


### Other functionalities

* `df.sort_index(axis)` sorts the dataframe based on the **values** of specified axis, `df.sort_values()` sorts by actual table values of column   


## pd.read_csv()

Sample input: `pd.read_csv(`hello there   i   am`, sep='\t', header=None, index_col=None)`
        
Sample output:
````
       0      1  2   3
0  hello  there  i  am
````

# R  