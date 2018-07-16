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


## --config variables

`snakemake -s test.Snakefile --config variable='one'`

`variable = config['variable']`

`rule all:`
	
	input:
		'one.done' if variable == 'one' else 'two.done'

`rule first:` 

	output:
		touch('one.done' if variable == 'one' else 'two.done')
	run:
		pass 

So it seems that we have to restate the config variable in the snakefile

**SUPER IMPT, every snakefile is initialized with an invisible `config` dictionary with no keys and values**
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

* `snakemake -s test.Snakefile --verbose -p`

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


Packages installed on atlas:

1. gffutils  
   




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
* `os.chdir('')` this is gonna be super sweet to guide the path for every snakefile. Even better, it does not screw up the pathing for `config` file.  
* `os.listdir( path )` to list every file in path 
 
### os.path

* `os.path.isdir(path)`: Returns boolean of whether directory specified by path exists.   
* `os.path.join(path, *paths)`: Returns the concatenation of 2 paths, both are strings. Ah and also note, every comma essentially == `/`    
* Also, if we use the comma, we dont need `/` at the end of each directory. 

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

## Misc
Consider using both for more advanced printing in the future  
`print(feature, "parenttype=" + ft,sep=";")`  
`"{} {}".format(feature, "p")`

Reading in files	sys.stdout = open('DMr6.21excludeUextra.fasta','w')
	lines = open('dmel-all-chromosome-r6.21.fasta').read().split("\n")
	for line in lines:
		if line.startswith('>2Cen_mapped_Scaffold_10_D1684'):
			break
		else:
			print(line)
	sys.stdout.close()

##Writing to files
* `f= open("guru99.txt","w+")`, `+` creates file if dosent already exist and opens it for writing 

# Pandas 


## Series

This is essentially just a horizontal array, initialized with `s = pd.Series(data,index)` , data can be a dictionary, index is a list

* `s[:3]`, slicing goes by rows and `s[row]` we can retrieve rows by index 


## Dataframes

columns (column labels) , index (row labels)  

### Creating a Dataframe 


`sales` =   `{'account': ['Jones LLC', 'Alpha Co', 'Blue Inc'],'Jan': [150, 200, 50],
         'Feb': [200, 210, 90],
         'Mar': [140, 215, 95]}`

`df = pd.DataFrame.from_dict(sales)`

So in this example, columns are the keys and our row values are kept in a list. 

### Column and row manipulation

* `df['one']`: indexing is done via columns  
* `df['new']`: adding a new column by calling a new header  
* `df['new'] = df['one'] + df['two']`: this works for text no problem. 


**`.loc[]` works on labels of your index, `.iloc[]` works on the positions in your index.**  
**Both methods return a series and try to use .loc[]**   
* Essentially its the same as working with columns, just throw a `.loc[]` at the back.   
* When adding/setting rows with `.loc[]`, added entry must fit current dimensions


* `df.drop(df.index[], axis = 1)`: For dropping entire rows/columns Optionally axis=1 ensures that we refer to a column and not row.  


* For the code snippet below, seems like indexing -> splitting a column returns a `Dataframe` of all the split data. Ah so it seems that **if we index columns, we get a `series` but after we `split` it we get a dataframe**. Also note, we can just split a series directly without prior indexing `(geneID_field.str.split(';',expand=True))`

def f(gff):

	gffdata = pd.read_csv(gff,sep='\t', header=None, index_col=None)
	split = (gffdata[8].str.split(';',expand=True))
	print(split)
	print(type(split))


* Whenever we `df[:6] index` , pandas sees it as row manipulation
* `bed.drop([6,7,8,9], axis =1, inplace = True)`, seems like we have no choice but to manually key in the columns.
* `bed[3] = bed[3].astype(str) + ";parenttype=mRNA"`, note for adding a string to every value in column. 


### Data entry manipulation 

* `df.loc[df[0] == 'goodbye']`: Removes rows where entry in target column is not 'goodbye'
* `df.loc[df[10].str.contains('LINE')]` checks for contains. seems like `Series.str.` can be quite helpful
* `df.loc[df[0].str.startswith('>')]`
* `df.drop_duplicates(0 /[0,1], keep = '' )`: Removes rows with duplicate entries in target column(s).
* `out = df[[4,5,6,10,7,8]]` **note, the resultant df keeps old header labels** 


### Other functionalities

* `df.sort_index(axis)` sorts the dataframe based on the **values** of specified axis, `df.sort_values()` sorts by actual table values of column  
* cool stuff, did abit of iterrows today, found out you can index the row object.  

			for index,row in tab_data.iterrows():
				data = row[:2]
				for line in data:
					print(line)

* `len(df or series)`
* `.astype('int')` for datatype converison and `series.sum()`
## pd.read_csv()

* `f = pd.read_csv(test , sep = '\s+', header = None, skiprows = 2)`

The regex is to deal with multiple variable whitespace delimiter and `skiprows` essentially just skips to row 2, so instead to configure paths and such properly, we use the Rstudio server @ `http://atlas.cbis.nus.edu.sg:8787` accessed through browser.



# R  

Ok, so what we learned is that if we open a file with our local Rstudio, when the script is run, it is run on the local system.

hmnn ask about global/ local packages later 

## General functionalities 


* `iris %>% head() %>% summary()` is summary(head(iris)), where `%>%` acts a pipe and if we want a line break, 
we need to end the line with `%>%`  
* `library(module)` is essentially your import
* `newCol <- arrange(newCol, petal.width)` seems like no issue re-assigning variables back to themselve
* `k <- c(0.5, 1)` returns `[0.5,1]`, **ALSO NOTE, R INDEXING IS 1-BASED** aka k[1]
* `paste(*args, sep = '')` converts its arguments to character strings, and concatenates them with the sep
* `print()` also a print function
* `rm(var)` to remove from global environment 
* only 2 main R data types to think about `numeric()` and `character()`
* `$` is good for general purpose indexing, but when we want to index with a predefined variable/ string, we have to use something like this `config_data[[paste(run_ID, "_index_list", sep="")]]`
* `-` are read as `.` in R 


##Datatypes

###Lists 
* So we learned today of this thing called **named lists**  
`lst = list(bob=c(2, 3, 5), john=c("aa", "bb"))`  and we access the 'keys' in here via `$bob`
* `names(list/df)` to get items/ column names
* `test$a = lst[test$a]` to use a named list for key matchings in order to modify all values in `test$a` but ...
* Note how, `typeof(lst[name]) > list` but `typeof(lst$name) > etc`, so since using variables require [], we have to use `lst[name] %>% unlist()`
* Wait a second, [[]] makes a huge ass difference when indexing from config_data, its kinda like extracting from an internal list. Example being:  

`var = count_df['B133',]` #slicing a row from df, typeof() = list  
`var['dmel-all-chromosome-r6.21']` #indexing column, **typeof() = list**  
`var[['dmel-all-chromosome-r6.21']]`#indexing column with [[]], **typeof() = integer**

* [[1]] means that R is showing the first element of a list.
* * you can use both [] and [[]] on lists, the [] construct will return a subset of  the list (which will be a list) while [[]] will return a single element of the list (which could be a list or a vector or whatever that element may be)


##Rscript (Running from command line and feeding in arguements) 
* `Rscript script.R` apparently to run r scripts from the command line and 
* `Rscript --vanilla script.R var1 var2`
* Ahh its because every R script is initialized with a `args` string vector


##here 
* `here(path_var,'file.txt')`, essentially a pjion that **STARTS AT the directory where the script is called from NOT where the script is located** 
* So, when running from Rstudio, the default here() path is `mnt/raid0/home/jinwee/`


##yaml
* `yaml.load_file("path_to_file/file.yaml")`
* This then reads the yaml file into a **named list** structure, where each dictionary set becomes a named list, accordingly the nested dictionaries as well. 

 


## Tidyverse and Dataframes 

Alright so apparently Tidyverse contains a whole host of subpackages `stringr`, `ggplot2`,
`dplyr`.... et..

###Dataframes
* For all questions so far, refer to `Dataframe_cheatsheet`
* `nrow(dataset)` and also `ncols(dataset)`
* `norm_count_df$Mapped_index =  config_data$index_dict[as.character(norm_count_df$Mapped_index)] %>% unlist()` 
Apparently, this is how we change every value in mapped_index to its corresponding value in the dictionary/ config data .
* Alright so apparently this works to divide every row by its spikein count

`spikein <- count_df$spikeinset01`  
`test<- count_df[,-ncol(count_df)] * 1000 / spikein`

A simpler illustration of how we can apply lists to matching rows in dataframes using operators  

`test <- data.frame('A'=c(1,2), 'B'=c(3,4), 'C' =c(5,6))`  
`colC = test$C`  
`test + colC`  
A  B  C  
6  8 10  
8 10 12

###rownames
* `has_rownames(df)`
* `remove_rownames(df)`
* `rownames_to_column(df, var = "rowname")`
* `rowid_to_column(df, var = "rowid")`
* `column_to_rownames(df, var = "rowname")`

###summary
calling summary on a list gives us

`> summary(spikein)`  
  ```Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   9076   11197   15072   15902   16342   39376```

###apply and lapply and sapply

`test<- apply(count_df,1, function(x) x * 1000 / x[["spikeinset01"]])`

* Ok kids, lapply always returns a list while `sapply` returns the **simplest** data structure, which is apparently a vector.
* apparently today we learnt that we can apply() a vector column/row wise to a dataframe if the vector vs column/row have the same length. 


###dplyr

* `filter(iris, species == "virginica", length > 6)` Filter is for getting rows based on column values
* `select(iris, length, width, plength)` Select is for columns and will somehow always carry over the goddamn column name. Note this only works with unique column names
* * `subset(iris,rownames(dataset) == names)` Subset can be for both row and columns, 
*  `mutate(iris, new_col = width > 0.5 * length)` Adding a new column of boolean values, note, try to see if we can add an entire series/ existing column.
*  `arrange(newCol, width)` default ascending, alternatively `arrange(newCol, desc(width))`


###readr

*`read_tsv(file, col_names = TRUE / col_names = c("col1","col2"))` 


###filtering rows based on condition 2 slightly different examples here
* `test <- norm_DE_df[apply(norm_DE_df>5, 1, all) & rowSums(norm_DE_df)>500,]` 	

###filtering vector based on conditional
* `test4 <- colnames(norm_DE_df)[grepl('Control', colnames(norm_DE_df))]`


###Group_by and aggregate,
seems like group_by is for rows,
###summarise and summarise_all



















