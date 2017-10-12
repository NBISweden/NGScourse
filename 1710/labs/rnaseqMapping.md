---
layout: default
title:  'RNAseq'
---

# RNA-seq reference-genome based: introduction
<font color="red">**Please read everything carefully!**</font>

A common problem in the analysis of RNA-seq data is to relate it to a known genome sequence and use that information to study the expression of genes: within a sample or across multiple conditions, such as tissues or time points. A popular pipeline to perform such an analysis is the **Tuxedo protocol**, which consists of set of programs that can be used to go from mapping short reads to the reference genomes all the way to the detection of differentially expressed genes. The two main programs included in the package are i) **Tophat**, a short read mapper and ii) **Cufflinks**, performing analysis of the mapped reads.

In this exercise you will learn how to use some of these tools to study gene expression differences between different human tissues.
**The main goal is to find genes that are differentially expressed between two different tissues.**

We will use data that come from Illumina Bodymap2.0. [Illumina Bodymap2.0](http://www.ebi.ac.uk/gxa/experiments/E-MTAB-513) data consists of 16 human tissues that were sequenced using both single-end (SE) and pair-end (PE) technologies. The mapped reads can be visualised in the genome browser, e.g. [Ensembl genome browser](http://www.ensembl.info/blog/2011/05/24/human-bodymap-2-0-data-from-illumina/) or [IGV](https://www.broadinstitute.org/igv/)

In this tutorial, due to time constraints, we will focus on a limited number of tissues and narrow down analyses to one chromosome, chromosme 1 of the human genome. For all included tissues there is one single-end library and one pair-end library available. In order to identify signficially different genes, differentially expressed bewteen tissues, more than one replicate from each tissue is needed. We will therefore use the two different library types as replicates. For instace, for comparing brain and kidney one should include in the analyses both single-end (ERR030890) and pair-end (ERR030882) libraries for brain and single-end (ERR030893) and pair-end (ERR030885) libraries for kidney.

Here is the summary of the data and tissues available:

* Single-end reads, 75bp
    * ERR030888: Female adipose
    * ERR030890: Female brain
    * ERR030892: Female colon
    * ERR030893: Female kidney
    * ERR030901: Female ovary
<br/>
* Pair-end reads, 2 x 50bp
    * ERR030880: Female adipose
    * ERR030882: Female brain
    * ERR030884: Female colon
    * ERR030885: Female kidney
    * ERR030874: Female ovary
* human reference genome (or in this lab, only chromosome 1 named: rm.chr.1.fa)
* genome index for aligning reads with Bowtie2
* reference genome annotation based on the [EnsEMBL](http://www.ensembl.org/index.html) database named: Homo_sapiens.GRCh38_Chr1.77.gtf

<font color="red">**Please start by testing Brain vs Kidney!**</font> 
Afterwards you can go back and test more tissues against each other.

Note: Do not simply copy the various unix commands quoted throughout this tutorial.
Many include placeholders, e.g. folder names, so make sure you alter
the command to actually reflect whatever file names you have created.

**Tophat short description**

Tophat is a script pipeline built on-top of the popular short-read aligner Bowtie.
It is used for aligning RNA-Seq reads to a reference genome and can identify exon-exon splice junctions.
More specifically, it produces data that we can use not only to study the expression of genes, but also the expression patterns across different isoforms.
You will have a bit of waiting time during the exercises as the more complex analyses are running, so please check out some of the details of [tophat](http://ccb.jhu.edu/software/tophat/index.shtml) when waiting.

**Cufflinks short description**

Cufflinks is a collection of programs that perform different steps in the analysis of aligned RNA-seq reads ([Details](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html)).
The output is usually a list of transcribed loci (primarily ‘genes’) and their expression levels within and/or between samples.
For the analysis of multiple data sets, the general workflow in cufflinks consists of the following steps:

* **Cufflinks:** Assemble the aligned reads of a given sample, identify transcribed loci and determine expression
* **Cuffmerge:** Reconcile data on transcribed loci across multiple
  samples to produce a consensus annotation of loci (Note that if a good annotation is available this can be stepped, e.g. here one can use the prepared
  Homo_sapiens.GRCh38_Chr1.77.gtf)
* **Cuffdiff:** Compare read data across samples, guided by consensus annotation, and determine differential expression of loci, test for significance. The main output we are interested in comes from the cuffdiff analysis and consists of differential expression estimates for a set of genes.

## Step-by-Step Tutorial

1. Prepare your data
1. Load software
1. Run Tophat on individual samples
1. Run Cufflinks on individual samples
1. Run Cuffmerge merge detected transcript over all samples
1. Run Cuffdiff

#### 1) Book a node

We have reserved half a node for each student during this course.
By now, you are probably already familiar with the procedure:

<font color="red">**NB! Do this only once and make sure you do not have multiple
reservations running at the same time, otherwise you will take away resources from the other course participants!**</font>


```bash
$ salloc -A g2016008 -t 08:00:00 -p core -n 8 --no-shell --reservation=g2016008_4 &

```

#### 2) Prepare your data

<font color="red">**Note: It is completely up to your how you organize your data - what follows below is merely a suggestion:**</font>

<font color="red">**Remember that you have to run the analysis on SE and PE for two tissues! Examples below shows how to do this for one tissue...**</font>

* create a folder for your project

```bash
$ cd ~/glob
$ mkdir transcriptome
$ cd transcriptome
$ mkdir results
```

* sym-link the required files and folders (this will create a symbolic link to the original folders/files and saves you the trouble of always typing the full path - BUT: Do not write into these linked folders, because that data is shared across everyone working with these folders...)

```bash
ln -s /sw/courses/ngsintro/transcriptome_map/reads/PE
ln -s /sw/courses/ngsintro/transcriptome_map/reads/SE
ln -s /sw/courses/ngsintro/transcriptome_map/reference
```

Your directory structure should look like this:

* working directory (your choice)
  * PE/
  * SE/
  * results/
  * reference/

We are skipping a few steps here, namely obtaining a reference annotation and genome sequence and preparing the latter for use with the Bowtie2 aligner.
We have taken care of that for you (located in the subfolder /reference).
This is due to two main factors.
First, it takes a lot of CPU hours to convert a genome sequence into a Bowtie index.
Second, finding the latest release of a genome sequence free of unmapped fragments and haplotype data as well as a fully Tophat/Cufflinks-compatible annotation of that sequence is an exercise in frustration for beginners.
One useful resource here is the FTP server of Illumina [here](http://support.illumina.com/sequencing/sequencing_software/igenome.html).
Finally, also note that we are providing you with the outputs of the different steps.
This is to make sure that if you run into some trouble, like software crashing half-way through analysis, you can still continue with your exercises.

#### 3) Load software

You have done this before, but here is a quick reminder:

```bash
$ module load bioinfo-tools samtools/0.1.19 bowtie2/2.2.3 tophat/2.0.12 cufflinks/2.2.1
```

If any of these packages does not load as expected, you can check that module names are correct using the command

```bash
$ module avail
```

It may be that you need to load a different version.

#### 4) Run Tophat

What goes in and what comes out? 

In:

* One or several [FastQ](http://en.wikipedia.org/wiki/FASTQ_format) files (one for single-end reads, two for paired-end)
* An indexed reference genome (and optional a reference annnotation (eg. a [GTF or GFF file](http://www.ensembl.org/info/website/upload/gff.html)))

Out:

* A read alignment (BAM)

NOTE: The /reads folder contains a small subset of an actual FASTQ, limited to a specific chromosome.
For the next step, please chose two tissues that you want to analyse and make sure to align both a pair-end and a single-end library from you tissues of choice.
We do this since the time needed to align an actual read file with data from all human chromosomes can take several hours.
With these sub-sampled data sets, it should be possible to align them with tophat in a reasonable time frame.
However, if this should still take too long (>20mins), you may wish to abort this step.
You already have the corresponding output in the subfolder results/tophat/.

Tophat will take one or multiple FASTQ files and align the reads therein to a genomic reference.
A common command may look like this:

```bash
$ tophat -o results/tophat_outputSE30888 --solexa-quals -p 8 --library-type=fr-unstranded reference/rm.chr.1 SE/ERR030888.fq.gz

$ tophat -o results/tophat_outputPE30880 --solexa-quals -p 8 -r 200 --mate-std-dev 90 --library-type=fr-unstranded reference/rm.chr.1 PE/ERR030880_1.fq.gz PE/ERR030880_2.fq.gz
```

We specify the output location (-o), the number of CPUs to use (-p),
which type of sequencing library was used to produce the data (here
‘fr-unstranded’), in which format the quality information was stored
(here ‘solexa’, pre 1.3), the location of the reference annotation,
the location of the Bowtie2-formatted index file for the genome sequence and finally a FASTQ file. For pair-end data we have two FASTQ files and also define the expected size and variation in size of the fragments sequenced.

While this is running, you may want to head over to the [tophat manual](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html) and have a look at the available options and technical details.

The aligned reads are found in the output directory **you have chosen**, and named accepted_hits.bam.
For convenience, you may want to sym-link this file into your main
project folder. This links the accepted_hits.bam in the folder
results/tophat_outputSE30888 to results/SE30888.bam making it easier
to refer to the correct bam file in later steps.

```bash
$ cd results
$ ln -s tophat_outputSE30888/accepted_hits.bam SE30888.bam
$ ln -s tophat_outputPE30880/accepted_hits.bam PE30880.bam
$ cd ..
```
#### 5) Cufflinks: Assembly and transcript calling

What goes in and what comes out? 

**In:** A read alignment in BAM format (SAM is also an option, but should not be used due to it being uncrompressed)

**Out:** A number of files, including a transcriptome annotation reconstructed from the read distribution

Please adjust the name of the bam file as well as the output folder name to the appropriate names for your analysis.
**Note: Remember the path to the output folder that you choose below, you will need it later!!

** Run cufflinks on all tophat results 

General command format: 

```bash
$ cufflinks -o my_output_folder -p 8 -g reference/Homo_sapiens.GRCh38_Chr1.77.gtf my_infile.bam
```

Note: The Cufflinks step can take a while - so this is now a good time to get a coffee, or read through the available documentation on the Cufflinks website.
While we have tried to cover the technical details of these tools to some degree in the lecture, there are a lot of details that will help you use these programs to their greatest effect.

Here we specify where to store the output, how many CPUs to use as well as where to find the reference files (genome sequence and annotation).
Depending on the size of the BAM/SAM file, this step may require several hours to complete.
To make this analysis feasible within the time limits of the course, we have created the chromosome-limited files you have been using.

The command line output will read something like:

```bash
> Processed 48858 loci. [* * * * *] 100%

> Map Properties:

> Total Map Mass: 1893657.20

> Fragment Length Distribution: Truncated Gaussian (default)

> Default Mean: 200

> Default Std Dev: 80
```

One important thing that can be noted here:

* Processed loci - these are the transcribed regions, or 'genes'.
How does this number compare to offical estimates of human gene content?
* Total Map Mass - a measure for the size of your read library
* Length distribution - this value measures the distance between mate-paired reads.
It can either be specified (if known) or will be determined by Cufflinks.
Since we are using single-end reads, this value should not matter.
The output of this run can then be found under my_output_folder/ and includes a total of 4 files:

genes.fpkm_tracking

isoforms.fpkm_tracking

skipped.gtf

transcripts.gtf

The first two files contain basic information about expressed genes and transcripts, respectively - those known from the annotation file as well as novel loci identified by cufflinks -and the strength of their expression, given in FPKM.
FPKM stands for ‘Fragments Per Kilobase of exon per Million fragments mapped’ and is a normalized measure of transcript abundance.
That is the short explanation.
The longer version for the more mathematically inclined among us can be found at [Mortazavi et al. 2008](http://www.nature.com/nmeth/journal/v5/n7/abs/nmeth.1226.html).

These output files are tab-delimited and can e.g. be opened in e.g. Microsoft Excel (or similiar) to be analyzed and/or visualized.

#### 6) Cuffmerge: Reconciling different transcript models

What goes in, what comes out:

**In:** An optional reference annotation and a list of transcript annotations to merge

**Out:** A consensus annotation, taking into account all input annotations

It is important to keep in mind that reference annotations are very likely incomplete.
This is because some genes or individual exons may be expressed at very low levels or under specific conditions, thus having evaded prior detection.
Moreover, many vertebrate genomes have only been annotated by reference to other genomes, which themselves may only be poorly characterized.
Using the expression data obtained through cufflinks may hence allow us to improve existing annotations.
Cuffmerge is a tool that takes cufflinks-derived annotation files (known & ‘novel’ loci) and reconciles them into a consensus annotation, discarding e.g. spuriously transcribed loci and merging overlapping loci into larger transcription units where possible.

Again, the commands below are **just examples**, your files and folder will be called differently.

```bash
$ cd ~/glob
$ mkdir cuffmerge
$ cd cuffmerge
$ ln -s ../cufflinks.brainSE/transcripts.gtf brainSE.gtf
$ ln -s ../cufflinks.brainPE/transcripts.gtf brainPE.gtf
```

*Note: If this didn't work (check that the linked files actually exist and have content), then you probably chose a different way of organizing your folders and will have to figure out where the cufflink files your generated earlier are located ;)*

Now that we have both transcript model files in one location, we can attempt to merge them.
For this, we first have to create a text file that contains a list of GTF files to merge (quite inconvenient, I know).
Use whichever tool you feel comfortable with and write the name of
each gtf file line by line, then save it as transcripts.txt.
If you would like to be fancy, you could try to use unix commands to do this (tip: use ls and ">")

```bash
$ cuffmerge -o merged -g reference/Homo_sapiens.GRCh38_Chr1.77.gtf -p 8 -s reference/rm.chr.1.fa transcripts.txt
```

This will save the reconciled annotation file as merged/merged.gtf.
Symlink this file into your main project folder.

```bash
$ cd ~/glob/transcriptome/
$ ln -s cuffmerge/merged/merged.gtf
```

Now we are ready to check for differential expression in our read data from chromosome 1.

#### 7) Cuffdiff: Differential expression analysis

What goes in, what comes out:

**In:** A consensus annotation (or just the reference annotation), read alignments for all samples that are to be compared and quantified

**Out:** A number of files, including tests for differential expression for all pairwise comparisons (gene_exp.diff)

Cuffdiff takes aligned reads from two or more samples, estimates comparable expression values and performs statistical analysis of the resulting data to determine which genes exhibit significantly different activity profiles between any two samples.
The nitty-gritty details of the underlying mathematics can be found here: [cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html).

For running Cuffdiff, we type something like this (being in the main directory of our project):

```bash
$ cuffdiff -o cuffdiff.brain_vs_kidney -L brain,kidney -p 8 -u merged.gtf brainPE.bam,brainSE.bam kidneyPE.bam,kidneySE.bam
```

Adopt this to your data - if uncertain, run cuffdiff -h to learn more about the options.
Note that the labels you give them are arbitrarily chosen, pick names that make sense to you.

This will write the output of the analysis into the subfolder cuffdiff.brain_vs_kidney (or whatever folder name you chose) and do a pairwise comparison of the samples (Note: the order in which you list the labels needs to be the same as the order of SAM/BAM files!).
NB! To reduce the computing time we do not use the flag -b that correct for genome sequence (see manual for details), but we resolve issues arising from reads mapping to multiple loci in the genome using the flag -u.

The main file of interest to us is gene_exp.diff.
It includes the analysis of differential expression.
A quick way to find the cases of interest is to filter the file for genes that show evidence of differential expression (identified by the tag ‘yes’ in the ‘significant’ column).

```bash
$ head -n1 gene_exp.diff > results.txt
$ grep yes gene_exp.diff >> results.txt
```

(This copies the header of the output file as well as all rows tagged as significant into a new text file - open this file in a text editor or spread sheet program).


### Where to go next: summarizing and visualizing the results
After obtaing the differential expression output, the data analysis part begins. Quite often, one is interested in summarizing the differential expression results, e.g. reporting number of down- and up-regulated genes given a significance and fold change threshold or showign the results in a graphical form. There are no standard and easy solutions here and the data analysis part tend to be guided by the specific biological questions. Here, quite often a custom-written R and Python scripts are used. 

To give you a feel of running one of such scripts for a very basic results summary and visualization, a very simple, dependency-free R script is available (summarizeDE.R). It can be run from the command line to summarize the Cuffdiff results and create a Volcano plot (what is Volcano plot? Google it).

```bash
$ cd ~glob/transcriptome
$ cp /sw/courses/ngsintro/transcriptome_map/extras/summarizeDE.R ~/glob/transcriptome
$ Rscript summarizeDE.R ~/glob/transcriptome/cuffdiff.brain_vs_kidney/gene_exp.diff 0.1 1
```
where,

'0.1' is the FDR threshold
and 
'1' is the log2FC

The script should print the number of differentially expressed, down- and up-regulated genes given the FDR and log2FC threshold. It also prints the entries for the top 10 down-up and up-regulated genes. Finally, it creates a .pdf file with a Volcano plot (google Volcano plot for explanation if not clear)

To open the Volcano plot from Uppmax

```bash
$ firefox Volcano_0.1_1.pdf 
```
You can experiment with summarizing your results using different cut-offs for FDR and log2FC. Does the number of differentially expressed genes change? Does the Volcano change?

After obtaining the differential expression results, it it common to transfer the files to the local computer. One can get far by just looking at the results in e.g. Excel, running R scripts or using online tools. 


### Closing remarks

This tutorial has introduced you to a very straight-forward, but somewhat simplified pipeline for the analysis of RNA-seq data by use of a reference genome to study transcription.
Both Cufflinks and Tophat come with additional parameters that we have not touched upon to avoid unnecessary confusion.
Likewise, the read data we have used was strand-unspecific.
This has some drawbacks, specifically with respect to accuracy in the isoform analysis.
Or perhaps you are not interested in comparing expression between pairs of samples but in a time series.
For this reason as well as others, you may need to adjust one or several parameters to get the best results - depending on the nature of your data.
<font color="red">We therefore highly recommend you to carefully read the manuals (and the original publications)  to familiarize yourself with these additional options.</font> 

### Other alternatives

There are several different tools available that have the ability to
do many of the steps described above with the Tuxedo pipeline. Here
are a few options for the different steps

#### Short read mappers that are suitable for RNA-seq data

* [Star](https://github.com/alexdobin/STAR)
* [Subread](http://subread.sourceforge.net)
* [Gsnap](http://research-pub.gene.com/gmap/)

#### Counting reads from mapped data (bam files) and a annotation file (GTF/GFF file)

* [HTseq](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html#count)
* [Featurecounts](http://bioinf.wehi.edu.au/featureCounts/)

#### Mapping & counting RNASeq
* [RSEM](http://deweylab.github.io/RSEM/) 

#### Detect differential gene expression

Most software suitable for detection of differential gene expression
are developed in R and available as add-on packages to R at
[bioconductor](http:/www.bioconductor.org)

* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* [Limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
* [DEXseq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html)

