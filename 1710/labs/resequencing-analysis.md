---
layout: default
title:  'Resequencing Analysis'
---


# Resequencing Analysis

The data we will work with comes from the 1000 Genomes Project.
Because whole human genomes are very difficult to work with, we will use only a small portion of the human genome, a little over a megabase from chromosome 17.
Samtools have been used to extract the data from the 1000 Genomes ftp site for just this region from all of the individuals from the CEU (CEPH Europeans from Utah) population who were low coverage (2-4x average) whole genome shotgun sequenced.
We have 81 low coverage Illumina sequences, plus 63 Illumina exomes, and 15 low coverage 454 samples.
There are 55 of these samples that were done both ways.

We will walk through alignment, alignment processing and cleanup, quality recalibration, variant calling, and variant filtering.

In order to do these exercises, you will need to know a few things.

## Book a node

By now, you are probably already familiar with the procedure:

```bash
salloc -A g2017019 -t 08:00:00 -p core -n 8 --no-shell --reservation=g2017019_wed &
```

Make sure you ony do this once, otherwise you will take away resources from the other course participants! 
Once your job allocation has been granted you can connect to the node you got using ssh, just like in the [Uppmax Introduction exercise](uppmax-intro) yesterday. I.e. use 
```
squeue -u username
```
to find out the name of your node, and then 
```
ssh -Y nodename
```
To connect to the node. 

## Accessing Programs

First, we're going to run several programs that are installed under the module system.
To access the bioinformatics modules you first need to load the bioinfo-tools module:

```bash
module load bioinfo-tools
```
This makes it possible to load the individual programs we need:

```bash
module load bwa
module load samtools
```

We will also use Picard and GATK.
These are java programs, which means that we need to explicitly invoke java each time we run them and we need to know where the code actually lives. For various parts of this exercise, you will need to know the location of the executable jarfiles for GATK and Picard:

```bash
/sw/apps/bioinfo/GATK/3.4-46/
/sw/apps/bioinfo/picard/1.69/kalkyl/
```

For the other programs, you can just type the name of the program and it will run.
You can even tab complete the name.
This is what the module system does for you.

## Accessing Data

You need to know where your input data are and where your output will go.

All input data for the first steps of these exercises is located in this folder:

```bash
/sw/courses/ngsintro/gatk
```

Since we're all sharing the same data, we've made these file read-only so that no one accidentally deletes them or writes over the raw data or someone else's output.

Instead, you are going to write your output to the glob directory in your home directory.
Remember that your home directory can be represented by the '~' character.
This may save you a lot of typing.
The glob space is not backed up and is occasionally deleted, and is meant to be used for temporary storage.
(You could also write these files to your regular home directory space, but you may run out of space, and it is not good practice to keep large amounts of data in your home directory, so please do not do that.)

This creates some complexity, because your input data and your output data are not in the same place.
This is a common data processing problem, and one you should get used to dealing with.
It does mean that you'll need to type a lot.
There are a few ways to deal with this.

1. Remember where you are (your current working directory, `pwd`) and use relative or absolute paths as necessary to type the least.
This is a quick but sloppy solution, and error prone, but if you are doing things one time by hand, it works.
We all do it sometimes.
1. Use the full paths to everything, regardless of where you are actually working.
This is the most time consuming, and requires that you remember where everything is, but it is also the safest, because you always know that you are telling the computer exactly where you want to read and write.
This method is not dependent on keeping track of your current directory, because there are no relative paths, and you are much less likely to write data out to the wrong place by mistake.
Any time you get to the point of writing code or batch scripts to automate your data processing, you should do this.
For purposes of these exercises, it does not really matter which of these you do.
This is part of learning to work on the command line.
For purposes of example, the full paths will be given, but there will be examples where only the general syntax will be given, and you will have to find your data.

Also, remember that tab completion can be very helpful in typing paths to files, not just because it saves keystrokes but also because it validates that you have typed a valid path (if the file is not there, tab completion will not work).

So that we don't clutter up the top level of our globs and get in the way of later exercises, we will make a subdirectory in there

```bash
mkdir ~/glob/gatk
```

## Running commands

Throughout the exercises, we will illustrate commands on the format:  
```bash
command <parameter1> <parameter2> ...  
```
The convention is that you should replace &lt;parameter&gt; with your specific parameter, for example your input file name, output file name, directory name, etc.  
If you don't know what you should specify, please ask.  
We do this for two reasons.  
First, as you all work, not everyone will create files with exactly the same names, so there is no way to make standard instructions for everyone.  
Second, you need to learn how to figure out what goes into these spaces.
Usually, if you type a command without input parameters the documentation of the tool including possible input parameters will be displayed on the screen.

That brings us to copying and pasting.
It is possible to copy some of the commands out of this wiki and paste them into your terminal and make them work.
This is not recommended.
First, there can be formatting differences (especially how return characters are handled) between the browser and the terminal that make these commands not work properly.
Second, and more important, when you are doing this on your own data, there will be no cutting and pasting.
You will learn more by typing.
Remember that tab completion will help you with this.

NOTE - one you have typed the command for a step in the exercise below and seen that it works with your specifications, we recommend you to save the command in a plain text file. Plain text format will not cause formatting errors. You will perform the same procedure for at least two samples, so storing the commands in a text file will make the second analysis much faster.    

We will align our data to the reference using BWA, a popular aligner based on the Burrows-Wheeler transform.

## Indexing the reference genome

Before we can run BWA at all, we need a reference genome, and we need to perform the Burrows-Wheeler transform on the reference and build the associated files.
For our exercises, we'll use only human chromosome 17.
You can copy this from the project directory to your workspace.
(Normally copying references is a bad thing, but this is so that everyone can see the full BWA process.)

```bash
cp /sw/courses/ngsintro/gatk/refs/human_17_v37.fasta ~/glob/gatk
```

Check to see that this worked.

```bash
ls -l ~/glob/gatk
```

should show you:

```bash
-rw-r--r-- 1 mczody uppmax 82548517 Sep 23 21:44 human_17_v37.fasta
```

except with your username. The size of the file in bytes is the number showing after your username. 

If your file is not there or if it's the wrong size, something went wrong with your copy and you need to figure out what before you move on.
Checking the existence and size of files from each step in a process before performing the next step is a good practice that save a lot of time.
A common mistake people make is to attempt to load input files that do not exist or create output files where they cannot write.

Now we need to build the Burrows-Wheeler transform

```bash
bwa index -a bwtsw ~/glob/gatk/human_17_v37.fasta
```

BWA is a single program that takes a series of different commands as the first argument.
That command says to index the specified reference and use the bwtsw algorithm (BWA also has another indexing method for small genomes that we will not use).

This command will take about 2 minutes to run and should create 5 new files in your gatk directory with the same base name as the reference and different extensions.

While we're doing this, we will also build two different sequence dictionaries for the reference, which just lists the names and lengths of all the chromosomes.
Other programs will need these as input later and they are used to make sure the headers are correct.

```bash
samtools faidx ~/glob/gatk/human_17_v37.fasta
```

```bash
java -Xmx16g -jar /sw/apps/bioinfo/picard/1.69/kalkyl/CreateSequenceDictionary.jar R=~/glob/gatk/human_17_v37.fasta O=~/glob/gatk/human_17_v37.dict
```

## Mapping - Making Single Read Alignments for Each of the Reads in the Paired End Data

Running BWA for paired end data is done in multiple steps.
First we align each set of reads, then we combine the paired alignments together (which also includes a realignment step using a more sensitive algorithm for unplaced mates).
Let's start with one chunk of whole genome shotgun data from individual NA06984.

```bash
bwa aln ~/glob/gatk/human_17_v37.fasta /sw/courses/ngsintro/gatk/fastq/wgs/NA06984.ILLUMINA.low_coverage.17q_1.fq > ~/glob/gatk/NA06984.ILLUMINA.low_coverage.17q_1.sai
```

Note that if you have to use a file redirect ( &gt;) for your output, otherwise bwa will print the output directly to stdout, i.e. your screen. Which means that forgetting the redirect can be very disappointing.

While that's running, take a minute to look at the input file path.
This is a fastq file, so I put it in a directory called fastq.
It is from whole genome shotgun sequencing, so it is in a subdirectory called wgs.
The file name has 6 parts, separated by . or \_:

1. NA06984 - this is the individual name
1. ILLUMINA - these reads came from the Illumina platform
1. low_coverage - these are low coverage whole genome shotgun reads
1. 17q - I have sampled these reads from one region of 17q
1. 1 - these are the first reads in their paired sets
1. fq - this is a fastq file

Now we need to do this again for the second read file.
Everything is that same except with 2s instead of 1s.
Don't forget to change your output file also!

Before we go on to the next step, take a minute and look at the fastq files and understand the format and contents of these files.
Use

```bash
less
```

to read one of those .fq files in the project directory.

## Merging Alignments and Making SAM Files

The sai files are a binary format internal to BWA.
We now need to process those into something we can use.
For paired ends, this is done with the sampe function of BWA.
(Note that if you ever forget the syntax for a function, you can just type

```bash
bwa <function>
```

and it will list the parameters and options.
Run it for your files:

```bash
bwa sampe <ref> <sai1> <sai2> <fq1> <fq2> > ~/glob/gatk/<sample>.sam
```

The sampe function takes a lot of arguments.
It needs the reference and the reads, because the sai files just have the definitions of the alignments, not the sequences.
It needs the sai files to get the alignments.
It outputs a SAM format file.
I would suggest that you give it the same name prefix as the others, but if you are getting tired of typing that, pick something shorter.
Retain the sample name and the fact that it is the 17q low coverage data.

## Creating a BAM File

SAM files are nice, but bulky, so there is a compressed binary format, BAM.
We want to convert our SAM into BAM for everything that comes downstream.

Typically the BAM has the same name as the SAM but with the .sam extension replaced with .bam.

We need to add something called read groups which defines information about the sequencing run to our BAM file, because GATK is going to need this information.
Normally, you would do this one sequencing run at a time, but because of the way I downloaded these data from 1000 Genomes, our data are pulled from multiple runs and merged.
We will pretend that we just have one run for each sample, but on real data, you should not do this.

Now, we use the Picard package to add read group information.
However, it turns out that Picard is a very smart program, and we can start with the sam file and ask it to simultaneously add read groups, sort the file, and spit it out as BAM.
(It does, however, have a very awkward calling syntax.)

```bash
java -Xmx16g -jar /sw/apps/bioinfo/picard/1.69/kalkyl/AddOrReplaceReadGroups.jar INPUT=<sam file> OUTPUT=<bam file> SORT_ORDER=coordinate RGID=<sample>-id RGLB=<sample>-lib RGPL=ILLUMINA RGPU=<sample>-01 RGSM=<sample>
```

Note that the arguments to Picard tools are parsed (read by the computer) as single words, so it is important that there is no whitespace between the upper case keyword, the equals, and the value specified, and that you quote ('write like this') any arguments that contain whitespace.

We specify the input, the output (assumed to be BAM), the SORT_ORDER, meaning we want Picard to sort the reads according to their genome coordinates, and a lot of sample information.
The sample names for each of these 1000 Genomes runs is the Coriell identifier, the two letters and five numbers at the start of the file names (e.g., NA11932).
We're going to use this for all our read group information.

* RGID is the group ID. This is usually derived from the combination of the sample id and run id, or the SRA/EBI id. We will just add -id to the sample name.
* RGLB is the group library. This will come from your library construction process. You may have multiple read groups per library if you did multiple sequencing runs, but you should only have one library per read group. We will add -lib the sample name.
* RGPL is the platform. It is a restricted vocabulary. These reads are ILLUMINA.
* RGPU is the run identifier. It would normally be the barcode of your flowcell. You may have multiple read groups per run, but only one run per read group. We will just fake it as &lt;sample&gt;-01.
* RGSM is the sample name. Use the actual sample name. You can have multiple read groups, libraries, runs, and even platforms per sample, but you can only have one sample per read group. (If you are pooling samples without barcoding, there is no way to separate them later, so you should just designate the pool itself as a sample, but downstream analyses like SNP calling will be blind to that knowledge.) 

Lastly, we need to index this BAM, so that programs can randomly access the sorted data without reading the whole file.
This creates a file called &lt;input bam&gt;.bai, which contains the index.
You do not have to specify this because the index file always has the exact same name as the BAM except that it has .bai instead of the .bam extension.
This is how programs know to find the index associated with a BAM file.
If you manually mix these things up (like you change a BAM without changing its name and do not reindex it), you can cause problems for programs that expect them to be in sync.

```bash
java -Xmx16g -jar /sw/apps/bioinfo/picard/1.69/kalkyl/BuildBamIndex.jar INPUT=<bam file>
```

## Processing the BAM file with GATK

Now, we want to use the Genome Analysis Toolkit (GATK) to perform a couple of alignment and quality improvement steps, although on our data, they may not actually do much, due to the nature of the data and some of the shortcuts we have taken in identifying our read groups.

First, we'll realign locally around potential indels.
This is done in two steps.
First, we identify possible sites to realign:

```bash
java -Xmx16g -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -I <input bam file> -R <ref file> -T RealignerTargetCreator -o <intervals file>
```

The &lt;bam file&gt; should be your sorted and indexed BAM with read groups added from before.
Note that the option flag preceding the input bam is a capital I (as in Input), not a lower case l.
The &lt;ref file&gt; is the reference you used for alignment, and the &lt;intervals file&gt; is an output text file that will contain the regions GATK thinks should be realigned.
Give it the extension ".intervals".
Note that there is an additional option we are not using, which is to specify a list of known indels that might be present in the data (i.e., are known from other sequencing experiments).
Using this speeds up the process of identifying potential realignment sites, but because our data set is so small, we won't use that.

Now we feed our intervals file back into GATK with a different argument to actually do the realignments:

```bash
java -Xmx16g -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -I <input bam> -R <ref file> -T IndelRealigner -o <realigned bam> -targetIntervals <intervals file>
```

Note that we need to give it the intervals file we just made, and also specify a new output bam (&lt;realigned bam&gt;).
GATK is also clever and automatically indexes that bam for us (you can type ls and look at the list of files to verify this).

Next, we're going to go back to Picard and mark duplicate reads:

```bash
java -Xmx16g -jar /sw/apps/bioinfo/picard/1.69/kalkyl/MarkDuplicates.jar INPUT=<input bam> OUTPUT=<marked bam> METRICS_FILE=<metrics file>
```

Note that you need to feed it an &lt;input bam&gt;, which should be your realigned BAM from before, and you need to specify an output, the &lt;marked bam&gt; which will be a new file used in the following steps.
There is also a &lt;metrics file&gt;, which is a output text file.
We will take a look at that now.

Picard do not automatically index the .bam file so you need to do that before proceeding.

```bash
java -Xmx16g -jar /sw/apps/bioinfo/picard/1.69/kalkyl/BuildBamIndex.jar INPUT=<bam file>
```

Now we can look at the duplicates we marked with Picard, using a filter on the bit flag.
The mark for duplicates is the bit for 1024, so we can look at only duplicate marked reads with that.

```bash
samtools view -f 1024 <bam file> | less
```

If we just want a count of the marked reads, we can use the -c option.

```bash
samtools view -f 1024 -c <bam file>
```

Finally, we want to perform quality recalibration with GATK.
We do this last, because we want all the data to be as clean as possible when we get here.
This also happens in two steps.
First, we compute all the covariation of quality with various other factors:

```bash
java -Xmx64g -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T BaseRecalibrator -I <input bam> -R <ref file> -knownSites /sw/courses/ngsintro/gatk/ALL.chr17.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf -o <calibration table>
```

We need to feed it our bam file and our ref file.
We also need a list of known sites.
Otherwise, GATK will think all the real SNPs in our data are errors.
We're using calls from 1000 Genomes, which is a good plan for human (although a bit circular in our case).
If you are sequencing an organism with few known sites, you could try calling once and then using the most confident variants as known sites (which should remove most of the non-erroneous bases).
Failure to remove real SNPs from the recalibration will result in globally lower quality scores.
We also give it the name of a table file we want it to write out containing the covariation data.
We will take a look at this.
It will be used in the next step:

```bash
java -Xmx16g -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T PrintReads -BQSR <calibration table> -I <input bam> -R <ref file> -o <output bam>
```

The &lt;input bam&gt; in this step is the same as the last step, because we haven't changed it yet, but the &lt;output bam&gt; is new and will have the recalibrated qualities.
The &lt;calibration table&gt; is the file we created in the previous step.

## Variant Calling

Now we'll run the GATK HaplotypeCaller on our bam and output a gVCF file that will later be used for joint genotyping.

```bash
java -Xmx16g -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T HaplotypeCaller -R <ref file> -I <input bam> --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o <output>
```

The &lt;ref file&gt; is our old reference fasta again.
The &lt;input bam&gt; is the output from the recalibration step.
The output file is &lt;filename.g.vcf&gt;.
It needs to have a .g.vcf extension because it is a gvcf file.
The beginning part should be identifiable as associated with your bam file name (like the name root you use before the .bam) so you can tell later which vcf file came from which BAM). The --variant_index_type LINEAR --variant_index_parameter 128000 to set the correct index strategy for the output gVCF.

Rerun the mapping and variant calling steps for at least one more sample before continuing with the next step. 

# Joint genotyping

Now you will call variants on all the gvcf-files produced in the previous step by using the GenotypeGVCFs. This takes the output from the Haplotypecaller that was run on each sample to create raw SNP and indel VCFs.

```bash
java -Xmx16g -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T GenotypeGVCFs -R <ref file> --variant <sample1>.g.vcf --variant <sample2>.g.vcf ... -o <output>.vcf
```

As an alternative try also to run the same thing but with all the gvcf for all low_coverage files in the course directory. A gvcf file where these have been merged can be found in the course directory, /sw/courses/ngsintro/gatk/vcfs/ILLUMINA.low_coverage.17q.g.vcf. In the next step when viewing the data in IGV, look at both and try to see if there is a difference for a your sample. 

```bash
java -Xmx16g -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T GenotypeGVCFs -R <ref file> --variant /sw/courses/ngsintro/gatk/vcfs/ILLUMINA.low_coverage.17q.g.vcf -o <output>
```

## Filtering Variants

The last thing we will do is filter variants.
We do not have enough data that the VQSR technique for training filter thresholds on our data are likely to work, so instead we're just going to use the "best practices" parameters suggested by the GATK team (http://www.broadinstitute.org/gatk/guide/topic?name=best-practices).

The parameters are slightly different for SNPs and indels, but we have called ours together.
I would suggest trying both and seeing what you get.
Why do you think that some of these parameters are different between the two types of variants?

An example command line is:

```bash
java -Xmx16g -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T VariantFiltration -R <reference> -V <input vcf> -o <output vcf> --filterExpression "QD<2.0" --filterName QDfilter --filterExpression "MQ<40.0" --filterName MQfilter --filterExpression "FS>60.0" --filterName FSfilter
```

Note two things about this.
First, each filterName option has to immediately follow the filterExpression it matches.
This is an exception to the rule that options can come in any order.
However, the order of these pairs, or their placement relative to other arguments, can vary.
Second, the arguments to filterExpression are in quotation marks (").
Why is that?

If you want to run the indel filtering version, you can look on the web page above and get those value and substitute them.

Once you have the filtered calls, open your filtered VCF with less and page through it.
It has all the variant lines, still, but one of the fields that was blank before is now filled in, indicating that the variant on that line either passed filtering or was filtered out, with a list of the filters it failed.
Note also that the filters that were run are described in the header section.

## Look at Your Data with IGV

Next, we want to know how to look at these data.
For that, we will use IGV (Integrative Genomics Viewer).
We will launch IGV from our desktops because it runs faster that way.
Go to your browser window and Google search for IGV.
Find the downloads page.
You will be prompted for an email address.
If you have not already downloaded IGV from that email address, it will prompt you to fill in some information and agree to a license.
When you go back to your own lab, you can just type in your email and download the software again without agreeing to the license.

Now launch the viewer through webstart.
The 1.2 Gb version should be sufficient for our data.
It will take a minute or two to download IGV and start it up.
While that's going on, we need to download some data to our local machines so the viewer can find it (IGV can also look at web hosted data, but we are not going to set that up for our course data).
When it prompts you to save the IGV program, if you are working on a Mac put it in the Applications folder, otherwise just save it in your home directory. 

Open a new terminal or xterm _on your local machine_ (i.e., do not log in to uppmax again).
You should be in your home directory.
Now we're going to use the command scp (secure copy) to get some data copied down:

We will start with the merged bam files.
We want to get both the bams and bais for the low coverage and exome data.

```bash
scp <username>@milou.uppmax.uu.se:/sw/courses/ngsintro/gatk/processed/MERGED.illumina.\* ./
```

Because your uppmax user name is different than the user name on the local machine, you have to put your uppmax user name in front of the @ in the scp so that it knows you want to log in as your uppmax user, not as macuser.
After the colon, we give the path to the files we want.
The wildcard (*) character indicates that we want all the files that start with "MERGED.illumina".
However, in this case, we need to add a backslash ('\') in front of the wildcard ('*').
This is known as "escaping", because ordinarily your local shell would try to expand the wildcard in your local directory, but we want it expanded on the remote machine.
The './' means copy the files to your current directory.

It will prompt you for your uppmax password, then it should download four files.

We will also want to load the vcfs into IGV, so you can look at what calls got made.

```bash
scp <username>@milou.uppmax.uu.se:/sw/courses/ngsintro/gatk/vcfs/MERGED.illumina.\* ./
```

Do the same thing for the vcf that you have created in your home directory. 

By now, IGV should be launching.
The first thing we want to do is make sure we have the right reference.
In IGV, go to the popup menu in the upper left and set it to "Human 1kg (b37+decoy)".
This is the latest build of the human genome (also known as GRCh37).

Now, go under the Tools menu and selection "Run igvtools..." Change the command to "Count" and then use the Browse button next to the Input File line to select the bams (not the bai) that you just downloaded.
It will autofill the output file.
Now hit the Run button.
This generates a .tdf file for each bam.
This allows us to see the coverage value for our BAM file even at zoomed at views.
(We could also do this offline using a standalone version of igvtools.)

Now close the igvtools window and go back to the File menu, select "Load from File..." and select your bams (not the .bai or the .tdf).
They should appear in the tracks window.
Click on chromosome 17 to zoom in to there.
You can now navigate with the browser to look at some actual read data.
If you want to jump directly to the region we looked at, you can type MAPT in the text box at the top and hit return.
This will jump you to one of the genes in the region.

Let's look at a few features of IGV.

Go under the View menu and select Preferences.
Click on the Alignments tab.
There are a number of things we can configure.
Feel free to play with them.
Two important ones for our data are near the top.
Because we have multiple samples and the exome coverage is very deep, we want to turn off downsampling (upper left).
However, this will cause us to load more reads, so we want to reduce the visible range threshold (top).
I would suggest 5 kb.

Next, we want to look at some of the track features.
If you control-click (or right click for PCs or multi-button mice on Macs) on the track name at the left, you will get a popup menu with several options.
For example, click the gene track and play with the view (collapsed, squished, expanded).
I would suggest squished for our purposes.

If you go to a read alignment track, you can control some useful features of the display.
One is how you color the reads (by sample is an interesting one here).
Another is the grouping.
Group by sample is again useful (having grouped by sample, we could then use color for something else).

You can look at just the calls you made, or you can look at the calls from the full set, where you may see more of a difference between different types and depths of sequencing and between the calls with and without filtering. (IGV displays the filtered variant site in lighter shades, so you only need to load the filtered file).
You can even load these data all together.
Are there calls that were made using only one or two samples that were not made in the full data set or vice versa?

Try to browse around in your data and get a feeling for the called variants. Can you find a variant that have an allele frequency of exactly 0.5? A variant that have calls for all individuals? Exonic variants?


## [Extra labs](resequencing-extra)

If you have more time there are a couple of extra exercises where you will perform downstream analysis of the called variants in your .vcf file. [Extra labs](resequencing-extra)

## More info Quality Scores  
Here is a technical documentation of Illumina Quality Scores: [technote_Q-Scores.pdf](technote_Q-Scores.pdf) 
