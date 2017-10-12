---
layout: default
title:  'Instructions for the ChIP-seq Data Analysis Class'
---

# Instructions for the ChIP-seq Data Analysis Class

## Introduction

REST (NRSF) is a transcriptional repressor that represses neuronal genes in non-neuronal tissues. It is a member of the Kruppel-type zinc finger transcription factor family. It represses transcription by binding a DNA sequence element called the neuron-restrictive silencer element (NRSE). The protein is also found in undifferentiated neuronal progenitor cells and it is thought that this repressor may act as a master negative regular of neurogenesis. In addition, REST has been implicated as tumour suppressor, as the function of REST is lost in breast, colon and small cell lung cancers.

ChIP sequencing (ChIP-seq) is a method of identifying genome-wide occupancy patterns of proteins of interest (transcription factors, chromatin binding proteins, histones, DNA / RNA polymerases, etc.). A pressing question in ChIP-seq experiments is “Did my ChIP work?”, i.e. did the antibody-treatment enrich sufficiently so that the ChIP signal can be separated from the background signal? (After all, around 90% of all DNA fragments in a ChIP experiment represent the genomic background). This question is impossible to answer by counting number of peaks or by visual inspection of mapped reads in a genome browser; however, several quality control methods have been developed to address this task, as will be demonstrated in the first part of this lab. The second part of the lab deals with identification of binding sites and signal visualisation. The third part is a collection of bonus exercises which expand on techniques used in parts 1 and 2, as well as an example workflow for peak annotation and differential occupancy. You can select one or more (or all) of these additional exercises. Combined they present a complete workflow for detection of binding sites in ChIP-seq experiments.


## Data and Methods

### Data

Data you will use originates from the [ENCODE](www.encodeproject.org) project. It consists of duplicates of ChIP-seq of a transcription factor REST in several human cell lines and in vitro differentialed neural cells. The ChIP data contains matching input chromatin samples.The accession numbers are listed in the Table 1. Each sample accession number used in this exercise is listed in the Table 2 (Appendix).

The reads were mapped to the human genome assembly version hg19.


| No |  Accession  | Cell line | Description                                            |
| --- | ----------- | --------- | ------------------------------------------------------ |
| 1  | ENCSR000BMN |  HeLa     | adenocarcinoma (Homo sapiens, 31 year female)          |
| 2  | ENCSR000BOT |  HepG2    | hepatocellular carcinoma (Homo sapiens, 15 year male)  |
| 3  | ENCSR000BOZ |  SK-N-SH  | neuroblastoma (Homo sapiens, 4 year female)            |
| 4  | ENCSR000BTV |  neural   | in vitro differentiated (Homo sapiens, embryonic male) |

Table 1. ENCODE accession numbers for data sets used in this class.
 

In the interest of time, all steps are performed using a scaled down data set (containing reads mapped to chromosomes 1 and 2). However, for the post peak-calling QC and differential occupancy, the peaks on chromosomes 1 and 2 are selected from the peaks called using the full data set.


### Methods

After mapping the reads to the genome using a preferred short read aligner performing **ungapped (global)** alignment (bowtie in this example), the bam files are preprocessed to remove reads which could confound the subsequent analysis. The files you are working with are mapped by the ENCODE consortium, and only the reads with **one best alignment** are reported (in some sources these are referred to as "unique alignments" or "uniquely aligned reads"). This means that alignments of reads mapped to multiple locations in the genome are not present in the data. If other mapping strategy is used, such multi-mapping reads need to be removed from the data prior to analyses.

Initial steps of data analysis in a ChIP-seq experiment are focused on addressing two questions:

1. How successful was the ChIP-seq experiment? (quality control steps, part 1)

2. Identification of binding sites for the factor of interest. (peak calling steps, part 2)


The typical workflow in a ChIP-seq experiment (after an initial quality control, read preprocessing and read mapping to reference genome - steps not covered in this exercise)  includes:

1. Post-alignment quality control (QC) to assess successful enrichment in ChIP by computing cross-correlation profile and associated metrics; the tool of choice is phantompeakqualtools (part 1 of this class).

2. Processing of bam files to remove alignments of reads which 
(i) are duplicates and 
(ii) map to blacklisted regions ("hyper-chippable" regions). 
These reads may result from experimental artefacts and their presence may bias downstream analyses.
Additionally, in  this step you will compute normalised coverage tracks; these will be later viewed in a genome browser
(part 1 of this class).

3. Other post-alignment QC steps (part 3 of this class, additional exercises) include calculation of
	(i) cumulative enrichment
	(ii) sample clustering in genomic bins mode.

4. Identification of binding sites. Enriched regions aka peaks are identified using MACS2. After peak calling, several checks can be performed on the peak regions, for example selection of peaks consistently identified in both biological replicates. Subsequently, to assess similarity of binding profiles between samples, another round of sample clustering can be performed, this time of signal resulting from reads mapped to regions identified as the binding sites in all libraries.
(part 2 of this class)


5. Signal visualisation in a genome browser. The read alignments, normalised coverage tracks and binding sites are visualised for general signal assessment (part 2 of this class). 
Integrative genome browser [IGV](https://www.broadinstitute.org/igv/) is an easy to use tool compatible with all operating systems - please install it in your local computer if you haven't done so already.

6. Additional exercises (part 3 of this class):

	6.1. Computing cumulative enrichment as another metric to assess enrichment on ChIP-seq (deepTools).

	6.2. Computing correlation between libraries in bins mode: in this method the genome is divided into bins of specified size and reads mapped to each bin are counted; the resulting signal profiles are used to cluster libraries to identify groups of similar signal profile (deepTools).

	6.3. Alternative method for post-peak calling QC (R script, library ChIPQC).

	6.5. Peak annotation and differential occupancy (R script, libraries ChIPpeakAnno and DiffBind).

	6.6. Visualisation of the ChIP signal associated with genomic regions (deepTools).


Please note that all methods used in this exercise perform significantly better when used on complete (i.e. non-subset) data sets. Their accuracy most often scales with the number of mapped reads in each library, but so does the run time. As a reference, for some of the steps, plots generated analysing the complete data set are also presented in the Appendix of this document.

Last but not least, we have prepared intermediate files in case some step won't work; these will allow you to progress through the analysis. You will find them in the `/results` directory.


## Before you start...

First, you need to book a node on milou. We have reserved half a node for each student during this course. By now, you are probably already familiar with the procedure:

```bash
salloc -A g2016008 -t 08:00:00 -p core -n 8 --no-shell --reservation=g2016008_3 &
```

<font color='red'>Please do the node booking only once and make sure you do not have multiple reservations running at the same time, otherwise you will take away resources from other course participants! If you are disconnected from Uppmax, check with `squeue -u username if you still have your reservation running. If it is running, skip this step and connect to that reservation.</font>

There are many files which are part of the data set, and some additional files with annotations are required for the analysis, therefore saving the files in a structured manner is essential to keep track of the analysis steps. For the same reason, the paths to files are already listed in each command - you are of course welcome to change the directory structure and adjust the paths to reflect these changes.

Many commands are quite long, as there are many input files, and several parameters to set up; as a consequence the entire command spans several lines of text. The backslash character ("\\") is used to indicate to the interpreter that the command input will continue in the following line, and the command is not executed prematurely.

To see all options for applications used throughout the class, you can type `command -h` to view usage help.

You will notice that you will load and unload modules practically before and after each command; this is done because there are often dependency conflicts between the modules used in this exercise, and if not unloaded, some modules will cause error messages from the module system on milou.

To save time and minimise issues with data access, you will use bash scripts which will set up data and environment for this exercise. 

* `chiseq_data.sh` this script sets up directory structure and creates symbolic links to data or copies smaller files; needs to be sourced only once, before the exercise is started;
* `chipseq_env.sh` this script sets several environmental variables you will use in the exercise; it needs to be sourced every time you open a new connection to milou (i.e. after you accidentally log out or the connection times out).

To copy the scripts to your home directory and execute them:

```bash
cp /sw/courses/ngsintro/chipseq/scripts/setup/chipseq_data.sh ./

source chipseq_data.sh

source chipseq_env.sh
```

You should see a directory named "chipseq":

```bash
cd chipseq
cd analysis
```

...and you are ready to start.

[//]: #
[//]: # (x correlation)
[//]: #

## Part 1: Quality control and alignment processing

### Assessment of the enrichment in the ChIP-seq data by computing cross correlation

You will calculate cross correlation for REST ChIP-seq in HeLa cells. This section is performed using data subset to chromosomes 1 and 2.

```bash
mkdir xcor
cd xcor

module load phantompeakqualtools/1.1

run_spp.R -c=../../data/ENCFF000PED.chr12.bam -savp=hela1_xcor.pdf \
-out=xcor_metrics_hela.txt

module unload phantompeakqualtools/1.1
```

This step takes a few minutes, and phantompeakqualtools prints messages as it progresses through different stages of the analysis.

Inspect the resulting file `xcor_metrics_hela.txt`. 

The metrics file is tabulated, and the fields are:

```bash
COL1: Filename

COL2: numReads: effective sequencing depth i.e. total number of mapped reads
 in input file

COL3: estFragLen: comma separated strand cross-correlation peak(s) in decreasing 
 order of correlation. In almost all cases, the top (first) value in the list 
 represents the predominant fragment length.

COL4: corr_estFragLen: comma separated strand (Pearson) cross-correlation value(s)
 in decreasing order (col3 follows the same order)

COL5: phantomPeak: Read length/phantom peak strand shift

COL6: corr_phantomPeak: Correlation value at phantom peak

COL7: argmin_corr: strand shift at which cross-correlation is lowest

COL8: min_corr: minimum value of cross-correlation

COL9: Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8

COL10: Relative strand cross-correlation coefficient 
 (RSC) = (COL4 - COL8) / (COL6 - COL8)

COL11: QualityTag: Quality tag based on thresholded RSC 
  (codes: -2:veryLow; -1:Low; 0:Medium; 1:High; 2:veryHigh)
```

The columns to pay attention to are:

```bash
COL3: gives the fragment length as estimated from the data;

COL9: NSC; NSC>1.1 (higher values indicate more enrichment; 1 = no enrichment)

COL10: RSC; RSC>0.8 (0 = no signal; <1 low quality ChIP; >1 high enrichment

COL11: Quality tag based on thresholded RSC
(codes: -2:veryLow,-1:Low,0:Medium,1:High; 2:veryHigh)
```

For comparison, the cross correlation metrics computed for the entire data set using non-subset data are available by:

```bash
cat ../../results/xcor/rest.xcor_metrics.txt
```

In addition to inspecting the metrics file, it is always recommended to view the plots, as the shape of the cross correlation is in some cases more informative than just numbers.

Compare the plot `hela1_xcor.pdf` (cross correlation of the first replicate of REST ChIP in HeLa cells, using subset data) with cross correlation computed using the non-subset data set presented on figures 1 - 3. Compare to the ChIP using the same antibody performed in HepG2 cells (figures 4 - 6). 

If you have enabled X-forwarding on your local computer, this command opens pdf files directly on the server:

```bash
evince hela1_xcor.pdf &
```

If this above command does not work for you (there may be problems, depending on the configuration of your local computer, and connection to milou), you can copy the pdf files to your local computer and open using your preferred pdf viewer (recommended). Using terminal, cd to the desired destination directory in your local system and substitute the text in CAPS with appropriate names:


```bash
scp <LOGIN>@milou.uppmax.uu.se:~/chipseq/analysis/<FOLDER_NAME>/*pdf ./
```

-----

|Figure 1. <br> HeLa, REST ChIP  <br>  replicate 1, QScore:2 | Figure 2. <br> HeLa, REST ChIP <br> replicate 2, QScore:2  | Figure 3. <br> HeLa, input <br> QScore:-1                                         |
| --- | ----------- | --------- |
|<img src="files/chipseq_fig/ENCFF000PEDxcorrelationplot.png" alt="xcor, hela chip" style="width: 200px;"/>| <img src="files/chipseq_fig/ENCFF000PEExcorrelationplot.png" alt="xcor, hela chip" style="width: 200px;"/>| <img src="files/chipseq_fig/ENCFF000PETxcorrelationplot.png" alt="xcor, hela input" style="width: 200px;"/>|


----


|Figure 4. <br> HepG2, REST ChIP  <br>  replicate 1, QScore:0 | Figure 5. <br> HepG2, REST ChIP <br> replicate 2, QScore:1  | Figure 6. <br> HepG2, input <br> QScore:0                                        |
| --- | ----------- | --------- |
|<img src="files/chipseq_fig/ENCFF000PMGppqtxcorrelationplot.png" alt="xcor, hepg2 chip" style="width: 200px;"/>| <img src="files/chipseq_fig/ENCFF000PMJppqtxcorrelationplot.png" alt="xcor, hepg2 chip" style="width: 200px;"/>| <img src="files/chipseq_fig/ENCFF000POMppqtxcorrelationplot.png" alt="xcor, hepg2 input" style="width: 200px;"/>|

-----

The usual questions to ask after performing this step are:
How would you rate these particular two data sets? Are all samples of good enough quality? Which data set would you rate higher in terms of how successful the ChIP was? Would any of the samples fail this QC step? Why? In short: did the ChIP(-seq) experiment work?

[//]: #
[//]: # (bam processing)
[//]: #


### Alignment Processing

This section is performed using data subset to chromosomes 1 and 2. First, duplicated reads are marked and removed using MarkDuplicates tool from Picard. Marking as "duplicates" is based on their alignment location, not sequence. 

```bash
cd ..
mkdir bam_preproc
cd bam_preproc

module load samtools/1.1
module load java/sun_jdk1.8.0_40
module load picard/1.141

java -Xmx64G -jar $PICARD_HOME/picard.jar MarkDuplicates \
I=../../data/ENCFF000PED.chr12.bam O=ENCFF000PED.chr12.rmdup.bam \
M=dedup_metrics.txt VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=true ASSUME_SORTED=true
```

Second, reads mapped to [ENCODE blacklisted regions](https://sites.google.com/site/anshulkundaje/projects/blacklists) are removed.

```bash
bamutils filter ENCFF000PED.chr12.rmdup.bam \
ENCFF000PED.chr12.rmdup.filt.bam \
-excludebed ../../hg19/wgEncodeDacMapabilityConsensusExcludable.bed nostrand
```

Finally, the processed bam files are sorted and indexed:

```bash
samtools sort -T sort_tempdir -o ENCFF000PED.chr12.rmdup.filt.sort.bam \
ENCFF000PED.chr12.rmdup.filt.bam

samtools index ENCFF000PED.chr12.rmdup.filt.sort.bam

module unload samtools/1.1
module unload java/sun_jdk1.8.0_40
module unload picard/1.141
```


[//]: #
[//]: # (coverage)
[//]: #


### Computing read coverage

In this section you will compute the read coverage normalised to 1x coverage using tool bamCoverage from deepTools. This procedure is useful for comparing libraries sequenced to a different depth when viewing them in a genome browser such as IGV. It will be used later in section Visualisation. This section is performed using data subset to chromosomes 1 and 2;  hence the effective genome size used is 492449994 (4.9e8) (for hg19 the effective genome size is 2.45e9 (see [publication](http://www.nature.com/nbt/journal/v27/n1/fig_tab/nbt.1518_T1.html)). The reads are extended to 110 nt (the fragment length obtained from the cross correlation computation). The reads are summarised in 50 bp bins. No smoothing is applied at this step. In this step and several subsequent steps you will use a tool available as part of [deepTools](http://deeptools.readthedocs.org/en/latest/) a suite of tools developed for ChIP-seq analysis and visualisation.

```bash
module load deepTools/2.0.1

bamCoverage --bam ENCFF000PED.chr12.rmdup.filt.sort.bam \
 --outFileName ENCFF000PED.chr12.cov.norm1x.bedgraph \
 --normalizeTo1x 492449994 --extendReads 110 --binSize 50 \
 --outFileFormat bedgraph

module unload deepTools/2.0.1
```

[//]: #
[//]: # (Peak calling)
[//]: #

## Part 2: Identification of binding sites
You will identify peaks in the ChIPseq data using MACS2. MACS2 is one of the most popular peak callers, and it performs very well on data sets with good enrichment of transcription factors ChIP-seq. Peaks should be called on each replicate separately (do not pool the replicates) and peaks identified consistently across replicates should be used in downstream analyses of differential occupancy, annotation etc.

As before, this section is performed using data subset to chromosomes 1 and 2; hence the effective genome size used is 492449994 (4.9e8). To avoid long paths in the command line, you will link the necessary bam files with chip and input data. You will call peaks only in one ChIP-seq library; the rest of the work is already done, and the peaks are in the directory `/results/peaks_macs` (for HeLa cells only) and in bed format `/results/peaks_bed`.

```bash
cd ..
mkdir peak_calling
cd peak_calling

ln -s /sw/courses/ngsintro/chipseq/data/bam/hela/ENCFF000PED.chr12.rmdup.sort.bam \
./ENCFF000PED.preproc.bam
ln -s /sw/courses/ngsintro/chipseq/data/bam/hela/ENCFF000PET.chr12.rmdup.sort.bam \
./ENCFF000PET.preproc.bam
```



### Peak Calling

There are several parameters which affect peak calling itself, as well as result reporting. It is important to understand them to modify the command to the needs of your data set.

Meaning of the parameters:

-t, -c, -f, -n denote treatment, control, file format and output file names, respectively

-g is the genome size, in this case it is already encoded in MACS:
-g hs   =  -g 2.7e9;
-g mm   =  -g 1.87e9;
-g ce   =  -g 9e7;
-g dm   =  -g 1.2e8.

-q 0.01 is the q value (FDR) cutoff for reporting peaks; this is recommended over reporting using raw (undajusted) p value.

In our case, you use -g  4.9e8 because for peak calling you use data from chromosomes 1 and 2 only (total size is 4.9e8 bp).

You will see the progress of the analysis printed in the terminal (MACS prints messages as it progresses through different stages of the process). This step may take more than 10 minutes.


```bash
module load MACS/2.1.0

macs2 callpeak -t ENCFF000PED.preproc.bam -c ENCFF000PET.preproc.bam \
-f BAM -g 4.9e8 -n hela_1_REST.chr12.macs2 -q 0.01

module unload MACS/2.1.0
module unload python/2.7.6
```

The output of a MACS2 run consists of several files. You can inspect their content using

```bash
head -n 50 filename
```

You will use the `narrowPeak` files in the subsequent parts. These files are in [bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format, which is one of the popular file format used in genomics. It is used to store information on genomic ranges (such as ChIP-seq peaks, gene models, transcription starts sites, etc). Bed files can be also used for visualisation in genome browsers, including the popular [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks) and [IGV](https://www.broadinstitute.org/igv), which you will use later in part "Visualisation".


Files used in this step are derived from the `*.narrowPeak` files by selecting relevant columns, for example:

```bash
cut -f 1-3 hela_1_REST.chr12.macs2_peaks.narrowPeak >hela_1_chr12_peaks.bed
```


Peaks detected on chromosomes 1 and 2 are present in directory `/results/peaks_bed`. These peaks were detected using the non-subset data, and therefore there may be differences between the peaks present in the prepared `hela_1_peaks.bed` file compared to the peaks you have just detected. Use these pre-made peak bed files instead of the file you have just created. You can check how many peaks were detected in each library by listing number of lines in each file:

```bash
wc -l ../../results/peaks_bed/*.bed
```

```bash
cp ../../results/peaks_bed/*.bed ./
```


By checking for overlaps in the peak lists from different libraries you can detect peaks present in both libraries. This will give you an idea which peaks are reproducible between replicates. 
[Bedtools](http://bedtools.readthedocs.org/en/latest/) is a suite of utilities developed for manipulation of bed files.

In the command used here, the arguments are

-a, -b : two files to be intersected

-f 0.50 : fraction of the overlap between features in each file to be reported as an overlap

-r : reciprocal overlap fraction required


Select two replicates of the same condition, for example:

```bash
module load BEDTools/2.25.0

bedtools intersect -a hela_1_peaks.chr12.bed -b hela_2_peaks.chr12.bed -f 0.50 -r \
> peaks_hela.chr12.bed

wc -l peaks_hela.chr12.bed
```

This way you can compare peaks from replicates of the same condition, and peaks present in different conditions. You will need to create files with peaks common to replicates for the cell types you wish to compare (as in the example above). You can inspect which peaks were reproducibly found in two different cell lines, for example:

```bash
bedtools intersect -a peaks_hepg2.chr12.bed -b peaks_hela.chr12.bed -f 0.50 -r \
> peaks_hepg2_hela.chr12.bed

wc -l peaks_hepg2_hela.chr12.bed
```

When you have done all intersections you were interested in, unload the Bedtools module:

```bash
module unload BEDTools/2.25.0
```


Now you will generate a merged list of all peaks detected in the experiment, which will be needed in the next step. In this list, all overlapping intervals are merged, and a list of non-overlapping intervals is created. To do this, you will use another suite of tools for manipulation of bed files, [BEDOPS](http://bedops.readthedocs.org/en/latest/). Files used in this step are derived from the `*.narrowPeak` files by selecting relevant columns, as before.

These files are already prepared and present in the directory `peak_calling` (your current directory, if you follow the directory structure in the exercise).

```bash
module load BEDOPS/2.4.3

bedops -m hela_1_peaks.chr12.bed hela_2_peaks.chr12.bed hepg2_1_peaks.chr12.bed hepg2_2_peaks.chr12.bed \
neural_1_peaks.chr12.bed neural_2_peaks.chr12.bed sknsh_1_peaks.chr12.bed sknsh_2_peaks.chr12.bed \
>REST_peaks.chr12.bed

module unload BEDOPS/2.4.3

wc -l REST_peaks.chr12.bed 
```

In case things go wrong at this stage, you can find the merged list of all peaks in the `/results` directory. Link the file to your current directory


```bash
ln -s ../../results/peaks_bed/rest_peaks.chr12.bed ./rest_peaks.chr12.bed
```

[//]: #
[//]: # (Clustering in peaks)
[//]: #

### Quality Control After Peak Calling: Are the peaks consistent between replicates?

This section is performed using data subset to chromosomes 1 and 2. For the subsequent parts, you will use peaks which were called in the non-subset data set, and peaks located on chromosomes 1 and 2 were selected.

You will use multiBamSummary and plotCorrelation in the BED-file mode (deepTools). You will provide a list of intervals where reads are counted; this list is a merged list of all peaks detected in all libraries you have created in the previous step. You will begin by making a new directory and creating links to the data there.


```bash
cd ..
mkdir plots
cd plots

mkdir hela
mkdir hepg2
mkdir sknsh
mkdir neural
ln -s /sw/courses/ngsintro/chipseq/data/bam/hela/* ./hela
ln -s /sw/courses/ngsintro/chipseq/data/bam/hepg2/* ./hepg2
ln -s /sw/courses/ngsintro/chipseq/data/bam/sknsh/* ./sknsh
ln -s /sw/courses/ngsintro/chipseq/data/bam/neural/* ./neural

module load deepTools/2.0.1

multiBamSummary BED-file --BED ../peak_calling/REST_peaks.chr12.bed --bamfiles \
hela/ENCFF000PED.chr12.rmdup.sort.bam \
hela/ENCFF000PEE.chr12.rmdup.sort.bam hela/ENCFF000PET.chr12.rmdup.sort.bam \
hepg2/ENCFF000PMG.chr12.rmdup.sort.bam hepg2/ENCFF000PMJ.chr12.rmdup.sort.bam \
hepg2/ENCFF000POM.chr12.rmdup.sort.bam hepg2/ENCFF000PON.chr12.rmdup.sort.bam \
neural/ENCFF000OWM.chr12.rmdup.sort.bam neural/ENCFF000OWQ.chr12.rmdup.sort.bam \
neural/ENCFF000OXB.chr12.rmdup.sort.bam neural/ENCFF000OXE.chr12.rmdup.sort.bam \
sknsh/ENCFF000RAG.chr12.rmdup.sort.bam sknsh/ENCFF000RAH.chr12.rmdup.sort.bam \
sknsh/ENCFF000RBT.chr12.rmdup.sort.bam sknsh/ENCFF000RBU.chr12.rmdup.sort.bam \
--outFileName multiBamArray_bed_bam_chr12.npz \
--extendReads=110 \
--labels hela_1 hela_2 hela_i hepg2_1 hepg2_2 hepg2_i1 hepg2_i2 neural_1 \
neural_2 neural_i1 neural_i2 sknsh_1 sknsh_2 sknsh_i1 sknsh_i2 

plotCorrelation --corData multiBamArray_bed_bam_chr12.npz \
--plotFile correlation_peaks.pdf --outFileCorMatrix correlation_peaks_matrix.txt \
--whatToPlot heatmap --corMethod pearson --plotNumbers --removeOutliers

module unload deepTools/2.0.1
```

In the pattern of library clustering consistent with biological differences between the samples and the cell lines? 

[//]: #
[//]: # (IGV)
[//]: #

### Visualisation of mapped reads, coverage profiles and peaks in a genome browser

This part of the exercise requires installation of Interactive Genome Browser [IGV](https://www.broadinstitute.org/igv/) and it should be performed on your own computer (IGV is much faster when run locally). In fact, you are already familiar with IGV, as you used it during the variant calling class. You will view the following files (copy them from Uppmax):

REST ChIP in HeLa cells, replicate 1 (the files you created so far, or):

```bash
../data/bam/hela/ENCFF000PED.chr12.rmdup.sort.bam
../data/bam/hela/ENCFF000PED.chr12.rmdup.sort.bam.bai

../results/coverage/ENCFF000PED.cov.norm1x.bedgraph

../results/peaks_macs/hela_1_REST.chr12.macs2_peaks.narrowPeak 
```

input in HeLa cells:

```bash
../data/bam/hela/ENCFF000PET.chr12.rmdup.sort.bam
../data/bam/hela/ENCFF000PET.chr12.rmdup.sort.bam.bai

../results/coverage/ENCFF000PEE.cov.norm1x.bedgraph
```

Remember setting the reference genome to hg19 - this is the assembly the reads were mapped to.

In IGV, in menu "File": "Load from file", select files you want to visualise. Load the selected tracks. Click on chromosome 1 or 2, zoom in and explore! Go to one of the locations you found interesting (for example, of the REST binding peaks detected in both HeLa samples), available in the file peaks_hela.bed, which you generated earlier.

You can change the signal display mode in the tracks in the left hand side panel. Right click in the bam file track, select from the menu "display" - squishy; "color by" - read strand and "group by" - read strand.

Is the read distribution in the peaks (bam file tracks) consistent with the expected bimodal distribution? Are the detected peaks associated with annotated genes?


[//]: #
[//]: # (Additional exercises)
[//]: #


## Part 3: Additional analyses
The following parts are optional, and can be preformed independently of one another.

The parts using R can be performed either by executing the scripts provided with the exercise from the command line, or by typing in all commands in an R console. The latter is recommended for people who have some previous exposure to the R environment. Post-peak calling QC is performed using the R / Bioconductor package ChIPQC, and some steps are redundant with the steps already performed; it is an alternative to the already presented workflow. Differential Occupancy and Peak Annotation sections present examples using R packages developed specifically for analysis of ChIP-seq data (there are many more useful packages in Bioconductor).

There are examples of two more QC plots generated using deepTools which allow you to inspect enrichment in the libraries (cumulative enrichment aka bamFingerprint) and assess overall signal similarities between libraries (correlation plot in the "bin" mode).

Finally, you can visualise the distribution of the ChIP signal in relation to genomic features using deepTools. 


### Quality Control: cumulative enrichment

The cumulative enrichment (aka "bam Fingerprint") will be computed for the HeLa REST ChIP and corresponding input samples.
This section is performed using data subset to chromosomes 1 and 2. You will use the processed file you have created earlier in `/analysis/bam_preproc`, and two other libraries from the same data set, already prepared.

```bash
module load deepTools/2.0.1

plotFingerprint --bamfiles ENCFF000PED.chr12.rmdup.filt.sort.bam \
../../data/bam/hela/ENCFF000PEE.chr12.rmdup.sort.bam  \
../../data/bam/hela/ENCFF000PET.chr12.rmdup.sort.bam  \
 --extendReads 110  --binSize=1000 --plotFile HeLa.fingerprint.pdf \
--labels HeLa_rep1 HeLa_rep2 HeLa_input

module unload deepTools/2.0.1
```

Inspect the plot; does it indicate good sample quality (i.e. enrichment in ChIP samples and no indications of enrichment in input samples)? 
For comparison, similar plots generated for other samples used in this exercise are presented in figures 7 and 8.


|Figure 7. <br> Cumulative enrichment for REST ChIP and corresponding inputs    <br> in HepG2 cells | Figure 8. <br> Cumulative enrichment for REST ChIP and corresponding inputs    <br> in SK-N-SH cells |
| --- | ----------- |
|<img src="files/chipseq_fig/hepg2fingerprint.png" alt="fingerprint, hepg2" style="width: 280px;"/>| <img src="files/chipseq_fig/sknshfingerprint.png" alt="fingerprint, sknsh" style="width: 280px;"/>|

Can you tell which samples are ChIP and which are input? Are the cumulative enrichment plots in agreement with the cross-correlation metrics computed earlier?

### Sample clustering

To assess overall similarity between libraries from different samples and data sets, you will compute sample clustering heatmaps using multiBamSummary and plotCorrelation in bins mode from deepTools. In this method the genome is divided into bins of specified size (--binSize parameter) and reads mapped to each bin are counted; the resulting signal profiles are used to cluster libraries to identify groups of similar signal profile.

This section is performed using data subset to chromosomes 1 and 2.

First, to avoid very long paths in the command line, create subdirectories and link preprocessed bam files:

```bash
mkdir hela
mkdir hepg2
mkdir sknsh
mkdir neural
ln -s /sw/courses/ngsintro/chipseq/data/bam/hela/* ./hela
ln -s /sw/courses/ngsintro/chipseq/data/bam/hepg2/* ./hepg2
ln -s /sw/courses/ngsintro/chipseq/data/bam/sknsh/* ./sknsh
ln -s /sw/courses/ngsintro/chipseq/data/bam/neural/* ./neural

```


```bash
module load deepTools/2.0.1

multiBamSummary bins --bamfiles hela/ENCFF000PED.chr12.rmdup.sort.bam \
hela/ENCFF000PEE.chr12.rmdup.sort.bam hela/ENCFF000PET.chr12.rmdup.sort.bam \
hepg2/ENCFF000PMG.chr12.rmdup.sort.bam hepg2/ENCFF000PMJ.chr12.rmdup.sort.bam \
hepg2/ENCFF000POM.chr12.rmdup.sort.bam hepg2/ENCFF000PON.chr12.rmdup.sort.bam \
neural/ENCFF000OWM.chr12.rmdup.sort.bam neural/ENCFF000OWQ.chr12.rmdup.sort.bam \
neural/ENCFF000OXB.chr12.rmdup.sort.bam neural/ENCFF000OXE.chr12.rmdup.sort.bam \
sknsh/ENCFF000RAG.chr12.rmdup.sort.bam sknsh/ENCFF000RAH.chr12.rmdup.sort.bam \
sknsh/ENCFF000RBT.chr12.rmdup.sort.bam sknsh/ENCFF000RBU.chr12.rmdup.sort.bam \
 --outFileName multiBamArray_dT201_preproc_bam_chr12.npz --binSize=5000 \
--extendReads=110 --labels hela_1 hela_2 hela_i hepg2_1 hepg2_2 hepg2_i1 hepg2_i2 \
neural_1 neural_2 neural_i1 neural_i2 sknsh_1 sknsh_2 sknsh_i1 sknsh_i2 

plotCorrelation --corData multiBamArray_dT201_preproc_bam_chr12.npz \
--plotFile REST_bam_correlation_bin.pdf --outFileCorMatrix corr_matrix_bin.txt \
--whatToPlot heatmap --corMethod spearman

module unload deepTools/2.0.1
```

Inspect the resulting pdf. How does it differ from the clustering performed on signal in peaks?


### Alternative Quality Control Workflow in R

Before you start using R, you need to set the environmental variable that holds information on where R libraries are installed. In this exercise you will use libraries installed in the class home directory (simply because it takes quite some time to install all dependencies). The path should already be set, inspect it by:

```bash
echo $R_LIBS
```

The result should be `/sw/courses/ngsintro/chipseq/software/R_lib`. If it is not, source the `chipseq_env.sh` script.


You will start in directory `/analysis/R`. The file `REST_samples.txt` contains information on files location, and the paths are given in relation to `/analysis/R`; if you choose to start in another directory, please modify the paths in `REST_samples.txt`. This script takes a while to run.

```bash
cd ../analysis/R
Rscript chipqc.R
```

Inspect the html output of this script.

### Differential Occupancy and Peak Annotation in R

Before you start using R, you need to set the environmental variable that holds information on where R libraries are installed. In this exercise you will use libraries installed in the class home directory (simply because it takes quite some time to install all dependencies). The path should already be set, inspect it by:

```bash
echo $R_LIBS
```

The result should be `/sw/courses/ngsintro/chipseq/software/R_lib`. If it is not, source the `chipseq_env.sh` script.

Please note that normally three biological replicates are required for statistical analysis of factor occupancy. There are only two replicates each in the ENCODE data sets used in this class - hence you use duplicates for demonstration sake. This script uses the same sample information file as the previous one; the paths are all relative to the `/analysis/R directory`.

```bash
cd ../analysis/R
Rscript diffbind_annot.R
```

People familiar with R can modify and execute the commands from within the R terminal, selecting different contrasts of interest, for example.

### Signal visualisation using deepTools

You will visualise ChIP signal in relation to annotated transcription start sites (TSS) on chromosomes 1 and 2. A description of all visualisation options is given at [deepTools](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html). Create a separate directory in `/analysis`; cd to it. Check if all the paths to create links are correct for the location of your directory. For details on the options of the applications used, please consult the documentation available at [computeMatrix](http://deeptools.readthedocs.org/en/latest/content/tools/computeMatrix.html) and [plotHeatmap](http://deeptools.readthedocs.org/en/latest/content/tools/plotHeatmap.html).

First you will compute the matrix of values using computeMatrix. This program takes [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) files as input; you will need to convert bedgraph to bigWig using UCSC utilities:

```bash
cd analysis/
mkdir vis
cd vis

cp ../../hg19/chrom.sizes.hg19 chrom.sizes.hg19
cp ../bam_preproc/ENCFF000PED.chr12.cov.norm1x.bedgraph ./

module load ucsc-utilities/v287

bedGraphToBigWig ENCFF000PED.chr12.cov.norm1x.bedgraph chrom.sizes.hg19 hela_1.bw

module unload ucsc-utilities/v287
```

You are now ready to compute the matrix of scores for visualisation. You will need a bed file with positions of TSS; you can copy it to your current directory.

```bash
cp ../../hg19/refGene_hg19_TSS_chr12.bed ./

module load deepTools/2.0.1

computeMatrix reference-point -S hela_1.bw \
-R refGene_hg19_TSS_chr12.bed -b 5000 -a 5000 \
--outFileName matrix.tss.dat --outFileNameMatrix matrix.tss.txt \
--referencePoint=TSS --numberOfProcessors=max
```

Having the matrix of scores ready, you can now plot the binding profile around TSS and the heatmap:

```bash
plotHeatmap --matrixFile matrix.tss.dat \
--outFileName tss.hela_1.pdf \
--sortRegions descend --sortUsing mean

module unload deepTools/2.0.1
```

## Concluding remarks

The workflow presented in this exercise is similar to a typical one used for analysis of ChIP-seq data. There are more types of analyses you can do, which were not discussed here. One typical task is to identify short sequence motifs enriched in the regions bound by the assayed factor (peaks). There are several tools available, and I recommend testing at least two tools for your data.

[Homer](http://homer.salk.edu/homer/)

[GEM](http://groups.csail.mit.edu/cgs/gem/)

[RSAT](http://floresta.eead.csic.es/rsat/peak-motifs_form.cgi)

[MEME](http://meme-suite.org/)

## Appendix

### Data used in class


| No |  Accession  | Cell line | Replicate |    Input    |
| --- | ----------- | --------- | --------- | ----------- |
| 1  | ENCFF000PED | HeLa      | 1         | ENCFF000PET |
| 2  | ENCFF000PEE | HeLa      | 2         | ENCFF000PET |
| 3  | ENCFF000PMG | HepG2     | 1         | ENCFF000POM |
| 4  | ENCFF000PMJ | HepG2     | 2         | ENCFF000PON |
| 5  | ENCFF000OWQ | neural    | 1         | ENCFF000OXB |
| 6  | ENCFF000OWM | neural    | 2         | ENCFF000OXE |
| 7  | ENCFF000RAG | SK-N-SH   | 1         | ENCFF000RBT |
| 8  | ENCFF000RAH | SK-N-SH   | 2         | ENCFF000RBU |

Table 2. ENCODE accession numbers for samples used in this exercise.


### Figures generated during class

<img src="files/chipseq_fig/resENCFF000PEDchr12xcor.png" alt="" style="width: 400px;"/><br>

Figure 9. Cross correlation plot for REST ChIP in Hela cells, replicate 1, chromosome 1 and 2


----

<img src="files/chipseq_fig/peaksbedchr12pears.png" alt="" style="width: 400px;"/><br>

Figure 10. Sample clustering (pearson) by reads mapped in merged peaks; only chromosomes 1 and 2 included


----

<img src="files/chipseq_fig/resHelaChr12Fingerprint.png" alt="" style="width: 400px;"/><br>

Figure 11. Fingerprint plot for REST ChIP in Hela cells, replicate 1, chromosome 1 and 2


----

<img src="files/chipseq_fig/bin5kchr12spear.png" alt="" style="width: 400px;"/><br>


Figure 12. Sample clustering (spearman) by reads mapped in bins genome-wide; only chromosomes 1 and 2 included


----

<img src="files/chipseq_fig/resHelaProfileTSS.png" alt="" style="height: 400px;"/><br>

Figure 13. Binding profile in HeLa replicate 1, centered on TSS; data subset to chromosome 1 and 2


----

### Figures generated using the full (i.e. not subset) data set


<img src="files/chipseq_fig/helaprocfingerprint.png" alt="" style="width: 400px;"/><br>

Figure 14. Cumulative enrichment in  HeLa replicate 1, aka bam fingerprint


----

<img src="files/chipseq_fig/bin5kspear.png" alt="" style="width: 400px;"/><br>

Figure 15. Sample clustering (spearman) by reads mapped in bins genome-wide


----

<img src="files/chipseq_fig/peaksbedpears.png" alt="" style="width: 400px;"/><br>

Figure 16. Sample clustering (pearson) by reads mapped in merged peaks


----




