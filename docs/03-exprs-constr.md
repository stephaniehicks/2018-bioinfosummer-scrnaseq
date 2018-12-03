---
output: html_document
---

# Construction of expression matrix



Many analyses of scRNA-seq data take as their starting point 
an __expression matrix__. By convention, the each row of the 
expression matrix represents a gene and each column represents
a cell (although some authors use the transpose). Each entry 
represents the expression level of a particular gene in a 
given cell. The units by which the expression is meassured 
depends on the protocol and the normalization strategy used.

In this section, we will describe how to go from a set of 
raw sequencing reads to an expression matrix. 

## Raw sequencing reads QC

### FastQ

The output from a scRNA-seq experiment is a large collection of
cDNA reads. FastQ is the most raw form of scRNA-seq data you 
will encounter. The first step is to ensure that the reads are of 
high quality. 

All scRNA-seq protocols are sequenced with paired-end sequencing. 
Barcode sequences may occur in one or both reads depending on the 
protocol employed. However, protocols using unique molecular
identifiers (UMIs) will generally contain one read with the cell
and UMI barcodes plus adapters but without any transcript sequence. 
Thus reads will be mapped as if they are single-end sequenced
despite actually being paired end. 

### Check the quality of the raw sequencing reads

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
is a good tool for this. This software produces a FastQC report than 
be used to evaluate questions such as 
_How good quality are the reads?_ 
_Is there anything we should be concerned about?_, etc.

### Trim sequencing adapaters from sequencing reads

For example, [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) 
is wrapper for trimming reads from FastQ files using
[cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html). Once 
you trim the adapaters, you can re-run FastQC to confirm successful removal. 


## Reads alignment

After trimming low quality bases from the reads, the remaining sequences can
be mapped to a reference genome. Again, there is no need for a special purpose
method for this, so we can use the
[STAR](https://github.com/alexdobin/STAR) or the [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml) aligner. For large full-transcript datasets from well annotated organisms (e.g. mouse, human) pseudo-alignment methods (e.g. [Kallisto](https://pachterlab.github.io/kallisto/), [Salmon](http://salmon.readthedocs.io/en/latest/salmon.html)) may out-perform conventional alignment. For drop-seq based datasets with tens- or hundreds of thousands of reads pseudoaligners become more appealing since their run-time can be several orders of magnitude less than traditional aligners.

### Genome (`FASTA`, `GTF`)

To map your reads you will need the reference genome and in many cases 
the genome annotation file (in either GTF or GFF format). These can be 
downloaded for model organisms from any of the main genomics databases: 
[Ensembl](http://www.ensembl.org/info/data/ftp/index.html), 
[NCBI](ftp://ftp.ncbi.nih.gov/genomes/), or [UCSC Genome Browser](http://hgdownload.soe.ucsc.edu/downloads.html). 

GTF files contain annotations of genes, transcripts, and exons. They must contain: 

1. seqname : chromosome/scaffold 
2. source : where this annotation came from
3. feature : what kind of feature is this? (e.g. gene, transcript, exon)
4. start : start position (bp)
5. end : end position (bp)
6. score : a number
7. strand : + (forward) or - (reverse)
8. frame : if CDS indicates which base is the first base of the first codon (0 = first base, 1 = second base, etc..)
9. attribute : semicolon-separated list of tag-value pairs of extra information (e.g. names/IDs, biotype)

Empty fields are marked with "."

In our experience: 

* Ensembl is the easiest of these to use, and has the largest 
set of annotations
* NCBI tends to be more strict in including only high 
confidence gene annotations
* Whereas UCSC contains multiple geneset annotations
that use different criteria

If your experimental system includes non-standard sequences these must be added 
to both the genome fasta and gtf to quantify their expression. Most commonly
this is done for the ERCC spike-ins, although the same must be done for CRISPR-
related sequences or other overexpression/reporter constructs. 

### Alignment using Star

An example of how to map `reads.fq` (the `.gz` means the file is zipped)
using STAR is

```
$<path_to_STAR>/STAR --runThreadN 1 --runMode alignReads
--readFilesIn reads1.fq.gz reads2.fq.gz --readFilesCommand zcat --genomeDir <path>
--parametersFiles FileOfMoreParameters.txt --outFileNamePrefix <outpath>/output
```

__Note__, if the _spike-ins_ are used, the reference sequence
should be augmented with the DNA sequence of the _spike-in_ 
molecules prior to mapping.

__Note__, when UMIs are used, their barcodes should be removed 
from the read sequence. A common practice is to add the barcode
to the read name.

Once the reads for each cell have been mapped to the reference genome,
we need to make sure that a sufficient number of reads from each cell
could be mapped to the reference genome. In our experience, the
fraction of mappable reads for mouse or human cells is 60-70%. However, 
this result may vary depending on protocol, read length and settings for 
the read alignment. As a general rule, we expect all cells to have a similar
fraction of mapped reads, so any outliers should be inspected and
possibly removed. A low proportion of mappable reads usually indicates contamination.

### Pseudo-alignment using Salmon

An example of how to quantify expression using Salmon is
```
$<path_to_Salmon>/salmon quant -i salmon_transcript_index -1 
    reads1.fq.gz -2 reads2.fq.gz -p #threads -l A -g genome.gtf 
    --seqBias --gcBias --posBias
```
__Note__ Salmon produces estimated read counts and estimated 
transcripts per million (tpm) in our experience the latter 
over corrects the expression of long genes for scRNA-seq, 
thus we recommend using read counts. 


### Mapped reads: `BAM` file format

`BAM` file format stores mapped reads in a standard and efficient manner. The 
human-readable version is called a SAM file, while the `BAM` file is the highly
compressed version. `BAM`/`SAM` files contain a header which typically includes  
information on the sample preparation, sequencing and mapping; and a
tab-separated row for each individual alignment of each read. 

`BAM`/`SAM` files can be converted to the other format using 
['samtools'](http://www.htslib.org):



#### Checking read quality from `BAM` file

Assuming that our reads are in `experiment.bam`, we can run FastQC as
```
<path_to_fastQC>/fastQC experiment.bam
```

Below is an example of the output from FastQC for a dataset
of 125 bp reads. The plot reveals a technical error which 
resulted in a couple of bases failing to be read correctly 
in the centre of the read. However, since the rest of the 
read was of high quality this error will most likely have 
a negligible effect on mapping efficiency.

<div class="figure" style="text-align: center">
<img src="figures/per_base_quality.png" alt="Example of FastQC output" width="90%" />
<p class="caption">(\#fig:exprs-constr-fastqc)Example of FastQC output</p>
</div>

Additionally, it is often helpful to visualize the data 
using the 
[Integrative Genomics Browser (IGV)](https://www.broadinstitute.org/igv/) or [SeqMonk](http://www.bioinformatics.babraham.ac.uk/projects/seqmonk/) tools.


### CRAM

[`CRAM`](https://www.ebi.ac.uk/ena/software/cram-usage) files are similar 
to `BAM` files only they contain information in the header 
to the reference genome used in the mapping in the header. This allow the bases
in each read that are identical to the reference to be further compressed. CRAM
also supports some data compression approaches to further optimize storage
compared to `BAM`s. `CRAM`s are mainly used by the Sanger/EBI sequencing facility.


### Reads Mapping QC

#### Total number of reads mapped to each cell 

The histogram below shows the total number of reads mapped to each
cell for an scRNA-seq experiment. Each bar represents one cell, and
they have been sorted in ascending order by the total number of reads
per cell. The three red arrows indicate cells that are outliers in
terms of their coverage and they should be removed from further
analysis. The two yellow arrows point to cells with a surprisingly
large number of unmapped reads. In this example we kept the cells 
during the alignment QC step, but they were later removed during 
cell QC due to a high proportion of ribosomal RNA reads. 

<div class="figure" style="text-align: center">
<img src="figures/Bergiers_exp1_mapping_by_cell.png" alt="Example of the total number of reads mapped to each cell." width="90%" />
<p class="caption">(\#fig:exprs-constr-total-num-cells)Example of the total number of reads mapped to each cell.</p>
</div>

#### Mapping quality

After mapping the raw sequencing to the genome we need to evaluate 
the quality of the mapping. There are many ways to measure the mapping
quality, including: amount of reads mapping to rRNA/tRNAs, proportion 
of uniquely mapping reads, reads mapping across splice junctions,
read depth along the transcripts. Methods developed for bulk RNA-seq, 
such as [RSeQC](http://rseqc.sourceforge.net/), are applicable to 
single-cell data. 

However the expected results will depend on the experimental protocol, 
e.g. many scRNA-seq methods use poly-A selection to avoid sequencing 
rRNAs which results in a 3' bias in the read coverage across the genes 
(aka gene body coverage). The figure below shows this 3' bias as well 
as three cells which were outliers and removed from the dataset:

<div class="figure" style="text-align: center">
<img src="figures/Exp1_RSEQC_geneBodyCoverage_plot_Combined.png" alt="Example of the 3' bias in the read coverage." width="90%" />
<p class="caption">(\#fig:exprs-constr-3-bias)Example of the 3' bias in the read coverage.</p>
</div>

## Reads quantification

The next step is to quantify the expression level of each gene for
each cell. For mRNA data, we can use one of the tools which has been
developed for bulk RNA-seq data, e.g. 
[HT-seq](http://www-huber.embl.de/users/anders/HTSeq/) or
[FeatureCounts](http://subread.sourceforge.net/) (can read about 
[parameters in user guide here](http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf)).

```
# include multimapping
<featureCounts_path>/featureCounts -O -M -Q 30 -p -a genome.gtf -o outputfile input.bam
# exclude multimapping
<featureCounts_path>/featureCounts -Q 30 -p -a genome.gtf -o outputfile input.bam
```

[Unique molecular identifiers (UMIs)](http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1778.html) 
make it possible to count the absolute number of molecules and they have proven popular for [scRNA-seq](http://www.nature.com/nmeth/journal/v11/n2/full/nmeth.2772.html). 

## Unique Molecular Identifiers (UMIs)

Unique Molecular Identifiers are short (4-10bp) random barcodes 
added to transcripts during reverse-transcription. They enable 
sequencing reads to be assigned to individual transcript molecules
and thus the removal of amplification noise and biases from scRNASeq data. 

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-protocol.png" alt="UMI sequencing protocol" width="90%" />
<p class="caption">(\#fig:intro-umi-protocol)UMI sequencing protocol</p>
</div>

When sequencing UMI containing data, techniques are used to
specifically sequence only the end of the transcript containing 
the UMI (usually the 3' end).

For more information processing reads from a UMI experiment, 
including mapping and counting barcodes, I refer you to 
the [UMI chapter](http://hemberg-lab.github.io/scRNA.seq.course/construction-of-expression-matrix.html#umichapter) 
in original course material. 

However, determining how to best process and use UMIs is currently 
an active area of research in the bioinformatics community. 
We are aware of several methods that have recently been developed, including:

* [UMI-tools](https://github.com/CGATOxford/UMI-tools)
* [PoissonUMIs](https://github.com/tallulandrews/PoissonUMIs)
* [zUMIs](https://github.com/sdparekh/zUMIs)
* [dropEst](https://github.com/hms-dbmi/dropEst)




