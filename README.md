# BioInfoSummer 2018 (scRNA-seq workshop)

The material for this work was kindly borrowed with permission and adapted
from the fantastic online course 
[Analysis of single cell RNA-seq data](http://hemberg-lab.github.io/scRNA.seq.course/index.html)
from Vladimir Kiselev (<a href = 'https://twitter.com/wikiselev'>wikiselev</a>), Tallulah Andrews (<a href = 'https://twitter.com/talandrews'>talandrews</a>), Jennifer Westoby (<a href = 'https://twitter.com/Jenni_Westoby'>Jenni_Westoby</a>), Davis McCarthy (<a href = 'https://twitter.com/davisjmcc'>davisjmcc</a>), Maren Büttner (<a href = 'https://twitter.com/marenbuettner'>marenbuettner</a>) and Martin Hemberg (<a href = 'https://twitter.com/m_hemberg'>m_hemberg</a>). 

The material in the course above covers about 1.5 days 
and we will be taking a subset of the material for our 
2-3 hour workshop for 2018 BioInfoSummer where 
we will be discussing the statistical analysis and comprehension 
of single cell RNA-sequencing data in R/Bioconductor. 

## R packages to install 

You will need to install the following R packages: 

```
install.packages("devtools")
install.packages("BiocManager")
install.packages("RColorBrewer", "reshape2", 
			  "matrixStats", "mclust", "pheatmap", "mvoutlier")
devtools::install_github("hemberg-lab/scRNA.seq.funcs")
devtools::install_github("theislab/kBET")
BiocManager::install("scater", "scran", "Rtsne", "sva", 
				"DESeq2", "edgeR", "SC3", "zinbwave")

```

## GitHub

The orginal and complete course material is available at: 

<a href="https://github.com/hemberg-lab/scRNA.seq.course" target="blank">https://github.com/hemberg-lab/scRNA.seq.course</a>

The adapted material for this course at BioInfoSummer 2018 is
available at: 

<a href="https://github.com/stephaniehicks/2018-bioinfosummer-scrnaseq" target="blank">https://github.com/stephaniehicks/2018-bioinfosummer-scrnaseq</a>

## Installation

The course material is available on the
[course GitHub repository](https://github.com/stephaniehicks/2018-bioinfosummer-scrnaseq) 
which can be cloned using 

```{bash, eval=FALSE}
git clone https://github.com/stephaniehicks/2018-bioinfosummer-scrnaseq
```

## License

The license from the [original course material](https://github.com/hemberg-lab/scRNA.seq.course)
is licensed under <b>GPL-3</b> and that license is maintained here. 
Anyone is welcome to go through the material in order to
learn about analysis of scRNA-seq data. If you plan to use 
the material for your own teaching, the original authors have 
requested that they would appreciate it if you tell them about
it in addition to providing a suitable citation. Please contact 
the original lead author 
<a href="mailto:vladimir.yu.kiselev@gmail.com">Vladimir Kiselev</a>. 

## Prerequisites

The course is intended for those who have basic familiarity 
with Unix and the R scripting language. We will also assume that 
you are familiar with mapping and analyzing bulk RNA-seq data as 
well as with the commonly available computational tools.


## Questions/Comments? 

If you have any __comments__, __questions__ or __suggestions__ about the 
original and complete course material, please contact 
<a href="mailto:vladimir.yu.kiselev@gmail.com">Vladimir Kiselev</a>.

If you have questions about the material presented in this course at 
BioInfoSummer 2018, you can reach me at 
<a href="mailto:shicks19@jhu.edu">Stephanie Hicks</a>
