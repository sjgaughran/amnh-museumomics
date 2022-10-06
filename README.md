# amnh-museumomics

- [Workshop Background](#Workshop Background)
  * [Sub-heading](#sub-heading)
    + [Sub-sub-heading](#sub-sub-heading)




## Workshop background ##

Welcome to the museum genomics workshop at SCCS-NY 2022, hosted by the Center for Biodiversity and Conservation at the American Museum of Natural History! This two day workshop provides an introduction to using museum collections for conservation-related genomics work. The first day covers background of the field, and a chance for breakout group discussion with organizers about project ideas and designs. The second day consists of hands-on data analysis using a publicly available genomics data set. 

This workshop was organized by:  
* Mary Blair (CBC/AMNH)  
* Luca Pozzi (UTSA)  
* Anna Penna (UTSA, NMNH)  
* Alexander Salis (AMNH)  
* Lauren Clark (AMNH)  
* Megan Wallace (AMNH)  
* Suzanne Macey (CBC/AMNH)  
* Melina Giakoumis (CUNY, AMNH)  
* Stephen Gaughran (Princeton U, AMNH)  

## Preparing for the workshop ##

### Required: Getting started with R ###

For the hands-on portion of this workshop, we will be analyzing our SNP data set in R. R is an open-source statistical software package. We use it in conjunction with RStudio, an easy way to manage R code. 

Please try to get R set up on your computer before the workshop. To download R, go to:  
* Windows: https://cran.r-project.org/bin/windows/base/   
* Mac: https://cran.rstudio.com/bin/macosx/   

To download RStudio, go to: https://www.rstudio.com/products/rstudio/download/ 

Install RStudio Desktop (Open Source) accepting the default pathways.  

Now, you will need to tell RStudio where R lives. Open RStudio, go to Tools > Global Options and change the R version by navigating to where you saved R. Unless you have a reason against it, you should use the 64 bit version of R (64 and 32 bit R are both downloaded). 

### Required: Installing R packages

There are several R packages we will be using for this workshop. To install them, open R studio and click on the “Packages” tab in the lower right hand pane. 

Click the “Install” button (all the way to the left in the “Packages” tab) and a window should pop up. It should automatically show that we are installing from the CRAN repository. Type “adegenet” and then click “Install.”

Once that is installed, do this again but type “vcfR” into the pop-up window. 

The third package is not on the CRAN repository, so it requires a different type of installation. In your R console (top left pane), type the following: 

`if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")`

`BiocManager::install("LEA")`

Then click “Run” in the top right of your console (or type command + enter). 

Alternatively, you can install the first two packages by typing the following into your R console: 

`install.packages("adegenet")`

`install.packages(“vcfR”)`

If you run into trouble, no worries! We will have time to sort it out at the beginning of the workshop. 

### Optional: getting started with anaconda and the command line ###

While we will only be doing some basic analyses in the SCCS-NY workshop, this GitHub repo provides a fairly complete workflow for generating the SNP data set and analyzing it. Most of the other steps require tools found not in R, but in the command line. We won't be providing an introduction to working on the command line, but there are many useful (and free) resources out there on how to get oriented with the basic command line interfaces.

All of the command line tools can be downloaded and accessed through anaconda. Anaconda is a way for python packages to be distributed, and the package `conda` allows for the efficient download, installation, and management of anaconda packages. You can find instructions on how to install anaconda here:

https://docs.anaconda.com/anaconda/install/

Once installed, you should use `conda` to create an environment for this workshop (or maybe one for analzying ancient DNA genomic data more generally). "Environment" just means a virtual confined space where you can install packages and run them without worrying about incompatibilities across your system. You can create the environment by running:

`conda create --name ancient-dna`

To activate this environment, run:

`conda activate ancient-dna` or `source activate ancient-dna`

*(different versions of conda have different ways of activating an environment, so you may need to try both)*

Now that you are in your `ancient-dna` environment, you can safely install all of the packages used in this tutorial. To do that, run each of these commands:

`conda install -c bioconda sra-tools`  
`conda install -c bioconda bwa`  
`conda install -c bioconda bcftools`  

To install the paleomix package, follow these instructions:

https://paleomix.readthedocs.io/en/stable/installation.html#conda-installation

Remember that every time you close out of your terminal, your environment will be deactivated. Next time you're ready to use it, reactivate it by running:

`conda activate ancient-dna`

## Study system: changing chipmunks ##

The hands-on portion of this workshop uses a publicly available data set from a study entitled *Temporal genomic contrasts reveal rapid evolutionary responses in an alpine mammal during recent climate change* by Bi *et al.* (2019). Briefly, that study performs exon capture on historic samples (collected from 1911-1916) and modern samples of three populations of chipmunks in the genus *Tamias*. For ease, we will be working only with the samples from *T. alpinus*, which is found in a limited range in the Sierra Mountains of California. You can read more about the species here:

https://en.wikipedia.org/wiki/Alpine_chipmunk

## Fetching publicly available data ##

### SRA ###

The list of SRR codes for historical samples was:

`["SRR3172000", "SRR3172001", "SRR3172003", "SRR3172004", "SRR3172005", "SRR3172006", "SRR3172007", "SRR3172008", "SRR3172009", "SRR3172010", "SRR3172011", "SRR3172012", "SRR3172014", "SRR3172015", "SRR3172016", "SRR3172018", "SRR3172020", "SRR3172021", "SRR3172022", "SRR3172023", "SRR3172024", "SRR3172025", "SRR3172027", "SRR3172028", "SRR3172029", "SRR3172030", "SRR3172031", "SRR3172032", "SRR3172033", "SRR3172035", "SRR3172036", "SRR3172037", "SRR3172039", "SRR3172040", "SRR3172041", "SRR3172042", "SRR3172043", "SRR3172044", "SRR3172045", "SRR3172047", "SRR3172048", "SRR3172049", "SRR3172051", "SRR3172052", "SRR3172053", "SRR3172054", "SRR3172055", "SRR3172056", "SRR3172057", "SRR3172059", "SRR3172060", "SRR3172061"]`

The list of SRR codes for modern samples was:

`["SRR3171968", "SRR3171969", "SRR3171970", "SRR3171971", "SRR3171972", "SRR3171973", "SRR3171974", "SRR3171975", "SRR3171976", "SRR3171977", "SRR3171978", "SRR3171979", "SRR3171980", "SRR3171981", "SRR3171982", "SRR3171983", "SRR3171984", "SRR3171985", "SRR3171986", "SRR3171987", "SRR3171988", "SRR3171989", "SRR3171990", "SRR3171991", "SRR3171992", "SRR3171993", "SRR3171994", "SRR3171995", "SRR3171996", "SRR3171997", "SRR3171998", "SRR3171999", "SRR3172002", "SRR3172013", "SRR3172026", "SRR3172038", "SRR3172050", "SRR3172062"]`

### Reference genome ###

Reference genomes are becoming an increasingly important part of genomics, as more and more become available for non-model organisms. These assemblies allow you to map raw sequencing reads to a single reference for the species (or a closely related species). This not only increases repeatability across studies, it also allows for more accurate mapping and variant calling compared to *de novo* methods. The reference genome for *T. minimus* was recently published in: 

https://onlinelibrary.wiley.com/doi/full/10.1111/evo.14546

And the reference assembly is available at:

https://figshare.com/articles/dataset/Tamias_minimus_de_novo_genome_assembly_fasta_file/19853902/1

After you download the reference, you will have to index it to make it usable for mapping. Indexing creates a number of files that make it easier for our read-mapper (bwa) to efficiently map reads. To do this, run:

`bwa index`


## Trimming and aligning data with Paleomix

https://paleomix.readthedocs.io/en/stable/

### Looking at Paleomix output ###

## Variant calling with BCFtools ##

## Variant filtering with SAMtools ##

*Variant filtering depends on many factors including your data set, sample quality, analyses you'll run, and questions you're asking. Some analyses (especially those that depend on individual sites, like GWAS) are very sensitive to genotyping error. Others, like those calculated from genome-wide averages, are often more foregiving. There are no solid heuristics for filter thresholds, and filtering should be viewed as an iterative process. We are going to do a few filtering steps here to illustrate the process, but these filters are not exhaustive and would likely not be considered sufficient by reviewers!*

Before we get to filtering, let's take a look at the samples included in our VCF. This is always good to check, as Sample IDs occassionally do not come out as you'd expect. You can take a look at them by running:

`bcftools query -l Tminimus_SS.vcf`

which should give you something like:

Sample_SRR3172000
Sample_SRR3172001  
Sample_SRR3172003  
Sample_SRR3172004  
...  
Sample_SRR3172050  
Sample_SRR3172062

That looks good for now, so we can move on to filtering!

Let's start by filtering this variant set by setting a minimum PHRED-scaled quality and a minimum depth for each site. This weeds out the majority of sites that were covered by a few stray reads. We can use bcftools to filter with the command:

`bcftools filter -i 'QUAL>20 && INFO/DP>100' Tminimus_SS.vcf > Tminimus_SS_minQ20minDP100.vcf`

Next let's filter on individual genotype calls. 

`bcftools filter -S . -i 'FMT/DP>2 | FMT/GQ>20' Tminimus_SS_minQ20minDP100.vcf > Tminimus_SS_minQ20minDP100_GenoDP3GQ20.vcf`

Finally, we want to make sure our SNPs are bi-allelic. 

`bcftools view -m2 -M2 Tminimus_SS_minQ20minDP100_GenoDP3GQ20.vcf > Tminimus_SS_minQ20minDP100_GenoDP3GQ20_bi.vcf`

Our last two steps will be removing sites with too much missing data (>20%), and removing individuals with too many missing genotypes (>20%). Before we do that, let's take a look at where these measures stand. 


Sample_SRR3171971 has a suspiciously high number of singletons. 

Looking at the output, some of our samples do have high amounts of missing data (>30%). However, some of these may be affected by sites that are poor quality across most samples. Let's first remove those sites with missing data:

`bcftools filter -i 'F_MISSING<0.2' Tminimus_SS_minQ20minDP100_GenoDP3GQ20_bi.vcf > Tminimus_SS_minQ20minDP100_GenoDP3GQ20_bi_lowmiss.vcf` 

Checking the missing % again, we're now in much better shape! It looks like all individuals now have less than 10% missing genotype calls, which is good news! This is a good example of why it

There are several other filtering steps we could take, and the Bi *et al.* 2019 paper covers a few more. For the sake of ease, though, we'll stop our filtering here.

However, there is one final genotype quality aspect that we should consider specifically because we are working with historical DNA. As we discussed in lecture, as DNA degrades over time, cytosines are deaminated (i.e. lose their amino group), which turns them into uracil. This gets prepared in our sequencing libraries as thymine ("T", the DNA version of uracil), which can produce a C/T heterozygote or T/T homozygote at a site that was actually C/C in the organisms genome. (**Note**: because sequencers read in both directions, the above is also true for G/A heterozygotes or A/A homozygotes at sites that were originally G/G.)

The easiest way for us to deal with this issue is to remove sites where the reference allele is C and the alternate allele is T ("C-to-T") or where the reference allele is G and the alternate allele is A ("G-to-A"). This is an overly conservative step, and will remove lots of perfectly valid SNPs from our data set. Filtering for deamination is an active topic of research, and more nuanced ways of filtering deaminated sites while retaining valid transitions are being developed. For now, though, we'll follow Bi and colleagues' lead, and remove C->T and G->A sites. We can filter for this with:

`bcftools filter -e 'REF="C" & ALT="T"' Tminimus_SS_minQ20minDP100_GenoDP3GQ20_bi_lowmiss.vcf | bcftools filter -e 'REF="G" & ALT="A"' - > Tminimus_SS_minQ20minDP100_GenoDP3GQ20_bi_lowmiss_noTransit.vcf`

While not perfect, this VCF will be the final set of SNPs we use in our population genetic analyses.


## Importing VCF and converting to ##

*Melina*

## Running DAPC ##

*Melina*

## Running structure ## 

*Melina*

## Running an outlier test ##

*Melina*

## Wrap up ##

As you'll see in the Bi *et al.* paper, we have only covered a subset of analyses in this workshop. Two other main analyses they do in that paper are demographic inference and outlier loci detection. Now that you have the data set, you can follow the methods that paper and try those analyeses yourselves!

Beyond that, there are many more possible analyses out there, and new methods are being developed every year. As aDNA sequencing methods and analytical resources become more advanced, more and more studies will become possible using museum specimens of non-model organisms. We hope this workshop has given you a place to start in planning one of those future innovative studies. 

If you have any questions about the workshop material, or want to chat more about project ideas, please feel free to reach out to Stephen Gaughran (sjgaughran@princeton.edu) or Melina Giakoumis (m.giakoumis1@gmail.com). 



