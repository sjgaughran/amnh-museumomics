# amnh-museumomics

## Workshop contents

- [Workshop Background](#Workshop-Background)
- [Preparing for the workshop](#Preparing-for-the-workshop)
  * [Required: Getting started with R](#Required-Getting-started-with-R)
  * [Required: Installing R packages](#Required-Installing-R-packages)
  * [Optional: Getting started with anaconda and the command line](#Optional-Getting-started-with-anaconda-and-the-command-line)
    + [Sub-sub-heading](#sub-sub-heading)
- [Study system: changing chipmunks](#Study-system-changing-chipmunks)
- [Fetching publicly available data](#Fetching-publicly-available-data)
  * [SRA](#SRA)
  * [Reference genome](#Reference-genome)
- [Trimming and aligning reads with Paleomix](#Trimming-and-aligning-reads-with-Paleomix)
  * [Looking at Paleomix output](#Looking-at-Paleomix-output)
- [Variant calling with BCFtools](#Variant-calling-with-BCFtools)
- [Variant filtering with SAMtools](#Variant-calling-with-SAMtools)
- [Importing a VCF into R](#Importing-a-VCF-into-R)
- [Running a PCA and DAPC](#Running-a-PCA-and-DAPC)
- [Assessing population structure with sNMF in LEA](#Assessing-population-structure-with-sNMF-in-LEA)
- [Running tests for outlier loci](#Running-tests-for-outlier-loci)
- [Wrap up](#Wrap-up)



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

### Required: Installing R packages ###

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

### Optional: Getting started with anaconda and the command line ###

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


## Trimming and aligning reads with Paleomix ##

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

`bcftools filter -i 'QUAL>20000 && INFO/DP>100' Tminimus_SS.vcf > Tminimus_SS_minQ20kminDP100.vcf`

Next let's filter on individual genotype calls. 

`bcftools filter -S . -i 'FMT/DP>2 | FMT/GQ>20 | FMT/DP<150' Tminimus_SS_minQ20kminDP100.vcf > Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150.vcf`

Finally, we want to make sure our SNPs are bi-allelic. 

`bcftools view -m2 -M2 Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150.vcf > Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150_bi.vcf`

Our last two steps will be removing sites with too much missing data (>20%), and removing individuals with too many missing genotypes (>20%). Before we do that, let's take a look at where these measures stand. Run:

`bcftools stats -s - Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150_bi.vcf`

Sample_SRR3171971 has a suspiciously high number of singletons. 

Looking at the output, some of our samples do have high amounts of missing data (>30%). However, some of these may be affected by sites that are poor quality across most samples. Let's first remove those sites with missing data:

`bcftools filter -i 'F_MISSING<0.2' Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150_bi.vcf > Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150_bi_lowmiss.vcf` 

Checking the missing % again, we're now in much better shape! It looks like all individuals now have less than 10% missing genotype calls, which is good news!

There are several other filtering steps we could take, and the Bi *et al.* 2019 paper covers a few more. For the sake of ease, though, we'll stop our filtering here.

However, there is one final genotype quality aspect that we should consider specifically because we are working with historical DNA. As we discussed in lecture, as DNA degrades over time, cytosines are deaminated (i.e. lose their amino group), which turns them into uracil. This gets prepared in our sequencing libraries as thymine ("T", the DNA version of uracil), which can produce a C/T heterozygote or T/T homozygote at a site that was actually C/C in the organisms genome. (**Note**: because sequencers read in both directions, the above is also true for G/A heterozygotes or A/A homozygotes at sites that were originally G/G.)

The easiest way for us to deal with this issue is to remove sites where the reference allele is C and the alternate allele is T ("C-to-T") or where the reference allele is G and the alternate allele is A ("G-to-A"). This is an overly conservative step, and will remove lots of perfectly valid SNPs from our data set. Filtering for deamination is an active topic of research, and more nuanced ways of filtering deaminated sites while retaining valid transitions are being developed. For now, though, we'll follow Bi and colleagues' lead, and remove C->T and G->A sites. We can filter for this with:

`bcftools filter -e 'REF="C" & ALT="T"' Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150_bi_lowmiss.vcf | bcftools filter -e 'REF="G" & ALT="A"' - > Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150_bi_lowmiss_noTransit.vcf`

While not perfect, this VCF will be the final set of SNPs we use in our population genetic analyses.


## Importing a VCF into R ##
First, download the files from Github and put them in a folder on your desktop labelled "SCCS_tutorial"
Then set this folder as your working directory by typing: 

`setwd("Desktop/SCCS_tutorial/")`

Into your R console and clicking Run.

Now let's call our libraries: 

`library(adegenet)`

`library(vcfR)`

And now we can finally read in the vcf file:  

`vcf <- read.vcfR("Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150_bi_lowmiss_noTransit_12klines.vcf")`

And then make it a genind object, which is the preferred format for adegenet: 

`data <- vcfR2genind(vcf)`

## Running a PCA and DAPC ##

This program uses an algorithm to find the likeliest number of clusters in your data. The first function we will run transforms the data using PCA, then runs a k-means algorithm (testing increasing numbers of k, or populations) and produces summary statistics to evaluate which is the likeliest number of clusters. 
You can read more about it in the DAPC tutorial here: https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
Since we haven't looked at our data yet, let's test a large number of possible clusters. Here we are testing up to 10 clusters: 

`grp<-find.clusters(data,max.n.clust=10)`

The program will ask you to make some choices about your data. 
For this algorithm, we should keep all PCs because it is not computationally intensive and will not suffer from overfitting.
We also generally want to pick the number of clusters with the smallest BIC- in this case, 2. 

Let's check out what has been output into the data frame: 

`names(grp)`

In this data frame, Kstat are the BIC values for each k value, stat is the selected k value, and grp shows you which samples have been assigned to which group. 

`head(grp$Kstat, 5)` 

`grp$stat`

`head(grp$grp, 10)`

Now let's use the DAPC algorithm to describe the clusters, and assign membership probabilities to the samples using a discriminant analysis on our inferred groups.

`dapc1 <- dapc(data, grp$grp)`

This algorithm is a little more prone to overfitting, so let's retain only ~80 PCs. 

`dapc1`

Let's plot it! 

`scatter(dapc1,scree.da=FALSE,bg="white",pch=20,cell=0,cstar=0,col=myCol,solid=.4, cex=3,clab=0,leg=TRUE,txt.leg=paste("Cluster",1:2))`

Not very interesting, but let's see which samples are associated with each population! 

`assignplot(dapc1, cex.lab = 0.4)`

In  this figure, red signifies assignment into a cluster.
The blue marks where these assignments match with our previous estimates of clustering (in the k-means algorithm).
Looks like the modern samples are clustering together and the and historical samples are clustering together!

We know from the paper that there is variation within these groups- North vs South and historical vs modern.
So let's take a look at another value of K, and see the variation within these two distinct groups. 

`dapc1 <- dapc(data, grp$grp)`

Now choose a K of 4, and retain ~80 PCs. 

`scatter(dapc1,scree.da=FALSE,bg="white",pch=20,cell=0,cstar=0,col=myCol,solid=.4, cex=3,clab=0,leg=TRUE,txt.leg=paste("Cluster",1:4))`

Interesting! Let's see which samples are in which group.

`assignplot(dapc1, cex.lab = 0.4)`

So there is a North/South and Historical/Modern divide! 

## Assessing population structure with sNMF in LEA ## 

Now let's take a look at more fine-scale structure in the data. 
First, let's load our libraries: 

`library(LEA)`

`library(vcfR)`

LEA likes things in geno format, so let's make our vcf into a geno file:

`geno1 = vcf2geno("Tminimus_SS_minQ20kminDP100_GenoDP3GQ20DP150_bi_lowmiss_noTransit_12klines.vcf", "Tminimus.geno", force = TRUE)`

Let's make a new project:

`project1 = NULL`

And start our snmf run! 
This program uses admixture analysis, similar to STRUCTURE (Pritchard et al. 2000) using sparse non negative matrix factorization (sNMF). In the end, we will end up with estimates of ancestry proportion for each sample. 
Much of this tutorial was taken from the sNMF vignette. There is also lots more information and analysis in that tutorial! Link to the vignette is here: http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf

`project1= snmf("Tminimus.geno",
                           K = 1:5,
                           entropy = TRUE,
                           repetitions = 10,
                           project = "new")`
                         
This script will test K values of 1 through 5 and record the cross-entropy scores so we can evaluate them.
It will repeat this process 10 times for each value of k, although this is just for the sake of time- you should have many more repetitions when running this on real data. Once it is done, we will plot cross-entropy criterion for all runs in the snmf project. This is a way to evaluate model fit- which number of ancestral populations best explains the genotype data.
Once it is done running, we can plot the cross-entropy: 

`plot(project1, col = "blue", pch = 19, cex = 1.2, main = "Cross Entropy")`

Usually we try to choose the "knee"- so, usually the lowest in this CE graph is probably best.

Then we select the "best" run for each value of K, and continue our downstream analysis with these values: 

`best1 = which.min(cross.entropy(project1, K = 1))`

`best2 = which.min(cross.entropy(project1, K = 2))`

`best3 = which.min(cross.entropy(project1, K = 3))`

`best4 = which.min(cross.entropy(project1, K = 4))`

`best5 = which.min(cross.entropy(project1, K = 5))`

Now, we give the program our list of our sample names for plotting: 

`id <- c("HS_SRR3172000", "HS_SRR3172001", "HS_SRR3172003", "HS_SRR3172004", "HS_SRR3172005",
        "HS_SRR3172006", "HS_SRR3172007", "HE_SRR3172008", "HE_SRR3172009", "HE_SRR3172010",
        "HE_SRR3172011", "HE_SRR3172012", "HS_SRR3172014", "HS_SRR3172015", "HS_SRR3172016",
        "HS_SRR3172018", "HS_SRR3172020", "HS_SRR3172021", "HS_SRR3172022", "HS_SRR3172023",
        "HS_SRR3172024", "HS_SRR3172025", "HS_SRR3172027", "HS_SRR3172028", "HS_SRR3172029",
        "HN_SRR3172030", "HN_SRR3172031", "HN_SRR3172032", "HN_SRR3172033", "HN_SRR3172035",
        "HN_SRR3172036", "HN_SRR3172037", "HN_SRR3172039", "HN_SRR3172040", "HN_SRR3172041", 
        "HN_SRR3172042", "HN_SRR3172043", "HN_SRR3172044", "HN_SRR3172045", "HN_SRR3172047",
        "HN_SRR3172048", "HN_SRR3172049", "HN_SRR3172051", "HN_SRR3172052", "HN_SRR3172053",
        "HN_SRR3172054", "HN_SRR3172055", "HN_SRR3172056", "HN_SRR3172057", "HN_SRR3172059",
        "HN_SRR3172060", "HN_SRR3172061", "MN_SRR3171968", "MN_SRR3171969", "MS_SRR3171970", 
        "MS_SRR3171971", "MN_SRR3171972", "MN_SRR3171973", "MN_SRR3171974", "MN_SRR3171975", 
        "MN_SRR3171976", "MS_SRR3171977", "MS_SRR3171978", "MS_SRR3171979", "MS_SRR3171980", 
        "MS_SRR3171981", "MS_SRR3171982", "MS_SRR3171983", "MS_SRR3171984", "MS_SRR3171985", 
        "MS_SRR3171986", "MS_SRR3171987", "MS_SRR3171988", "MS_SRR3171989", "MS_SRR3171990", 
        "MS_SRR3171991", "MS_SRR3171992", "MS_SRR3171993", "MS_SRR3171994", "MS_SRR3171995", 
        "MS_SRR3171996", "MS_SRR3171997", "MS_SRR3171998", "MS_SRR3171999", "MS_SRR3172002", 
        "MS_SRR3172013", "MS_SRR3172026", "MS_SRR3172038", "MS_SRR3172050", "MS_SRR3172062")`
        
For the sake of figure readability, I have shortened our sample names- H stands for historical, M for modern, N for North and S for South.
 
We also need to choose our colors: 

`my.colors <- c("tomato","lightblue", "olivedrab", "gold")`

And now let's take a look at the data! First for K=2: 

`barchart(project1, K = 2, run = best2, lab = TRUE, 
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix K=2", sort.by.Q = FALSE) -> bp
axis(1, at = 1:length(bp$order),
     labels = id, las=2,
     cex.axis = .3)`
     
Now let's do K=3 and K=4. See if you can figure out how to change this script to plot those K values! 

Also FYI, You don't need to run this from scratch each time! You can load old projects with this command

`project = load.snmfProject("Tminimus.snmfProject")`
 
## Running tests for outlier loci ##
Finally, let's take a look at some outlier loci between the groups. sNMF can produce population differentiation statistics computed from the ancestry coefficients that we just calculated. p-values are then returned for all loci. 
If you're going to do this for real, you will want to take some extra quality control steps, like imputing missing data. We won't be doing that here, but be sure to take a look at the LEA manual to find out more details: http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf

So let's do this first using a k of 2. This script is using a genomic control, which helps reduce the effetcs of population structure. The lambda value is used as an inflation factor to rescale the cho-squared statistics in the computation of p-values.

`p2 = snmf.pvalues(project1, entropy = TRUE, ploidy = 2, K = 2, genomic.control = T, lambda = 2.5)`

`pvalues2 = p2$pvalues`

`par(mfrow = c(2,1))` 

`hist(pvalues2, col = "orange")` 

`plot(-log10(pvalues2), pch = 19, col = "blue", cex = .5)`

Looks like we have a lot of highly differentiated SNPs! 
Let's try with K=3.

`p3 = snmf.pvalues(project1, entropy = TRUE, ploidy = 2, K = 3, genomic.control = T, lambda = 2.5)` 

`pvalues3 = p3$pvalues`

`par(mfrow = c(2,1))` 

`hist(pvalues3, col = "orange")`

`plot(-log10(pvalues3), pch = 19, col = "blue", cex = .5)`

These are some highly differentiated populations! 

## Wrap up ##

As you'll see in the Bi *et al.* paper, we have only covered a subset of analyses in this workshop. Two other main analyses they do in that paper are demographic inference and outlier loci detection. Now that you have the data set, you can follow the methods that paper and try those analyeses yourselves!

Beyond that, there are many more possible analyses out there, and new methods are being developed every year. As aDNA sequencing methods and analytical resources become more advanced, more and more studies will become possible using museum specimens of non-model organisms. We hope this workshop has given you a place to start in planning one of those future innovative studies. 

If you have any questions about the workshop material, or want to chat more about project ideas, please feel free to reach out to Stephen Gaughran (sjgaughran@princeton.edu) or Melina Giakoumis (m.giakoumis1@gmail.com). 



