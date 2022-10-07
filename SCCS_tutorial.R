#First, download the files from Github and put them in a folder on your desktop labelled "SCCS_tutorial"
#Then set this folder as your working directory
setwd("Desktop/SCCS_tutorial/")
#Making a PCA with DAPC
library(adegenet)
#read in the vcf
vcf <- read.vcfR("Tminimus_SS_minQ20minDP100_GenoDP3GQ20_bi_lowmiss_noTransit.vcf.gz")
#make it a genind object for adegenet
data <- vcfR2genind(vcf)
#This program uses an algorithm to find the likeliest number of clusters in your data. 
#This first function transforms the data using PCA, then runs a k-means algorithm (testing increasing numbers of k, or populations) and produces summary statistics. 
#You can read more about it in the DAPC tutorial here: https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
grp<-find.clusters(data,max.n.clust=5)
#for this algorithm, we should keep all PCs because it is not computationally intensive and will not suffer from overfitting.
#We generally want to pick the number of clusters with the smallest BIC- in this case, 2. 
names(grp)
#check out the values in this dataframe
head(grp$Kstat, 5)
#check out the actual BIC numbers for each k value 
grp$stat
#And the selected k value by the algorithm
head(grp$grp, 10)
#You can see here which samples are in which group
#Now let's use the DAPC algorithm to describe the clusters and assign membership probabilities to the samples using a discriminant analysis on our inferred groups.
dapc1 <- dapc(data, grp$grp)
#this algorithm is a little more prone to overfitting, so let's retain ~80 PCs.
myCol <- c("darkblue", "purple", "green", "orange", "red", "blue")
scatter(dapc1,scree.da=FALSE,bg="white",pch=20,cell=0,cstar=0,col=myCol,solid=.4, cex=3,clab=0,leg=TRUE,txt.leg=paste("Cluster",1:2))
#Not very interesting, but let's see which samples are associated with each population! 
assignplot(dapc1, cex.lab = 0.4)
#In  this figure, red signifies assignment into a cluster.
#The blue marks where these assignments match with our previous estimates of clustering (k-means).
#looks like the modern samples are clustering together and the and historical samples are clustering together!
#This isn't super interesting, and we know that there might be more variation within the groups based on sampling locality (North vs South).
#So let's take a look at the PCA when we have a k of 4! 
dapc1 <- dapc(data, grp$grp)
#let's retain ~80 PCs. 
scatter(dapc1,scree.da=FALSE,bg="white",pch=20,cell=0,cstar=0,col=myCol,solid=.4, cex=3,clab=0,leg=TRUE,txt.leg=paste("Cluster",1:4))
assignplot(dapc1, cex.lab = 0.4)
#Now let's take a look at more fine-scale structure in the data. 
#Population Structure Analysis with sNMF in LEA 
library(LEA)
library(vcfR)
#LEA likes things in geno format, so let's make our vcf into a geno file
#First you will need to decompress the vcf with gunzip on your terminal. 
geno1 = vcf2geno("Tminimus_SS_minQ20minDP100_GenoDP3GQ20_bi_lowmiss_noTransit.vcf", "Tminimus.geno", force = TRUE)
#Let's make a new project
project1 = NULL
#And start our snmf run! 
#This program uses admixture analysis, similar to STRUCTURE (Pritchard et al. 2000) using sparse non negative matrix factorization (snmf).
#In the end, we will end up with estimates of ancestry proportion for each sample. 
#Much of this tutorial was taken from the sNMF vignette.
#Link to the vignette is here: http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf
project1= snmf("Tminimus.geno",
                           K = 1:5,
                           entropy = TRUE,
                           repetitions = 10,
                           project = "new")
#This script will test K values of 1 through 5 and record the cross-entropy scores so we can evaluate them
#It will repeat this process 10 times for each value of k, although this is just for the sake of time- you should have many more repetitions when running this on real data.
#Once it is done, we will plot cross-entropy criterion for all runs in the snmf project.
#This is a way to evaluate model fit- which number of ancestral populations best explains the genotype data.
#Usually we try to choose the "knee"- so, usually the lowest in this CE graph is probably best.
plot(project1, col = "blue", pch = 19, cex = 1.2, main = "Cross Entropy")
# select the best run for each K 
best1 = which.min(cross.entropy(project1, K = 1))
best2 = which.min(cross.entropy(project1, K = 2))
best3 = which.min(cross.entropy(project1, K = 3))
best4 = which.min(cross.entropy(project1, K = 4))
best5 = which.min(cross.entropy(project1, K = 5))
#Now, we give the program our list of sample names 
id <- c("HS_SRR3172000", "HS_SRR3172001", "HS_SRR3172003", "HS_SRR3172004", "HS_SRR3172005",
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
        "MS_SRR3172013", "MS_SRR3172026", "MS_SRR3172038", "MS_SRR3172050", "MS_SRR3172062")
#And our colors
my.colors <- c("tomato","lightblue", "olivedrab", "gold")
#And now let's take a look at the data! First for K=2
barchart(project1, K = 2, run = best2, lab = TRUE, 
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix K=2", sort.by.Q = FALSE) -> bp
axis(1, at = 1:length(bp$order),
     labels = id, las=2,
     cex.axis = .3)
#Then K=3
barchart(project1, K = 3, run = best3, lab = TRUE, 
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix K=3", sort.by.Q = FALSE) -> bp
axis(1, at = 1:length(bp$order),
     labels = id, las=2,
     cex.axis = .3)
#Then K=4
barchart(project1, K = 4, run = best4, lab = TRUE, 
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix K=4", sort.by.Q = FALSE) -> bp
axis(1, at = 1:length(bp$order),
     labels = id, las=2,
     cex.axis = .3)

#Also FYI, You don't need to run this from scratch each time! You can load old projects with this command
project = load.snmfProject("Tminimus.snmfProject")

#Finally, let's take a look at some outlier loci between the groups. 
#sNMF can produce population differentiation statistics computed from the ancestry coefficients that we just calculated. 
#p-values are then returned for all loci. 
#If you're going to do this for real, you will want to take some extra quality control steps, like imputing missing data. 
#We won't be doing that here, but be sure to take a look at the LEA manual to find out more details
#So let's do this first using a k of 2. 
#This script is using a genomic control, which helps reduce the effetcs of population structure.
#The lambda value is used as an inflation factor to rescale the cho-squared statistics in the computation of p-values. 
p2 = snmf.pvalues(project1, entropy = TRUE, ploidy = 2, K = 2, genomic.control = T, lambda = 2.5) 
pvalues2 = p2$pvalues
par(mfrow = c(2,1)) 
hist(pvalues2, col = "orange") 
plot(-log10(pvalues2), pch = 19, col = "blue", cex = .5)
#Looks like we have a lot of highly differentiated SNPs! 
#Let's try with K=3.
p3 = snmf.pvalues(project1, entropy = TRUE, ploidy = 2, K = 3, genomic.control = T, lambda = 2.5) 
pvalues3 = p3$pvalues
par(mfrow = c(2,1)) 
hist(pvalues3, col = "orange") 
plot(-log10(pvalues3), pch = 19, col = "blue", cex = .5)
#These are some highly differentiated populations! 
