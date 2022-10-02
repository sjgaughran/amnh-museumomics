# amnh-museumomics
Workshop on museum genomics for SCCS-NY 2022 at AMNH

## Workshop background ##

## Study background: changing chipmunks ##

Bi et al. 2015, Bi et al. 2019, Herrera et al. 2022

## Fetching publicly available data ##

### SRA ###

The list of SRR codes for historical samples was:

`["SRR3172000", "SRR3172001", "SRR3172003", "SRR3172004", "SRR3172005", "SRR3172006", "SRR3172007", "SRR3172008", "SRR3172009", "SRR3172010", "SRR3172011", "SRR3172012", "SRR3172014", "SRR3172015", "SRR3172016", "SRR3172018", "SRR3172020", "SRR3172021", "SRR3172022", "SRR3172023", "SRR3172024", "SRR3172025", "SRR3172027", "SRR3172028", "SRR3172029", "SRR3172030", "SRR3172031", "SRR3172032", "SRR3172033", "SRR3172035", "SRR3172036", "SRR3172037", "SRR3172039", "SRR3172040", "SRR3172041", "SRR3172042", "SRR3172043", "SRR3172044", "SRR3172045", "SRR3172047", "SRR3172048", "SRR3172049", "SRR3172051", "SRR3172052", "SRR3172053", "SRR3172054", "SRR3172055", "SRR3172056", "SRR3172057", "SRR3172059", "SRR3172060", "SRR3172061"]`

The list of SRR codes for modern samples was:

`["SRR3171968", "SRR3171969", "SRR3171970", "SRR3171971", "SRR3171972", "SRR3171973", "SRR3171974", "SRR3171975", "SRR3171976", "SRR3171977", "SRR3171978", "SRR3171979", "SRR3171980", "SRR3171981", "SRR3171982", "SRR3171983", "SRR3171984", "SRR3171985", "SRR3171986", "SRR3171987", "SRR3171988", "SRR3171989", "SRR3171990", "SRR3171991", "SRR3171992", "SRR3171993", "SRR3171994", "SRR3171995", "SRR3171996", "SRR3171997", "SRR3171998", "SRR3171999", "SRR3172002", "SRR3172013", "SRR3172026", "SRR3172038", "SRR3172050", "SRR3172062"]`

### Reference genome ###

The reference genome was recently published in: 

https://onlinelibrary.wiley.com/doi/full/10.1111/evo.14546

And the reference assembly is available at:

https://figshare.com/articles/dataset/Tamias_minimus_de_novo_genome_assembly_fasta_file/19853902/1




## Trimming and aligning data with Paleomix

https://paleomix.readthedocs.io/en/stable/

### Looking at Paleomix output ###

## Variant calling with BCFtools ##

## Variant filtering with SAMtools ##

Remove individuals with low coverage/calls
Remove calls lower than 5X coverage
Remove sites with more than 20% missing data
Remove sites with very low or high coverage
(From Bi et al) Remove sites with biases associated with reference and alternative allele Phred quality (1e100), mapping quality (0), and distance of alleles from the ends of reads (0.0001). Also
remove sites that show a bias towards sequencing reads coming from the forward or reverse
strand (0.0001).
Remove min map q <20
Remove min base q <20
Remove transitions (CT or GA)

## Getting started with R ##

## Importing VCF and converting to ##

## Running DAPC ##

## Running structure ## 

## Running an outlier test ##

## Warp up ##
