# This is a python3 script for fetching SRA files.
# You will need to have the SRA toolkit installed. You can either
# install it on your own from the NCBI website, or use anaconda to 
# install: https://anaconda.org/bioconda/sra-tools

# This fetching workflow was adapter from: https://github.com/erilu/python-fastq-downloader

import subprocess

# list of SRA IDs for modern and historical samples
hist_sra_numbers = ["SRR3172000", "SRR3172001", "SRR3172003", "SRR3172004", "SRR3172005", "SRR3172006", "SRR3172007", "SRR3172008", "SRR3172009", "SRR3172010", "SRR3172011", "SRR3172012", "SRR3172014", "SRR3172015", "SRR3172016", "SRR3172018", "SRR3172020", "SRR3172021", "SRR3172022", "SRR3172023", "SRR3172024", "SRR3172025", "SRR3172027", "SRR3172028", "SRR3172029", "SRR3172030", "SRR3172031", "SRR3172032", "SRR3172033", "SRR3172035", "SRR3172036", "SRR3172037", "SRR3172039", "SRR3172040", "SRR3172041", "SRR3172042", "SRR3172043", "SRR3172044", "SRR3172045", "SRR3172047", "SRR3172048", "SRR3172049", "SRR3172051", "SRR3172052", "SRR3172053", "SRR3172054", "SRR3172055", "SRR3172056", "SRR3172057", "SRR3172059", "SRR3172060", "SRR3172061"]
mod_sra_numbers = ["SRR3171968", "SRR3171969", "SRR3171970", "SRR3171971", "SRR3171972", "SRR3171973", "SRR3171974", "SRR3171975", "SRR3171976", "SRR3171977", "SRR3171978", "SRR3171979", "SRR3171980", "SRR3171981", "SRR3171982", "SRR3171983", "SRR3171984", "SRR3171985", "SRR3171986", "SRR3171987", "SRR3171988", "SRR3171989", "SRR3171990", "SRR3171991", "SRR3171992", "SRR3171993", "SRR3171994", "SRR3171995", "SRR3171996", "SRR3171997", "SRR3171998", "SRR3171999", "SRR3172002", "SRR3172013", "SRR3172026", "SRR3172038", "SRR3172050", "SRR3172062"]


# use a simple wget command to fetch .sra files from the historical list
for sra_id in hist_sra_numbers:
    prefetch = "wget https://sra-pub-run-odp.s3.amazonaws.com/sra/" + sra_id + "/" + sra_id
    subprocess.call(prefetch, shell=True)

# use a simple wget command to fetch .sra files from the modern list
for sra_id in mod_sra_numbers:
    prefetch = "wget https://sra-pub-run-odp.s3.amazonaws.com/sra/" + sra_id + "/" + sra_id
    subprocess.call(prefetch, shell=True)

# run the fastq-dump command to extra the fastq files from the historical .sra files
for sra_id in hist_sra_numbers:
    fastq_dump = "fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip " + sra_id
    subprocess.call(fastq_dump, shell=True)

# run the fastq-dump command to extra the fastq files from the modern .sra files
for sra_id in mod_sra_numbers:
    fastq_dump = "fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip " + sra_id
    subprocess.call(fastq_dump, shell=True)
