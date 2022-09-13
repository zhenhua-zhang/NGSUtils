#!/bin/bash
#SBATCH --time 8:59:0
#SBATCH --mem 70G
#SBATCH --cpus-per-task 8
#SBATCH --output %j-%u-mtppl.log
#SBATCH --job-name mtppl

#
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : 2020 May 24
# Version   : v0.1.0
# License   : MIT

set -Eeux -o pipefail

# Metatranscriptome analysis pipeline

fq_1=
fq_2=

time_stamp=$(date +%Y%m%d_%H%M%S)

pjdir=~/Documents/projects/wp_liyu
ipdir=${pjdir}/inputs
opdir=${pjdir}/outputs/mtppl_${time_stamp}_${fq_1%%_1.fq.gz}
hostnm=$(hostname)


# Step 0. Setup workspace
mkdir -p \
    ${opdir}/{fastqc/{raw,trimmed},trim_galore,vsearch,sortmerna,humann2,metaphlan2}

#
## Step 1. QC and filtering
#

fastqc=~/tools/bin/fastqc
trim_galore=~/tools/bin/trim_galore
vsearch=~/tools/bin/vsearch
if [[ ${hostnm} =~ "peregrine" ]]; then
    module load FastQC cutadapt
    fastqc=${fastqc}
fi


## QC
${fastqc} \
    --outdir ${opdir}/fastqc/raw \
    --threads 8 \
    ${ipdir}/fastqs/${fq_1} ${ipdir}/fastqs/${fq_2}


## trimming
# Not adapter specified, then trim_galore will detect them
${trim_galore} \
    --paired \
    --retain_unpaired \
    --trim-n \
    --polyA \
    --implicon \
    --output_dir ${opdir}/trim_galore \
    --fastqc_args "--outdir ${opdir}/fastqc/trimmed" \
    --cores 8 \
    ${ipdir}/fastqs/${fq_1} ${ipdir}/fastqs/${fq_2}

## QC
fq_tm_1=${fq_1/.fq.gz/_8bp_UMI_R1.fastq.gz}
fq_tm_2=${fq_2/.fq.gz/_8bp_UMI_R2.fastq.gz}
${fastqc} \
    --outdir ${opdir}/fastqc/trimmed \
    --threads 8 \
    ${opdir}/trim_galore/${fq_tm_1} ${opdir}/trim_galore/${fq_tm_2}

exit

## Dereplication
fa_tm_dr_1=${fq_tm_1/_R1.fastq.gz/.dr.R1.fasta}
fa_tm_dr_2=${fq_tm_2/_R2.fastq.gz/.dr.R2.fasta}
${vsearch} \
    --derep_fulllength ${opdir}/trim_galore/${fq_tm_1} \
    --output ${opdir}/vsearch/${fa_tm_dr_1}

${vsearch} \
    --derep_fulllength ${opdir}/trim_galore/${fq_tm_2} \
    --output ${opdir}/vsearch/${fa_tm_dr_2}

rm -fr ${opdir}/trim_galore/*
#
## Step 2. Removing rRNA sequence
#
# Database: SILVA and Rfam
silva_db=${ipdir}/references/rRNA_databases/silva-
rfam_db=${ipdir}/references/rRNA_databases/rfam-
aligned_opt=${fq_1/_1.fq.gz/.aligned}
other_opt=${fq_1/_1.fq.gz/.other}
~/tools/bin/sortmerna \
    --reads ${opdir}/vsearch/${fa_tm_dr_1} \
    --reads ${opdir}/vsearch/${fa_tm_dr_2} \
    --ref ${silva_db}arc-16s-id95.fasta \
    --ref ${silva_db}arc-23s-id98.fasta \
    --ref ${silva_db}bac-16s-id90.fasta \
    --ref ${silva_db}bac-23s-id98.fasta \
    --ref ${silva_db}euk-18s-id95.fasta \
    --ref ${silva_db}euk-28s-id98.fasta \
    --ref ${rfam_db}5.8s-database-id98.fasta \
    --ref ${rfam_db}5s-database-id98.fasta   \
    --workdir ${opdir}/sortmerna \
    --aligned ${opdir}/sortmerna/out/${aligned_opt} \
    --other ${opdir}/sortmerna/out/${other_opt} \
    --threads 8 \
    --fastx \
    --otu_map \

rm -fr ${opdir}/vsearch/*


#
## Step 3. Taxonmy
#

# fasta, trimmed, dereplication, remove rRNA
fa_tm_dr_rr=${other_opt}.fasta

module purge
module load MetaPhlAn2
module load Bowtie2
# module load Biopython/1.71-foss-2018a-Python-3.6.4

source ~/Documents/projects/wp_liyu/scripts/.env/bin/activate
mp2_profile=${other_opt}.txt
metaphlan2.py ${opdir}/sortmerna/out/${fa_tm_dr_rr} \
    --input_type fasta \
    --nproc 8 > ${opdir}/metaphlan2/${mp2_profile}


#
## Step 4. Metabolic assignment
#
module purge
module load MetaPhlAn2
module load Bowtie2
module list

source ~/Documents/projects/wp_liyu/scripts/.env/bin/activate

humann2 \
    --bypass-prescreen \
    --bypass-nucleotide-index \
    --threads 8 \
    --identity-threshold 35.0 \
    --translated-subject-coverage-threshold 40.0 \
    --translated-query-coverage-threshold 40.0 \
    --input ${opdir}/sortmerna/out/${fa_tm_dr_rr} \
    --output ${opdir}/humann2/

mv ${opdir}/sortmerna/out/{*.log,*_otus.txt} ${opdir}/sortmerna
rm -fr ${opdir}/humann2/*_humann2_temp

mv ${opdir}/sortmerna/out/${fa_tm_dr_rr}.bowtie2out.txt ${opdir}/metaphlan2
