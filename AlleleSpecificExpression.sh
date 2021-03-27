#!/bin/bash
# Author:  Zhenhua Zhang
# E-mail:  zhenhua.zhang217@gmail.com
# Version: 0.2.0
# License: MIT

#
## ASE analysis pipeline by fastp, STAR, WASP, SAMtools, and GATK/ASEReadCounter
#

#
## WASP (v0.3.4) dependencies:
#
# Python 3.7.x
# Python packages:
#   - numpy==1.19.4
#   - pandas==1.2.3
#   - scipy==1.5.4
#   - pysam==0.16.0.1
#   - pyfaidx==0.5.9.5
#   - PyVCF==0.6.8

#
## Meta config {
#

# Load Easy-build modules
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
# } // Meta config

# ENVs {
export ASEPL_VERSION
ASEPL_VERSION="0.2.0"
ASEPL_FILENAME=$(basename "$0")
# } // ENVs

set -Ee -o pipefail

# Text color and format {
if [[ -x /usr/bin/tput ]] && tput setaf 1 &> /dev/null; then
    tput sgr0
    fmtBold=$(tput bold)
    fmtReset=$(tput sgr0)

    fmtRed=$(tput setaf 1)
    fmtBlue=$(tput setaf 4)
    fmtCyan=$(tput setaf 6)
    fmtGrey=$(tput setaf 8)
    fmtBlack=$(tput setaf 0)
    fmtGreen=$(tput setaf 2)
    fmtWhite=$(tput setaf 7)
    fmtYellow=$(tput setaf 3)
    fmtPurple=$(tput setaf 5)
else
    fmtBold="\e[1m"
    fmtReset="\e[0m"

    fmtRed="\e[1;31m"
    fmtBlue="\e[1;34m"
    fmtCyan="\e[1;36m"
    fmtGrey="\e[1;37m"
    fmtBlack="\e[1;30m"
    fmtGreen="\e[1;32m"
    fmtWhite="\e[1;37m"
    fmtYellow="\e[1;33m"
    fmtPurple="\e[1;35m"
fi
# } // Text color and format


#
## Utilities {
#
# Error information, exit -1
echoErro() {
    echo -e "$fmtBold$fmtRed[E]: $1$fmtReset" >&2 && exit -1
}

# Warning information
echoWarn() {
    echo -e "$fmtBold$fmtYellow[W]:$fmtReset $*" >&2
}

# General information
echoInfo() {
    echo -e "$fmtBold$fmtWhite[I]:$fmtReset $*"
}

echoVersion() {
    cat << EOF

$ASEPL_FILENAME, Version ${ASEPL_VERSION:=UNKNOWN}

EOF
}

# Echo help for the script
echoHelp() {
    cat <<EOF

$ASEPL_FILENAME, Version ${ASEPL_VERSION:=UNKNOWN}

Help:
  -c, --conf     Required.
    The confituration file to supply values for the variables
  -h, --help     Optional.
    Print this help context and exit.
  -v, --version  Optional. Action: store_true
    Print version of current script and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
}
# } // Utilities

#
## CLI options {
#
opt=$(getopt -o "c:hv" -l "conf:,help,version" -- "$@")
eval set -- $opt
while true; do
    case $1 in
        -c|--conf) shift && config=$1 && break ;;
        -h|--help) echoHelp && exit 0 ;;
        -v|--version) echoVersion && exit 0 ;;
        --) shift && break ;;
    esac
    shift
done

# Check the config file
[[ -f ${config:?-c/--conf is required} ]] && source $config || echoErro "Not found $config to read."
# } // CLI options

# This is an example config. The user is allowed to use a config file via CLI -c/--conf option.
# #
# ## Config { # This is a template
# #
# cohortId=CODAM
# fastqId=AD10W1ACXX-1-11
# sampleId=188_2233
# 
# # Master directory
# projDir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp
# workDir=$projDir/outputs/aseQuan_v2
# 
# # SLURM logs
# logDir=$workDir/$cohortId/logDir 
# 
# # Temporary directory and files
# tmpDir=$workDir/$cohortId/tmpDir
# fastpTmpDir=$tmpDir/$fastqId/fastpTmpDir
# starTmpDir=$tmpDir/$fastqId/starTmpDir
# waspTmpDir=$tmpDir/$fastqId/waspTmpDir
# gatkTmpDir=$tmpDir/$fastqId/gatkTmpDir
# aseqTmpDir=$tmpDir/$fastqId/aseqTmpDir
# 
# # The final output results
# optDir=$workDir/$cohortId/optDir
# fastpOptDir=$optDir/$fastqId/fastpOptDir
# starOptDir=$optDir/$fastqId/starOptDir
# waspOptDir=$optDir/$fastqId/waspOptDir
# gatkOptDir=$optDir/$fastqId/gatkOptDir
# aseqOptDir=$optDir/$fastqId/aseqOptDir
# 
# # Genome sequence and gene structure
# genomeDir=$workDir/genomeDir
# genomeFastaFile=$projDir/inputs/GRCh37_reference/human_g1k_v37.fasta
# genomeAnnotationFile=$projDir/inputs/Ensembl_references/Homo_sapiens.GRCh37.75.gtf
# 
# # Genotypes (GT field is required)
# snpH5dbDir=$workDir/snpH5dbDir/$cohortId
# vcfFile=$projDir/outputs/phasing/all-$cohortId-singleAlt.vcf.gz
# chromInfoFile=$projDir/inputs/Ensembl_references/human_g1k_v37_chrom_info.txt
# 
# # FASTQ files
# fastqDir=$workDir/$cohortId/tmpDir
# fastqPrefix=
# fastqSuffix=.fq.gz
# 
# # Gene ids
# geneIdFile=$projDir/inputs/Ensembl_references/protein_coding_gene_id.txt
# 
# # Tools versions, Python virtual env
# pyEnv=~/Documents/projects/wp_ase_dlp/scripts/.env
# GCCVer=GCC/7.3.0-2.30
# GATKVer=GATK/4.1.4.1-Java-8-LTS
# HDF5Ver=HDF5/1.8.14-foss-2018b
# STARVer=STAR/2.6.1c-foss-2018b
# PythonVer=Python/3.7.4-GCCcore-7.3.0-bare
# BCFtoolsVer=BCFtools/1.11-GCCcore-7.3.0
# SAMtoolsVer=SAMtools/1.9-foss-2018b

# # Job dependency
afterok="--dependency=afterok"  # XXX: the colon after the afterok is maully added
# # } // Config

#
## Create the working directory tree
#
# FIXME: Some directories are not necessary
mkdir -p $logDir \
    $tmpDir/$fastqId/{fastp,star,wasp,gatk,aseq}TmpDir \
    $optDir/$fastqId/{fastp,star,wasp,gatk,aseq}OptDir \

#
## STAR: Generate STAR genome index {
#
stepName=STARBuildGenomeIndex
if [ -d $genomeDir ]; then
    echoWarn "Found $genomeDir, skip $stepName"
else
    [[ $genomeFastaFile"xxx" == "xxx" ]] \
        && echoErro "No genome index found, please give genome fasta file"

    [[ $genomeAnnotationFile"xxx" == "xxx" ]] \
        && echoErro "No genome index found, please give genome annotation file"

    mem=150G
    cpus=15
    time=0:59:00
    STARGenomeIndexJobId=$(sbatch \
        --mem=$mem \
        --time=$time \
        --cpus-per-task=$cpus \
        --job-name=$stepName \
        --output=$logDir/%j-$stepName.log \
        <<EOF | cut -d' ' -f4
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
set -Ee -o pipefail

mkdir -p $workDir/genomeDir

module load $STARVer
module list

STAR \
    --runMode genomeGenerate \
    --genomeDir $genomeDir \
    --genomeFastaFiles $genomeFastaFile \
    --sjdbGTFfile $genomeAnnotationFile \
    --runThreadN $cpus \
    --outFileNamePrefix $genomeDir/starGenomeIndex
EOF
)
echoInfo "$stepName was submitted: $STARGenomeIndexJobId ..."
fi
# } // STAR: Generate STAR genome index

#
## WASP: Generate SNP HDF5 database for WASP {
#
stepName=WASPCreateSnpHdf5Database
if [ -d $snpH5dbDir ]; then
    echoWarn "Found $snpH5dbDir, skip $stepName"
else
    mem=5G
    cpus=1
    time=3:59:00
    WaspCreateSnpHdf5DatabaseJobId=$(sbatch \
        --mem=$mem \
        --time=$time \
        --cpus-per-task=$cpus \
        --job-name=$stepName \
        --output=$logDir/%j-%u-$stepName.log \
        <<EOF | cut -f4 -d ' '
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
set -Ee -o pipefail

mkdir -p $snpH5dbDir
mkdir -p $tmpDir/snpH5dbTmpDir

# Split the VCF file by chromosome
module load $BCFtoolsVer
module list

## Index the VCF file, if the default one does not exist. Make SURE you have write permission.
if [[ ! -e $vcfFile.bai ]]; then
    bcftools index \
        --tbi \
        --force \
        --threads $cpus \
        $vcfFile
fi

## Split
seq 1 22 | xargs \
    -n 1 \
    -I '{}' bcftools view -Oz -o $tmpDir/snpH5dbTmpDir/$cohortId-chr{}.vcf.gz $vcfFile {}

# Generate SNP HDF5 database
module purge
module load $HDF5Ver
module list
~/tools/bin/snp2h5 \
    --chrom $chromInfoFile \
    --format vcf \
    --snp_tab $snpH5dbDir/snps_tab.h5 \
    --snp_index $snpH5dbDir/snps_index.h5 \
    --haplotype $snpH5dbDir/haplotype.h5 \
    $tmpDir/snpH5dbTmpDir/$cohortId-chr*.vcf.gz

# Clean up
rm -fr $tmpDir/snpH5dbTmpDir
EOF
)
echoInfo "$stepName was submitted: $WaspCreateSnpHdf5DatabaseJobId ..."
fi
# } // Generate SNP HDF5 database for WASP

#
## FASTP: Preprocessing fastq files. {
#
mem=2G
cpus=1
time=0:19:00
stepName=FastpPreproc
dependency=$afterok
if [[ $STARGenomeIndexJobId != '' ]]; then
    dependency=$dependency:$STARGenomeIndexJobId
fi

if [[ $WaspCreateSnpHdf5DatabaseJobId != '' ]]; then
    dependency=$dependency:$WaspCreateSnpHdf5DatabaseJobId
fi

if [[ $dependency == "--dependency=afterok" ]]; then dependency=""; fi
fastpPreprocJobId=$(sbatch $dependency \
    --mem=$mem \
    --time=$time \
    --cpus-per-task=$cpus \
    --job-name=$fastqId-$stepName \
    --output=$logDir/%j-$fastqId-$stepName.log \
    <<EOF | cut -d' ' -f4
#!/bin/bash
set -Ee -o pipefail
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc

module load $GCCVer # fastp was compiled using GCC v7.3.0
module list

# add a CLI option to fastq files
~/tools/bin/fastp \
    --in1 $fastqDir/$fastqPrefix$fastqId"_R1"$fastqSuffix \
    --in2 $fastqDir/$fastqPrefix$fastqId"_R2"$fastqSuffix \
    --out1 $fastpTmpDir/${fastqId}_paired_R1.fq.gz \
    --out2 $fastpTmpDir/${fastqId}_paired_R2.fq.gz \
    --unpaired1 $fastpTmpDir/${fastqId}_unpaired_R1.fq.gz \
    --unpaired2 $fastpTmpDir/${fastqId}_unpaired_R2.fq.gz  \
    --failed_out $fastpTmpDir/${fastqId}_failed.fq.gz \
    --html $fastpTmpDir/${fastqId}_report.html \
    --json $fastpTmpDir/${fastqId}_report.json \
    --thread $cpus \
    --overrepresentation_sampling 100 \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --correction \
    --trim_poly_g \
    --trim_poly_x
EOF
)
echoInfo $stepName was submitted: $fastpPreprocJobId ...
# } // FASTP: Preprocessing fastq files.

#
## STAR: Mapping reads to genome {
#
mem=35G  # XXX: For human genome, it requires at least 40G memory.
cpus=4
time=0:29:00
stepName=STARMapping
dependency=$afterok:$fastpPreprocJobId
if [[ $dependency == "--dependency=afterok:" ]]; then dependency=""; fi
STARMappingJobId=$(sbatch $dependency \
    --qos=priority \
    --mem=$mem \
    --time=$time \
    --cpus-per-task=$cpus \
    --job-name=$fastqId-$stepName \
    --output=$logDir/%j-$fastqId-$stepName.log \
    <<EOF | cut -d' ' -f4
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
set -Ee -o pipefail

# Mapping
module load $STARVer
module list

STAR --runMode alignReads \
    --genomeDir $genomeDir \
    --runThreadN $cpus \
    --readFilesIn $fastpTmpDir/${fastqId}_paired_R{1,2}.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode "Basic" \
    --outFileNamePrefix $starTmpDir/$fastqId.

# Index the output BAM file
module load $SAMtoolsVer
module list

# Add read group information
samtools addreplacerg \
    --threads $cpus \
    --output-fmt BAM \
    -r "ID:$fastqId\tSM:$fastqId\tPL:illumina" \
    -o $starTmpDir/$fastqId.newrg.bam \
    $starTmpDir/$fastqId".Aligned.sortedByCoord.out.bam"

mv -f $starTmpDir/$fastqId.newrg.bam $starTmpDir/$fastqId.bam

samtools index -@ $cpus $starTmpDir/$fastqId.bam
EOF
)
echoInfo "$stepName was submitted: $STARMappingJobId ..."
# } // Mapping reads to genome 

#
## WASP: Find intersected SNPs {
#
mem=10G
cpus=1
time=1:29:00
stepName=WASPRMBFindIntersectedSnps
dependency=$afterok:$STARMappingJobId
if [[ $dependency == "--dependency=afterok:" ]]; then dependency=""; fi
WaspRMBFindIntersectedSnpsJobId=$(sbatch $dependency \
    --mem=$mem \
    --time=$time \
    --cpus-per-task=$cpus \
    --job-name=$fastqId-$stepName \
    --output=$logDir/%j-$fastqId-$stepName.log \
    <<EOF | cut -f4 -d' '
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
set -Ee -o pipefail

module load $HDF5Ver $PythonVer
module list

if [ $pyEnv"xxx" != "xxx" ]; then source $pyEnv/bin/activate; fi

# Find intersecting SNPs
python ~/tools/WASP/mapping/find_intersecting_snps.py \
    --is_sorted \
    --is_paired_end \
    --output_dir $waspTmpDir \
    --snp_tab $snpH5dbDir/snps_tab.h5 \
    --snp_index $snpH5dbDir/snps_index.h5 \
    --haplotype $snpH5dbDir/haplotype.h5 \
    --samples $sampleId \
    $starTmpDir/$fastqId.bam

mv -f $waspTmpDir/$fastqId.keep.bam $waspTmpDir/$fastqId.keep.sam
mv -f $waspTmpDir/$fastqId.to.remap.bam $waspTmpDir/$fastqId.to.remap.sam

# Convert the SAM into BAM. Note: find_intersecting_snps.py output SAM not BAM
module purge
module load $SAMtoolsVer
module list

samtools view -hbo $waspTmpDir/$fastqId.keep.bam $waspTmpDir/$fastqId.keep.sam
samtools view -hbo $waspTmpDir/$fastqId.to.remap.bam $waspTmpDir/$fastqId.to.remap.sam
rm -f $waspTmpDir/$fastqId.keep.sam $waspTmpDir/$fastqId.to.remap.sam
EOF
)
echoInfo "$stepName was submitted: $WaspRMBFindIntersectedSnpsJobId ..."
# } // Find intersected SNPs

#
## STAR: Remapping by STAR for WASP {
#
mem=35G
cpus=2
time=0:19:00
stepName=WASPRMBRemapping
dependency=$afterok:$WaspRMBFindIntersectedSnpsJobId
if [[ $dependency == "--dependency=afterok:" ]]; then dependency=""; fi
WaspRMBRemappingJobId=$(sbatch $dependency \
    --qos=priority \
    --mem=$mem \
    --time=$time \
    --cpus-per-task=$cpus \
    --job-name=$fastqId-$stepName \
    --output=$logDir/%j-$fastqId-$stepName.log \
    <<EOF | cut -f4 -d' '
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
set -Ee -o pipefail

module load $STARVer
module list

STAR \
    --runMode alignReads \
    --genomeDir $genomeDir \
    --runThreadN $cpus \
    --readFilesIn $waspTmpDir/$fastqId.remap.fq{1,2}.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode "Basic" \
    --outFileNamePrefix $waspTmpDir/$fastqId.remap

# Create index files for BAMs
module load $SAMtoolsVer
module list

# FIXME: if you change the STAR options, the output suffix will change.
mv -f $waspTmpDir/$fastqId.remapAligned.sortedByCoord.out.bam $waspTmpDir/$fastqId.remap.bam

samtools index -@ $cpus $waspTmpDir/$fastqId.remap.bam
EOF
)
echoInfo "$stepName was submitted: $WaspRMBRemappingJobId ..."
# } // STAR: Remapping by STAR for WASP

#
## WASP: RMB filtering and removing duplication {
#
mem=4G
cpus=1
time=1:29:00
stepName=WASPRMBFilterAndRmDup
dependency=$afterok:$WaspRMBRemappingJobId
if [[ $dependency == "--dependency=afterok:" ]]; then dependency=""; fi
WaspRMBFilterAndRmDupJobId=$(sbatch $dependency \
    --mem=$mem \
    --time=$time \
    --cpus-per-task=$cpus \
    --job-name=$fastqId-$stepName \
    --output=$logDir/%j-$fastqId-$stepName.log \
        <<EOF | cut -f4 -d' '
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
set -Ee -o pipefail

# Filter remapped reads
module load $HDF5Ver $PythonVer
module list

if [ $pyEnv"xxx" != "xxx" ]; then source $pyEnv/bin/activate; fi
python ~/tools/WASP/mapping/filter_remapped_reads.py \
    $waspTmpDir/${fastqId}.to.remap.bam \
    $waspTmpDir/${fastqId}.remap.bam \
    $waspTmpDir/${fastqId}.remap.keep.bam

# Merge original and remapped BAM, then sort and index it
module purge
module load $SAMtoolsVer
module list

## Merge
samtools merge \
    -f \
    $waspTmpDir/${fastqId}.keep.merged.bam \
    $waspTmpDir/${fastqId}.remap.keep.bam \
    $waspTmpDir/${fastqId}.keep.bam

## Sort
samtools sort \
    -@ $cpus \
    -o $waspTmpDir/${fastqId}.keep.merged.sorted.bam \
    $waspTmpDir/${fastqId}.keep.merged.bam

## Index
samtools index \
    -@ $cpus \
    $waspTmpDir/${fastqId}.keep.merged.sorted.bam

# # The allelic read counts by bam2h5.py of WASP are in messy, port to GTAK/ASEReadCounter instead.
# # Count allelic reads (with duplications)
# module purge
# module load $HDF5Ver $PythonVer
# module list
# 
# if [ $pyEnv"xxx" != "xxx" ]; then source $pyEnv/bin/activate; fi
# python ~/tools/WASP/CHT/bam2h5.py \
#     --chrom $chromInfoFile \
#     --individual $sampleId \
#     --snp_index $snpH5dbDir/snps_index.h5 \
#     --snp_tab $snpH5dbDir/snps_tab.h5 \
#     --haplotype $snpH5dbDir/haplotype.h5 \
#     --txt_counts $waspTmpDir/${fastqId}.wasp.allAlleleReadCounts.txt \
#     --read_counts $waspTmpDir/${fastqId}.wasp.allAlleleReadCounts.h5 \
#     --ref_as_counts $waspTmpDir/${fastqId}.wasp.refAlleleReadCounts.h5 \
#     --alt_as_counts $waspTmpDir/${fastqId}.wasp.altAlleleReadCounts.h5 \
#     --other_as_counts $waspTmpDir/${fastqId}.wasp.otherAlleleReadCounts.h5 \
#     $waspTmpDir/${fastqId}.keep.merged.sorted.bam \
#     2>&1 | grep -v WARNING


# Remove duplications
module purge
module load $HDF5Ver $PythonVer
module list

if [ $pyEnv"xxx" != "xxx" ]; then source $pyEnv/bin/activate; fi
python ~/tools/WASP/mapping/rmdup_pe.py \
    $waspTmpDir/${fastqId}.keep.merged.sorted.bam \
    $waspTmpDir/${fastqId}.keep.merged.sorted.rmdup.bam

# Sort and index the non-duplicated bam file
module purge
module load $SAMtoolsVer
module list

## Sort
samtools sort \
    -@ $cpus \
    -o $waspTmpDir/${fastqId}.keep.merged.sorted.rmdup.sorted.bam \
    $waspTmpDir/${fastqId}.keep.merged.sorted.rmdup.bam

## Index
samtools index \
    -@ $cpus \
    $waspTmpDir/${fastqId}.keep.merged.sorted.rmdup.sorted.bam

# The allelic read counts by bam2h5.py of WASP are in messy, port to GTAK/ASEReadCounter instead.
# # Count allelic reads (without duplications)
# module purge
# module load $HDF5Ver $PythonVer
# module list
# 
# if [ $pyEnv"xxx" != "xxx" ]; then source $pyEnv/bin/activate; fi
# python ~/tools/WASP/CHT/bam2h5.py \
#     --chrom $chromInfoFile \
#     --individual $sampleId \
#     --snp_index $snpH5dbDir/snps_index.h5 \
#     --snp_tab $snpH5dbDir/snps_tab.h5 \
#     --haplotype $snpH5dbDir/haplotype.h5 \
#     --txt_counts $waspTmpDir/${fastqId}.nondup.wasp.allAlleleReadCounts.txt \
#     --read_counts $waspTmpDir/${fastqId}.nondup.wasp.allAlleleReadCounts.h5 \
#     --ref_as_counts $waspTmpDir/${fastqId}.nondup.wasp.refAlleleReadCounts.h5 \
#     --alt_as_counts $waspTmpDir/${fastqId}.nondup.wasp.altAlleleReadCounts.h5 \
#     --other_as_counts $waspTmpDir/${fastqId}.nondup.wasp.otherAlleleReadCounts.h5 \
#     $waspTmpDir/${fastqId}.keep.merged.sorted.rmdup.sorted.bam
#     2>&1 | grep -v WARNING
EOF
)
echoInfo "$stepName was submitted: $WaspRMBFilterAndRmDupJobId ..."
# } // WASP: RMB filtering and removing duplication

#
## GATK: ASEReadCounter counts reads allelicly {
#
mem=5G
cpus=1
time=1:29:00
stepName=GATKASEReadCounter
dependency=$afterok:$WaspRMBFilterAndRmDupJobId
if [[ $dependency == "--dependency=afterok:" ]]; then dependency=""; fi
GatkASEReadCounterJobId=$(sbatch $dependency \
    --mem=$mem \
    --time=$time \
    --cpus-per-task=$cpus \
    --job-name=$fastqId-$stepName \
    --output=$logDir/%j-$fastqId-$stepName.log \
    <<EOF | cut -d' ' -f4
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
set -Ee -o pipefail

module load $GATKVer
module list

gatk ASEReadCounter \
    --variant $vcfFile \
    --output-format CSV \
    --reference $genomeFastaFile \
    --input $waspTmpDir/${fastqId}.keep.merged.sorted.bam \
    --output $gatkTmpDir/$fastqId.gatkAlleleReadCounts.csv \
    2>&1 | grep -v WARN

gatk ASEReadCounter \
    --variant $vcfFile \
    --output-format CSV \
    --reference $genomeFastaFile \
    --input $waspTmpDir/${fastqId}.keep.merged.sorted.rmdup.sorted.bam \
    --output $gatkTmpDir/$fastqId.nondup.gatkAlleleReadCounts.csv \
    2>&1 | grep -v WARN
EOF
)
echoInfo "$stepName was submitted: $GatkASEReadCounterJobId ..."
# } // GATK: ASEReadCounter counts reads allelicly

#
## aseq estimates ASE effects {
#
mem=3G
cpus=2
time=1:59:00
stepName=AseqASEEstimation
dependency=$afterok:$GatkASEReadCounterJobId
if [[ $dependency == "--dependency=afterok:" ]]; then dependency=""; fi
AseqASEEstimationJobId=$(sbatch $dependency \
    --mem=$mem \
    --time=$time \
    --cpus-per-task=$cpus \
    --job-name=$fastqId-$stepName \
    --output=$logDir/%j-$fastqId-$stepName.log \
    <<EOF | cut -d' ' -f4
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc
set -Ee -o pipefail

module load $PythonVer
if [ $pyEnv"xxx" != "xxx" ]; then source $pyEnv/bin/activate; fi

python $projDir/scripts/asedlp quantify \
    -v $vcfFile \
    -i $sampleId \
    -G $geneIdFile \
    -s $genomeFastaFile \
    -f $genomeAnnotationFile \
    -r $gatkTmpDir/$fastqId.gatkAlleleReadCounts.csv \
    -T $aseqTmpDir/$fastqId-train-set.fa.gz \
    -R $aseqTmpDir/$fastqId-train-report.txt &

python $projDir/scripts/asedlp quantify \
    -v $vcfFile \
    -i $sampleId \
    -G $geneIdFile \
    -s $genomeFastaFile \
    -f $genomeAnnotationFile \
    -r $gatkTmpDir/$fastqId.nondup.gatkAlleleReadCounts.csv \
    -T $aseqTmpDir/$fastqId-nondup-train-set.fa.gz \
    -R $aseqTmpDir/$fastqId-nondup-train-report.txt &

wait
EOF
)
echoInfo "$stepName was submitted: $AseqASEEstimationJobId ..."
# } // aseq estimates ASE effects

#
## Collect output and clean up {
#
mem=500M
cpus=1
time=0:5:00
stepName=CollectOutputAndCleanup
dependency=$afterok:$AseqASEEstimationJobId
if [[ $dependency == "--dependency=afterok:" ]]; then dependency=""; fi
CollectOutputAndCleanupJobId=$(sbatch $dependency \
    --mem=$mem \
    --time=$time \
    --cpus-per-task=$cpus \
    --job-name=$fastqId-$stepName \
    --output=$logDir/%j-$fastqId-$stepName.log \
    <<EOF | cut -d' ' -f4
#!/bin/bash
[[ -f /apps/modules/modules.bashrc ]] && source /apps/modules/modules.bashrc

set -Ee -o pipefail

# Copy fastp reports to output directory.
cp -fr $fastpTmpDir/$fastqId*.{html,json} $fastpOptDir

# Copy GATK/ASEReadCounter outputs to output directory
cp -fr $gatkTmpDir/$fastqId*.csv $gatkOptDir

# Copy asedlp/quantify outputs to output directory
cp -fr $aseqTmpDir/$fastqId*.{fa.gz,txt} $aseqOptDir

# Remove fastq files and temporary directory (tmpDir/$fastqId)
rm -fr $tmpDir/$fastqId
rm -fr $tmpDir/$fastqId"_"R{1,2}.fq.gz
EOF
)
echoInfo "$stepName was submitted. Job id: $CollectOutputAndCleanupJobId ..."
# } // Collect output and clean up

# vim: set ai nowrap nospell ft=sh:
