#!/bin/bash
# FileName : metatranscriptome-pipeline.sh
# Created  : 2021-01-31
# E-mail   : zhenhua.zhang217@gmail.com
# Author   : Zhenhua Zhang
# Version  : V0.1.0
# License  : MIT

# Scripts settings
# set -E -o pipefail

# Debug
DEBUG=0
if [[ $DEBUG -eq 1 ]]; then
    extraParameters=--test-only
else
    extraParameters=""
fi

# Text color and format {
if [[ -x /usr/bin/tput ]] && tput setaf 1 &> /dev/null; then
    fmtBold=$(tput bold)
    fmtReset=$(tput sgr0)
    fmtItalic=$(tput sitm)

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
    fmtItalic="\e[3m"

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

Not implemented yet.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
}

echoHelp() {
    cat <<EOF

Not implemented yet.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
}

#
## Basic setting up
#
prjdir=~/Documents/projects/RNA-seq
hostnm=$(hostname)

#
## Step 0: Setup a workspace
#
speciesName=$1
analysisName=transcriptome

iptdir=${prjdir}/inputs
optdir=${prjdir}/buffer/$analysisName/$speciesName

# Workspace tree
fastqId=$2
mkdir -p $optdir/{logs,genomeIndex,$fastqId/{fastp,samtools,hisat2,featureCounts}}

#
## fastp: Filter FASTQ files.
#
# Run settings
stepName=fastpPreproc
stepPath=$optdir/$fastqId/fastp
stepLogs=$optdir/logs/%j-%u-$fastqId-$stepName.log

stepTool=~/tools/bin/fastp
stepTime=0:19:00
stepMem=5G
stepCpu=5  # Current only for paired-end FASTQs. It can only use up to 16 CPUs.
stepCmd=$stepTool

# Input
fastqRawR1=$iptdir/fastqs/raw.split.${fastqId}.1.fq
fastqRawR2=$iptdir/fastqs/raw.split.${fastqId}.2.fq

# Output
fastqCleanR1=$stepPath/${fastqId}_R1_paired.fq.gz
fastqCleanR2=$stepPath/${fastqId}_R2_paired.fq.gz

stepJobId=$(sbatch \
    --mem=$stepMem \
    --time=$stepTime \
    --cpus-per-task=$stepCpu \
    --output=$stepLogs \
    --job-name=$stepName \
    $extraParameters \
    <<EOF | cut -d' ' -f4
#!/bin/bash
set -Ee -o pipefail
if [[ ! -x $stepTool ]]; then module load $stepTool; fi

$stepCmd \
    --in1 $fastqRawR1 \
    --in2 $fastqRawR2 \
    --out1 $stepPath/${fastqId}_R1_paired.fq.gz \
    --out2 $stepPath/${fastqId}_R2_paired.fq.gz \
    --unpaired1 $stepPath/${fastqId}_R1_unpaired.fq.gz \
    --unpaired2 $stepPath/${fastqId}_R2_unpaired.fq.gz \
    --failed_out $stepPath/${fastqId}_failed.fq.gz \
    --html $stepPath/${fastqId}_report.html \
    --json $stepPath/${fastqId}_report.json \
    --thread $stepCpu \
    --overrepresentation_sampling 100 \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --correction \
    --trim_poly_g \
    --trim_poly_x
    # --dont_overwrite

echo fastp preprocessing done.
EOF
)
echoInfo "$stepName was submitted and job could be tracked by jobid $stepJobId"
runDependency="--dependency=afterok:$stepJobId"

#
## Hisat2: Make Hisat2 genomic index.
#
stepName=Hisat2GenomeIndex
stepPath=$optdir/genomeIndex
stepLogs=$optdir/logs/%j-%u-$fastqId-$stepName.log

stepTool=HISAT2
stepTime=0:9:00
stepMem=5G
stepCpu=1
stepCmd=hisat2-build  # Exact command tool will be used.

# Genomic index
genomeFastaFile=$iptdir/bathy-genome/$speciesName/*-newChrId.fna
genomeIndexPref=$optdir/genomeIndex/$speciesName-hisat2Idx

if [[ -d $stepPath && -e $genomeIndexPref.1.ht2 ]]; then
    echoWarn "Found $stepName outputs, ingnore current job!"
else

    stepJobId=$(sbatch \
        --mem=$stepMem \
        --time=$stepTime \
        --cpus-per-task=$stepCpu \
        --output=$stepLogs \
        --job-name=$stepName \
        $extraParameters \
        <<EOF | cut -d' ' -f4
#!/bin/bash
set -Ee -o pipefail
if [[ ! -x $stepTool ]]; then module load $stepTool; fi

$stepCmd \
    -p $stepCpu \
    $genomeFastaFile \
    $genomeIndexPref

echo Hisat2 index build done.
EOF
)
    runDependency=$runDependency:$stepJobId
    echoInfo "$stepName was submitted and job could be tracked by jobid $stepJobId"
fi

#
## Hisat2: Align processed reads against target genome.
#
# Step variables
stepName=Hisat2Alignment
stepPath=$optdir/$fastqId/hisat2
stepLogs=$optdir/logs/%j-%u-$fastqId-$stepName.log

stepTool=HISAT2
stepTime=0:59:00
stepMem=25G
stepCpu=5
stepCmd=hisat2  # Exact command tool will be used.

# Genome index
genomeIndex=$genomeIndexPref

# Alignment SAM files
alignmentSamFile=$stepPath/$fastqId-alignment.sam

stepJobId=$(sbatch \
    --mem=$stepMem \
    --time=$stepTime \
    --cpus-per-task=$stepCpu \
    --output=$stepLogs \
    --job-name=$stepName \
    $runDependency \
    $extraParameters \
    <<EOF | cut -d' ' -f4
#!/bin/bash
set -Ee -o pipefail
if [[ ! -x $stepTool ]]; then module load $stepTool; fi

$stepCmd \
    -1 $fastqCleanR1 \
    -2 $fastqCleanR2 \
    -x $genomeIndex \
    -S $alignmentSamFile \
    --rg-id $fastqId \
    --threads $stepCpu \
    --very-sensitive \
    --no-spliced-alignment

echo Hisat2 alignment done.
EOF
)
echoInfo "$stepName was submitted and job could be tracked by jobid $stepJobId"
runDependency=$runDependency:$stepJobId

#
## SAMtools: compression, sort, and index.
#
stepName=SAMtoolsCompSortIndex
stepPath=$optdir/$fastqId/samtools
stepLogs=$optdir/logs/%j-%u-$fastqId-$stepName.log

stepTool=SAMtools
stepTime=0:59:00
stepMem=5G
stepCpu=5
stepCmd=samtools  # Exact command tool will be used.

# Input and output
inputSamFile=$alignmentSamFile
outputBamFile=$(basename $alignmentSamFile)
outputBamFile=$stepPath/${outputBamFile/.sam/-sort.bam}

stepJobId=$(sbatch \
    --mem=$stepMem \
    --time=$stepTime \
    --cpus-per-task=$stepCpu \
    --output=$stepLogs \
    --job-name=$stepName \
    $runDependency \
    $extraParameters \
    <<EOF | cut -d' ' -f4
#!/bin/bash
set -Ee -o pipefail
if [[ ! -x $stepTool ]]; then module load $stepTool; fi

$stepCmd view -@ $stepCpu -hb $inputSamFile \
    | $stepCmd sort -@ $stepCpu -T $stepPath/$fastqId -o $outputBamFile

$stepCmd index -@ $stepCpu $outputBamFile

echo SAMtools compression, sort, and index done
EOF
)
echoInfo "$stepName was submitted and job could be tracked by jobid $stepJobId"
runDependency=$runDependency:$stepJobId

#
## featureCounts: Count reads per gene.
#
stepName=featureCounts
stepPath=$optdir/$fastqId/featureCounts
stepLogs=$optdir/logs/%j-%u-$fastqId-$stepName.log

stepTool=~/tools/bin/featureCounts
stepTime=0:19:00
stepMem=5G
stepCpu=5
stepCmd=$stepTool  # Exact command tool will be used.

# Input and output
inputBamFile=$outputBamFile
inputAnnotationFile=$iptdir/bathy-gene-structure/$speciesName/*_genomic.gtf
outputReadCountsFile=$stepPath/$fastqId-readcounts.txt
targetFeature=gene
targetFeatureId=gene_id

stepJobId=$(sbatch \
    --mem=$stepMem \
    --time=$stepTime \
    --cpus-per-task=$stepCpu \
    --output=$stepLogs \
    --job-name=$stepName \
    $runDependency \
    $extraParameters \
    <<EOF | cut -d' ' -f4
#!/bin/bash
set -Ee -o pipefail
if [[ ! -x $stepTool ]]; then module load $stepTool; fi

$stepCmd \
    -t $targetFeature \
    -g $targetFeatureId \
    -T $stepCpu \
    -a $inputAnnotationFile \
    -o $outputReadCountsFile \
    $inputBamFile

echo FeatureCounts done.
EOF
)
echoInfo "$stepName was submitted and job could be tracked by jobid $stepJobId"
runDependency=$runDependency:$stepJobId

#
## Cleanup: Count reads per gene.
#
stepName=Cleanup
stepPath=$optdir/$fastqId
stepLogs=$optdir/logs/%j-%u-$fastqId-$stepName.log

stepTool=/usr/bin/rm
stepTime=0:9:00
stepMem=4G
stepCpu=1
stepCmd=$stepTool  # Exact command tool will be used.

stepJobId=$(sbatch \
    --mem=$stepMem \
    --time=$stepTime \
    --cpus-per-task=$stepCpu \
    --output=$stepLogs \
    --job-name=$stepName \
    $runDependency \
    $extraParameters \
    <<EOF | cut -d' ' -f4
#!/bin/bash
set -Ee -o pipefail
if [[ ! -x $stepTool ]]; then module load $stepTool; fi

mv $stepPath/fastp/*.{html,json} $stepPath
mv $stepPath/{samtools,featureCounts}/* $stepPath
rm -fr $stepPath/{fastp,hisat2,samtools,featureCounts}

echo Clean up done.
EOF
)
echoInfo "$stepName was submitted and job could be tracked by jobid $stepJobId"
