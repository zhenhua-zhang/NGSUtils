#
## Parameters
#


rule all:
  input:
    "GoShifter_enrichment.csv", "gchromVAR_enrichment.csv"


rule Enrich_gchromVAR:
  input:
    "GRCh38.vcf.gz"
  output:
    "gchromVAR_enrichment.csv"
  shell:
    """
      touch {output}
    """


rule Enrich_GoShifter:
  input:
    snp_map = "snp_map.tsv"
    ld_matrix = "LD_matrix.txt",
    genomic_annotations = "genomic_annotations.bed"
  output:
    "GoShifter_enrichment.csv"
  params:
    python_exe = "$(/usr/bin/env python3)",
    goshifter_path = "~/tools/bin/goshifter.py"
  shell:
    """
    {params.python_exe} {params.goshifter_path} \
        -s {input.snp_map}
    """


rule Prepare_LDMatrix:
  input:
    reference_panel = "GRCh38.vcf.gz",
    target_snps = "target_snps.tsv"
  output:
    "LD_matrix.txt"
  params:
    pf_jar = "",
    jv_exe = "/usr/bin/java"
  shell:
    """
    {params.jv_exe} -jar {params.pf_jar} \
        -i {input.target_snps} \
        -o {output}
    """


"""
$> java -jar ProxyFinder-0.2.jar --help

Provide reference
Provide output
usage:
 -h,--hwep <arg>           HWE-P threshold [default: 0.0001]
 -i,--snps <arg>           SNP path (format: 3 or 6 columns, tab
                           separated, one or two snps per line: chr pos
                           rsid)
    --locusld              Perform Pairwise LD calculation within a region
                           (provide region with --regions)
 -m,--maf <arg>            MAF threshold [default: 0.005]
    --matchrsid            Match variants on RS id
 -o,--out <arg>            Output path
    --pairwise             Perform Pairwise LD calculation
    --proxy
 -r,--regions <arg>        Region bed file path
    --samplefilter <arg>   Limit samples to individuals in this list (one
                           sample per line)
 -t,--threshold <arg>      R-squared threshold [default: 0.8]
    --tabix <arg>          Prefix for tabix path [format
                           /path/to/chrCHR.vcf.gz]. Replace the chromosome
                           number with CHR (will be replaced by chr number
                           depending on input SNP or region).
    --threads <arg>        Nr of threads [default: 1]
    --vcf <arg>            Use non-indexed VCF as input
 -w,--windowsize <arg>     Window size [default 1000000]
"""
