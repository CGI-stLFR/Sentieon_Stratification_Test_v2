
import os.path

configfile: "run_binning.config"

stratification_map = {}
with open(config["stratification_region_tsv"]) as fh:
  for line in fh:
    name, strat_bed = line.rstrip().split('\t')
    stratification_map[name] = config["stratification_bed_dir"] + '/' + strat_bed

def get_all_output_files():
  '''
  Find all of the expected outputs
  '''
  expected_files = []
  for sample in config["samples"]:
    for truth_hap in ["hapA", "hapB"]:
      expected_files.append("synthetic_binned_fasta/" + sample + "/" + truth_hap + "/binned.fasta")
    for ref_name in ("hs37d5", ):
      expected_files.append("synthetic_stratified_counts/" + sample + '/' + ref_name + "/all/counts.txt")
      expected_files.append("synthetic_binned_truth_vcf_annotated/" + sample + '/' + ref_name + "/calls.vcf.gz")
      expected_files.append("stratified_counts_table/" + sample + '/' + ref_name + "/counts_table.txt")
      for truth_hap in ['hapA', 'hapB']:
        expected_files.append("synthetic_ref_align/" + sample + '/' + truth_hap + '/' + ref_name + "/align.paf.gz")
  print(expected_files)
  return expected_files

rule all:
  input:
    all_outputs = get_all_output_files()

#############################################
# Align each contig to each truth haplotype #
#############################################

rule bin_synthetic_align:
  input:
    input_fa = config["input_dir"] + "/{sample}/LinkedReads.contig.phase_{hap}.fasta",
    truth_fa = config["truth_dir"] + "mergeScaftig_normalized_{truth_hap}.fa",
    minimap_binary = config["minimap2"]
  output:
    aligned = "synthetic_bin_align/{sample}/{hap}/{truth_hap}/align.paf.gz"
  threads:
    32
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.aligned})
    mkdir -p $outdir
    {input.minimap_binary} -x asm5 -t {threads} --paf-no-hit --secondary=no --cs -y {input.truth_fa} <(sed 's/=/:i:/g' {input.input_fa}) | gzip -c > {output.aligned}
    """

##################################
# Align each contig to reference #
##################################

rule ref_synthetic_align:
  input:
    binned_fasta = "synthetic_binned_fasta/{sample}/{truth_hap}/binned.fasta",
    ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref_name]),
    minimap_binary = config["minimap2"]
  output:
    aligned = "synthetic_ref_align/{sample}/{truth_hap}/{ref_name}/align.paf.gz"
  threads:
    32
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.aligned})
    mkdir -p $outdir
    {input.minimap_binary} -x asm5 -t {threads} --paf-no-hit --secondary=no --cs -y {input.ref} <(sed 's/=/:i:/g' {input.binned_fasta}) | gzip -c > {output.aligned}
    """

################################
# Actually perform the binning #
################################

rule bin_by_pid:
  input:
    input_fa = expand(config["input_dir"] + "/{sample}/LinkedReads.contig.phase_{hap}.fasta", sample="{sample}", hap=[0, 1]),
    aligned = expand("synthetic_bin_align/{sample}/{hap}/{truth_hap}/align.paf.gz", sample="{sample}", hap=[0, 1], truth_hap=["hapA", "hapB"]),
    bin_by_pid_script = config["bin_by_pid_script"]
  output:
    binned_results = "synthetic_bin/{sample}/truth_bin_results.txt"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    source /home/eanderson/.virtualenvs/General3/bin/activate

    outdir=$(dirname {output.binned_results})
    mkdir -p $outdir
    python3 {input.bin_by_pid_script} {input.input_fa} {input.aligned} {output.binned_results}
    """

##########################
# Perform the subsetting #
##########################

def get_truth_hap_id(wildcards):
  '''
  Given a truth hap name, get a 0/1 ID
  '''
  truth_hap = wildcards.truth_hap
  if truth_hap == "hapA":
    return '0'
  elif truth_hap == "hapB":
    return '1'
  else:
    raise ValueError("Unknown 'truth_hap': " + wildcards.truth_hap)


rule subset_binned:
  input:
    bin_results = "synthetic_bin/{sample}/truth_bin_results.txt",
    input_0 = config["input_dir"] + "/{sample}/LinkedReads.contig.phase_0.fasta",
    input_1 = config["input_dir"] + "/{sample}/LinkedReads.contig.phase_1.fasta",
    truth_fa = config["truth_dir"] + "mergeScaftig_normalized_{truth_hap}.fa",
    seqtk = config["seqtk"]
  output:
    binned_fasta = "synthetic_binned_fasta/{sample}/{truth_hap}/binned.fasta"
  params:
    truth_hap_id = get_truth_hap_id
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.binned_fasta})
    mkdir -p $outdir
    rm {output.binned_fasta} || $(exit 0)
    contigs=({input.input_0} {input.input_1})
    for hap_id in 0 1; do
      grep "^${{hap_id}}	{params.truth_hap_id}" {input.bin_results} | cut -f 3 > $outdir/tmp_contigs_${{hap_id}}.txt
      {input.seqtk} subseq ${{contigs[$hap_id]}} $outdir/tmp_contigs_${{hap_id}}.txt >> {output.binned_fasta}
    done
    """

##########################################
# Align everything to a reference genome #
##########################################

def get_input_fa(wildcards):
  """
  Get the input file depending on the sample
  """
  if wildcards.sample == "truth":
    return config["truth_dir"] + "mergeScaftig_normalized_" + wildcards.truth_hap + ".fa"
  else:
    return "synthetic_binned_fasta/" + wildcards.sample + '/' + wildcards.truth_hap + "/binned.fasta"

rule align_all:
  input:
    input_fa = get_input_fa,
    ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref_name]),
    minimap_binary = config["minimap2"],
    k8_binary = config["k8"],
    sam_flt = config["sam_flt"],
    samtools = config["samtools"]
  output:
    aligned = "synthetic_binned_aligned/{sample}/{truth_hap}/{ref_name}/aligned.bam"
  threads:
    32
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.aligned})
    mkdir -p $outdir
    {input.minimap_binary} -a -x asm5 -t 16 --paf-no-hit --secondary=no -r2k {input.ref} {input.input_fa} | {input.k8_binary} {input.sam_flt} /dev/stdin | {input.samtools} sort -m4G -@4 -o {output.aligned}
    """

##################################################################################
# Call variants from the aligned BAMs. Merge multiple samples in to a single VCF #
##################################################################################
rule call_all_against_ref:
  input:
    hapA_sample = "synthetic_binned_aligned/{sample}/hapA/{ref_name}/aligned.bam",
    hapB_sample = "synthetic_binned_aligned/{sample}/hapB/{ref_name}/aligned.bam",
    hapA_truth = "synthetic_binned_aligned/truth/hapA/{ref_name}/aligned.bam",
    hapB_truth = "synthetic_binned_aligned/truth/hapB/{ref_name}/aligned.bam",
    ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref_name]),
    htsbox = config["htsbox"],
    sentieon = config["sentieon"]
  output:
    sample_truth_vcf = "synthetic_binned_truth_vcf/{sample}/{ref_name}/calls.vcf.gz"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.sample_truth_vcf})
    mkdir -p $outdir
    {input.htsbox} pileup -q5 -evcf {input.ref} {input.hapA_sample} {input.hapB_sample} {input.hapA_truth} {input.hapB_truth} | {input.sentieon} util vcfconvert - {output.sample_truth_vcf}
    """

################################
# Annotate the called variants #
################################
rule annotate_against_ref:
  input:
    sample_truth_vcf = "synthetic_binned_truth_vcf/{sample}/{ref_name}/calls.vcf.gz",
    label_calls_script = config["label_calls"]
  output:
    annotated_truth_vcf = "synthetic_binned_truth_vcf_annotated/{sample}/{ref_name}/calls.vcf.gz"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.annotated_truth_vcf})
    mkdir -p $outdir
    export PYTHONPATH=/opt/sentieon-genomics-201911/lib/python/sentieon
    python3 {input.label_calls_script} {input.sample_truth_vcf} {output.annotated_truth_vcf}
    """

#################################################
# Get variant counts for each stratified region #
#################################################

rule get_filter_list:
  input:
    annotated_truth_vcf = "synthetic_binned_truth_vcf_annotated/{sample}/{ref_name}/calls.vcf.gz",
    bcftools = config["bcftools"]
  output:
    filter_list = "synthetic_filters/{sample}/{ref_name}/filters.txt"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.filter_list})
    mkdir -p $outdir
    {input.bcftools} view -h {input.annotated_truth_vcf} | grep "^##FILTER" | sed -e 's/^##FILTER=<ID=//' -e 's/,Description=.*$//' > {output.filter_list}
    """

rule filter_stratified:
  input:
    stratification_bed = lambda wildcards: stratification_map[wildcards.stratificaiton_name],
    annotated_truth_vcf = "synthetic_binned_truth_vcf_annotated/{sample}/{ref_name}/calls.vcf.gz",
    bcftools = config["bcftools"],
    sentieon = config["sentieon"]
  output:
    stratified_truth_vcf = "synthetic_stratified_vcf/{sample}/{ref_name}/{stratificaiton_name}/stratified.vcf.gz"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.stratified_truth_vcf})
    mkdir -p $outdir
    {input.bcftools} view -T {input.stratification_bed} {input.annotated_truth_vcf} | {input.sentieon} util vcfconvert - $outdir/stratified.vcf.gz
    """

def get_count_vcf(wildcards):
  '''
  Find the stratified VCF from the wildcards
  '''
  if wildcards.stratificaiton_name == "all":
    return "synthetic_binned_truth_vcf_annotated/" + wildcards.sample + '/' + wildcards.ref_name + "/calls.vcf.gz"
  else:
    return "synthetic_stratified_vcf/" + wildcards.sample + '/' + wildcards.ref_name + '/' + wildcards.stratificaiton_name + "/stratified.vcf.gz"

rule count_filters:
  input:
    filter_list = "synthetic_filters/{sample}/{ref_name}/filters.txt",
    vcf = get_count_vcf,
    bcftools = config["bcftools"]
  output:
    stratified_truth_counts = "synthetic_stratified_counts/{sample}/{ref_name}/{stratificaiton_name}/counts.txt"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.stratified_truth_counts})
    mkdir -p $outdir
    echo "Total variants:" > {output.stratified_truth_counts}
    {input.bcftools} view -H {input.vcf} | wc -l >> {output.stratified_truth_counts}
    while IFS= read -r filter_name; do
      echo "Counting filter $filter_name" >> {output.stratified_truth_counts}
      {input.bcftools} view -H -f "$filter_name" {input.vcf} | wc -l >> {output.stratified_truth_counts}
    done < "{input.filter_list}"
    """

###########################################
# Turn the counts data in to a nice table #
###########################################
rule counts_to_table:
  input:
    "synthetic_stratified_counts/{sample}/{ref_name}/all/counts.txt",
    all_counts_files = expand("synthetic_stratified_counts/{sample}/{ref_name}/{stratificaiton_name}/counts.txt", sample="{sample}", ref_name="{ref_name}", stratificaiton_name=stratification_map.keys()),
    stratification_to_table = config["stratification_to_table"]
  output:
    stratified_counts_table = "stratified_counts_table/{sample}/{ref_name}/counts_table.txt"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.stratified_counts_table})
    mkdir -p $outdir
    python3 {input.stratification_to_table} $(dirname $(dirname {input.all_counts_files[0]})) > {output.stratified_counts_table}
    """
