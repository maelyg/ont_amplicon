#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage () {
    log.info """
    ont_amplicon
    Marie-Emilie Gauthier

    Usage:
    Run the command
    nextflow run ont_amplicon {arguments}...

    Required arguments:
      --analysis_mode                 clustering, map2ref
                                      Default: '' [required]

    Optional arguments:
      --help                          Will print this usage document
      -resume                         Resume a failed run
      --outdir                        Path to save the output file
                                      'results'
      --samplesheet '[path/to/file]'  Path to the csv file that contains the list of
                                      samples to be analysed by this pipeline.
                                      Default:  'index.csv'
 

    Contents of index.csv:
      sampleid,sample_files
      SAMPLE01,/user/folder/sample.fastq.gz
      SAMPLE02,/user/folder/*.fastq.gz

      #### Pre-processing and QC options ####
      --merge                         Merge fastq files with the same sample name
      --qc_only                       Only perform preliminary QC step using Nanoplot
                                      Default: false
      --preprocessing_only            Only perform preprocessing steps specied
                                      Default: false

      --adapter_trimming              Run porechop step
                                      Default: false
      --porechop_options              Porechop_ABI options
                                      Default: ''
      --porechop_custom_primers       Limit porechop search to custom adapters specified under porechop_custom_primers_path
                                      Default: ''
      --porechop_custom_primers_path  Path to custom adpaters for porechop
                                      Default: ''

      --qual_filt                     Run quality filtering step using chopper
                                      [False]
      --chopper_options               Chopper options
                                      Default: ''

      --host_filtering                Run host filtering step using Minimap2
                                      Default: false
      --host_fasta                    Path to the fasta file of nucleotide sequences to filter
                                      Default: ''

      #### Analysis mode and associated parameters ####
      ### Clustering (clustering) ###
      --rattle_clustering_options     Rattle clustering options
                                      Default: ''
      --rattle_polishing_options      Rattle polishing options
                                      Default: ''

      ### Map to reference (map2ref) ###
      --reference                     Path to the reference fasta file to map reads to
                                      Default: ''
      --medaka_consensus_options      Medaka options
                                      Default: ''
      --bcftools_min_coverage         Minimum coverage required by bcftools for annotation
                                      Default: '20'

      #### Blast options ####
      --blast_mode                    Specify whether megablast search is against NCBI or a custom database
                                      Default: ''. Select from 'ncbi' or 'localdb'
      --blast_threads                 Number of threads for megablast
                                      Default: '4'
      --blastn_db                     Path to blast database
                                      Default: '4'
      --blast_vs_ref                  blast versus reference specified in fasta file.
                                      Default: 'false'

      ### Reporting ###
      --contamination_flag_threshold  Percentage of maximum FPKM value to use as threshold for flagging detection as potential contamination
                                      Default: '0.01'

    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
if (params.blastn_db != null) {
    blastn_db_name = file(params.blastn_db).name
    blastn_db_dir = file(params.blastn_db).parent
}
if (params.blastn_COI != null) {
    blastn_COI_name = file(params.blastn_COI).name
    blastn_COI_dir = file(params.blastn_COI).parent
}

//if (params.taxdump != null) {
//    taxdump_dir = file(params.taxdump).parent
//}

if (params.reference != null) {
    reference_name = file(params.reference).name
    reference_dir = file(params.reference).parent
}
if (params.host_fasta != null) {
    host_fasta_dir = file(params.host_fasta).parent
}

if (params.porechop_custom_primers == true) {
    porechop_custom_primers_dir = file(params.porechop_custom_primers_path).parent
}

switch (workflow.containerEngine) {
  case "singularity":
    bindbuild = "";
    if (params.blastn_db != null) {
      bindbuild = (bindbuild + "-B ${blastn_db_dir} ")
    }
    if (params.blastn_COI != null) {
      bindbuild = (bindbuild + "-B ${blastn_COI_dir} ")
    }
    if (params.taxdump != null) {
      bindbuild = (bindbuild + "-B ${taxdump} ")
    }
    if (params.reference != null) {
      bindbuild = (bindbuild + "-B ${reference_dir} ")
    }
    if (params.host_fasta != null) {
      bindbuild = (bindbuild + "-B ${host_fasta_dir} ")
    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}

process BLASTN {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '*_megablast*.txt'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly)
  output:
    path("${sampleid}*_megablast*.txt")
    tuple val(sampleid), path("${sampleid}*_megablast_top_10_hits.txt"), emit: blast_results

  script:
  def tmp_blast_output = assembly.getBaseName() + "_megablast_top_10_hits_temp.txt"
  def blast_output = assembly.getBaseName() + "_megablast_top_10_hits.txt"
  
  if (params.blast_mode == "ncbi") {
    """
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${tmp_blast_output} \
      -evalue 1e-3 \
      -word_size 28 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length nident pident mismatch gaps gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames sskingdoms' \
      -max_target_seqs 10
    
    cat <(printf "qseqid\tsgi\tsacc\tlength\tnident\tpident\tmismatch\tgaps\tgapopen\tqstart\tqlen\tqend\tsstart\tsend\tslen\tsstrand\tevalue\tbitscore\tqcovhsp\tstitle\tstaxids\tqseq\tsseq\tsseqid\tqcovs\tqframe\tsframe\n") ${tmp_blast_output} > ${blast_output}
    """
  }
}

process BLASTN_COI {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '*_megablast*.txt'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly), val(gene_target)
  output:
    path("${sampleid}*_megablast_COI_top_hit.txt")
    tuple val(sampleid), path("${sampleid}_ids_to_reverse_complement.txt"), emit: coi_blast_results, optional: true

  script:
  def blast_output_COI = assembly.getBaseName() + "_megablast_COI_top_hit.txt"
    """
    blastn -query ${assembly} \
      -db ${params.blastn_COI} \
      -out ${blast_output_COI} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send evalue bitscore sstrand' \
      -max_target_seqs 1 \
      -max_hsps 1
      grep minus ${blast_output_COI} | cut -f1 > ${sampleid}_ids_to_reverse_complement.txt
    """
}

process BLASTN2 {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '*_megablast*.txt'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly), val(gene_target)
  output:
    path("${sampleid}*_megablast*.txt")
    tuple val(sampleid), path("${sampleid}*_megablast_top_10_hits.txt"), emit: blast_results

  script:
  def tmp_blast_output = assembly.getBaseName() + "_megablast_top_10_hits_temp.txt"
  def blast_output = assembly.getBaseName() + "_megablast_top_10_hits.txt"
  
  if (params.blast_mode == "ncbi") {
    """
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${tmp_blast_output} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length nident pident mismatch gaps gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
      -max_target_seqs 10
      
    cat <(printf "qseqid\tsgi\tsacc\tlength\tnident\tpident\tmismatch\tgaps\tgapopen\tqstart\tqlen\tqend\tsstart\tsend\tslen\tsstrand\tevalue\tbitscore\tqcovhsp\tstitle\tstaxids\tqseq\tsseq\tsseqid\tqcovs\tqframe\tsframe\n") ${tmp_blast_output} > ${blast_output}
    """
  }
}

process CHOPPER {
  publishDir "${params.outdir}/${sampleid}/preprocessing/chopper", pattern: '*_chopper.log', mode: 'link'
  tag "${sampleid}"
  label 'setting_3'

  input:
    tuple val(sampleid), path(sample)

  output:
    path("${sampleid}_chopper.log")
    tuple val(sampleid), path("${sampleid}_filtered.fastq.gz"), emit: chopper_filtered_fq

  script:
  def chopper_options = (params.chopper_options) ? " ${params.chopper_options}" : ''
    """
    gunzip -c ${sample} | chopper ${chopper_options} --threads ${task.cpus} 2> ${sampleid}_chopper.log | gzip > ${sampleid}_filtered.fastq.gz
    """
}
/*
process FASTPLONG {
  publishDir "${params.outdir}/${sampleid}/preprocessing/chopper", pattern: '*_chopper.log', mode: 'link'
  tag "${sampleid}"
  label 'setting_3'

  input:
    tuple val(sampleid), path(sample), val(fwd_primer), val(rev_primer)

  output:
    path("${sampleid}_fastp.log")
    path("${sampleid}_fastp.report.html")
    tuple val(sampleid), path("${sampleid}_filtered.fastq.gz"), emit: fastp_filtered_fq

  script:
// if (fwd_primer != null & rev_primer != null) {
///    adapter_trimming_options = '-s ${fwd_primer} -e ${rev_primer}'
 // }
// else {
 //   adapter_trimming_options = '-A'
// }
    """
    if [[ -z ${fwd_primer} && -z ${rev_primer} ]]; then 
      fastplong -i ${sample} -o ${sampleid}_filtered.fastq.gz -V -x --cut_front --cut_tail --cut_mean_quality ${params.fastplong_min_quality} -m ${params.fastplong_min_quality} -l ${params.fastplong_min_length} -y -h ${sampleid}_fastp.report.html -R '${sampleid}_fastp_report.html' -w ${task.cpus} -s ${fwd_primer} -e ${rev_primer} > ${sampleid}_fastp.log
    else 
      fastplong -i ${sample} -o ${sampleid}_filtered.fastq.gz -V -x --cut_front --cut_tail --cut_mean_quality ${params.fastplong_min_quality} -m ${params.fastplong_min_quality} -l ${params.fastplong_min_length} -y -h ${sampleid}_fastp.report.html -R '${sampleid}_fastp_report.html' -w ${task.cpus} -A > ${sampleid}_fastp.log
    fi
    """
}
*/
process COVSTATS {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/mapping_to_consensus", mode: 'copy'

  input:
    tuple val(sampleid), path(bed), path(consensus), path(coverage), path(top_hits), path(nanostats), val(target_size)
  output:
    path("*top_blast_with_cov_stats.txt"), optional: true
    tuple val(sampleid), path("*top_blast_with_cov_stats.txt"), emit: detections_summary, optional: true

  script:
    """
    if compgen -G "*_coverage.txt" > /dev/null;
      then
        derive_coverage_stats.py --sample ${sampleid} --blastn_results ${top_hits} --nanostat ${nanostats} --coverage ${coverage} --bed ${bed} --target_size ${target_size}
    fi
    """
}

process EXTRACT_READS {
  tag "${sampleid}"
  label "setting_11"
  publishDir "${params.outdir}/${sampleid}/host_filtering", mode: 'copy', pattern: '{*.fastq.gz,*reads_count.txt}'

  input:
  tuple val(sampleid), path(fastq), path(unaligned_ids)
  output:
  path("*reads_count.txt"), emit: read_counts
  file("${sampleid}_unaligned_reads_count.txt")
  file("${sampleid}_unaligned.fastq.gz")
  tuple val(sampleid), path("*_unaligned.fastq"), emit: unaligned_fq

  script:
  """
  seqtk subseq ${fastq} ${unaligned_ids} > ${sampleid}_unaligned.fastq
  gzip -c ${sampleid}_unaligned.fastq > ${sampleid}_unaligned.fastq.gz
  
  n_lines=\$(expr \$(cat ${sampleid}_unaligned.fastq | wc -l) / 4)
  echo \$n_lines > ${sampleid}_unaligned_reads_count.txt
  """
}

process CUTADAPT {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/polishing", pattern: '{*.fasta,*_cutadapt.log}', mode: 'copy'
  tag "${sampleid}"
  label 'setting_1'
  

  input:
    tuple val(sampleid), path(consensus), val(fwd_primer), val(rev_primer)

  output:
    file("${sampleid}_cutadapt.log")
    file("${sampleid}_final_polished_consensus.fasta")
    tuple val(sampleid), path("${sampleid}_final_polished_consensus.fasta"), emit: trimmed

  script:
//  String fwd_primer_trimmed = fwd_primer[-10..-1]
//  String rev_primer_trimmed = rev_primer[0..9]
    """
    if fwd_primer != null && rev_primer != null;
      then
        cutadapt -n 2 -j ${task.cpus} -g "${fwd_primer};max_error_rate=0.1" -a "${rev_primer};max_error_rate=0.1" --trim-n -o ${sampleid}_final_polished_consensus.fasta ${consensus} > ${sampleid}_cutadapt.log
    else
        cutadapt -n 2 -j ${task.cpus} --trim-n -o ${sampleid}_final_polished_consensus.fasta ${consensus} > ${sampleid}_cutadapt.log
    fi
    """
}

//cutadapt -n 2 -j ${task.cpus} -g "${fwd_primer_trimmed};max_error_rate=0.1;min_overlap=10" -a "${rev_primer_trimmed};max_error_rate=0.1;min_overlap=10" --trim-n -o ${sampleid}_final_polished_consensus.fasta ${consensus} > ${sampleid}_cutadapt.log

process EXTRACT_BLAST_HITS {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '{*fasta}'
  tag "${sampleid}"
  label "setting_1"

  input:
    tuple val(sampleid), path(blast_results), val(spp_targets), val(gene_targets), val(target_size)

  output:
    file "${sampleid}_final_polished_consensus_match.fasta"
    file "${sampleid}_reference_match.fasta"
    path("${sampleid}*_megablast_top_hits_tmp.txt")

    tuple val(sampleid), path("${sampleid}*_megablast_top_hits_tmp.txt"), emit: topblast, optional: true
//    tuple val(sampleid), file("${sampleid}*_megablast_top_hits.txt"), emit: blast_results, optional: true
    tuple val(sampleid), path("${sampleid}_reference_match.fasta"), emit: reference_fasta_files, optional: true
    tuple val(sampleid), path("${sampleid}_final_polished_consensus_match.fasta"), emit: consensus_fasta_files, optional: true

  script:
    """
    select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}*_top_10_hits.txt --mode ${params.blast_mode} --spp_targets ${spp_targets}

    # extract segment of consensus sequence that align to reference
    awk  -F  '\\t' 'NR>1 { printf ">%s\\n%s\\n",\$2,\$24 }' ${sampleid}*_top_hits_tmp.txt | sed 's/-//g' > ${sampleid}_final_polished_consensus_match.fasta

    # extract segment of reference that align to consensus sequence
    awk  -F  '\\t' 'NR>1 { printf ">%s_%s\\n%s\\n",\$2,\$4,\$25 }' ${sampleid}*_top_hits_tmp.txt | sed 's/-//g' > ${sampleid}_reference_match.fasta
    """
}

/*
process EXTRACT_REF_FASTA {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy', pattern: '*fasta'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    path("*fasta"), optional: true
    tuple val(sampleid), path("*fasta"), emit: fasta_files, optional: true
  
  script:
    """
    cut -f 3 ${blast_results} | grep -v sacc | sed 's/^/>/' > seq_name.txt
    cut -f21  ${blast_results} | grep -v 'qseq' | sed 's/-//g' > sequence.txt
    cat seq_name.txt sequence.txt > reference.fasta
    """
}
*/

process FASTCAT {
  publishDir "${params.outdir}/${sampleid}/qc/fastcat", mode: 'copy'
  tag "${sampleid}"
  label "setting_1"

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_stats.tsv")
    path("histograms/*")
    tuple val(sampleid), path("${sampleid}.fastq.gz"), emit: merged

  script:
    """
    fastcat \
        -s ${sampleid} \
        -f ${sampleid}_stats.tsv \
        --histograms histograms \
        ${fastq} \
        | bgzip > ${sampleid}.fastq.gz
    """
}

process FASTQ2FASTA {
  publishDir "${params.outdir}/${sampleid}/clustering", mode: 'copy', pattern: '*_rattle.fasta'
  tag "${sampleid}"
  label "setting_1"

  input:
  tuple val(sampleid), path(fastq), path(assembly)
  output:
  tuple val(sampleid), path(fastq), path("${sampleid}_rattle.fasta"), emit: fasta
  tuple val(sampleid), path("${sampleid}_rattle.fasta"), emit: fasta2

  script:
  """
  cut -f1,3 -d ' ' ${assembly} | sed 's/ total_reads=/_RC/' > ${sampleid}_tmp.fastq
  seqtk seq -A -C ${sampleid}_tmp.fastq > ${sampleid}_rattle.fasta
  """
}

process FASTA2TABLE {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy'

  input:
    tuple val(sampleid), path(tophits), path(fasta)
  output:
    file("${sampleid}*_megablast_top_hits.txt")
    tuple val(sampleid), file("${sampleid}*_megablast_top_hits.txt"), emit: blast_results, optional: true

  script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${tophits}
    """
}
/*
process DERIVE_MASKED_CONSENSUS {
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai), path(vcf)

  output:
    tuple val(sampleid), path("${sampleid}_masked.fasta"), path(bam), path(bai), path(vcf), emit: masked_fasta

  script:
    """
    # Get consensus fasta file
    bedtools genomecov -ibam ${bam} -bga > ${sampleid}_genome_cov.bed
    # Assign N to nucleotide positions that have zero coverage
    awk '\$4==0 {print}' ${sampleid}_genome_cov.bed > ${sampleid}_zero_cov_genome_cov.bed
    bedtools maskfasta -fi ${ref} -bed ${sampleid}_zero_cov_genome_cov.bed -fo ${sampleid}_masked.fasta
    """
}

process FILTER_VCF {
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai), path(vcf)

  output:
    path("${sampleid}_medaka.consensus.fasta")
    path("${sampleid}_medaka.annotated.vcf.gz")

  script:
    """
    bcftools reheader ${vcf} -s <(echo '${sampleid}') \
    | bcftools filter \
        -e 'INFO/DP < ${params.bcftools_min_coverage}' \
        -s LOW_DEPTH \
        -Oz -o ${sampleid}_medaka.annotated.vcf.gz

    # create consensus
    bcftools index ${sampleid}_medaka.annotated.vcf.gz
    bcftools consensus -f ${ref} -i 'FILTER="PASS"' ${sampleid}_medaka.annotated.vcf.gz --iupac-codes -H I -o ${sampleid}_medaka.consensus.fasta
    """
}

process MAPPING_BACK_TO_REF {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy', pattern: '*sorted.bam*'
  //publishDir "${params.outdir}/01_VirReport/${sampleid}/alignments/NT", mode: 'link', overwrite: true, pattern: "*{.fa*,.fasta,metrics.txt,scores.txt,targets.txt,stats.txt,log.txt,.bcf*,.vcf.gz*,.bam*}"

  input:
    tuple val(sampleid), path(results)

  output:
    path("*bam"), optional: true
    path("*bam.bai"), optional: true
    tuple val(sampleid), path("*sorted.bam"), emit: bam_files, optional: true
    tuple val(sampleid), path("*sorted.bam.bai"), emit: bai_files, optional: true

  script:
    """
    if compgen -G "*.fasta" > /dev/null;
      then
        mapping_back_to_ref.py --fastq ${sampleid}_preprocessed.fastq.gz
    fi
    """
}

process MEDAKA1 {
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai)

  output:
    tuple val(sampleid), path(ref), path(bam), path(bai), path("${sampleid}_medaka.annotated.unfiltered_sorted.vcf"), emit: unfilt_vcf

  script:
  def medaka_consensus_options = (params.medaka_consensus_options) ? " ${params.medaka_consensus_options}" : ''
    """
    medaka consensus ${bam} ${sampleid}_medaka_consensus_probs.hdf \
      ${medaka_consensus_options} --threads ${task.cpus}

    medaka variant ${ref} ${sampleid}_medaka_consensus_probs.hdf ${sampleid}_medaka.vcf
    medaka tools annotate --dpsp ${sampleid}_medaka.vcf ${ref} ${bam} \
          ${sampleid}_medaka.annotated.unfiltered.vcf
    cat ${sampleid}_medaka.annotated.unfiltered.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > ${sampleid}_medaka.annotated.unfiltered_sorted.vcf
    """
}
*/

process MEDAKA2 {
  publishDir "${params.outdir}/${sampleid}/polishing", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(fastq), path(assembly)

  output:
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), path("${sampleid}_medaka_consensus.bam"), path("${sampleid}_medaka_consensus.bam.bai"), path("${sampleid}_samtools_consensus.fasta")
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), path("${sampleid}_medaka_consensus.bam"), path("${sampleid}_medaka_consensus.bam.bai"), emit: consensus1
   tuple val(sampleid), path("${sampleid}_samtools_consensus.fasta"), emit: consensus2

  script:
  def medaka_consensus_options = (params.medaka_consensus_options) ? " ${params.medaka_consensus_options}" : ''
    """
    medaka_consensus -i ${fastq} -d ${assembly} -t ${task.cpus} -o ${sampleid}
    
    cp ${sampleid}/calls_to_draft.bam ${sampleid}_medaka_consensus.bam
    cp ${sampleid}/calls_to_draft.bam.bai ${sampleid}_medaka_consensus.bam.bai
    cp ${sampleid}/consensus.fasta ${sampleid}_medaka_consensus.fasta
    samtools consensus -f fasta -a -A ${sampleid}_medaka_consensus.bam --call-fract 0.5 -H 0.5 -o ${sampleid}_samtools_consensus.fasta
    """
}

process MINIMAP2_CONSENSUS {
  tag "${sampleid}"
  label 'setting_2'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(consensus), path(fastq)

  output:
    tuple val(sampleid), path(consensus), file("${sampleid}_aln.sam"), emit: aligned_sample

  script:
    """
    minimap2 -ax map-ont -t ${task.cpus} --MD --sam-hit-only ${consensus} ${fastq} > ${sampleid}_aln.sam
    """
}

process MINIMAP2_RACON {
  tag "${sampleid}"
  label "setting_2"

  input:
  tuple val(sampleid), path(fastq), path(assembly)

  output:
  tuple val(sampleid), path(fastq), path(assembly), path("${sampleid}_pre-racon.paf"), emit: draft_mapping

  script:
    """
    minimap2 -L -x ava-ont -t ${task.cpus} ${assembly} ${fastq} > ${sampleid}_pre-racon.paf
    """
}

process MINIMAP2_REF {
  tag "${sampleid}"
  label 'setting_2'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(ref), path(fastq)

  output:
    tuple val(sampleid), path(ref), file("${sampleid}_aln.sam"), emit: aligned_sample

  script:
    """
    minimap2 -ax map-ont --MD -t ${task.cpus} --sam-hit-only ${ref} ${fastq} > ${sampleid}_aln.sam
    """
}

process MOSDEPTH {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/mapping_to_consensus", mode: 'copy'

  input:
    tuple val(sampleid), path(consensus), path(bam), path(bai), path(bed)

  output:
    path("*.mosdepth.global.dist.txt"), optional: true
    path("*.per-base.bed"), optional: true
    path("*regions.bed"), optional: true
    tuple val(sampleid), path("${sampleid}.thresholds.bed"), emit: mosdepth_results, optional: true

  script:
    """
    mosdepth --by ${bed} --thresholds 30 -t ${task.cpus} ${sampleid} ${bam}
    gunzip *regions.bed.gz
    gunzip *.per-base.bed.gz
    gunzip *.thresholds.bed.gz
    """
}

process PYFAIDX {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/mapping_to_consensus", mode: 'copy'

  input:
    tuple val(sampleid), path(fasta)

  output:
    tuple val(sampleid), path("${sampleid}.bed"), emit: bed, optional: true

  script:
    """
    faidx --transform bed ${fasta} > ${sampleid}.bed
    """
}

process PORECHOP_ABI {
  tag "${sampleid}"
  publishDir "$params.outdir/${sampleid}/preprocessing/porechop",  mode: 'copy', pattern: '*_porechop.log'
  label "setting_2"

  input:
    tuple val(sampleid), path(sample)

  output:
    file("${sampleid}_porechop_trimmed.fastq.gz")
    file("${sampleid}_porechop.log")
    tuple val(sampleid), file("${sampleid}_porechop_trimmed.fastq.gz"), emit: porechopabi_trimmed_fq

  script:
  def porechop_options = (params.porechop_options) ? " ${params.porechop_options}" : ''
    """
    if [[ ${params.porechop_custom_primers} == true ]]; then
      porechop_abi -i ${sample} -t ${task.cpus} -o ${sampleid}_porechop_trimmed.fastq.gz --custom_adapters ${params.porechop_custom_primers_path} ${porechop_options}  > ${sampleid}_porechop.log
    else
      porechop_abi -i ${sample} -t ${task.cpus} -o ${sampleid}_porechop_trimmed.fastq.gz ${porechop_options}  > ${sampleid}_porechop.log
    fi
    """
}

process QCREPORT {
  publishDir "${params.outdir}/qc_report", mode: 'copy', overwrite: true
  containerOptions "${bindOptions}"

  input:
    path multiqc_files

  output:
    path("run_qc_report_*txt")
    path("run_qc_report_*html")
    path("run_qc_report_*html"), emit: qc_report_html
    path("run_qc_report_*txt"), emit: qc_report_txt

  script:
    """
    seq_run_qc_report.py --host_filtering ${params.host_filtering} --adapter_trimming ${params.adapter_trimming} --quality_trimming ${params.qual_filt}
    """
}

process RACON {
  publishDir "${params.outdir}/${sampleid}/polishing", mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
   tuple val(sampleid), path(fastq), path(assembly), path(paf)

  output:
   tuple val(sampleid), path(fastq), path("${sampleid}_racon_polished.fasta")
   tuple val(sampleid), path(fastq), path("${sampleid}_racon_polished.fasta"), emit: polished

  script:
    """
    racon -m 8 -x -6 -g -8 -w 500 -t ${task.cpus} -q -1 --no-trimming -u \
        ${fastq} ${paf} ${assembly} \
        > ${sampleid}_racon_polished_tmp.fasta
    cut -f1 -d ' ' ${sampleid}_racon_polished_tmp.fasta > ${sampleid}_racon_polished.fasta
    """
}

process RATTLE {
  tag "${sampleid}"
  label 'setting_10'

  input:
    tuple val(sampleid), path(fastq), val(target_size)

  output:
    file("transcriptome.fq")
    tuple val(sampleid), path("transcriptome.fq"), emit: clusters
    tuple val(sampleid), path("${fastq}"), path("transcriptome.fq"), emit: clusters2

  script:
  def rattle_clustering_options = params.rattle_clustering_options ?: ''
  def rattle_polishing_options = params.rattle_polishing_options ?: ''
  if (params.rattle_clustering_min_length != null) {
    rattle_clustering_min_length_set = params.rattle_clustering_min_length
  }
  else {
    if (target_size != null & target_size.toInteger() <= 300) {
      rattle_clustering_min_length_set = '100'}
    else if (target_size != null & target_size.toInteger() > 300) {
      rattle_clustering_min_length_set = '150'}
    else { 
      rattle_clustering_min_length_set = '150'}
  }
    """
    (
      set +eo pipefail
      rattle cluster -i ${fastq} -t ${task.cpus} --lower-length ${rattle_clustering_min_length_set} ${rattle_clustering_options} -v ${params.rattle_clustering_max_variance} -o .
      rattle cluster_summary -i ${fastq} -c clusters.out > ${sampleid}_cluster_summary.txt
      mkdir clusters
      rattle extract_clusters -i ${fastq} -c clusters.out -l ${sampleid} -o clusters --fastq
      rattle correct -t ${task.cpus} -i ${fastq} -c clusters.out -t ${task.cpus} -l ${sampleid}
      rattle polish -i consensi.fq -t ${task.cpus} --summary ${rattle_polishing_options}
      trap 'if [[ \$? -eq 139 ]]; then echo "segfault !"; fi' CHLD
    ) 2>&1 | tee ${sampleid}_rattle.log

    
    if [[ ! -s transcriptome.fq ]]
    then
        touch transcriptome.fq
        echo "Rattle clustering and polishing failed."
    else
      echo "Rattle clustering and polishing completed successfully."
    fi
    """
}

//grep '@cluster' transcriptome.fq  | cut -f1,3 -d ' '  | sed 's/total_reads=//' | sort -k2,2 -rn | sed 's/ /_RC/' | sed 's/@//' | head -n 10 > ${sampleid}_most_abundant_clusters_ids.txt

process REFORMAT {
  tag "${sampleid}"
  label "setting_3"
  publishDir "$params.outdir/${sampleid}/preprocessing", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_preprocessed.fastq.gz")
    path("${sampleid}_basecalling_model_inference.txt")
    tuple val(sampleid), path("${sampleid}_preprocessed.fastq.gz"), emit: reformatted_fq
    tuple val(sampleid), path("${sampleid}_preprocessed.fastq.gz"), emit: cov_derivation_ch

  script:
    """
    reformat.sh in=${fastq} out=${sampleid}_preprocessed.fastq.gz trd qin=33
    zcat ${fastq} | head -n1 | sed 's/^.*basecall_model_version_id=//' > ${sampleid}_basecalling_model_inference.txt
    """
}

process REVCOMP {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '{*fasta}'
  tag "${sampleid}"
  label "setting_1"
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(contigs), path(ids_to_revcomp)


  output:
    file "${sampleid}_final_polished_consensus_rc.fasta"
    tuple val(sampleid), path("${sampleid}_final_polished_consensus_rc.fasta"), emit: revcomp, optional: true

  script:
    """
    reverse_complement.py --sample ${sampleid} --ids_to_rc ${ids_to_revcomp} --fasta ${contigs}
    
    """
}

process SAMTOOLS {
  publishDir "${params.outdir}/${sampleid}/mapping_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(ref), path(sample)

  output:
    path "${sampleid}_aln.sorted.bam"
    path "${sampleid}_aln.sorted.bam.bai"
    path "${sampleid}_coverage.txt"
    path "${sampleid}_histogram.txt"
    path "${sampleid}_samtools_consensus_from_ref.fasta"
    tuple val(sampleid), path(ref), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_sample

  script:
    """
    samtools view -Sb -F 4 ${sample} | samtools sort -o ${sampleid}_aln.sorted.bam
    samtools index ${sampleid}_aln.sorted.bam
    samtools coverage ${sampleid}_aln.sorted.bam > ${sampleid}_coverage.txt
    samtools coverage -A -w 50 ${sampleid}_aln.sorted.bam > ${sampleid}_histogram.txt
    samtools consensus -f fasta -a -A ${sampleid}_aln.sorted.bam --call-fract 0.5 -H 0.5 -o ${sampleid}_samtools_consensus_from_ref.fasta
    """
}

process SAMTOOLS_CONSENSUS {
  publishDir "${params.outdir}/${sampleid}/mapping_to_consensus", mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(consensus), path(sample)

  output:
    path "${sampleid}_aln.sorted.bam"
    path "${sampleid}_aln.sorted.bam.bai"
    path "${sampleid}_coverage.txt"
    path "${sampleid}_histogram.txt"
    tuple val(sampleid), path(consensus), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_bams
    tuple val(sampleid), path("${sampleid}_coverage.txt"), emit: coverage
    tuple val(sampleid), path("${sampleid}_contigs_reads_ids.txt"), emit: contig_seqids
  script:
    """
    samtools view -S -F 4 ${sample} | cut -f1,3 | sort | uniq > ${sampleid}_contigs_reads_ids.txt
    samtools view -Sb -F 4 ${sample} | samtools sort -o ${sampleid}_aln.sorted.bam
    samtools index ${sampleid}_aln.sorted.bam
    samtools coverage ${sampleid}_aln.sorted.bam  > ${sampleid}_coverage.txt
    samtools coverage -A -w 50 ${sampleid}_aln.sorted.bam > ${sampleid}_histogram.txt
    """
}

/*
samtools view -S -F 4 ${sample} | cut -f3 | sort | uniq > contigs_id.txt
    for id in `cat contigs_id.txt`;
      do 
        samtools view  -S -F 4 ${sample} | grep ${id} | cut -f1 >  ${sampleid}_${id}_aligning_ids.txt; 
      done
*/
process TIMESTAMP_START {
    publishDir "${params.outdir}/logs", mode: 'copy', overwrite: true
    cache false
    output:
    path "*nextflow_start_timestamp.txt"
    path("*nextflow_start_timestamp.txt"), emit: timestamp

    script:
    """
    START_TIMESTAMP=\$(date "+%Y%m%d%H%M%S")
    echo "\$START_TIMESTAMP" > "\${START_TIMESTAMP}_nextflow_start_timestamp.txt"
    """
}

process HTML_REPORT {
  publishDir "${params.outdir}/${sampleid}/html_report", mode: 'copy', overwrite: true
  containerOptions "${bindOptions}"
  label 'setting_3'

  input:
    tuple val(sampleid), path(consensus_fasta), path(consensus_match_fasta), path(aln_sorted_bam), path(aln_sorted_bam_bai), path(raw_nanoplot), path(filtered_nanoplot), path(top_blast_hits), path(blast_with_cov_stats)
//   path(timestamp), path(report), path(index)
    path("*")

  output:
    path("*")

  script:
    """
    build_report.py .
    """
}

process SEQTK {
  tag "${sampleid}"
  label "setting_2"

  input:
  tuple val(sampleid), path(contig_seqids), path(fastq)
  output:
  tuple val(sampleid), path("${sampleid}.fasta"), emit: fasta
  tuple val(sampleid), path(contig_seqids), emit: contig_seqids

  script:
  """
  seqtk seq -A -C ${fastq} > ${sampleid}_all_reads.fasta
  cut -f1 ${contig_seqids} | sort | uniq > reads_ids.txt
  seqtk subseq ${sampleid}_all_reads.fasta reads_ids.txt > ${sampleid}.fasta
  """
}

/*
cut -f2 ${contig_seqids} | sort | uniq > contigs.txt
  
  for id in `cat contigs.txt`;
    do
      grep ${id} ${contig_seqids} | cut -f1 >  ${sampleid}_${id}_aligning_ids.txt;
      seqtk subseq ${sampleid}.fasta ${sampleid}_${id}_aligning_ids.txt > ${sampleid}_${id}.fasta
    done



process EXTRACT_READ_LENGTHS {
  tag "${sampleid}"
  label "setting_2"

  input:
  tuple val(sampleid), path(fasta), path(contig_seqids), path(results_table)
  output:
  path "${sampleid}.fasta"

  script:
  """
  seq_length.py --sample ${sampleid} --contig_seqids ${contig_seqids} --fasta ${fasta} --results_table ${results_table}
  """
}
*/

include { NANOPLOT as QC_PRE_DATA_PROCESSING } from './modules.nf'
include { NANOPLOT as QC_POST_DATA_PROCESSING } from './modules.nf'
//include { BLASTN as BLASTN } from './modules.nf'
//include { BLASTN as BLASTN2 } from './modules.nf'

workflow {
  TIMESTAMP_START ()
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), file(row.sample_files)) }
      .set{ ch_sample }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.spp_targets), (row.gene_targets), (row.target_size)) }
      .set{ ch_targets }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.target_size)) }
      .set{ ch_target_size }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.gene_targets)) }
      .set{ ch_gene_targets }

    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.fwd_primer), (row.rev_primer)) }
      .set{ ch_primers }
    
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.gene_targets)) }
      .filter { sampleid, gene_targets -> gene_targets.contains("COI") }
      .set{ ch_coi }
    
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.gene_targets)) }
      .filter { sampleid, gene_targets -> !gene_targets.contains("COI") }
      .set{ ch_other }
  
  } else { exit 1, "Input samplesheet file not specified!" }


  if ( params.analysis_mode == 'clustering') {
    if (!params.blast_vs_ref) {
      if ( params.blastn_db == null) {
        error("Please provide the path to a blast database using the parameter --blastn_db.")
      }
    }
    else if (params.blast_vs_ref ) {
      if ( params.reference == null) {
      error("Please provide the path to a reference fasta file with the parameter --reference.")
      }
    }
  }
  else if ( params.analysis_mode == 'map2ref' ) {
    if ( params.reference == null) {
      error("Please provide the path to a reference fasta file with the parameter --reference.")
      }
  }
  
  if (params.merge) {
    //Merge split fastq.gz files
    FASTCAT ( ch_sample )
    //Run Nanoplot on merged raw fastq files before data processing
    QC_PRE_DATA_PROCESSING ( FASTCAT.out.merged )
    fq = FASTCAT.out.merged
  }
  else {
    fq = ch_sample
    QC_PRE_DATA_PROCESSING ( fq )
  }

  // Data pre-processing
  if (!params.qc_only) {
    // Remove adapters using PORECHOP_ABI
    if (params.adapter_trimming) {
      PORECHOP_ABI ( fq )
      trimmed_fq = PORECHOP_ABI.out.porechopabi_trimmed_fq
    }
    else {
      trimmed_fq = fq
    }

    // Perform quality filtering of reads using chopper
    if (params.qual_filt) {
      CHOPPER ( trimmed_fq)
      filtered_fq = CHOPPER.out.chopper_filtered_fq
//      FASTPLONG ( trimmed_fq.join(ch_primers))
 //     filtered_fq = FASTPLONG.out.fastp_filtered_fq
    }
    else { filtered_fq = trimmed_fq
    }

    //Reformat fastq read names after the first whitespace
    REFORMAT( filtered_fq )

    //Run Nanoplot on merged raw fastq files after data processing
    if ( params.qual_filt & params.adapter_trimming | !params.qual_filt & params.adapter_trimming | params.qual_filt & !params.adapter_trimming) {
      QC_POST_DATA_PROCESSING ( filtered_fq )
    }
    
/*
    //Legacy code from ontvisc to filter host sequences, consider removing if not needed
    if (params.host_filtering) {
      if ( params.host_fasta == null) {
        error("Please provide the path to a fasta file of host sequences that need to be filtered with the parameter --host_fasta.")
      }
      else {
        MINIMAP2_ALIGN_RNA ( REFORMAT.out.reformatted_fq, params.host_fasta )
        EXTRACT_READS ( MINIMAP2_ALIGN_RNA.out.sequencing_ids )
        final_fq = EXTRACT_READS.out.unaligned_fq
      }
    }
    else {
      final_fq = REFORMAT.out.reformatted_fq
    }
*/

    final_fq = REFORMAT.out.reformatted_fq
    //Derive QC report if any preprocessing steps were performed
    //if ( params.qual_filt & params.host_filtering | params.adapter_trimming & params.host_filtering ) {
    if ( params.qual_filt | params.adapter_trimming ) {
      ch_multiqc_files = Channel.empty()
      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      QCREPORT(ch_multiqc_files.collect())
    }


    if (!params.preprocessing_only) {
      //Currently only one analysis mode in ont_amplicon, this is legacy from ontvisc, consider removing if no other mode is added to this pipeline
      //We have had talks about including an option to just a map to a reference of interest, but this is not implemented yet
      if ( params.analysis_mode == 'clustering' ) {
        //Perform clustering using Rattle and convert to fasta file
        ch_fq_target_size = (final_fq.join(ch_target_size))
        RATTLE ( ch_fq_target_size )
        FASTQ2FASTA( RATTLE.out.clusters2 )

        //Polish consensus sequence using Racon followed by Medaka
        if (params.polishing) {
          MINIMAP2_RACON ( FASTQ2FASTA.out.fasta )
          RACON ( MINIMAP2_RACON.out.draft_mapping)
          MEDAKA2 ( RACON.out.polished )
          //Remove trailing Ns and primer sequences from consensus sequence
          CUTADAPT ( MEDAKA2.out.consensus2.join(ch_primers) )
          consensus = CUTADAPT.out.trimmed
        }
        //If no polishing is skipped, directly use the clusters generated by Rattle for blast search
        else {
          consensus = FASTQ2FASTA.out.fasta2
        }

        //Limit blast homology search to a reference (legacy from ontvisc, placeholder at this stage, not tested in ont_amplicon)
        //if (params.blast_vs_ref) {
        //  BLASTN2REF ( consensus )
        //}

        //else {
        //Blast steps for samples targetting COI
        ch_coi_for_blast = (consensus.join(ch_coi))
        //Blast to COI database
        BLASTN_COI(ch_coi_for_blast)
        //Identify consensus that are in the wrong orientation and reverse complement them
        ch_revcomp = (consensus.join(BLASTN_COI.out.coi_blast_results))
        REVCOMP ( ch_revcomp )
        //Blast to NCBI nt database
        BLASTN ( REVCOMP.out.revcomp )

        //Directly blast to NCBI nt database all other samples
        ch_other_for_blast = (consensus.join(ch_other))
        BLASTN2 ( ch_other_for_blast )

        //Merge blast results from all samples
        ch_blast_merged = BLASTN.out.blast_results.mix(BLASTN2.out.blast_results.ifEmpty([]))

        //Extract top blast hit, assign taxonomy information to identify consensus that match target organism
        EXTRACT_BLAST_HITS ( ch_blast_merged.join(ch_targets) )
        //Add consensus sequence to blast results summary table
        FASTA2TABLE ( EXTRACT_BLAST_HITS.out.topblast.join(consensus) )

        //MAPPING BACK TO REFERENCE
        mapping_ch = (EXTRACT_BLAST_HITS.out.reference_fasta_files.join(REFORMAT.out.cov_derivation_ch))
        //Map filtered reads back to the reference sequence which was retrieved from blast search
        MINIMAP2_REF ( mapping_ch )
        //Derive bam file and consensus fasta file
        SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )

        //MAPPING BACK TO CONSENSUS
        mapping2consensus_ch = (EXTRACT_BLAST_HITS.out.consensus_fasta_files.join(REFORMAT.out.cov_derivation_ch))
        //Map filtered reads back to the portion of sequence which returned a blast hit
        MINIMAP2_CONSENSUS ( mapping2consensus_ch )
        //Derive bam file and coverage statistics
        SAMTOOLS_CONSENSUS ( MINIMAP2_CONSENSUS.out.aligned_sample )
        //Derive bed file for mosdepth to run coverage statistics
        PYFAIDX ( EXTRACT_BLAST_HITS.out.consensus_fasta_files )
        MOSDEPTH (SAMTOOLS_CONSENSUS.out.sorted_bams.join(PYFAIDX.out.bed))

        

        //Derive summary file presenting coverage statistics alongside blast results
        cov_stats_summary_ch = MOSDEPTH.out.mosdepth_results.join(EXTRACT_BLAST_HITS.out.consensus_fasta_files)
                                                             .join(SAMTOOLS_CONSENSUS.out.coverage)
                                                             .join(FASTA2TABLE.out.blast_results)
                                                             .join(QC_POST_DATA_PROCESSING.out.filtstats)
                                                             .join(ch_target_size)

        COVSTATS(cov_stats_summary_ch)

        SEQTK (SAMTOOLS_CONSENSUS.out.contig_seqids.join(final_fq))
        //EXTRACT_READ_LENGTHS ((SEQTK.out.fasta).join(SEQTK.out.contig_seqids).join(COVSTATS.out.detections_summary))

        files_for_report_ind_samples_ch = SAMTOOLS_CONSENSUS.out.sorted_bams.join(consensus)
                                                                            .join(QC_PRE_DATA_PROCESSING.out.rawnanoplot)
                                                                            .join(QC_POST_DATA_PROCESSING.out.filtnanoplot)
                                                                            .join(ch_blast_merged)
                                                                            .join(COVSTATS.out.detections_summary)

        files_for_report_global_ch = TIMESTAMP_START.out.timestamp
            .concat(QCREPORT.out.qc_report_html)
            .concat(QCREPORT.out.qc_report_txt)
            .concat(Channel.from(params.samplesheet).map { file(it) }).toList()
        HTML_REPORT(files_for_report_ind_samples_ch, files_for_report_global_ch)

      //DETECTION_REPORT(COVSTATS.out.detections_summary.collect().ifEmpty([]))
//        }
      }
/*
      //Perform direct alignment to a reference
      else if ( params.analysis_mode == 'map2ref') {
        MINIMAP2_REF ( final_fq )
        SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )
        MEDAKA ( SAMTOOLS.out.sorted_sample )
        FILTER_VCF ( MEDAKA.out.unfilt_vcf )
      }
*/    
      else {
        error("Analysis mode (clustering) not specified with e.g. '--analysis_mode clustering' or via a detectable config file.")
      }
      
    }
  }
//  TIMESTAMP_END ()
}
