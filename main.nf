#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage () {
    log.info """
    ont_amplicon
    Marie-Emilie Gauthier
    Cameron Hyde

    Usage:
    Run the command
    nextflow run main.nf -profile singularity -params-file {params.yml}

    Required arguments:
      --analysis_mode                 clustering, map2ref
                                      Default: 'clustering' [required]
      --analyst_name                  Name of the analyst
                                      Default: null [required]
      --facility                      Name of the facility where the analyst is performing the analysis
                                      Default: null [required]

    Optional arguments:
      --help                          Will print this usage document
      -resume                         Resume a failed run
      --outdir                        Path to save the output file
                                      'results'
      --samplesheet '[path/to/file]'  Path to the csv file that contains the list of
                                      samples to be analysed by this pipeline.
                                      Default:  'index.csv'. Needs to have .csv suffix


    Contents of index.csv:
      sampleid,fastq_path,target_organisms,target_gene,target_size,fwd_primer,rev_primer
      VE24-1279_COI,/work/tests/mtdt_data/barcode01_VE24-1279_COI/*fastq.gz,drosophilidae,COI,711,GGTCAACAAATCATAAAGATATTGG,ATTTTTTGGTCACCCTGAAGTTTA

      #### Pre-processing and QC options ####
      --merge                         Merge fastq files with the same sample name
                                      Default: true
      --qc_only                       Only perform preliminary QC step using Nanoplot
                                      Default: false
      --preprocessing_only            Only perform preprocessing steps specied
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

      #### Polishing ###
      --polishing                     Run polishing step
                                      Default: true

      #### Analysis mode and associated parameters ####
      ### Clustering (clustering) ###
      --rattle_clustering_min_length   Filter out reads shorter than this value
                                       Default: ''
      --rattle_clustering_max_length   Filter out reads longer than this value
                                       Default: ''
      --rattle_raw                     Use all the reads without any length filtering
                                       Default: false
      --rattle_clustering_max_variance Max allowed variance for two reads to be in the same gene cluster
                                       Default: '10000'
      --rattle_clustering_max_variance  Use all the reads without any length filtering
                                       Default: false
      --rattle_clustering_options      Rattle clustering options
                                       Default: ''
      --rattle_polishing_options       Rattle polishing options
                                       Default: ''

      #### Blast options ####
      --blast_threads                 Number of threads for megablast
                                      Default: '2'
      --blastn_db                     Path to blast database [required if not performing qc_only or preprocessing_only]
                                      Default: ''
      --blastn_COI                    Path to blast database for COI [required if performing analysis on COI gene]
                                      Default: ''
      --taxdump                       Path to taxonomykit database directory [required if not performing qc_only or preprocessing_only]
                                      Default: ''

      #### Mapping back to ref options ####
      --mapping_back_to_ref           Mapped back to reference blast match
                                      Default: 'true'

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

//if (params.reference != null) {
//    reference_name = file(params.reference).name
 //   reference_dir = file(params.reference).parent
//}
//if (params.host_fasta != null) {
//   host_fasta_dir = file(params.host_fasta).parent
//}

if (params.porechop_custom_primers == true) {
    porechop_custom_primers_dir = file(params.porechop_custom_primers_path).parent
}

def isNonEmptyFile(file) {
    return file.exists() && file.size() > 0
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
      bindbuild = (bindbuild + "-B ${params.taxdump} ")
    }
//    if (params.reference != null) {
//      bindbuild = (bindbuild + "-B ${reference_dir} ")
//    }
//   if (params.host_fasta != null) {
//      bindbuild = (bindbuild + "-B ${host_fasta_dir} ")
//    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}

process BLASTN {
  publishDir "${params.outdir}/${sampleid}/04_megablast", mode: 'copy', pattern: '{*_megablast_top_10_hits.txt,*_blast_status.txt}'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly)
  output:
    path("${sampleid}*_megablast_top_10_hits.txt")
    path("${sampleid}_blast_status.txt")
    tuple val(sampleid), path("${sampleid}*_megablast_top_10_hits.txt"), path("${sampleid}_blast_status.txt"), emit: blast_results


  script:
  def tmp_blast_output = assembly.getBaseName() + "_megablast_top_10_hits_temp.txt"
  def blast_output = assembly.getBaseName() + "_megablast_top_10_hits.txt"
  def status_file = sampleid + "_blast_status.txt"
    """
    STATUS="failed"
    echo "failed" > "${status_file}"
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${tmp_blast_output} \
      -evalue 1e-3 \
      -word_size 28 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length nident pident mismatch gaps gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
      -max_target_seqs 10

    cat <(printf "qseqid\tsgi\tsacc\tlength\tnident\tpident\tmismatch\tgaps\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tsstrand\tevalue\tbitscore\tqcovhsp\tstitle\tstaxids\tqseq\tsseq\tsseqid\tqcovs\tqframe\tsframe\n") ${tmp_blast_output} > ${blast_output}
    if [[ \$(wc -l < *_megablast_top_10_hits.txt) -ge 2 ]]
      then
        STATUS="passed"
        echo "passed" > "${status_file}"
    fi
    """

}

process BLASTN_COI {
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly), val(target_gene)
  output:
    tuple val(sampleid), path("${sampleid}_ids_to_reverse_complement.txt"), emit: coi_blast_results

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

    if [[ ! -s ${blast_output_COI} ]];
      then
        touch ${sampleid}_ids_to_reverse_complement.txt
    else
        grep minus ${blast_output_COI} | cut -f1 > ${sampleid}_ids_to_reverse_complement.txt
    fi
    """
}

process BLASTN2 {
  publishDir "${params.outdir}/${sampleid}/04_megablast", mode: 'copy', pattern: '{*_megablast_top_10_hits.txt,*_blast_status.txt}'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly), val(target_gene)
  output:
    path("${sampleid}*_megablast_top_10_hits.txt")
    path("${sampleid}_blast_status.txt")
    tuple val(sampleid), path("${sampleid}*_megablast_top_10_hits.txt"), path("${sampleid}_blast_status.txt"), emit: blast_results

  script:
  def tmp_blast_output = assembly.getBaseName() + "_megablast_top_10_hits_temp.txt"
  def blast_output = assembly.getBaseName() + "_megablast_top_10_hits.txt"
  def status_file = sampleid + "_blast_status.txt"
    """
    STATUS="failed"
    echo "failed" > "${status_file}"
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${tmp_blast_output} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length nident pident mismatch gaps gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
      -max_target_seqs 10

    cat <(printf "qseqid\tsgi\tsacc\tlength\tnident\tpident\tmismatch\tgaps\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tsstrand\tevalue\tbitscore\tqcovhsp\tstitle\tstaxids\tqseq\tsseq\tsseqid\tqcovs\tqframe\tsframe\n") ${tmp_blast_output} > ${blast_output}
    if [[ \$(wc -l < *_megablast_top_10_hits.txt) -ge 2 ]]
      then
        STATUS="passed"
        echo "passed" > "${status_file}"
    fi
    """
}

process CHOPPER {
  publishDir "${params.outdir}/${sampleid}/00_preprocessing/chopper", pattern: '*_chopper.log', mode: 'link'
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

process COVSTATS {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/05_mapping_to_consensus", mode: 'copy'

  input:
    tuple val(sampleid), path(bed), path(consensus), path(coverage), path(mapping_qual), path(top_hits), path(nanostats), val(target_size), path(reads_fasta), path(contig_seqids)
  output:
    path("*top_blast_with_cov_stats.txt")
    tuple val(sampleid), path("*top_blast_with_cov_stats.txt"), emit: detections_summary
    path("*top_blast_with_cov_stats.txt"), emit: detections_summary2

  script:
    """
    derive_coverage_stats.py --sample ${sampleid} --blastn_results ${top_hits} --nanostat ${nanostats} --coverage ${coverage} --bed ${bed} --target_size ${target_size} --contig_seqids ${contig_seqids} --reads_fasta ${reads_fasta} --consensus ${consensus} --mapping_quality ${mapping_qual}
    """
}
/*
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
*/
process CUTADAPT {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/03_polishing", pattern: '{*.fasta,*_cutadapt.log}', mode: 'copy'
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
//cutadapt -n 2 -j ${task.cpus} -g "${fwd_primer_trimmed};max_error_rate=0.1;min_overlap=10" -a "${rev_primer_trimmed};max_error_rate=0.1;min_overlap=10" --trim-n -o ${sampleid}_final_polished_consensus.fasta ${consensus} > ${sampleid}_cutadapt.log

    """
    if ${fwd_primer} != null && ${rev_primer} != null;
      then
        cutadapt -n 2 -j ${task.cpus} -g "${fwd_primer};max_error_rate=0.1" -a "${rev_primer};max_error_rate=0.1" --trim-n -o ${sampleid}_final_polished_consensus.fasta ${consensus} > ${sampleid}_cutadapt.log
    else
        cutadapt -n 2 -j ${task.cpus} --trim-n -o ${sampleid}_final_polished_consensus.fasta ${consensus} > ${sampleid}_cutadapt.log
    fi
    """
}


process EXTRACT_BLAST_HITS {
  tag "${sampleid}"
  label "setting_1"
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results), path(status), val(target_organism), val(target_gene), val(target_size)

  output:
    tuple val(sampleid), path("${sampleid}*_megablast_top_hits_tmp.txt"), emit: topblast
    tuple val(sampleid), path("${sampleid}_reference_match.fasta"), emit: reference_fasta_files
    tuple val(sampleid), path("${sampleid}_final_polished_consensus_match.fasta"), emit: consensus_fasta_files

  script:
    target_organism_str = (target_organism instanceof List)
    ? "\"${target_organism.join('|')}\""
    : "\"${target_organism}\""
    """
    if [[ \$(wc -l < *_megablast_top_10_hits.txt) -ge 2 ]]
      then
        select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}*_top_10_hits.txt --target_organism ${target_organism_str} --taxonkit_database_dir ${params.taxdump}

        # extract segment of consensus sequence that align to reference
        awk  -F  '\\t' 'NR>1 { printf ">%s\\n%s\\n",\$2,\$23 }' ${sampleid}*_top_hits_tmp.txt | sed 's/-//g' > ${sampleid}_final_polished_consensus_match.fasta

        # extract segment of reference that align to consensus sequence
        awk  -F  '\\t' 'NR>1 { printf ">%s_%s\\n%s\\n",\$2,\$4,\$24 }' ${sampleid}*_top_hits_tmp.txt | sed 's/-//g' > ${sampleid}_reference_match.fasta
    else
        echo "No hits found for ${sampleid} in the blast results. Skipping the extraction of consensus and reference fasta files." >&2
        touch ${sampleid}_final_polished_consensus_match.fasta
        touch ${sampleid}_reference_match.fasta
        touch ${sampleid}_megablast_top_hits_tmp.txt
    fi
    """
}

process FASTCAT {
  publishDir "${params.outdir}/${sampleid}/01_QC/fastcat", mode: 'copy'
  tag "${sampleid}"
  label "setting_1"

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_basecalling_model_inference.txt")

    tuple val(sampleid), path("${sampleid}.fastq.gz"), emit: merged

  script:
    """
    fastcat \
        -s ${sampleid} \
        -f ${sampleid}_stats.tsv \
        --histograms histograms \
        ${fastq} \
        | bgzip > ${sampleid}.fastq.gz
    zcat ${sampleid}.fastq.gz | head -n1 | sed 's/^.*basecall_model_version_id=//' > ${sampleid}_basecalling_model_inference.txt
    if (grep -q 'fast' ${sampleid}_basecalling_model_inference.txt) ; then
      echo "A Fast basecalling model was used for this sample! Please rerun the basecalling step using High accuracy." >&2
      exit 1
    fi
    """
}

process FASTQ2FASTA {
  publishDir "${params.outdir}/${sampleid}/02_clustering", mode: 'copy', pattern: '*_rattle.fasta'
  tag "${sampleid}"
  label "setting_1"

  input:
    tuple val(sampleid), path(fastq)

  output:
    tuple val(sampleid), path("${sampleid}.fasta"), emit: fasta

  script:
    """
    seqtk seq -A -C ${fastq} > ${sampleid}.fasta
    """
}

process CLUSTER2FASTA {
  publishDir "${params.outdir}/${sampleid}/02_clustering", mode: 'copy', pattern: '*_rattle.fasta'
  tag "${sampleid}"
  label "setting_1"

  input:
    tuple val(sampleid), path(fastq), path(assembly), path(status)

  output:
    tuple val(sampleid), path(fastq), path("${sampleid}_rattle.fasta"), emit: fasta
    tuple val(sampleid), path("${sampleid}_rattle.fasta"), emit: fasta2

  script:
    """
    if [[ ! -s ${assembly} ]];
      then
        echo "Clustering file is empty. Skipping conversion to FASTA." >&2
        touch ${sampleid}_rattle.fasta
    else
        echo "Converting assembly to FASTA format."
        cut -f1,3 -d ' ' ${assembly} | sed 's/ total_reads=/_RC/' > ${sampleid}_tmp.fastq
        seqtk seq -A -C ${sampleid}_tmp.fastq > ${sampleid}_rattle.fasta
    fi
    """
}

process FASTA2TABLE {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/04_megablast", mode: 'copy'

  input:
    tuple val(sampleid), path(tophits), path(fasta)
  output:
    file("${sampleid}*_megablast_top_hits.txt")
    tuple val(sampleid), file("${sampleid}*_megablast_top_hits.txt"), emit: blast_results

  script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${tophits}
    """
}

//If Racon polishing failed or produced an empty file, then Medaka will be run on the Rattle assembly.
//If Medaka polishing fails, then the script will run samtools consensus to produce a consensus sequence from the Rattle assembly.
process MEDAKA2 {
  publishDir "${params.outdir}/${sampleid}/03_polishing", mode: 'copy', pattern: '{*_consensus.fasta,*_consensus.fastq,*log}'
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(fastq), path(rattle_assembly), path(assembly)

  output:
   path("${sampleid}_medaka_consensus.fasta")
   path("${sampleid}_samtools_consensus.fasta")
   path("${sampleid}_samtools_consensus.fastq")
   path("${sampleid}_medaka.log"), optional: true
   path("${sampleid}_samtools_consensus.log"), optional: true
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), path("${sampleid}_medaka_consensus.bam"), path("${sampleid}_medaka_consensus.bam.bai"), path("${sampleid}_samtools_consensus.fasta")
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), path("${sampleid}_medaka_consensus.bam"), path("${sampleid}_medaka_consensus.bam.bai"), emit: consensus1
   tuple val(sampleid), path(rattle_assembly), path("${sampleid}_samtools_consensus.fasta"), env(STATUS), emit: consensus2

  script:
    """
    STATUS=failed
    if [[ ! -s ${rattle_assembly} ]];
      then
        touch ${sampleid}_medaka_consensus.fasta
        touch ${sampleid}_medaka_consensus.bam
        touch ${sampleid}_medaka_consensus.bam.bai
        touch ${sampleid}_samtools_consensus.fasta
        touch ${sampleid}_samtools_consensus.fastq
    else
      if [[ ! -s ${assembly} ]];
        then
          (
            set +eo pipefail
            medaka_consensus -i ${fastq} -d ${rattle_assembly} -t ${task.cpus} -o ${sampleid}
          ) 2>&1 | tee ${sampleid}_medaka.log
      else
          (
            set +eo pipefail
            medaka_consensus -i ${fastq} -d ${assembly} -t ${task.cpus} -o ${sampleid}
          ) 2>&1 | tee ${sampleid}_medaka.log
      fi

      if [[ ! -s ${sampleid}/consensus.fasta ]];
        then
          touch ${sampleid}_medaka_consensus.fasta
          touch ${sampleid}_medaka_consensus.bam
          touch ${sampleid}_medaka_consensus.bam.bai
          samtools consensus -f fasta -a -A -X r10.4_sup -o ${sampleid}_samtools_consensus.fasta ${rattle_assembly}
          samtools consensus -f fastq -a -A -X r10.4_sup -o ${sampleid}_samtools_consensus.fastq ${rattle_assembly}
      else
        cp ${sampleid}/calls_to_draft.bam ${sampleid}_medaka_consensus.bam
        cp ${sampleid}/calls_to_draft.bam.bai ${sampleid}_medaka_consensus.bam.bai
        cp ${sampleid}/consensus.fasta ${sampleid}_medaka_consensus.fasta
        samtools consensus -f fasta -a -A -X r10.4_sup -o ${sampleid}_samtools_consensus.fasta ${sampleid}_medaka_consensus.bam 2>&1 | tee ${sampleid}_samtools_consensus.log
        samtools consensus -f fastq -a -A -X r10.4_sup -o ${sampleid}_samtools_consensus.fastq ${sampleid}_medaka_consensus.bam
        if [[ ! -s ${sampleid}_samtools_consensus.fasta ]];
          then
            echo "Samtools consensus failed or produced an empty file." >> ${sampleid}_samtools_consensus.log
            touch ${sampleid}_samtools_consensus.fasta
            touch ${sampleid}_samtools_consensus.fastq
        else
          echo "Samtools consensus ran succesfully." >> ${sampleid}_samtools_consensus.log
          STATUS=passed
        fi
      fi
    fi
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
    if [[ ! -s ${consensus} ]]; then
      echo "Consensus file is empty or does not exist. Skipping minimap2 alignment." >&2
      touch ${sampleid}_aln.sam

    else
      minimap2 -ax map-ont -t ${task.cpus} --MD --sam-hit-only ${consensus} ${fastq} > ${sampleid}_aln.sam
    fi
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
    if [[ ! -s ${assembly} ]];
      then
        touch ${sampleid}_pre-racon.paf
    else
      minimap2 -L -x ava-ont -t ${task.cpus} ${assembly} ${fastq} > ${sampleid}_pre-racon.paf
    fi
    """
}

process MINIMAP2_REF {
  tag "${sampleid}"
  label 'setting_2'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(ref), path(fastq)

  output:
    tuple val(sampleid), path(ref), file("${sampleid}_ref_aln.sam"), emit: aligned_sample

  script:
    """
    minimap2 -ax map-ont --MD -t ${task.cpus} --sam-hit-only ${ref} ${fastq} > ${sampleid}_ref_aln.sam
    """
}

process MOSDEPTH {
  tag "$sampleid"
  label "setting_3"

  input:
    tuple val(sampleid), path(consensus), path(bam), path(bai), path(bed)

  output:
    tuple val(sampleid), path("${sampleid}.thresholds.bed"), emit: mosdepth_results

  script:
    """
    if [[ ! -s ${consensus} ]]; then
      touch ${sampleid}.thresholds.bed

    else
      mosdepth --by ${bed} --thresholds 30 -t ${task.cpus} ${sampleid} ${bam}
      gunzip *.per-base.bed.gz
      gunzip *.thresholds.bed.gz
    fi
    """
}

process PYFAIDX {
  tag "$sampleid"
  label "setting_3"

  input:
    tuple val(sampleid), path(fasta)

  output:
    tuple val(sampleid), path("${sampleid}.bed"), emit: bed

  script:
    """
    if [[ ! -s ${fasta} ]]; then
      touch ${sampleid}.bed
    else
      faidx --transform bed ${fasta} > ${sampleid}.bed
    fi
    """
}

process PORECHOP_ABI {
  tag "${sampleid}"
  publishDir "$params.outdir/${sampleid}/00_preprocessing/porechop",  mode: 'copy', pattern: '*_porechop.log'
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
  publishDir "${params.outdir}/00_QC_report", mode: 'copy', overwrite: true
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
    seq_run_qc_report.py
    """
}

//Polish draft consensus with `racon` (one round only). `racon` params are taken from epi2me-labs/wf-amplicon and pomoxis `mini_assemble`. 
//`--no-trimming` is added as it sometimes trims the consensus too aggressively. We trim the sequences downstream instead.

process RACON {
  publishDir "${params.outdir}/${sampleid}/03_polishing", mode: 'copy', pattern: '{*_racon_consensus.fasta,*.log}'
  tag "${sampleid}"
  label 'setting_2'

  input:
   tuple val(sampleid), path(fastq), path(rattle_assembly), path(paf)

  output:
   path("${sampleid}_racon_consensus.fasta")
   path("${sampleid}_racon.log"), optional: true
   tuple val(sampleid), path(fastq), path(rattle_assembly), path("${sampleid}_racon_consensus.fasta"), emit: polished

  script:
    """
    (
      set +eo pipefail
      if [[ ! -s ${rattle_assembly} ]];
        then
          touch ${sampleid}_racon_consensus.fasta
      else
          racon -m 8 -x -6 -g -8 -w 500 -t ${task.cpus} -q -1 --no-trimming -u \
              ${fastq} ${paf} ${rattle_assembly} \
              > ${sampleid}_racon_consensus_tmp.fasta
      fi
      
      ) 2>&1 | tee ${sampleid}_racon.log
      if [[ ! -s ${sampleid}_racon_consensus_tmp.fasta ]]; then
        echo "Racon polishing failed or produced an empty file." >> ${sampleid}_racon.log
        touch ${sampleid}_racon_consensus.fasta
      else
        cut -f1 -d ' ' ${sampleid}_racon_consensus_tmp.fasta > ${sampleid}_racon_consensus.fasta
        echo "Racon polishing ran succesfully." >> ${sampleid}_racon.log
      fi
    """
}

process RATTLE {
  tag "${sampleid}"
  label 'setting_10'
  publishDir "${params.outdir}/${sampleid}/02_clustering", mode: 'copy', pattern: '{*_rattle.log,*_rattle_status.txt}'

  input:
    tuple val(sampleid), path(fastq), val(target_size)

  output:
    file("${sampleid}_rattle.log")
    path("${sampleid}_rattle_status.txt")
    tuple val(sampleid), path("${sampleid}_rattle_status.txt"), emit: status
    tuple val(sampleid), path(fastq), path("transcriptome.fq"), path("${sampleid}_rattle_status.txt"), emit: clusters

  script:
  def status_file = sampleid + "_rattle_status.txt"
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
    echo "failed" > "${status_file}"
    (
      set +eo pipefail
      if [[ ${params.rattle_raw} == true ]]; then
        rattle cluster -i ${fastq} -t ${task.cpus} --raw ${rattle_clustering_options} -v ${params.rattle_clustering_max_variance} -o .
      else
        rattle cluster -i ${fastq} -t ${task.cpus} --lower-length ${rattle_clustering_min_length_set} --upper-length ${params.rattle_clustering_max_length} ${rattle_clustering_options} -v ${params.rattle_clustering_max_variance} -o .
      fi
      rattle cluster_summary -i ${fastq} -c clusters.out > ${sampleid}_cluster_summary.txt
      mkdir clusters
      rattle extract_clusters -i ${fastq} -c clusters.out -l ${sampleid} -o clusters --fastq
      rattle correct -t ${task.cpus} -i ${fastq} -c clusters.out -t ${task.cpus} -l ${sampleid}
      rattle polish -i consensi.fq -t ${task.cpus} --summary ${rattle_polishing_options}
      trap 'if [[ \$? -eq 139 ]]; then echo "segfault !"; fi' CHLD
    ) 2>&1 | tee ${sampleid}_rattle.log


    if [[ ! -s transcriptome.fq ]]; then
      touch transcriptome.fq
      echo "Rattle clustering and polishing failed." >> ${sampleid}_rattle.log
    else
      echo "Rattle clustering and polishing completed successfully." >> ${sampleid}_rattle.log
      echo "passed" > "${status_file}"
    fi
    """
}

process REFORMAT {
  tag "${sampleid}"
  label "setting_3"
  publishDir "$params.outdir/${sampleid}/00_preprocessing", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_preprocessed.fastq.gz")
    tuple val(sampleid), path("${sampleid}_preprocessed.fastq.gz"), emit: reformatted_fq
    tuple val(sampleid), path("${sampleid}_preprocessed.fastq.gz"), emit: cov_derivation_ch

  script:
    """
    reformat.sh in=${fastq} out=${sampleid}_preprocessed.fastq.gz trd qin=33
    """
}

process REVCOMP {
  publishDir "${params.outdir}/${sampleid}/04_megablast", mode: 'copy', pattern: '{*fasta}'
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
  publishDir "${params.outdir}/${sampleid}/06_mapping_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(ref), path(sample)

  output:
    path "${sampleid}_ref_aln.sorted.bam"
    path "${sampleid}_ref_aln.sorted.bam.bai"
    path "${sampleid}_samtools_consensus_from_ref.fasta"
    path "${sampleid}_coverage.txt"
    tuple val(sampleid), path(ref), path("${sampleid}_ref_aln.sorted.bam"), path("${sampleid}_ref_aln.sorted.bam.bai"), emit: sorted_sample

  script:
    """
    samtools view -Sb -F 4 ${sample} | samtools sort -o ${sampleid}_ref_aln.sorted.bam
    samtools index ${sampleid}_ref_aln.sorted.bam
    samtools coverage ${sampleid}_ref_aln.sorted.bam > ${sampleid}_coverage.txt
    samtools coverage -A -w 50 ${sampleid}_ref_aln.sorted.bam > ${sampleid}_histogram.txt
    samtools consensus -f fasta -a -A ${sampleid}_ref_aln.sorted.bam -X r10.4_sup -o ${sampleid}_samtools_consensus_from_ref.fasta
    """
}

process SAMTOOLS_CONSENSUS {
  publishDir "${params.outdir}/${sampleid}/05_mapping_to_consensus", mode: 'copy', pattern: '{*.bam,*.bai,*final_polished_consensus_match.*}'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(consensus), path(sample)

  output:
    path "${sampleid}_final_polished_consensus_match.fasta"
    path "${sampleid}_aln.sorted.bam"
    path "${sampleid}_aln.sorted.bam.bai"
    path "${sampleid}_final_polished_consensus_match.fastq"
    tuple val(sampleid), path(consensus), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_bams
    tuple val(sampleid), path("${sampleid}_coverage.txt"), emit: coverage
    tuple val(sampleid), path("${sampleid}_mapq.txt"), emit: mapping_quality
    tuple val(sampleid), path("${sampleid}_contigs_reads_ids.txt"), emit: contig_seqids
  script:
    """
    if [[ ! -s ${consensus} ]]; then
      touch ${sampleid}_contigs_reads_ids.txt
      touch ${sampleid}_aln.sorted.bam
      touch ${sampleid}_aln.sorted.bam.bai
      touch ${sampleid}_coverage.txt
      touch ${sampleid}_mapq.txt
      touch ${sampleid}_final_polished_consensus_match.fastq
    else
      samtools view -S -F 4 ${sample} | cut -f1,3 | sort | uniq > ${sampleid}_contigs_reads_ids.txt
      samtools view -Sb -F 4 ${sample} | samtools sort -o ${sampleid}_aln.sorted.bam
      samtools index ${sampleid}_aln.sorted.bam
      samtools coverage ${sampleid}_aln.sorted.bam  > ${sampleid}_coverage.txt
      samtools coverage -A -w 50 ${sampleid}_aln.sorted.bam > ${sampleid}_histogram.txt
      samtools view ${sampleid}_aln.sorted.bam | awk '{mapq[\$3]+=\$5; count[\$3]++} END {for (chr in mapq) printf "%s\\t%.2f\\n", chr, mapq[chr]/count[chr]}' > ${sampleid}_mapq.txt
      samtools consensus -f fastq -a -A -X r10.4_sup -o ${sampleid}_final_polished_consensus_match.fastq ${sampleid}_aln.sorted.bam
      samtools consensus -f pileup -a -A -X r10.4_sup -o ${sampleid}_final_polished_consensus_match.pileup ${sampleid}_aln.sorted.bam
    fi
    """
}

process TIMESTAMP_START {
  publishDir "${params.outdir}/01_pipeline_info", mode: 'copy', overwrite: true
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
  publishDir "${params.outdir}/${sampleid}/07_html_report", mode: 'copy', overwrite: true
  containerOptions "${bindOptions}"
  label 'setting_3'

  input:
    tuple val(sampleid), path(raw_nanoplot), path(filtered_nanoplot), path (rattle_status), path(consensus_fasta), path(top_blast_hits), path(blast_status), path(consensus_match_fasta), path(aln_sorted_bam), path(aln_sorted_bam_bai), path(blast_with_cov_stats),
    path(timestamp),
    path(qcreport_html),
    path(qcreport_txt),
    path(configyaml),
    path(samplesheet)

  output:
    path("*"), optional: true
    path("run_qc_report.html"), optional: true

  script:
  analyst_name = params.analyst_name.replaceAll(/ /, '_')
  facility = params.facility.replaceAll(/ /, '_')
    """
    cp ${qcreport_html} run_qc_report.html
    cp ${params.tool_versions} versions.yml
    cp ${params.default_params} default_params.yml

    build_report.py --samplesheet ${samplesheet} --result_dir . --params_file ${configyaml} --analyst ${analyst_name} --facility ${facility} --versions versions.yml --default_params_file default_params.yml
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
    if [[ ! -s ${contig_seqids} ]]; then
      touch ${sampleid}.fasta
    else
      seqtk seq -A -C ${fastq} > ${sampleid}_all_reads.fasta
      cut -f1 ${contig_seqids} | sort | uniq > reads_ids.txt
      seqtk subseq ${sampleid}_all_reads.fasta reads_ids.txt > ${sampleid}.fasta
    fi
    """
}

process SUBSAMPLE {
  tag "${sampleid}"
  label "setting_2"

  input:
    tuple val(sampleid), path(fastq)

  output:
    tuple val(sampleid), path("${sampleid}_downsampled.fastq.gz"), emit: subsampled_fq

  script:
    """
    seqkit sample -2 ${fastq} -n ${params.reads_downsampling_size} \
    | gzip  > ${sampleid}_downsampled.fastq.gz
    """
}

// Main workflow definition
include { NANOPLOT as QC_PRE_DATA_PROCESSING } from './modules.nf'
include { NANOPLOT as QC_POST_DATA_PROCESSING } from './modules.nf'
include { COPY_INPUTS } from './modules.nf'


workflow {
  TIMESTAMP_START ()
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      //.map{ row-> tuple((row.sampleid), file(row.fastq_path)) }
      .map { row ->
        // Check required fields
        if (!row.sampleid || !row.fastq_path)  {
          exit 1, "ERROR: samplesheet is missing required fields for sample_id."
        }
        else if (!row.fastq_path)  {
          exit 1, "ERROR: samplesheet is missing required field for fastq_path."
        }
        // Return parsed row
        tuple((row.sampleid), file(row.fastq_path)) }
      .set{ ch_sample }

    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row->
      // Loop through required fields and check if any are null or empty
        def requiredFields = ['sampleid', 'target_organism', 'target_gene', 'target_size']
        for (field in requiredFields) {
            def value = row[field]
            if (value == null || value.toString().trim() == '') {
                exit 1, "ERROR: samplesheet is missing or has empty value for required field '${field}'."
            }
        }

        def intFields = ['target_size']
        for (field in intFields) {
            try {
                row[field] = row[field].toInteger()
            } catch (Exception e) {
                throw new IllegalArgumentException("Invalid integer in field '${field}': '${row[field]}' in row: ${row}")
            }
        }
        def raw_val = row['target_organism']
          if (raw_val.contains(';')) {
              row['target_organism'] = raw_val.split(';')*.trim()
          } else {
              row['target_organism'] = [raw_val.trim()]
          }
        tuple((row.sampleid), (row.target_organism), (row.target_gene).replaceAll(' ', '_'), (row.target_size)) }
      .set{ ch_targets }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.target_size)) }
      .set{ ch_target_size }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.fwd_primer), (row.rev_primer)) }
      .set{ ch_primers }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map { row ->
          def sampleid = row.sampleid
          def gene = row.target_gene?.trim()?.toUpperCase()
          tuple(sampleid, gene)
      }
      .filter { sampleid, gene ->
          ['COI', 'CO1'].contains(gene)
      }
      .set { ch_coi }

    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.target_gene)) }
      .filter { sampleid, target_gene -> !['COI', 'CO1'].contains(target_gene?.toUpperCase()) }
      .set{ ch_other }

        
    // Workaround for cloudgene deleting input files after workflow execution
    COPY_INPUTS(params.samplesheet)

  } else { exit 1, "Input samplesheet file not specified!" }

  configyaml = Channel.fromPath(workflow.commandLine.split(" -params-file ")[1].split(" ")[0])

// Conditional logic to require database if needed
  ch_coi
    .collect()
    .subscribe { items ->
        if (items && !items.isEmpty()) {
            if (params.blastn_COI == null) {
                error("Please provide path to the MetaCOXI database.")
            }
        } else {
            println "[Info] No COI/CO1 samples detected — database not required."
        }
    }

  if ( params.analyst_name == null) {
    error("Please provide the name of the analyst who is performing the analysis.")
  }

  if ( params.facility == null) {
    error("Please provide the name of the facility where the analysis was performed.")
      }

  if ( params.analysis_mode == 'clustering') {
    if (!params.blast_vs_ref & !params.qc_only & !params.preprocessing_only) {
      if ( params.blastn_db == null) {
        error("Please provide the path to a blast database using the parameter --blastn_db.")
      }
      if ( params.taxdump == null) {
        error("Please provide the path to a taxonkit database using the parameter --taxdump.")
      }
    }

  }
  /*
  else if ( params.analysis_mode == 'map2ref' ) {
    if ( params.reference == null) {
      error("Please provide the path to a reference fasta file with the parameter --reference.")
      }
  }
  */
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
    PORECHOP_ABI ( fq )
    trimmed_fq = PORECHOP_ABI.out.porechopabi_trimmed_fq

    // Perform quality filtering of reads using chopper
    if (params.qual_filt) {
      CHOPPER ( trimmed_fq)
      filtered_fq = CHOPPER.out.chopper_filtered_fq
    }
    else { filtered_fq = trimmed_fq
    }

    //Reformat fastq read names after the first whitespace
    REFORMAT( filtered_fq )

    //Run Nanoplot on merged raw fastq files after data processing
    QC_POST_DATA_PROCESSING ( filtered_fq )

    if (params.subsample) {
      SUBSAMPLE ( REFORMAT.out.reformatted_fq )
      final_fq = SUBSAMPLE.out.subsampled_fq
    }
    else {
      final_fq = REFORMAT.out.reformatted_fq
    }

    //Derive QC report
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
    QCREPORT(ch_multiqc_files.collect())

    if (!params.preprocessing_only) {
      //Currently only one analysis mode in ont_amplicon, consider removing if no other mode is added to this pipeline
      //We have had talks about including an option to just a map to a reference of interest, but this is not implemented yet
      if ( params.analysis_mode == 'clustering' ) {
        //Perform clustering using Rattle and convert to fasta file
        //Branch outputs of Rattle into passed and failed
        //If the clustering step succeeds, it will proceed to the polishing step
        ch_fq_target_size = (final_fq.join(ch_target_size))
        RATTLE ( ch_fq_target_size )
        CLUSTER2FASTA ( RATTLE.out.clusters )

        //Polish consensus sequence using Racon followed by Medaka and samtools consensus
        if (params.polishing) {
          MINIMAP2_RACON ( CLUSTER2FASTA.out.fasta )
          RACON ( MINIMAP2_RACON.out.draft_mapping)
          MEDAKA2 ( RACON.out.polished )
          ch_branched = MEDAKA2.out.consensus2
            | branch { sampleid, clusters, consensus, status ->
              passed: status == "passed"
              failed: status == "failed"
            }

          ch_passed = ch_branched.passed
            | map { sampleid, clusters, consensus, status -> [sampleid, consensus] }
          ch_failed = ch_branched.failed
            | map { sampleid, clusters, consensus, status -> [sampleid, clusters] }

          consensus = ch_passed.concat(ch_failed.ifEmpty([]))
        }
        //If polishing is skipped, directly use the clusters generated by Rattle for blast search
        else {
          consensus = CLUSTER2FASTA.out.fasta2
        }

        //Remove trailing Ns and primer sequences from consensus sequence
        CUTADAPT ( consensus.join(ch_primers) )

        //Blast steps for samples targetting COI
        ch_coi_for_blast = (CUTADAPT.out.trimmed.join(ch_coi))
        //Blast to COI database
        BLASTN_COI(ch_coi_for_blast)
        //Identify consensus that are in the wrong orientation and reverse complement them
        ch_revcomp = (CUTADAPT.out.trimmed.join(BLASTN_COI.out.coi_blast_results))
        REVCOMP ( ch_revcomp )
        //Blast to NCBI nt database
        BLASTN ( REVCOMP.out.revcomp )

        //Directly blast to NCBI nt database all other samples
        ch_other_for_blast = (CUTADAPT.out.trimmed.join(ch_other))
        BLASTN2 ( ch_other_for_blast )

        //Merge blast results from all samples
        ch_blast_merged = BLASTN.out.blast_results.mix(BLASTN2.out.blast_results.ifEmpty([]))

        //ch_blast_merged2 = ch_blast_merged.map { sampleid, blast_results, status -> [sampleid, blast_results] }

        //Extract top blast hit, assign taxonomy information to identify consensus that match target organism
        EXTRACT_BLAST_HITS ( ch_blast_merged.join(ch_targets) )
        //Add consensus sequence to blast results summary table
        FASTA2TABLE ( EXTRACT_BLAST_HITS.out.topblast.join(consensus) )

        //MAPPING BACK TO CONSENSUS
        mapping2consensus_ch = (EXTRACT_BLAST_HITS.out.consensus_fasta_files.join(REFORMAT.out.cov_derivation_ch))
        //Map filtered reads back to the portion of sequence which returned a blast hit
        MINIMAP2_CONSENSUS ( mapping2consensus_ch )
        //Derive bam file and coverage statistics
        SAMTOOLS_CONSENSUS ( MINIMAP2_CONSENSUS.out.aligned_sample )
        //Derive bed file for mosdepth to run coverage statistics
        PYFAIDX ( EXTRACT_BLAST_HITS.out.consensus_fasta_files )
        MOSDEPTH (SAMTOOLS_CONSENSUS.out.sorted_bams.join(PYFAIDX.out.bed))
        SEQTK (SAMTOOLS_CONSENSUS.out.contig_seqids.join(final_fq))
        //Derive summary file presenting coverage statistics alongside blast results
        cov_stats_summary_ch = MOSDEPTH.out.mosdepth_results.join(EXTRACT_BLAST_HITS.out.consensus_fasta_files)
                                                             .join(SAMTOOLS_CONSENSUS.out.coverage)
                                                             .join(SAMTOOLS_CONSENSUS.out.mapping_quality)
                                                             .join(FASTA2TABLE.out.blast_results)
                                                             .join(QC_POST_DATA_PROCESSING.out.filtstats)
                                                             .join(ch_target_size)
                                                             .join(SEQTK.out.fasta)
                                                             .join(SEQTK.out.contig_seqids)

        COVSTATS(cov_stats_summary_ch)

        files_for_report_ind_samples_ch = QC_PRE_DATA_PROCESSING.out.rawnanoplot.join((QC_POST_DATA_PROCESSING.out.filtnanoplot)
                                                                                .join(RATTLE.out.status)
                                                                                .join(CUTADAPT.out.trimmed)
                                                                                .join(ch_blast_merged)
                                                                                .join(SAMTOOLS_CONSENSUS.out.sorted_bams)
                                                                                .join(COVSTATS.out.detections_summary))
        files_for_report_global_ch = TIMESTAMP_START.out.timestamp
            .concat(QCREPORT.out.qc_report_html)
            .concat(QCREPORT.out.qc_report_txt)
            .concat(configyaml)
            .concat(Channel.from(params.samplesheet).map { file(it) }).toList()

        HTML_REPORT(files_for_report_ind_samples_ch
            .combine(files_for_report_global_ch))

        //MAPPING BACK TO REFERENCE
        if (params.mapping_back_to_ref) {
          mapping_ch = (EXTRACT_BLAST_HITS.out.reference_fasta_files.join(REFORMAT.out.cov_derivation_ch))
          //Map filtered reads back to the reference sequence which was retrieved from blast search
          MINIMAP2_REF ( mapping_ch )
          //Derive bam file and consensus fasta file
          SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )
        }
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
}
