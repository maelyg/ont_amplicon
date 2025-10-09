if (params.blastn_db != null) {
    blastn_db_name = file(params.blastn_db).name
    blastn_db_dir = file(params.blastn_db).parent
}

//if (params.host_fasta != null) {
//    host_fasta_dir = file(params.host_fasta).parent
//}

switch (workflow.containerEngine) {
  case "singularity":
    bindbuild = "";
    if (params.blastn_db != null) {
      bindbuild = (bindbuild + "-B ${blastn_db_dir} ")
    }
//    if (params.host_fasta != null) {
//      bindbuild = (bindbuild + "-B ${host_fasta_dir} ")
//    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}

process NANOPLOT {
  publishDir "${params.outdir}/${sampleid}/01_QC/nanoplot",  pattern: '{*NanoPlot-report.html}', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/01_QC/nanoplot",  pattern: '{*NanoStats.txt}', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/01_QC/nanoplot",  pattern: '{*LengthvsQualityScatterPlot_dot.html}', mode: 'link'
  tag "${sampleid}"
  label "setting_2"

  input:
    tuple val(sampleid), path(sample)
  output:
    path("*NanoPlot-report.html"), optional: true
    path("*NanoStats.txt"), optional: true
    path("*LengthvsQualityScatterPlot_dot.html"), optional: true
    path("*NanoStats.txt"), emit: read_counts
    tuple val(sampleid), path("${sampleid}_filtered_NanoStats.txt"), emit: filtstats, optional: true
    tuple val(sampleid), path("${sampleid}_raw_NanoPlot-report.html"), emit: rawnanoplot, optional: true
    tuple val(sampleid), path("${sampleid}_filtered_NanoPlot-report.html"), emit: filtnanoplot, optional: true

  
  script:
  def fastq = sample.getBaseName() + ".fastq.gz"
  """
  
  if [[ ${sample} == *trimmed.fastq.gz ]] || [[ ${sample} == *filtered.fastq.gz ]] ;
  then
    if [ -n "\$(gunzip < ${sample} | head -n 1 | tr '\0\n' __)" ];
    then
        NanoPlot -t ${task.cpus} --fastq ${sample} --prefix ${sampleid}_filtered_ --plots dot --N50 --tsv_stats
    else
        echo "Metrics dataset\nnumber_of_reads\t0" > ${sampleid}_filtered_NanoStats.txt
        touch ${sampleid}_filtered_LengthvsQualityScatterPlot_dot.html
        touch ${sampleid}_filtered_NanoPlot-report.html
    fi
  else
    NanoPlot -t ${task.cpus} --fastq ${sample} --prefix ${sampleid}_raw_ --plots dot --N50 --tsv_stats
  fi
  """
}

/*
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
    cp ${blastn_db_dir}/taxdb.btd .
    cp ${blastn_db_dir}/taxdb.bti .
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${tmp_blast_output} \
      -evalue 1e-3 \
      -word_size 28 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length nident pident mismatch gaps gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames sskingdoms' \
      -max_target_seqs 10
    
    cat <(printf "qseqid\tsgi\tsacc\tlength\tnident\tpident\tmismatch\tgaps\tgapopen\tqstart\tqlen\tqend\tsstart\tsend\tslen\tsstrand\tevalue\tbitscore\tqcovhsp\tstitle\tstaxids\tqseq\tsseq\tsseqid\tqcovs\tqframe\tsframe\tspecies\tsskingdoms\n") ${tmp_blast_output} > ${blast_output}
    """
  }
}
*/

process COPY_INPUTS {
    tag "$samplesheet_file"
    
    input:
    path samplesheet_file
    
    output:
    path "samplesheet.csv", emit: samplesheet
    
    script:
    """
    # Copy input files to work directory to ensure they remain available
    # after the workflow execution (workaround for cloudgene deleting input files after workflow execution)
    # Use different temporary names first to avoid same-file copy errors
    
    if [ "${samplesheet_file}" != "samplesheet.csv" ]; then
        cp "${samplesheet_file}" samplesheet.csv
    else
        # File already has correct name, just ensure it exists
        if [ ! -f samplesheet.csv ]; then
            echo "ERROR: samplesheet file not found" >&2
            exit 1
        fi
    fi
    
    echo "Successfully prepared input files in work directory"
    echo "samplesheet file: \$(wc -l < samplesheet.csv) lines"
    """
}