# ont_amplicon

## Introduction

ont_amplicon is a Nextflow-based bioinformatics pipeline designed to derive consensus sequences from:
 amplicon sequencing data that was generated using rapid library preparation kit from Oxford nanopore Technologies. 

It takes compressed fastq files as input.


## Pipeline overview
- Data quality check (QC) and preprocessing
  - Merge fastq files (Fascat, optional)
  - Raw fastq file QC (Nanoplot)
  - Trim adaptors (PoreChop ABI - optional)
  - Filter reads based on length and/or quality (Chopper - optional)
  - Reformat fastq files so read names are trimmed after the first whitespace (bbmap)
  - Processed fastq file QC (if PoreChop and/or Chopper is run) (Nanoplot)
- QC report
  - Derive read counts recovered pre and post data processing and post host filtering
- Clustering mode
  - Read clustering (Rattle)
  - Convert fastq to fasta format (seqtk)
  - Polishing (Minimap, Racon, Medaka2, samtools - optional)
  - Remove adapters if provided (Cutadapt)
  - Megablast homology search against COI database (if COI is targetted) and reverse complement where required
  - Megablast homology search against NCBI database
  - Derive top candidate hits, assign preliminary taxonomy and target organism flag (pytaxonkit)
  - Map reads back to segment of consensus sequence that align to reference and derive BAM file and alignment statistics (minimap, samtools and mosdepth)
  - Map reads to segment of NCBI reference sequence that align to consensus and derive BAM file and consensus (minimap, samtools)


## Installation
### Requirements  

If the pipeline is run on a local machine, it will require between 300-800Gb of space for the installation of required containers and databases alone. This includes:  
- 80 Mb ont_amplicon pipeline  
- ~ 3.8Gb for containers  
- 600Mb for taxonkit  
- 280Gb/760Gb for the blast NCBI database coreNT/NT
- 3.4Gb for the MetaCOXI database  

To run the pipeline will also require at least 2 cores and ~40Gb of memory per sample.  
The pipeline will generate ~5-100Mb of files per sample, depending on the number of consensuses recovered per sample and if mapping back to reference is required. Make sure you have enough space available on your local machine before running several samples at the same time.  

1. Install Java if not already on your system. Follow the java download instructions provided on this page [`page`](https://www.nextflow.io/docs/latest/getstarted.html#installation).

2. Install Nextflow [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

3. Install [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) to suit your environment. The pipeline has been validated using singularity version 3.10.2-1 and apptainer version 1.3.6-1.el9 but has not yet been tested with singularity version 4.

3. Install taxonkit using the script install_taxonkit.sh or follow the steps described on this page [`page`](https://bioinf.shenwei.me/taxonkit/download/).


4. Install NCBI NT or coreNT.  
Download a local copy of the NCBI database of interest, following the detailed steps available at https://www.ncbi.nlm.nih.gov/books/NBK569850/. Create a folder where you will store your NCBI databases. It is good practice to include the date of download. For instance:
  ```
  mkdir blastDB/20230930
  ```
  You will need to use a current update_blastdb.pl script from the blast+ version used with the pipeline (ie 2.16.0).
  For example:
  ```
  singularity exec -B /scratch https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4  update_blastdb.pl --decompress nt
  singularity exec -B /scratch https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4  update_blastdb.pl --decompress taxdb
  tar -xzf taxdb.tar.gz
  ```
  
  Specify the path of your local NCBI blast directories in your nextflow command using ```--blastn_db = '/full/path/to/blastDB/20230930/nt'``` or specify the following lines in a user config file.
  For instance:
  ```
  params {
    --blastn_db = '/full/path/to/blastDB/20230930/nt'
  }
  ```

5. Download the Cytochrome oxydase 1 (COI1) database if you are planning to analyse COI samples.
  ```
  git clone https://github.com/bachob5/MetaCOXI.git
  #extract MetaCOXI_Seqs.fasta from the MetaCOXI_Seqs.tar.gz file
  tar -xvf MetaCOXI_Seqs.tar.gz
  #make a blast database from the fasta file
  singularity exec -B /scratch https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4 makeblastdb -in MetaCOXI_Seqs.fasta -parse_seqids -dbtype prot
  ```
Specify the path of your COI database in your nextflow command using ```--blastn_COI = '/full/path/to/MetaCOXI_Seqs.fasta'``` or specify the following lines in a user config file.
  For instance:
  ```
  params {
    --blastn_COI = '/full/path/to/MetaCOXI_Seqs.fasta'
  }
  ```

## Running the pipeline  

### Run the pipeline for the first time
- Run the command:
  ```
  nextflow run maelyg/ont_amplicon -profile singularity --samplesheet index.csv
  ```
  The first time the command runs, it will download the pipeline into your assets.  

  The source code can also be downloaded directly from GitHub using the git command:
  ```
  git clone https://github.com/maelyg/ont_amplicon
  ```

- Provide an index.csv file.  
  Create a **comma separated file (csv)** that will be the input for the workflow. 
  
  **Please note**: it is best to edit the csv file with an editor that does not add special characters/symbols (e.g. VSCode or Atom). If using other editors, check your files and if necessary, run dos2unix[`dos2unix`](https://www.linuxfromscratch.org/blfs/view/git/general/dos2unix.html) on the file to get rid of these unwanted characters/symbols as they will cause the pipeline to fail.  
  
  By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the ```--samplesheet [filename]``` in the nextflow run command, as long as it has a **.csv** suffix. This text file requires the following columns (which need to be included as a header): ```sampleid,sample_files,spp_targets,gene_targets,target_size,fwd_primer,rev_primer``` 

   - **sampleid** will be the sample name that will be given to the files created by the pipeline (required).  
   - **sample_path** is the full path to the fastq files that the pipeline requires as starting input (required).  
   - **spp_targets** is the organism targetted by the PCR (required).  
   - **gene_targets** is the gene targetted by the PCR (optional).  
   - **target_size** is the expected size of the amplicon (required).  
   - **fwd_primer** is the nucleotide sequence of the FWD primer (optional).  
   - **rev_primer** is the nucleotide sequence of the REV primer (optional).  



  For the fastq files path, the pipeline is currently expecting either 1) multiple fastq.gz files per sample located within one folder or 2) a single fastq.gz file per sample.  
  If there are **multiple fastq.gz files per sample**, their full path can be specified on one line using **an asterisk (*fastq.gz)** and you will need to specify the parameter ```--merge``` either on the command line or a config file.  
  See an example of an index.csv file for 2 MTDT samples:  
  ```
  sampleid,sample_files,spp_targets,gene_targets,target_size,fwd_primer,rev_primer
  VE24-1279_COI,tests/mtdt_data/barcode01_VE24-1279_COI/*fastq.gz,drosophilidae,COI,711,GGTCAACAAATCATAAAGATATTGG,ATTTTTTGGTCACCCTGAAGTTTA
  MP24-1051A_16S,tests/mtdt_data/barcode06_MP24-1051A_16S/*fastq.gz,bacteria,16s,1509,AGAGTTTGATCATGGCTCAG,AAGTCGTAACAAGGTAACCGT
  MP24-1096B_gyrB,tests/mtdt_data/barcode19_MP24-1096B_gyrB/*fastq.gz,bacteria,gyrB,1258,GAAGTCATCATGACCGTTCTGCAYGCNGGNGGNAARTTYGA,ATGACNGAYGCNGAYGTNGAYGGCTCGCACATCCGTACCCTGCT
  ```
  For samples with a single fastq.gz file, specify **the full path to the fastq.gz file.**


- Specify a profile:
  ```
  nextflow run maelyg/ont_amplicon -profile singularity --samplesheet index_example.csv
  ```
  setting the profile parameter to one of ```docker``` or ```singularity``` to suit your environment.
  
- Specify the analysis mode: ```--analysis_mode clustering``` (this is set to clustering by default; haivng this parameter in place will enable us to add other analysis modes like mapping to ref down the track if required).  

- Specify the ``--analyst_name`` and the ``--facility`` either on the nextflow command or in a config file.  The analysis cannot proceed without these being set.  

- To set additional parameters, you can either include these in your nextflow run command:
  ```
  nextflow run maelyg/ont_amplicon -profile singularity --samplesheet index_example.csv --adapter_trimming
  ```
  or set them to true in the nextflow.config file.
  ```
  params {
    adapter_trimming = true
  }
  ```
### Run tests
Two tests are currently provided to check if the pipeline was successfully installed and to demonstrate to users the current outputs generated by the pipeline prototype so they can provide feedback. The mtdt_test runs three samples provided by  MTDT (barcode01_VE24-1279_COI, barcode06_MP24-1051A_16S and barcode19_MP24-1096B_gyrB). The peq_test runs two samples provided by PEQ (ONT141 and ONT142). A small NCBI blast database and COI database have been derived to speed up the analysis run in test mode.  

To use the tests, change directory to ont_amplicon and run the following command for the MTDT test:
  ```
  nextflow run main.nf -profile mtdt_test,singularity
  ```
  and this command for the PEQ test:  
  ```
  nextflow run main.nf -profile peq_test,singularity
  ```
The tests should take less than 5 minutes to run to completion.

If the installation is successful, your screen should output something similar to this, with a completion and success status at the bottom:
```
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [berserk_gutenberg] DSL2 - revision: 14779efa3e

executor >  pbspro (76)
[8f/37f88f] process > TIMESTAMP_START                                     [100%] 1 of 1 ✔
[9a/4ae987] process > FASTCAT (barcode19_MP24-1096B_gyrB)                 [100%] 3 of 3 ✔
[4b/c28438] process > QC_PRE_DATA_PROCESSING (barcode19_MP24-1096B_gyrB)  [100%] 3 of 3 ✔
[3a/3245aa] process > PORECHOP_ABI (barcode19_MP24-1096B_gyrB)            [100%] 3 of 3 ✔
[be/1a6a24] process > CHOPPER (barcode19_MP24-1096B_gyrB)                 [100%] 3 of 3 ✔
[e1/178003] process > REFORMAT (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
[27/dfdee1] process > QC_POST_DATA_PROCESSING (barcode19_MP24-1096B_gyrB) [100%] 3 of 3 ✔
[ed/547789] process > QCREPORT                                            [100%] 1 of 1 ✔
[bf/0b1014] process > RATTLE (barcode19_MP24-1096B_gyrB)                  [100%] 3 of 3 ✔
[71/803823] process > FASTQ2FASTA (barcode19_MP24-1096B_gyrB)             [100%] 3 of 3 ✔
[6c/8c02c0] process > MINIMAP2_RACON (barcode19_MP24-1096B_gyrB)          [100%] 3 of 3 ✔
[42/2b5a12] process > RACON (barcode19_MP24-1096B_gyrB)                   [100%] 3 of 3 ✔
[13/02e7b3] process > MEDAKA2 (barcode19_MP24-1096B_gyrB)                 [100%] 3 of 3 ✔
[9f/10d248] process > CUTADAPT (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
[8d/636107] process > BLASTN_COI (barcode01_VE24-1279_COI)                [100%] 1 of 1 ✔
[77/5b10fc] process > REVCOMP (barcode01_VE24-1279_COI)                   [100%] 1 of 1 ✔
[c4/119cf8] process > BLASTN (barcode01_VE24-1279_COI)                    [100%] 1 of 1 ✔
[fe/284a69] process > BLASTN2 (barcode19_MP24-1096B_gyrB)                 [100%] 2 of 2 ✔
[8e/7623d2] process > EXTRACT_BLAST_HITS (barcode19_MP24-1096B_gyrB)      [100%] 3 of 3 ✔
[ee/d92377] process > FASTA2TABLE (barcode19_MP24-1096B_gyrB)             [100%] 3 of 3 ✔
[25/32f878] process > MINIMAP2_REF (barcode19_MP24-1096B_gyrB)            [100%] 3 of 3 ✔
[9b/21ada5] process > SAMTOOLS (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
[01/69286e] process > MINIMAP2_CONSENSUS (barcode19_MP24-1096B_gyrB)      [100%] 3 of 3 ✔
[bf/351d07] process > SAMTOOLS_CONSENSUS (barcode19_MP24-1096B_gyrB)      [100%] 3 of 3 ✔
[0b/783c9e] process > PYFAIDX (barcode19_MP24-1096B_gyrB)                 [100%] 3 of 3 ✔
[bb/c64067] process > MOSDEPTH (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
[4e/20a073] process > COVSTATS (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
[6b/b98138] process > SEQTK (barcode19_MP24-1096B_gyrB)                   [100%] 3 of 3 ✔
[73/fbfcfe] process > HTML_REPORT (3)                                     [100%] 3 of 3 ✔
Completed at: 03-Apr-2025 09:47:49
Duration    : 5m 49s
CPU hours   : 0.2
Succeeded   : 76
```

And the pipeline will have generated a results folder with the following structure:
```
results/
├── barcode01_VE24-1279_COI
│   ├── clustering
│   │   └── barcode01_VE24-1279_COI_rattle.fasta
│   ├── html_report
│   │   ├── bam-alignment.html
│   │   ├── example_report_context.json
│   │   └── report.html
│   ├── mapping_to_consensus
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam.bai
│   │   ├── barcode01_VE24-1279_COI_coverage.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_match.fasta
│   │   ├── barcode01_VE24-1279_COI_histogram.txt
│   │   ├── barcode01_VE24-1279_COI.per-base.bed
│   │   └── barcode01_VE24-1279_COI_top_blast_with_cov_stats.txt
│   ├── mapping_to_ref
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam.bai
│   │   ├── barcode01_VE24-1279_COI_coverage.txt
│   │   ├── barcode01_VE24-1279_COI_histogram.txt
│   │   ├── barcode01_VE24-1279_COI_reference_match.fasta
│   │   └── barcode01_VE24-1279_COI_samtools_consensus_from_ref.fasta
│   ├── megablast
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_match.fasta
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_megablast_COI_top_hit.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc.fasta
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc_megablast_top_10_hits_temp.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc_megablast_top_10_hits.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc_megablast_top_hits.txt
│   │   └── barcode01_VE24-1279_COI_reference_match.fasta
│   ├── polishing
│   │   ├── barcode01_VE24-1279_COI_cutadapt.log
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus.fasta
│   │   ├── barcode01_VE24-1279_COI_medaka_consensus.bam
│   │   ├── barcode01_VE24-1279_COI_medaka_consensus.bam.bai
│   │   ├── barcode01_VE24-1279_COI_medaka_consensus.fasta
│   │   ├── barcode01_VE24-1279_COI_preprocessed.fastq.gz
│   │   ├── barcode01_VE24-1279_COI_racon_polished.fasta
│   │   └── barcode01_VE24-1279_COI_samtools_consensus.fasta
│   ├── preprocessing
│   │   ├── barcode01_VE24-1279_COI_basecalling_model_inference.txt
│   │   ├── barcode01_VE24-1279_COI_preprocessed.fastq.gz
│   │   ├── chopper
│   │   │   └── barcode01_VE24-1279_COI_chopper.log
│   │   └── porechop
│   │       └── barcode01_VE24-1279_COI_porechop.log
│   └── qc
│       ├── fastcat
│       │   ├── barcode01_VE24-1279_COI.fastq.gz
│       │   ├── barcode01_VE24-1279_COI_stats.tsv
│       │   └── histograms
│       │       ├── length.hist
│       │       └── quality.hist
│       └── nanoplot
│           ├── barcode01_VE24-1279_COI_filtered_LengthvsQualityScatterPlot_dot.html
│           ├── barcode01_VE24-1279_COI_filtered_NanoPlot-report.html
│           ├── barcode01_VE24-1279_COI_filtered_NanoStats.txt
│           ├── barcode01_VE24-1279_COI_raw_LengthvsQualityScatterPlot_dot.html
│           ├── barcode01_VE24-1279_COI_raw_NanoPlot-report.html
│           └── barcode01_VE24-1279_COI_raw_NanoStats.txt
└── qc_report
    ├── run_qc_report_20250401-210340.html
    └── run_qc_report_20250401-210340.txt
```

### QC step
By default the pipeline will run a quality control check of the raw reads using NanoPlot.

- It is recommended to first run only the quality control step to have a preliminary look at the data before proceeding with downstream analyses by specifying the ```--qc_only``` parameter.

The command you would run would look like this:
```
nextflow run maelyg/ont_amplicon -profile singularity \
                            --merge \
                            --qc_only
```

### Preprocessing reads
If multiple fastq files exist for a single sample, they will first need to be merged using the `--merge` option using [`Fascat`](https://github.com/epi2me-labs/fastcat).
Then the read names of the fastq file created will be trimmed after the first whitespace, for compatiblity purposes with all downstream tools.  

Reads can also be optionally trimmed of adapters and/or quality filtered:  
- Search for presence of sequencing adapters in sequences reads using [`Porechop ABI`](https://github.com/rrwick/Porechop) by specifying the ``--adapter_trimming`` parameter. Porechop ABI parameters can be specified using ```--porechop_options '{options} '```, making sure you leave a space at the end before the closing quote. Please refer to the Porechop manual.  

  **Special usage:**  
  To limit the search to known adapters listed in [`adapter.py`](https://github.com/bonsai-team/Porechop_ABI/blob/master/porechop_abi/adapters.py), just specify the ```--adapter_trimming``` option.  
  To search ab initio for adapters on top of known adapters, specify ```--adapter_trimming --porechop_options '-abi '```.  
T  o limit the search to custom adapters, specify ```--adapter_trimming --porechop_custom_primers --porechop_options '-ddb '``` and list the custom adapters in the text file located under bin/adapters.txt following the format:  
    ```
     line 1: Adapter name
     line 2: Start adapter sequence
     line 3: End adapter sequence
     --- repeat for each adapter pair---
    ```

- Perform a quality filtering step using [`Chopper`](https://github.com/wdecoster/chopper) by specifying the ```--qual_filt``` parameter. The following parameters can be specified using the ```--chopper_options '{options}'```. Please refer to the Chopper manual.  
For instance to filter reads shorter than 1000 bp and longer than 20000 bp, and reads with a minimum Phred average quality score of 10, you would specify: ```--qual_filt --chopper_options '-q 10 -l 1000 --maxlength 20000'```.  **Based on our benchmarking, we recommend using the following parameters ```--chopper_options '-q 8 -l 100'``` as a first pass**.  

  If you are analysing samples that are of poor quality (i.e. failed the QC_FLAG) or amplifying a very short amplicon (e.g. <150 bp), then we recommend using the following setting ```--chopper_options '-q 8 -l 25'``` to retain reads of all lengths.  

A zipped copy of the resulting **preprocessed** and/or **quality filtered fastq file** will be saved in the preprocessing folder.  

If you trim raw read of adapters and/or quality filter the raw reads, an additional quality control step will be performed and a qc report will be generated summarising the read counts recovered before and after preprocessing for all samples listed in the index.csv file.

A qc report will be generated in text and html formats summarising the read counts recovered after the pre-processing step.  
It will include 3 flags:  
1) For the raw_reads_flag, if raw_reads < 5000, the column will display: "Less than 5000 raw reads".  
2) For the qfiltered_flag, if quality_filtered_reads < 1000 , the column will display: "Less than 1000 processed reads".  
3) QC_FLAG:
- GREEN = > 5000 starting reads, > 1000 quality filtered reads.
- ORANGE = < 5000 starting reads, > 1000 quality filtered reads.
- RED = < 5000 starting reads, < 1000 quality filtered reads.

If the user wants to check the data after preprocessing before performing downstream analysis, they can apply the parameter ``--preprocessing_only``.

### Clustering step (RATTLE)

In the clustering mode, the tool [`RATTLE`](https://github.com/comprna/RATTLE#Description-of-clustering-parameters) will be run. 

- The ont_amplicon pipeline will automatically set a **lower read length** of **100** bp during the RATTLE clustering step if the amplicon target_size specified in the csv file is **<=300 bp**.  
- If the amplicon target_size specified in the csv file is **>300 bp**, the default lower read length of **150 bp** will be applied at the RATTLE clustering step instead.  
- For poor quality samples (i.e. failed the QC_FLAG) or if your amplicon is known to be shorter than 150 bp, use the parameter ```--rattle_raw``` to use all the reads without any length filtering during the RATTLE clustering step.  
- Finally, the ``rattle_clustering_max_variance`` is set by default to 10000. It is recommended to drop it to 10 if analysing fastq files that were generated using a **fast** basecalling model.  

  **Special usage:**
  The parameters ``--rattle_clustering_min_length [number]``` (by default: 150) and ```--rattle_clustering_max_length [number]``` (by default: 100,000) can also be specified on the command line to restrict more strictly read size.  
  Additional parameters (other than raw, lower-length, upper-length and max-variance) can be set using the parameter ```--rattle_clustering_options '[additional paramater]'```.  

Example in which all reads will be retained during the clustering step:  
```
nextflow run maelyg/ont_amplicon -resume -profile singularity \
                            --analysis_mode clustering \
                            --adapter_trimming \
                            --qual_filt \
                            --chopper_options '-q 8 -l 25' \
                            --rattle_raw \
                            --blast_threads 2 \
                            --blastn_db /path/to/ncbi_blast_db/nt
```

Example in which reads are first quality filtered using the tool chopper (only reads with a Phread average quality score above 10 are retained). Then for the clustering step, only reads ranging between 500 and 2000 bp will be retained:  
```
nextflow run maelyg/ont_amplicon -resume -profile singularity \
                            --qual_filt \
                            --chopper_options chopper_options = '-q 8 -l 100' \
                            --analysis_mode clustering \
                            --rattle_clustering_min_length 200 \
                            --rattle_clustering_max_length 2000 \
                            --blast_threads 2 \
                            --blastn_db /path/to/ncbi_blast_db/nt
```

### Polishing step (optional)
The clusters derived using RATTLE can be polished. The reads are first mapped back to the clusters using Minimap2 and then the clusters are polished using Racon, Medaka2 and Samtools consensus. 
This step is performed by default by the pipeline but can be skipped by specifying the paramater ``--polishing false``.  

### Primer search
If the fwd_primer and the rev_primer have been provided in the csv file, clusters are then searched for primers using Cutadapt.  

### Blast homology search against NCBI
If the gene targetted is Cytochrome oxidase I (COI), a preliminary megablast homology search against a COI database will be performed; then based on the strandedness of the consensus in the blast results, some will be reverse complemented where required.  

Blast homology search of the consensuses against NCBI is then performed and the top 10 hits are returned.
A separate blast output is then derived using pytaxonkit, which outputs preliminary taxonomic assignment to the top blast hit for each consensus. The nucleotide sequence of qseq (ie consensus match) and sseq (ie reference match) are extracted to use when mapping reads back to consensus and reference respectively (see steps below).  

### Mapping back to consensus
The quality filtered reads derived during the pre-processing step are mapped back to the consensus matches using Mimimap2. Samtools and Mosdepth are then used to derive bam files and coverage statistics. A summary of the blast results, preliminary taxonimic assignment, coverage statistics and associated flags are then derived for each consensus using python.  
Currently applied flags include:  
1) 30X DEPTH FLAG:
  - GREEN = when mapping back to consensus match (ie qseq), the percentage of bases that attained at least 30X sequence coverage > 90
  - ORANGE = when mapping back to consensus match (ie qseq), the percentage of bases that attained at least 30X sequence coverage is between 75 and 90
  - RED = when mapping back to consensus match (ie qseq), the percentage of bases that attained at least 30X sequence coverage is < 75
  - GREY = consensus returned no blast hits
 
2) TARGET ORGANISM FLAG
  - GREEN = target organism detected and % blast identity > 90%
  - ORANGE = target organism detected and % blast identity is < 90%
  - RED = target organism not detected
  - GREY = consensus returned no blast hits
 
3) TARGET SIZE FLAG
  - GREEN = the consensus match length is within ±20% of the target_size
  - ORANGE = the consensus match length is ±20% to ±40% of target size
  - RED = the consensus match length is outside the range of ±40% of the target_size.
 
4) MAPPED READ COUNT FLAG
  - GREEN = when mapping back to the consensus match (ie qseq), read count is >=1000
  - ORANGE = when mapping back to the  consensus match (ie qseq), read count is between 200 and 1000
  - RED = when mapping back to the consensus match (ie qseq), read count is <200
  - GREY = consensus returned no blast hits
 
 5) MEAN COVERAGE FLAG
  - GREEN = when mapping back to the consensus match (ie qseq), the mean coverage >=500
  - ORANGE = when mapping back to the consensus match (ie qseq), the mean coverage is between 100 and 500
  - RED = when mapping back to the consensus match (ie qseq), the mean coverage is < 100
  
### Mapping back to reference (optional)
By default the quality filtered reads derived during the pre-processing step are also mapped back to the
reference blast match and samtools consensus is used to derive independent guided-reference consensuses. Their nucleotide sequences can be compared to that of the original consensuses to resolve ambiguities (ie low complexity and repetitive regions).  

### HTML report
An html summary report is generated for each sample, incorporating sample metadata, QC before and after 
preprocessing, blast results and coverage statistics. It also provides a link to the bam files generated when ampping back to consensus.  
The current proposed structure of the report can be found at: https://miro.com/app/board/uXjVLghknb4=/.  

## Output files

### Preprocessing and host read filtering outputs
If a merge step is required, fastcat will create a summary text file showing the read length distribution.  
Quality check will be performed on the raw fastq file using [NanoPlot](https://github.com/wdecoster/NanoPlot) which is a tool that can be used to produce general quality metrics e.g. quality score distribution, read lengths and other general stats. A NanoPlot-report.html file will be saved under the **SampleName/qc/nanoplot** folder with the prefix **raw**. This report displays 6 plots as well as a table of summary statistics.  

<p align="center"><img src="docs/images/Example_Statistics.png" width="1000"></p>

Example of output plots:
<p align="center"><img src="docs/images/Example_raw_WeightedHistogramReadlength.png" width="750"></p>
<p align="center"><img src="docs/images/Example_LengthvsQualityScatterPlot.png" width="750"></p>

A preprocessed fastq file will be saved in the **SampleName/preprocessing** output directory which will minimally have its read names trimmed after the first whitespace, for compatiblity purposes with all downstream tools. This fastq file might be additionally trimmed of adapters and/or filtered based on quality and length (if PoreChopABI and/or Chopper were run).  

After quality/length trimming, a PoreChopABI log will be saved under the **SampleName/preprocessing/porechop** folder.  

After adapter trimming, a Chopper log file will be saved under the **SampleName/preprocessing/chopper** folder.  

If adapter trimming and/or quality/length trimming is performed, a second quality check will be performed on the processsed fastq file and a NanoPlot-report.html file will be saved under the **SampleName/qc/nanoplot** folder with the prefix **filtered**.  

If the adapter trimming and/or the quality filtering options have been run, a QC report will be saved both in text and html format (i.e. **run_qc_report_YYYYMMDD-HHMMSS.txt** and **run_qc_report_YYYYMMDD-HHMMSS.html**) under the **qc_report** folder.  

Example of report:

| Sample| raw_reads | quality_filtered_reads | percent_quality_filtered | raw_reads_flag | qfiltered_flag | QC_FLAG |
| --- | --- | --- | --- | --- | --- | --- |
| ONT141 | 10929 | 2338 | 21.39 | | | GREEN |
| ONT142| 21849 | 4232 | 9.37 | | | GREEN |

### Clustering step outputs  
In this mode, the output from Rattle will be saved under **SampleName/clustering/rattle/SampleName_rattle.fasta**. The number of reads contributing to each clusters is listed in the header. The amplicon of interest is usually amongst the most abundant clusters (i.e. the ones represented by the most reads).  

### Polishing step outputs 
(in progress)  

### Blast search outputs  
All the top hits derived for each contig are listed in the file **SampleName/megablast/SampleName_final_polished_consensus_megablast_top_10_hits.txt**. This file contains the following 26 columns:
```
- qseqid
- sgi
- sacc
- length
- pident
- mismatch
- gaps
- gapopen
- qstart
- qend
- qlen
- sstart
- send
- slen
- sstrand
- evalue
- bitscore
- qcovhsp
- stitle
- staxids
- qseq
- sseq
- sseqid
- qcovs
- qframe
- sframe
```

A separate blast output called **SampleName/megablast/SampleName_final_polished_consensus_megablast_top_hit.txt** is then derived using pytaxonkit, which outputs preliminary taxonomic assignment to the top blast hit for each consensus. The nucleotide sequence of qseq (ie consensus match) and sseq (ie reference match) are extracted to use when mapping reads back to consensus and reference respectively (see steps below). These are called **SampleName/megablast/SampleName_final_polished_consensus_match.fasta** and **SampleName/megablast/SampleName_reference_match.fasta** respectively.  

### Outputs from the mapping reads back to consensus matches step
(in progress)

### Outputs from mapping reads back to reference matches step
By default the quality filtered reads derived during the pre-processing step are mapped back to the
reference blast match. A bam file is generated and Samtools consensus is used to derive independent guided-reference consensuses that are stored in a file called **SampleName/mapping_back_to_ref/samtools_consensus_from_ref.fasta** file. Their nucleotide sequences can be compared to that of the original consensuses to resolve ambiguities (ie low complexity and repetitive regions). 

### HTML report output
(in progress)  

## Authors
Marie-Emilie Gauthier gauthiem@qut.edu.au
Cameron Hyde c.hyde@qcif.edu.au

## To do :
Provide a quick start up  
Provide a config file example  
Add an image depicting the current flow of the pipeline  
Finish output section  

Improve reporting errors when RATTLE fails to produce clusters  
Add additional flags (basecalling model, contamination flag, % long reads)  
Force specification of COI database if COI gene specified  

Fix bug in reporting of contigs returning a blast hit vs total in html report.  
Incorporate basecalling model, analyst and facility in html report  
Incorporate cluster match fasta file in html report.  
Display colour for each flag  


Prevent pipeline from proceeding if fast basecalling model is detected?  
Generate a QC report even if preprocessing is not run to capture the raw read counts?  
Provide option to run only map to ref  
List current version of tools  
