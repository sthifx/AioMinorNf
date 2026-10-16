#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process hisat_index {
  publishDir "${params.out_dir}/reference/"
  tag "${fasta}"
  cpus 64
  conda "${params.envs_dir}/preprocess.env"

  input:
    path fasta
  output:
    path 'genome.*'

  shell:
  '''
  hisat2-build \
    -p !{task.cpus} \
    !{fasta} \
    genome
  '''
}

process hisat_and_remove_unmapped {
  publishDir "${params.out_dir}/host_removed_fastq/"
  tag "${sid}"
  conda "${params.envs_dir}/preprocess.env"
  cpus 8

  input:
    tuple val(sid), path(r1), path(r2)
    path (genome)
  output:
    tuple val(sid), path("${sid}.mapping_metrics")
    tuple val(sid), path("${sid}.unmapped_R1.fastq.gz"), path("${sid}.unmapped_R2.fastq.gz")

  shell:
  '''
  hisat2 \
    -x genome \
    -1 !{r1} \
    -2 !{r2} \
    --fr \
    -p !{task.cpus} \
    2> !{sid}.mapping_metrics | \
    samtools view -bSh -f 1 -f 12 -F 256 -F 2048 | \
    samtools sort -n | \
    bedtools bamtofastq -i - \
          -fq  !{sid}.unmapped_R1.fastq \
	  -fq2 !{sid}.unmapped_R2.fastq

    gzip -f !{sid}.unmapped_R1.fastq !{sid}.unmapped_R2.fastq
  '''
  }  

process run_aiominor{
  publishDir "${params.out_dir}/aio", pattern: "*.aio"
  tag "${sid}"
  conda "${params.envs_dir}/aiominor.env"
  cpus 16

  input:
    tuple val(sid), path(r1), path(r2)
  output:
    tuple val(sid), path("${sid}.aio/")
    tuple val(sid), path("${sid}.aio/alignment/alignment.support_oem/consensus.txt")

  shell:
  '''
  perl !{params.aiominor_dir}/AioMinor.pl \
    -t !{task.cpus} \
    -platform illumina \
    -method amplicon \
    -maxins 500 \
    -ref !{params.aiominor_dir}/References/genome_NC_045512.2.fasta \
    -codingRegion !{params.aiominor_dir}/References/CodingRegion_NC_045512.2.txt \
    -fq1 !{r1} \
    -fq2 !{r2} \
    -o !{sid}.aio \
    -primerbed !{params.primer_bed} \
    -rotation 2
  '''
}

process run_aiominor_vs_day1_consensus{
  publishDir "${params.out_dir}/aio_vs_d1/${pt_id}/"
  tag "${sample_id}"
  conda "${params.envs_dir}/aiominor.env"
  cpus 16

  input:
    tuple val(sample_id), val(pt_id), val(study_day), path(r1), path(r2), path(consensus_fa)
  output:
    tuple val(sample_id), path("${sample_id}.day_${study_day}/")

  shell:
  '''
  perl !{params.aiominor_dir}/AioMinor_vs_specified_consensus.pl \
    -t !{task.cpus} \
    -platform illumina \
    -method amplicon \
    -maxins 500 \
    -ref !{params.aiominor_dir}/References/genome_NC_045512.2.fasta \
    -codingRegion !{params.aiominor_dir}/References/CodingRegion_NC_045512.2.txt \
    -fq1 !{r1} \
    -fq2 !{r2} \
    -o !{sample_id}.day_!{study_day} \
    -primerbed !{params.primer_bed} \
    -rotation 2 \
    -variantRef !{consensus_fa}
  '''
}


workflow {
  main:
    // Check params
    def usage = """
    Required:
      --fastq_dir        Directory containing quality- and adapter-trimmed FASTQ files.
                         Each sample must have its own subdirectory named Sample_X-<sampleid>,
                         where X is one or more digits and <sampleid> matches the metadata CSV.
                         Each subdirectory must contain a paired set of FASTQ files matching
                         *_R{1,2}.fastq.gz.
      --host_fasta      Path to host (human) reference FASTA file
      --aiominor_dir    Full path to the AioMinor downloaded github directory 
                        (git clone https://github.com/Hiscox-lab/AioMinor.git). 
                        Must contain AioMinor_vs_specified_consensus.pl alongside AioMinor.pl
      --primer_bed      BED formatted primers file 
                        (e.g. 'NC_045512.2     30      54      nCoV-2019_1_LEFT        2       +       .')
      --envs_dir        Directory containing both the aiominor.env and preprocess.env conda environment directories
      --out_dir         Output directory
    
    Optional:
      --mdata_csv       Metadata CSV. Must contain columns called 
                          * 'sampleid'  : the ID of the sample, which must match the naming of the FASTQ files subdirectories (see --fastq_dir)
                          * 'pt_id'     : the anonymised ID of the patient to which the sample was derived from
                          * 'study_day' : the day on which the sample was taken (e.g. 1,2,3)
                        These metadata are used to identify day 1 samples for patients. When running the aiominor vs_day1 part of the analysis, 
                        it will call variants relative to the day1 sample for each patient. As such:
                          * There must only be one study_day 1 per pt_id
                          * There should be at least one study_day 1 per pt_id. If there isn't one, no variants relative to day 1 will be called
                        If not provided, this part of the analysis will be skipped
    """

    // Check variables
    if( !params.fastq_dir || !params.host_fasta || !params.aiominor_dir || !params.primer_bed || !params.envs_dir || !params.out_dir )      error "Missing parameter\n${usage}"

    // Gather FASTQs
    fastqs = Channel.fromFilePairs("${params.fastq_dir}/*/*_R{1,2}.fastq.gz", flat : true).ifEmpty { exit 1, "Read pairs could not be found:"}.map{x -> tuple(x[1].parent.name, x[1],x[2])}      

    // Index host and remove host reads
    ref_idx = hisat_index(params.host_fasta)
    (bam_metrics,filt_fastqs) = hisat_and_remove_unmapped(fastqs, ref_idx)

    ////
    // Run AioMinor against each samples own consensus
    (aio_out, consensus) = run_aiominor(filt_fastqs)

    ////
    // Run AioMinor against each samples day 1, if metadata is provided

    if (params.mdata_csv){
      // Load metadata
      mdata = Channel.fromPath(params.mdata_csv)
                     .splitCsv(header : true)
                     .map{row -> tuple(row.sampleid, row.pt_id, row.study_day)}

     // Merge consensus sequences with metadata
     d1_consensus = consensus.map{sid, consensus_fa -> tuple(sid.split("-")[1], consensus_fa)}
              .combine(mdata, by: 0)
              .filter{sampleid, consensus, pt_id, study_day -> study_day == "1"}
              .map{sampleid, consensus, pt_id, study_day -> tuple(pt_id, consensus)}
    
     // Merge filtered fastqs with mdata
     // Get sample ID as 2nd element
     filt_fastqs_mdata = filt_fastqs.map{sid, r1, r2 -> tuple(sid.split("-")[1], sid, r1, r2)}
                 .combine(mdata, by: 0)
                 .map{sampleid, sid, r1, r2, pt_id, study_day -> tuple(pt_id, sampleid, study_day, r1, r2)} 
  
     // Merge day 1 consensus sequences with samples from the same patient
     reads_and_consensus = d1_consensus.combine(filt_fastqs_mdata, by: 0)
                 .map{pt_id, consensus_fa, sample_id, study_day, r1, r2 -> tuple(sample_id, pt_id, study_day, r1, r2, consensus_fa)}

     // Run AioMinor vs Day 1
     run_aiominor_vs_day1_consensus(reads_and_consensus)
  }
}
