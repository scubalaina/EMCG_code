#!/user/bin/nextflow
nextflow.enable.dsl=2

/// INPUTS
DIR_input = "/mnt/stepanauskas_nfs/aweinheimer/microcapsules_rd3/seqs/spc3_untrim_contigs_min1k/"

// OUTPUT
CHECKV_output = "/mnt/stepanauskas_nfs/aweinheimer/microcapsules_rd3/detect_vir/spc3_virdetect/spc3_checkv_out/"


/// DEV MODE STUFF
params.dev = false // Lets user testrun nextflow command (by adding flag '--dev') which will have this pipeline run on JUST ONE SAG (again, as a test)
params.num_inputs = 2
//SUFFIX_of_input_files = ".fasta"
//PREFIX_of_input_files = "final_contigs"



workflow {
  CH_sample_id_AND_fasta = channel
    .fromPath("$DIR_input/*_1kb_all.fasta", checkIfExists:true)
    .map { it -> tuple( it.getBaseName(), it) } 
  CHECKV(CH_sample_id_AND_fasta.take(params.dev ? params.num_inputs : -1))
  //DVF(CH_sample_id_AND_fasta.take( params.dev ? params.num_inputs : -1))
}



process CHECKV {
  conda = '/home/aweinheimer/miniconda3/envs/checkv'
  errorStrategy = 'ignore' 
  //errorStrategy = 'retry' ; maxRetries = 4; memory={100.GB * task.attempt}
  tag "${sample_id}"
  publishDir CHECKV_output, mode: 'copy'
  cpu=4
  input: tuple val(sample_id), path(fasta)
  output: tuple val(sample_id), path("${sample_id}_checkv_out")
  """
  checkv end_to_end ${fasta} ${sample_id}_checkv_out/ -t 4 -d /home/aweinheimer/Tools/checkv-db-v1.5
  """
}

