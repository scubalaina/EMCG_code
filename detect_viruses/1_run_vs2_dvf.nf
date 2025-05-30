#!/user/bin/nextflow
nextflow.enable.dsl=2


/// INPUTS
DIR_input = "test_seq/"
SUFFIX_of_input_files = "_1kb_all.fasta"

// OUTPUTS
VS2_output = "test_seq_vs2_out/"
DFV_output = "test_seq_dvf_out/"

/// DEV MODE STUFF
params.dev = false // Lets user testrun nextflow command (by adding flag '--dev') which will have this pipeline run on JUST ONE SAG (again, as a test)
params.num_inputs = 2




workflow {
  CH_sample_id_AND_fasta = channel
    .fromPath("$DIR_input/*$SUFFIX_of_input_files", checkIfExists:true)
    .map { it -> tuple( it.getBaseName().split("_1kb_all")[0], it) } 
  VIRSORTER2(CH_sample_id_AND_fasta.take(params.dev ? params.num_inputs : -1))
  DVF(CH_sample_id_AND_fasta.take( params.dev ? params.num_inputs : -1))
}



process VIRSORTER2 {
  conda='/mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/vs2/'
  errorStrategy = 'ignore' 
  //errorStrategy = 'retry' ; maxRetries = 4; memory={100.GB * task.attempt}
  tag "${sample_id}"
  publishDir VS2_output, mode: 'copy'
  cpu=4
  input: tuple val(sample_id), path(fasta)
  output: tuple val(sample_id), path("${sample_id}_vs2_out")
  """
  virsorter run -w ${sample_id}_vs2_out -l ${sample_id}_vs2 -i ${fasta} --min-length 1000 -j 4 all --rm-tmpdir
  """
}

process DVF {
  conda='/mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/dvf/'
  errorStrategy = 'ignore' 
  //errorStrategy = 'retry' ; maxRetries = 4; memory={100.GB * task.attempt}
  tag "${sample_id}"
  publishDir DFV_output, mode: 'copy'
  cpu=4
  input: tuple val(sample_id), path(fasta)
  output: tuple val(sample_id), path("${sample_id}_dvf_out")
  """
  python /home/aweinheimer/Tools/DeepVirFinder/dvf.py -i ${fasta} -o ${sample_id}_dvf_out -l 1000 -c 4 -m /home/aweinheimer/Tools/DeepVirFinder/models/
  """
}
