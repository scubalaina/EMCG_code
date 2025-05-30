#!/user/bin/nextflow
nextflow.enable.dsl=2
 
 
/// INPUTS
DIR_input = "test_proteins/"
SUFFIX_of_input_files = ".faa"
 
// OUTPUT
DIR_output = "test_nog_out/"
 
///  HMM to search against
PATH_hmm = "nog.hmm"
 
 
 
/// DEV MODE STUFF
params.dev = false
params.num_inputs = 2
 
 
workflow {
	CH_sample_id_AND_faa = channel
		.fromPath("$DIR_input/*$SUFFIX_of_input_files", checkIfExists:true) // --> path(faa)
		.map { it -> tuple(it.getBaseName(), it) } // --> tuple( val(sample_id), path(faa)) 
	HMMER_v3_3_2( CH_sample_id_AND_faa.take( params.dev ? params.num_inputs : -1) ) //.take(params.dev ? params.num_inputs : -1)
}
 
process HMMER_v3_3_2 {
	container='docker://quay.io/biocontainers/hmmer:3.3.2--h87f3376_2'
	publishDir DIR_output, mode: "copy"
	cpu=6
	tag "${sample_id}"
 
	input: tuple val(sample_id), path(faa)
 
	output: tuple val(sample_id), path("${sample_id}_hmm_tbl.txt")
 
	"hmmsearch -E 0.00001 --tblout ${sample_id}_hmm_tbl.txt ${PATH_hmm} ${faa} "
}