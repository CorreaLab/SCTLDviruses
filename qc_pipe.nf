#!/usr/bin/env nextflow

/*
*/

def helpMessage() {
    log.info """

    This is a workflow to process raw reads and prepare for downstream RNA seq analyses

    One mandatory argument:
        --RunThat           run the whole thing

    """.stripIndent()
}
// Show help message
if (params.help) {
helpMessage()
exit 0
} else if (params.fullHelp) {
fullHelpMessage()
exit 0
}

if (params.readsTest) {
    println("\n\tRunning vAMPirus with TEST dataset\n")
    Channel
        .fromFilePairs(params.readsTest)
        .ifEmpty{ exit 1, "params.readTest was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch; reads_processing }
} else if (params.single) {
    println("\n\tLocating single read files...\n")
    Channel
        .fromFilePairs("${params.reads}", size: -1, checkIfExists: true)
        .ifEmpty{ exit 1, "params.reads was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch; reads_processing }
} else {
    println("\n\tLocating paired read files...\n")
    Channel
        .fromFilePairs("${params.reads}", checkIfExists: true)
        .ifEmpty{ exit 1, "params.reads was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch}
}

if ( ${params.Clean} )  {

    process QualityCheck_1 {

                label 'low_cpus'

                tag "${sample_id}"

                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/FastQC/PreClean", mode: "copy", overwrite: true

                conda (params.condaActivate ? "-c conda-forge bioconda::fastqc=0.11.9=0" : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" : "quay.io/biocontainers/fastqc:0.11.9--0")

                input:
                    tuple sample_id, file(reads) from reads_qc_ch

                output:
                    tuple sample_id, file("*_fastqc.{zip,html}") into fastqc_results

                script:
                    """
                    fastqc --quiet --threads ${task.cpus} ${reads}
                    """
            }


    process Adapter_Removal {

            label 'norm_cpus'

            tag "${sample_id}"

            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval", mode: "copy", overwrite: true, pattern: "*.filter.fq"
            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval/fastpOut", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

            conda (params.condaActivate ? "bioconda::fastp=0.23.2=hb7a2d85_2" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastp:0.23.2--hb7a2d85_2" : "quay.io/biocontainers/fastp:0.23.2--hb7a2d85_2")

            input:
                tuple sample_id, file(reads) from reads_ch

            output:
                tuple sample_id, file("*.fastp.{json,html}") into fastp_results
                tuple sample_id, file("*.fastp.{json}") into fastp_json
                tuple sample_id, file("*.filter.fq") into readsforqc2

            script:
            """
                echo ${sample_id}
                fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                --average_qual ${params.avQ} -q ${params.trimq} -l ${params.len} -y -Y ${params.comp} -g -x -n ${params.mN} -c --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json --thread ${task.cpus} \
                --report_title ${sample_id}
                """
            }

        process QualityCheck_2 {

                label 'low_cpus'

                tag "${sample_id}"

                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/FastQC/PostClean", mode: "copy", overwrite: true

                conda (params.condaActivate ? "-c conda-forge bioconda::fastqc=0.11.9=0" : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" : "quay.io/biocontainers/fastqc:0.11.9--0")

                input:
                    tuple sample_id, file(reads) from readsforqc2

                output:
                    tuple sample_id, file("*_fastqc.{zip,html}") into fastqc2_results

                script:
                    """
                    fastqc --quiet --threads ${task.cpus} ${reads}
                    """
            }

        }