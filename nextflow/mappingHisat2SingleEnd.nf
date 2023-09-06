#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
* params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
* params.transcriptome = "$projectDir/data/ggal/transcriptome.fa"
*/
baseDir="/home/akram/gitrepos/training/nf-training"

params.reads = "$baseDir/data/ggal/*_1.fq"
/*params.genome = "/home/akram/genomes/human/GRCh38.primary_assembly.genome.fa"*/
params.genome = "$baseDir/data/chr20_GRCh38.fa"
params.index_basename= "/home/akram/gitrepos/training/myNextflowScripts/results_hisat2/hisat2_idx/chr20_GRCh38"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results_hisat2"


log.info """\
         ----------------------------------
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         projectDir: ${projectDir}
         baseDir: ${baseDir}
         genome: ${params.genome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent(true)


process FASTQC {
    tag "FASTQC"
    publishDir "${params.outdir}/quality_control", mode: 'copy'

    input:
    path reads

    output:
    path "fastqc_${reads.baseName}_logs"

    script:
    """
    mkdir fastqc_${reads.baseName}_logs
    fastqc -o fastqc_${reads.baseName}_logs -f fastq -q ${reads}

    """
}


process INDEX {

    tag "Building Hisat2 index on $genome"
    publishDir "${params.outdir}/hisat2_idx", pattern: "*.ht2"

    input:
    path genome
    
    output:
    path "*"

    script:
    """
    
    hisat2-build ${genome} ${genome.baseName} -p $task.cpus

    """
}

process MAPPING {

    tag "Mapping reads to index genome using HISAT2"
    publishDir "${params.outdir}/aln", mode: 'copy'

    input:
    path reads

    output:
    path "${reads.baseName}.hisat2.sorted.bam", emit: bam
    path "hisat2_summary_${reads.baseName}.txt", emit: summary_hisat2
    
    script:
    """
    hisat2 -p $task.cpus -x ${params.index_basename} -U ${reads} --summary-file hisat2_summary_${reads.baseName}.txt | samtools view -@ $task.cpus -bS - | samtools sort -o ${reads.baseName}.hisat2.sorted.bam -@ $task.cpus

    """
}

process MULTIQC {
    tag "MULTIQC on FASTQC reports"
    publishDir "${params.outdir}/quality_control", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

workflow {

    Channel
    .fromPath(params.reads, checkIfExists: true )
    .set {reads_ch}

    reads_ch.view()

    fastqc_ch = FASTQC(reads_ch)
    
    index_ch = INDEX(params.genome)
    /*index_ch.view()*/
    
    mappingHisat2_ch = MAPPING(reads_ch)
    
    MULTIQC(mappingHisat2_ch.summary_hisat2.mix(fastqc_ch).collect())

}


workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n"
        """
        .stripIndent()
       
    log.info (workflow.success ? msg : "Oops .. something went wrong")
}
