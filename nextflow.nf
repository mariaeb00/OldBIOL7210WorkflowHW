#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.reads  = "$projectDir/data/*_{1,2}.fastq.gz"
params.outdir = "$projectDir/results"

// Process: Trimming reads using fastp
process fastp {
    tag "${pair_id}"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    // Expect a tuple consisting of a sample identifier and two files (R1 and R2)
    tuple val(pair_id), path(read1), path(read2)

    output:
    // Emit a tuple of the sample ID with the two trimmed files
    tuple val(pair_id), path("trimmed_${pair_id}_R1.fastq.gz"), path("trimmed_${pair_id}_R2.fastq.gz")

    script:
    """
    fastp --in1 ${read1} --in2 ${read2} \\
          --out1 trimmed_${pair_id}_R1.fastq.gz \\
          --out2 trimmed_${pair_id}_R2.fastq.gz \\
          --qualified_quality_phred 30 --length_required 20
    """
}

// Process: Assembly with SKESA
process skesa {
    tag "${pair_id}"
    publishDir "${params.outdir}/skesa", mode: 'copy'

    input:
    // Accept a tuple of sample ID and the two trimmed read files
    tuple val(pair_id), path(trimmed1), path(trimmed2)

    output:
    // Emit a tuple with the sample ID and the assembly FASTA file
    tuple val(pair_id), path("${pair_id}_assembly.fasta")

    script:
    """
    skesa --fastq ${trimmed1},${trimmed2} \\
          --contigs_out ${pair_id}_assembly.fasta
    """
}

// Process: Quality assessment with QUAST
process quast {
    tag "${pair_id}"
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    // Accept a tuple with sample ID and the assembly file
    tuple val(pair_id), path(assembly)

    output:
    // Emit a tuple with the sample ID and the QUAST report directory
    tuple val(pair_id), path("${pair_id}_quast_report")

    script:
    """
    mkdir -p ${pair_id}_quast_report
    quast -o ${pair_id}_quast_report ${assembly}
    """
}

// Process: Multilocus sequence typing with MLST
process mlst {
    tag "${pair_id}"
    publishDir "${params.outdir}/mlst", mode: 'copy'

    input:
    // Accept a tuple with sample ID and the assembly file
    tuple val(pair_id), path(assembly)

    output:
    // Emit a tuple with sample ID and the MLST result file
    tuple val(pair_id), path("${pair_id}_mlst.txt")

    script:
    """
    mlst ${assembly} > ${pair_id}_mlst.txt
    """
}

workflow {
    // Create input channel from paired-end reads. This uses Nextflow's built-in
    // fromFilePairs function which groups files by common prefix.
    read_pairs = Channel.fromFilePairs(params.reads, size: 2)
                       .map { pair_id, files -> tuple(pair_id, files[0], files[1]) }

    // Sequential execution: fastp -> skesa
    trimmed = fastp(read_pairs)
    assembly = skesa(trimmed)

    // Branching: run QUAST and MLST in parallel on the assembly output
    quast(assembly)
    mlst(assembly)
}
