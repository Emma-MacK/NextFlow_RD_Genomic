process baseRecalibrator {

    label 'process_low'
    container 'broadinstitute/gatk:4.1.4.0'

    tag "$bamFile"

    // Publish BQSR BAM files to the specified directory
    publishDir("$params.outdir/BAM", mode: "copy")

    input:
    tuple val(sample_id), file(bamFile), file(baiFile)
    path knownSites
    path indexFiles

    output:
    tuple val(sample_id), file("${bamFile.baseName}_recalibrated.bam"), file("${bamFile.baseName}_recalibrated.bai")

    script:
    """
    echo "Running BQSR"

    if [[ -n params.genome_file ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"

    referenceVCF=\$(find -L -name '*.vcf.gz')

    # get fasta dict, Homo_sapiens_assembly38.fasta.dict, rename to Homo_sapiens_assembly38.dict

    mv Homo_sapiens_assembly38.fasta.dict Homo_sapiens_assembly38.dict

    # Generate recalibration table for the input BAM file
    gatk --java-options "-Xmx8G" BaseRecalibrator \
        -R "\${genomeFasta}" \
        -I ${bamFile} \
        --known-sites "\${referenceVCF}"\
        -O ${bamFile.baseName}.recal_data.table

    # Apply BQSR to the input BAM file
    gatk --java-options "-Xmx8G" ApplyBQSR \
        -R "\${genomeFasta}" \
        -I ${bamFile} \
        --bqsr-recal-file ${bamFile.baseName}.recal_data.table \
        -O ${bamFile.baseName}_recalibrated.bam

    # Index the recalibrated BAM file
    samtools index ${bamFile.baseName}_recalibrated.bam ${bamFile.baseName}_recalibrated.bai

    echo "BQSR Complete"
    """
}