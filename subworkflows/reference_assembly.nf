process map_reads{
    /*
    Map reads to reference using minimap2.
    Filter reads by mapping quality.
    Filter internally-primed reads.
    */
    label "nf_bambu"
    cpus params.threads

    input:
       tuple val(sample_id), path (fastq_reads), path(index), path(reference)

    output:
       tuple val(sample_id), path("${sample_id}_reads_aln_sorted.bam"), emit: bam
       tuple val(sample_id), path("${sample_id}_read_aln_stats.tsv"), emit: stats
    script:
        def ContextFilter = """AlnContext: { Ref: "${reference}", LeftShift: -${params.poly_context},
        RightShift: ${params.poly_context}, RegexEnd: "[Aa]{${params.max_poly_run},}",
        Stranded: True,Invert: True, Tsv: "internal_priming_fail.tsv"} """
    """
    minimap2 -t ${params.threads} -ax splice ${params.minimap2_opts} ${index} ${fastq_reads}\
        | samtools view -q ${params.minimum_mapping_quality} -F 2304 -Sb -\
        | seqkit bam -j ${params.threads} -x -T '${ContextFilter}' -\
        | samtools sort -@ ${params.threads} -o "${sample_id}_reads_aln_sorted.bam" - ;
    ((cat "${sample_id}_reads_aln_sorted.bam" | seqkit bam -s -j ${params.threads} - 2>&1)  | tee ${sample_id}_read_aln_stats.tsv ) || true

    if [[ -s "internal_priming_fail.tsv" ]];
        then
            tail -n +2 "internal_priming_fail.tsv" | awk '{print ">" \$1 "\\n" \$4 }' - > "context_internal_priming_fail_start.fasta"
            tail -n +2 "internal_priming_fail.tsv" | awk '{print ">" \$1 "\\n" \$6 }' - > "context_internal_priming_fail_end.fasta"
    fi
    """
}

process output2 {
    debug true
    
    label "nf_bambu"

    // publish bam to a bam directory located in the output directory
    bam_dir = params.out_dir + "bam/"

    publishDir (
        bam_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}

workflow reference_assembly {
    take:
       index
       reference
       fastq_reads
    main:
        map_sample_ids_cls2 = {it ->
        /* Harmonize tuples
        output:
            tuple val(sample_id), path('*.gff')
        When there are multiple paths, will emit:
            [sample_id, [path, path ..]]
        when there's a single path, this:
            [sample_id, path]
        This closure makes both cases:
            [[sample_id, path][sample_id, path]].
        */
            if (it[1].getClass() != java.util.ArrayList){
                // If only one path, `it` will be [sample_id, path]
                return [it]
            }
            l = [];
            for (x in it[1]){
                l.add(tuple(it[0], x))
            }
            return l
        }


        bam_results = Channel.empty()

        map_reads(fastq_reads.combine(index).combine(reference))

        // map_reads.out.bam.view()
        // map_reads.out.bam.flatMap(map_sample_ids_cls2).view()
        // map_reads.out.bam.flatMap(map_sample_ids_cls2).map {it -> it[1]}.view()

        bam_results = map_reads.out.bam.flatMap(map_sample_ids_cls2).map {it -> [it[1], null]}.concat(bam_results)
        output2(bam_results)

    emit:
       bam = map_reads.out.bam
       stats = map_reads.out.stats
}
