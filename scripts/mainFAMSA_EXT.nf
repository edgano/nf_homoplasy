#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'homoplasy-analysis'.
 *
 *   homoplasy-analysis is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   homoplasy-analysis is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with homoplasy-analysis.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main homoplasy-analysis pipeline script
 *
 * @authors
 * Edgar Garriga
 * Cedric Notredame 
 */

/*
 * defaults parameter definitions
 */

// input sequences to align in fasta format
params.seqs ="/users/cn/egarriga/datasets/extHomfam/combinedSeqs/*.fa"
//params.seqs ="${baseDir}/test/*.fa"

// input reference sequences aligned in
params.refs = "/users/cn/egarriga/datasets/extHomfam/refs/*"
//params.refs ="${baseDir}/test/*.ref"

// input guide trees in Newick format. Or `false` to generate trees
params.trees ="/users/cn/egarriga/datasets/extHomfam/trees/*.{FAMSA,CLUSTALO,MAFFT_PARTTREE,dpparttree0}.dnd"

// which tree methods to run if `trees` == `false`
params.tree_method = "dpparttreednd0"
//codnd,dpparttreednd1,dpparttreednd2,dpparttreednd2size,fastaparttreednd,fftns1dnd,fftns1dndmem,fftns2dnd,fftns2dndmem,mafftdnd,parttreednd0,parttreednd1,parttreednd2,parttreednd2size

// which alignment methods to run
params.align_method = "FAMSA" //CLUSTALO,MAFFT-FFTNS1,MAFFT-SPARSECORE,UPP,MAFFT-GINSI"

// generate regressive alignments ?
params.regressive_align = true

// create standard alignments ?
params.progressive_align = false

// evaluate alignments ?
params.evaluate = true

params.homoplasy = true

params.metrics = true

// bucket sizes for regressive algorithm
params.buckets= '1000,5000,10000,20000'

// output directory
params.outdir = "$baseDir/resultsFAMSAcompletEXT"

log.info """\
         H O M O P L A S Y   A n a l y s i s  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA)               : ${params.refs}
         Input trees (NEWICK)                           : ${params.trees}
         Alignment methods                              : ${params.align_method}
         Tree methods                                   : ${params.tree_method}
         Generate progressive alignments                : ${params.progressive_align}
         Generate regressive alignments (DPA)           : ${params.regressive_align}
         Bucket Sizes for regressive alignments         : ${params.buckets}
         Perform evaluation? Requires reference         : ${params.evaluate}
         Capture Homoplasy metrics?                     : ${params.homoplasy}
         Capture runtime metrics?                       : ${params.metrics}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// Channels containing sequences
if ( params.seqs ) {
  Channel
  .fromPath(params.seqs)
  .map { item -> [ item.baseName, item] }
  .into { seqsCh; seqs2 }
}

if ( params.refs ) {
  Channel
  .fromPath(params.refs)
  .map { item -> [ item.baseName, item] }
  .set { refs }
}

// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  Channel
    .fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
    .set { trees }
}
else { 
  Channel
    .empty()
    .set { trees }
}

tree_methods = params.tree_method
align_methods = params.align_method

/*
 * GENERATE GUIDE TREES USING MEHTODS DEFINED WITH "--tree_method"
 *
 * NOTE: THIS IS ONLY IF GUIDE TREES ARE NOT PROVIDED BY THE USER
 * BY USING THE `--trees` PARAMETER
 */

process generate_trees {
    tag "${id}.${tree_method}"
    publishDir "${params.outdir}/guide_trees", mode: 'copy', overwrite: true
   
    input:
    set val(id), \
         file(seqs) \
         from seqsCh
    each tree_method from tree_methods.tokenize(',') 

   output:
     set val(id), \
       val(tree_method), \
       file("${id}.${tree_method}.dnd") \
       into treesGenerated

   when:
     !params.trees

   script:
    template "tree/generate_tree_${tree_method}.sh"

}

treesGenerated
  .mix ( trees )
  .combine ( seqs2, by:0 )
  .into {seqsAndTreesForRegressiveAlignment; seqsAndTreesForProgressiveAlignment }

process regressive_alignment {
    tag "${id}-${align_method}-${bucket_size}-${tree_method}"
    publishDir "${params.outdir}/alignments", pattern: '*.aln', mode: 'copy', overwrite: true

    input:
        set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForRegressiveAlignment

      each bucket_size from params.buckets.tokenize(',')
      each align_method from align_methods.tokenize(',')   

    when:
      params.regressive_align

    output:
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("reg_align"), \
        val(bucket_size), \
        file("${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.aln") \
        into regressiveOut

      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val(bucket_size), \
        file("${id}.homoplasy") \
        into homoReg

      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val(bucket_size), \
        val("reg_align"), \
        file(".command.trace") \
        into metricsReg

    script:
       template "reg_align/reg_align_${align_method}.sh"
}

process progressive_alignment {
    tag "${id}-${align_method}-${bucket_size}-${tree_method}"
    publishDir "${params.outdir}/alignments", mode: 'copy', overwrite: true

    input:
        set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForProgressiveAlignment
    each align_method from align_methods.tokenize(',') 

    when:
      params.progressive_align

    output:
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("prog_align"), \
        val("NA"), \
        file("${id}.prog_align.NA.${align_method}.with.${tree_method}.tree.aln") \
        into progressiveOut

      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("NA"), \
        val("prog_align"), \
        file(".command.trace") \
        into metricsProg

    script:
      template "prog_align/prog_align_${align_method}.sh"
}

metricsReg
  .mix ( metricsProg )
  .set { all_metrics }

process metrics{
    tag "${id}"
    publishDir "${params.outdir}/${workflow.runName}", mode: 'copy', overwrite: true

    input:
    set val(id), \
      val(align_method), \
      val(tree_method), \
      val(bucket_size), \
      val(mode), \
      val(metricsFile) \
      from all_metrics

    when:
      params.metrics

    output:
    set file("${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.metrics"), \
      file("*.realtime"), \
      file("*.rss"), \
      file("*.peakRss"), \
      file("*.vmem"), \
      file("*.peakVmem"), \
      file("*.metrics") \
        into metricsOut

    script:
    """    
    ## realtime > Task execution time i.e. delta between completion and start timestamp i.e. compute wall-time
    awk -F = '{ if (\$1=="") print \$2}' ${metricsFile} > ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.realtime

    ## rss > Real memory (resident set) size of the process
    awk -F = '{ if (\$1=="rss") print \$2}' ${metricsFile}> ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.rss

    ## peakRss > Peak of real memory
    awk -F = '{ if (\$1=="peakRss") print \$2}' ${metricsFile} > ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.peakRss

    ## vmem > Virtual memory size of the process
    awk -F = '{ if (\$1=="vmem") print \$2}' ${metricsFile} > ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.vmem

    ## peakVmem > Peak of virtual memory
    awk -F = '{ if (\$1=="peakVmem") print \$2}' ${metricsFile} > ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.peakVmem

    mv ${metricsFile} ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.metrics
    """
}

process homoplasy{
    tag "${id}"
    publishDir "${params.outdir}/homoplasy", mode: 'copy', overwrite: true

    input:
    set val(id), \
      val(align_method), \
      val(tree_method), \
      val(bucket_size), \
      file(homoplasy) \
      from homoReg

    when:
      params.homoplasy

    output:
    set file("*.homo"), \
        file("*.w_homo"), \
        file("*.w_homo2"), \
        file("*.len"), \
        file("*.ngap"), \
        file("*.ngap2") \
        into homoplasyOut

    script:
    """    
    ## homo
    awk -F : '{ if (\$1=="HOMOPLASY") print \$2}' ${homoplasy} > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.homo
    ## w_homo
    awk -F : '{ if (\$1=="WEIGHTED_HOMOPLASY") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo
    ## w_homo2
    awk -F : '{ if (\$1=="WEIGHTED_HOMOPLASY2") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo2
    ## len
    awk -F : '{ if (\$1=="LEN") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.len
    ## ngap
    awk -F : '{ if (\$1=="NGAP") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.ngap
    ## ngap2
    awk -F : '{ if (\$1=="NGAP2") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.ngap2 
    """
}

progressiveOut
  .mix ( regressiveOut )
  .set { all_alignments }

refs
  .cross (all_alignments )
  .map { it -> [it[0][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5], it[0][1]] }
  .set { toEvaluate }

process evaluation {
    tag "${id}.${align_method}.${tree_method}.${align_type}.${bucket_size}"
    publishDir "${params.outdir}/individual_scores", mode: 'copy', overwrite: true
    label 'process_low'

    input:
      set val(id), \
          val(align_method), \
          val(tree_method), \
          val(align_type), \
          val(bucket_size), \
          file(test_alignment), \
          file(ref_alignment) \
          from toEvaluate

    output:
      set val(id), \
          val(tree_method), \
          val(align_method), \
          val(align_type), \
          val(bucket_size), \
          file("*.sp"), \
          file("*.tc"), \
          file("*.col") \
          into scores

    when:
      params.evaluate

     script:
     """
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"

       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"

       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode column \
            | grep -v "seq1" | grep -v '*' | \
              awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"

    """
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
