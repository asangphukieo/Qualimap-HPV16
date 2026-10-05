#!/usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


params.help   = null
params.multiqc_config = 'NO_FILE'
params.feature_file   = 'NO_FILE'
params.cpu = 1
params.mem = 40
params.input_folder = null
params.output_folder = "."
params.output_format = "html"

log.info ""
log.info "----------------------------------------------------------------"
log.info "           Quality control with Qualimap  and MultiQC           "
log.info "----------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info "=========== Pre-installed required =============="
log.info "qualimap, samtools, multiqc, rkmh (~/bin), datamash "
log.info "--------------------------------------------------------"
if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "                     USAGE : version 1.1                            "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "-------------------QC-------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/Qualimap.nf   --qualimap /path/to/qualimap  --multiqc /path/to/multiqc --samtools /path/to/samtools --input_folder /path/to/bam  --output_folder /path/to/output"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input_folder         FOLDER               Folder containing bam files"
    log.info ""
    log.info "Optional arguments:"
    log.info "--feature_file         FILE                 Qualimap feature file for coverage analysis"
    log.info "--output_folder        PATH                 Output directory for html and zip files (default=.)"
    log.info "--output_format        STRING               Output format for individual qualimap files, html or pdf (default=html)"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=1)"
    log.info '--multiqc_config       STRING               Config yaml file for multiqc (default : none)'
    log.info "--mem                  INTEGER              Size of memory used. Default 4Gb"
    log.info "--script_folder        PATH                 Script directory (default=.)"
    log.info "--ref_folder        PATH                    Reference genome directory (default=.)" 
    log.info "--rkmh_data_folder        PATH              Reference genome directory used by rkmh (default=.)"    
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    exit 0
}else{
/* Software information */
  log.info "input_folder   = ${params.input_folder}"
  log.info "cpu            = ${params.cpu}"
  log.info "mem            = ${params.mem}"
  log.info "feature_file   = ${params.feature_file}"
  log.info "output_folder  = ${params.output_folder}"
  log.info "output_format  = ${params.output_format}"
  log.info "multiqc_config = ${params.multiqc_config}"
  log.info "script_folder = ${params.script_folder}"  
  log.info "ref_folder = ${params.ref_folder}" 
  log.info "rkmh_data_folder = ${params.rkmh_data_folder}"    
  log.info "help           = ${params.help}"
}

assert (params.input_folder != null) : "please provide the --input_folder option"

qualimap_ff = file(params.feature_file)
multiqc_config = file(params.multiqc_config)

bams_init = Channel.fromPath( params.input_folder+'*.bam')
              .ifEmpty { error "Cannot find any bam file in: ${params.input_folder}" }
data_folder = file(params.rkmh_data_folder)

bams_init.into{ bams }


process extractHPVandUNMAP {
    cpus params.cpu
    memory params.mem+'G'
    queue 'low_p'
    tag { bam_tag }

    publishDir "${params.output_folder}/individual_reports", mode: 'copy'

    input:
    file bam from bams

    output:
    file ("${bam_tag}_hpv_unmap.bam") into hpv_bam,hpv_bam2,hpv_bam3

    shell:
    bam_tag=bam.baseName

    '''
    samtools index !{bam}
    hpv_read=\$(samtools idxstats !{bam} | grep gi| cut -f1 | awk '{printf $0" "}')
    samtools view -F 0x4 -o !{bam_tag}_hpv.bam !{bam} \${hpv_read}
    samtools view -f 0x4 -hb !{bam} > !{bam_tag}_unmap.bam
    samtools merge -f -cp !{bam_tag}_hpv_unmap.bam !{bam_tag}_hpv.bam !{bam_tag}_unmap.bam

    '''
}

process qualimap {
    cpus params.cpu
    memory params.mem+'G'
    queue 'low_p'
    tag { bam_tag }

    publishDir "${params.output_folder}/individual_reports", mode: 'copy'

    input:
    file bam from hpv_bam
    file qff from qualimap_ff

    output:
    file ("${bam_tag}") into qualimap_results

    shell:
    bam_tag=bam.baseName
    feature = qff.name != 'NO_FILE' ? "--feature-file $qff" : ''
    mem_qm = params.mem -2 //params.mem.intdiv(2)
    '''
    qualimap bamqc -nt !{params.cpu} !{feature} --skip-duplicated -bam !{bam} -outdir !{bam_tag} -outformat !{params.output_format}
    '''
}

process flagstat {
    cpus params.cpu
    memory params.mem+'G'
    queue 'low_p'
    tag { bam_tag }

    publishDir "${params.output_folder}/individual_reports", mode: 'copy'

    input:
    file bam from hpv_bam2

    output:
    file ("${bam_tag}.stats.txt") into flagstat_results

    shell:
    bam_tag=bam.baseName
    '''
    samtools flagstat --threads !{params.cpu} !{bam} > !{bam_tag}.stats.txt
    '''
}

process hpv16check {
    errorStrategy 'ignore'
    cpus params.cpu
    memory params.mem+'G'
    queue 'low_p'
    tag { bam_tag }

    publishDir "${params.output_folder}/individual_reports", mode: 'copy'

    input:
    file bam from hpv_bam3

    output:
    file ("${bam_tag}.qc_bam_hpv16check.txt") into hpv16_results
    //file ("${bam_tag}.gi|333031|lcl|HPV16REF.1|.cov.pdf") into hpv16_results_pdf
    //file ("${bam_tag}_prim_NM_1.bam") 
    //file ("${bam_tag}_prim_NM_2.bam") 
    file ("${bam_tag}_prim_hpv16.bam") into hpv16_pri_bam,hpv16_pri_bam2
    file ("${bam_tag}_unmap.bam")
    file ("${bam_tag}_prim16_unmap.bam")

    shell:
    bam_tag=bam.baseName

    """
    samtools sort !{bam} -o !{bam_tag}_sorted.bam
    samtools index !{bam_tag}_sorted.bam
    
    total_read_original=\$(samtools view !{bam_tag}_sorted.bam |wc -l)
    total_read_qc_pri=\$(samtools view -q 30 -f 0x2 !{bam_tag}_sorted.bam| wc -l)
    most_read_map=\$(samtools idxstats !{bam_tag}_sorted.bam | grep gi | sort -nrk 4 | head -1| cut -f1)
    samtools view -f 0x2 -F 256 -q 30 -hb !{bam_tag}_sorted.bam "gi|333031|lcl|HPV16REF.1|:1-7906" > !{bam_tag}_prim_hpv16.bam
    samtools view -f 0x4 -hb !{bam_tag}_sorted.bam > !{bam_tag}_unmap.bam

    total_read=\$(samtools view !{bam_tag}_prim_hpv16.bam | wc -l)
    total_read_unmap=\$(samtools view !{bam_tag}_unmap.bam | wc -l)
    E6=\$(samtools view -q 30 -f 0x2 !{bam_tag}_sorted.bam 'gi|333031|lcl|HPV16REF.1|:104-559'| wc -l)
    E7=\$(samtools view -q 30 -f 0x2 !{bam_tag}_sorted.bam 'gi|333031|lcl|HPV16REF.1|:562-858'| wc -l)
    E1=\$(samtools view -q 30 -f 0x2 !{bam_tag}_sorted.bam 'gi|333031|lcl|HPV16REF.1|:865-2814'| wc -l)
    E2=\$(samtools view -q 30 -f 0x2 !{bam_tag}_sorted.bam 'gi|333031|lcl|HPV16REF.1|:2756-3853'| wc -l)
    E5_a=\$(samtools view -q 30 -f 0x2 !{bam_tag}_sorted.bam 'gi|333031|lcl|HPV16REF.1|:3850-4101'| wc -l)
    L2=\$(samtools view -q 30 -f 0x2 !{bam_tag}_sorted.bam 'gi|333031|lcl|HPV16REF.1|:4237-5658'| wc -l)
    L1=\$(samtools view -q 30 -f 0x2 !{bam_tag}_sorted.bam 'gi|333031|lcl|HPV16REF.1|:5639-7156'| wc -l)
    E6_f=\$(echo "(\${E6} * 1000000 * 1000) / (\${total_read} * 455)" | bc )
    E7_f=\$(echo "(\${E7} * 1000000 * 1000) / (\${total_read} * 296)" | bc )
    E1_f=\$(echo "(\${E1} * 1000000 * 1000) / (\${total_read} * 1949)" | bc )
    E2_f=\$(echo "(\${E2} * 1000000 * 1000) / (\${total_read} * 1097)" | bc )
    E5_f=\$(echo "(\${E5_a} * 1000000 * 1000) / (\${total_read} * 251)" | bc )
    L2_f=\$(echo "(\${L2} * 1000000 * 1000) / (\${total_read} * 1421)" | bc )
    L1_f=\$(echo "(\${L1} * 1000000 * 1000) / (\${total_read} * 1517)" | bc )

    sample_id=\$(echo !{bam_tag} | cut -f1 -d".")
    out=\$(echo -e "\${sample_id}\t\${total_read_original}\t\${total_read_qc_pri}\t\${most_read_map}\t\${total_read}\t\${total_read_unmap}\t\${E6}\t\${E7}\t\${E1}\t\${E2}\t\${E5_a}\t\${L2}\t\${L1}\t\${E6_f}\t\${E7_f}\t\${E1_f}\t\${E2_f}\t\${E5_f}\t\${L2_f}\t\${L1_f}" )
    echo \${out} > !{bam_tag}.qc_bam_hpv16check.txt

    samtools merge -f -cp !{bam_tag}_prim16_unmap.bam !{bam_tag}_prim_hpv16.bam !{bam_tag}_unmap.bam
    
    #samtools calmd !{bam_tag}_prim_hpv16.bam ${params.ref_folder}hpv16_ref.fa > !{bam_tag}_prim_NM.bam
    #samtools view -hbf 64 !{bam_tag}_prim_NM.bam > !{bam_tag}_prim_NM_1.bam
    #samtools view -hbf 128 !{bam_tag}_prim_NM.bam > !{bam_tag}_prim_NM_2.bam
    #python ${params.script_folder}MappedReadPlot.py !{bam_tag}_prim_NM_1.bam -r ${params.ref_folder}hpv16_ref.fa -a ${params.ref_folder}hpv16_gene_annot.tsv -o !{bam_tag}    
    """
}


process multiqc {
    cpus params.cpu
    memory params.mem+'G'

    publishDir "${params.output_folder}", mode: 'copy'

    input:
    file qualimap_results from qualimap_results.collect()
    file flagstat_results from flagstat_results.collect()
    file config_file from multiqc_config

    output:
    file("multiqc_report.html") into final_output
    file("multiqc_data/") into final_output_data
    file("multiqc_data/multiqc_qualimap_bamqc_genome_results.txt") into final_table
    file("multiqc_data/multiqc_samtools_flagstat.txt") into final_table2
    file("multiqc_data/multiqc_general_stats.txt") into final_table3
    

    shell:
    config = config_file.name != 'NO_FILE' ? "--config $config_file" : ''
    '''
    multiqc !{config} .
    '''
}

//to make input for rkmh
process bam2fq {
    cpus params.cpu
    memory params.mem+'G'
    queue 'low_p'
    tag { bam_tag }

    publishDir "${params.output_folder}", mode: 'copy'

    input:
    file bam from hpv16_pri_bam

    output:
    file("${bam_tag}.fastq") into fastq_out

    shell:
    bam_tag=bam.baseName
    
    """
    samtools bam2fq !{bam} > !{bam_tag}.fastq
    """
}

process rkmh {
    cpus params.cpu
    memory params.mem
    queue 'low_p'
    tag { bam_tag }

    publishDir "${params.output_folder}", mode: 'copy'

    input:
    file fastq from fastq_out
    path data_folder

    output:
    file("${bam_tag}.rk") 
    file("${bam_tag}.cls") 
    file("${bam_tag}_rkmh.out") into rkmh_out

    shell:
    bam_tag=fastq.baseName
    """
    rkmh hpv16 -t ${params.cpu} -f !{fastq} > !{bam_tag}.rk
    python ${params.script_folder}rkmh/scripts/score_real_classification.py < !{bam_tag}.rk > !{bam_tag}.cls
    sample_id=\$(echo !{bam_tag} |cut -f1 -d "." )
    paste <(echo \${sample_id}) <(cut -f1 -d ' ' !{bam_tag}.cls|cut -f1 -d';') <(cut -f2 -d ' ' !{bam_tag}.cls|cut -f1 -d';') <(cut -f4 -d ' ' !{bam_tag}.cls|cut -f3 -d':') <(cut -f5 -d ' ' !{bam_tag}.cls|cut -f3 -d':') > !{bam_tag}_rkmh.out
    """

}


process hpv16ReadDepth {
    errorStrategy 'ignore'
    cpus params.cpu
    memory params.mem+'G'
    queue 'high_p'
    tag { bam_tag }

    publishDir "${params.output_folder}/individual_reports", mode: 'copy'

    input:
    file bam from hpv16_pri_bam2

    output:
    file ("${bam_tag}.hpv16_depth.txt") into hpv16_depth_results

    shell:
    bam_tag=bam.baseName
    """
    paste <(cat <(echo "#IID") <(echo !{bam_tag} |cut -f1 -d ".")) \\
    <(cat <(echo -e "Average_Depth\tQ1_Depth\tQ3_Depth\tMedian_Depth\tStd_Depth") <(samtools depth !{bam} | cut -f3 | datamash mean 1 q1 1 q3 1 median 1 sstdev 1)) \\
    <(samtools mpileup !{bam} | perl ${params.script_folder}bam_coverage.pl) \\
    <(cat <(echo -e "average_read_length\tmaximum_read_length\taverage_quality\tinsert_size_average\tinsert_size_standard_deviation") <(samtools stats !{bam} | grep "average length\\|maximum length\\|average quality\\|insert size average\\|insert size standard deviation" |cut -f3 |tr '\n' '\t')) \\
     > !{bam_tag}.hpv16_depth.txt
    """
    
}

process summaryTable {
    cpus params.cpu
    memory params.mem+'G'

    publishDir "${params.output_folder}", mode: 'copy'

    input:
    file "multiqc_qualimap_bamqc_genome_results.txt" from final_table
    file "multiqc_samtools_flagstat.txt" from final_table2
    file "multiqc_general_stats.txt" from final_table3
    file hpv16_results from hpv16_results.collect()
    file rkmh_out from rkmh_out.collect()
    file hpv16_depth_results from hpv16_depth_results.collect()

    output:
    file("summary_table.txt") 
    file("summary_table.tmp1") 
    file("summary_table.tmp2") 
    file("summary_table.tmp3") 
    file("summary_table.tmp4") 
    file("summary_table.tmp5") 
    file("rkmh_table.tmp") 
    file("hpv16_depth_table.tmp")     

    script:
    """
    join -t \$'\t' -a1 -a2 -j 1 <(tail -n+2 multiqc_samtools_flagstat.txt| sort -k1 ) <(tail -n+2 multiqc_qualimap_bamqc_genome_results.txt |sort -k1 ) > summary_table.tmp
    join -t \$'\t' -a1 -a2 -j 1 <(sort -k1 summary_table.tmp) <(tail -n+2 multiqc_general_stats.txt |sort -k1 ) > summary_table.tmp1

    cat ${hpv16_results} | sed s/" "/"\t"/g > summary_table.tmp2
    cat ${rkmh_out}  | sed s/" "/"\t"/g |sed s/"_prim_hpv16"//g > rkmh_table.tmp
    grep -h -v "#" ${hpv16_depth_results}  | sed s/" "/"\t"/g > hpv16_depth_table.tmp

    join -t \$'\t' -a1 -a2 -j 1 -e "NA" -o auto <(sort -k1 summary_table.tmp1) <(sort -k1 summary_table.tmp2) | sed s/" "/"\t"/g > summary_table.tmp3
    join -t \$'\t' -a1 -a2 -j 1 -e "NA" -o auto <(sort -k1 summary_table.tmp3) <(sort -k1 rkmh_table.tmp ) > summary_table.tmp4
    join -t \$'\t' -a1 -a2 -j 1 -e "NA" -o auto <(sort -k1 summary_table.tmp4) <(sort -k1 hpv16_depth_table.tmp ) > summary_table.tmp5
    
    paste <(head -n1 multiqc_samtools_flagstat.txt) <(head -n1 multiqc_qualimap_bamqc_genome_results.txt| sed s/"^Sample\t"/""/g) <(head -n1 multiqc_general_stats.txt| sed s/"^Sample\t"/""/g) <(echo -e "total_read_original\ttotal_read_qc_pri\tmost_read_map\ttotal_read_hpv16\ttotal_read_unmap\tE6\tE7\tE1\tE2\tE5_a\tL2\tL1\tE6_f\tE7_f\tE1_f\tE2_f\tE5_f\tL2_f\tL1_f") <(echo -e "Lineage\tSublineage\tLineage_counts\tSublineage_counts") \
    <(echo -e "Average_Depth\tQ1_Depth\tQ3_Depth\tMedian_Depth\tStd_Depth\tBase_count\tGenome_length\tMean_coverage\tMin_coverage\tMax_coverage\taverage_length\tmaximum_length\taverage_quality\tinsert_size_average\tinsert_size_standard_deviation") > header.tmp
    cat header.tmp summary_table.tmp5 > summary_table.txt

    """
    //Note: 
    //1. Samtools flagstat column (Sample	total_passed	total_failed	secondary_passed	secondary_failed	supplementary_passed	supplementary_failed	duplicates_passed	duplicates_failed	mapped_passed	mapped_failed	mapped_passed_pct	mapped_failed_pct	paired in sequencing_passed	paired in sequencing_failed	read1_passed	read1_failed	read2_passed	read2_failed	properly paired_passed	properly paired_failed	properly paired_passed_pct	properly paired_failed_pct	with itself and mate mapped_passed	with itself and mate mapped_failed	singletons_passed	singletons_failed	singletons_passed_pct	singletons_failed_pct	with mate mapped to a different chr_passed	with mate mapped to a different chr_failed	with mate mapped to a different chr (mapQ >= 5)_passed	with mate mapped to a different chr (mapQ >= 5)_failed	flagstat_total)
    //   is the quality checking of starting file (HPV reads without human reads)
    //2. Qualimap columns (Sample	bam_file	total_reads	mapped_reads	mapped_bases	sequenced_bases	mean_insert_size	median_insert_size	mean_mapping_quality	general_error_rate	mean_coverage	percentage_aligned)
    //   is the quality checking of starting file (HPV reads without human reads) mapping to HPV16 genome
    //3. HPV16 check column (flagstat_total	total_read_original	total_read_qc_pri	most_read_map	total_read_hpv16	total_read_unmap	E6	E7	E1	E2	E5_a	L2	L1	E6_f	E7_f	E1_f	E2_f	E5_f	L2_f	L1_f    Lineage	Sublineage	Lineage_counts	Sublineage_counts)
    //   is the quality checking of high-quality reads aligned to HPV16 reference
    //
}