#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2020 IRB Barcelona

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

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  SimulateLcWES-nf: nextflow pipeline to simulate low coverage WES calling "
log.info "          with strelka2 based on existing BAM files          "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

params.help = null

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run main.nf --learning --trainingTable table.txt --features AO,DP,QVAL'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bam_folder                   FOLDER         Input folder containing BAM files to downsample (should be indexed).'
    log.info '    --downsampling_prop            INTEGER        Proportion of reads that will be randomly selected from the BAM (between 1 and 100).'
    log.info '    --ref                          FILE           Genome reference file.'
    log.info '    --strelka2                     PATH           Strelka2 installation dir.'
    log.info '    --tn_pairs                     FILE           Text file containing 2 columns: tumor=file name of tumor bams and normal=file name of normal bam (with colnames).'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --output_folder                FOLDER         Output folder (default: vf_output).'
    log.info '    --cpu                          INTEGER        Number of cpu to use with strelka2 (default=2).'
    log.info '    --mem                          INTEGER        Memory to be used in GB (default=8).'
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.bam_folder = null
params.downsampling_prop = null
params.ref = null
params.strelka2 = null
params.tn_pairs = null
params.output_folder = "calling_lowcovWES"
params.cpu = 2
params.mem = 8

if(params.bam_folder == null | params.downsampling_prop == null |  params.ref == null |  params.strelka2 == null | params.tn_pairs == null ){
  exit 1, "Please specify each of the following parameters: --bam_folder, --downsampling_prop, --ref, --strelka2, --tn_pairs"
}

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )

pairs = Channel.fromPath(params.tn_pairs).splitCsv(header: true, sep: '\t', strip: true)
  .map{ row -> [ file(params.bam_folder + "/" + row.tumor), file(params.bam_folder + "/" + row.tumor+'.bai'),
                 file(params.bam_folder + "/" + row.normal), file(params.bam_folder + "/" + row.normal+'.bai') ] }

workflow = params.strelka2 + '/bin/configureStrelkaSomaticWorkflow.py'

process samtoolsDownsampling {

  tag {bam_tag_t}

  input:
  file pair from pairs

  output:
  file '*bam*' into ds_bambai

  shell:
  bam_tag_t = pair[0].baseName
  bam_tag_n = pair[2].baseName
  '''
  #samtools view -s 3.!{params.downsampling_prop} -b !{pair[0]} -o !{bam_tag_t}_DS.bam
  #samtools index !{bam_tag_t}_DS.bam

  #samtools view -s 3.!{params.downsampling_prop} -b !{pair[2]} -o !{bam_tag_n}_DS.bam
  #samtools index !{bam_tag_n}_DS.bam

  touch !{bam_tag_t}_DS.bam
  touch !{bam_tag_t}_DS.bam.bai
  touch !{bam_tag_n}_DS.bam
  touch !{bam_tag_n}_DS.bam.bai
  '''
}

process strelka2Somatic {

     cpus params.cpu
     memory params.mem+'GB'

     tag {bam_tag_t}

     publishDir params.output_folder, mode: 'copy'

     input:
     file pair from ds_bambai
     file fasta_ref
     file fasta_ref_fai

     output:
     file '*vcf.gz' into vcffiles
     file '*bed.gz' optional true into regionfiles

     shell:
     bam_tag_t = pair[0].baseName
     '''
     !{workflow} --tumorBam !{pair[0]} --normalBam !{pair[2]} --referenceFasta !{fasta_ref} --exome --runDir strelkaAnalysis
     cd strelkaAnalysis
     ./runWorkflow.py -m local -j !{params.cpu} -g !{params.mem}
     cd ..
     mv strelkaAnalysis/results/variants/* .
     mv somatic.indels.vcf.gz !{pair[0]}_vs_!{pair[2]}.somatic.indels.vcf.gz
     mv somatic.snvs.vcf.gz !{pair[0]}_vs_!{pair[2]}.somatic.snvs.vcf.gz
     mv somatic.indels.vcf.gz.tbi !{pair[0]}_vs_!{pair[2]}.somatic.indels.vcf.gz.tbi
     mv somatic.snvs.vcf.gz.tbi !{pair[0]}_vs_!{pair[2]}.somatic.snvs.vcf.gz.tbi
     fixStrelkaOutput.sh *.vcf.gz
     if [ -d "strelkaAnalysis/results/regions" ]; then
          mv strelkaAnalysis/results/regions/* .
          mv somatic.callable.regions.bed.gz !{pair[0]}_vs_!{pair[2]}.somatic.callable.regions.bed.gz
          mv somatic.callable.regions.bed.gz.tbi !{pair[0]}_vs_!{pair[2]}.somatic.callable.regions.bed.gz.tbi
     fi
     '''
}
