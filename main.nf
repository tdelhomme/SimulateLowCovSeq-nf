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
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --output_folder                FOLDER         Output folder (default: vf_output).'
    log.info '    --cpu                          INTEGER        Number of cpu to use with strelka2 (default=2).'
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.bam_folder = null
params.downsampling_prop = null
params.ref = null
params.strelka2 = null

if(params.bam_folder == null | params.downsampling_prop == null |  params.ref == null |  params.strelka2 == null){
  exit 1, "Please specify each of the following parameters: --bam_folder, --downsampling_prop, --ref, --strelka2 "
}

bams = Channel.fromPath( params.bam_folder+'/*.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.input_folder}" }

bais = Channel.fromPath( params.bam_folder+'/*.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${params.input_folder}" }

process samtoolsDownsampling {

  tag {bam_tag}

  input:
  file bam from bams
  file bai from bais

  output:
  file("${bam_tag}*bam"), file("${bam_tag}*bai")  into ds_bambai

  shell:
  bam_tag = bam.baseName
  '''
  mkdir BAM_downsampled
  samtools view -s 3.!{params.downsampling_prop} -b !{bam} -o !{bam_tag}_DS.bam
  samtools index !{bam_tag}_DS.bam
  '''

}
