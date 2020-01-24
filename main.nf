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
    log.info '    --learning AND/OR --scoring        STRING               Indicate which mode (learning or scoring) to run.'
    log.info ''
    log.info '    if --learning:'
    log.info '      --trainingTable                  TXT                  File containing variant calls used to train the model. Must contain a column "status" if'
    log.info '                                                            options --duplicatedSeqVCF and --duplicatedSeqCov are not provided.'
    log.info '      --features                       LIST                 List of features used to train the model separated by commas, e.g. "--features AF,DP,RVSB".'
    log.info '                                                            Features are predifined when running with option --needlestack.'
    log.info '      if --replication'
    log.info '        --duplicatedSeqVCF             VCF                  Variant calls from an other sequencing, used to assign status to trainingTable.'
    log.info '        --duplicatedSeqCov             TXT                  File containing for each sequenced position in toTrainVCF the coverage in duplicatedSeq data'
    log.info '                                                            for quality check. Use iarcbioinfo/mpileup-nf pipeline on BAM/BED used to generate --toTrainVCF.'
    log.info ''
    log.info '    if --scoring:'
    log.info '      --targetTable                    TXT                  File containing variant calls on which models would be apply.'
    log.info '    if --scoring without --learning:'
    log.info '      --modelSNV                       RDATA                Variant calls from an other sequencing, used to assign status to toTrainTable.'
    log.info '      --modelINDEL                     RDATA                File containing for each sequenced position in toTrainVCF the coverage in duplicatedSeq data'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --output_folder                    FOLDER               Output folder (default: vf_output).'
    log.info '    --nsplit                           INTEGER              Split the input toTrainVCF to transform it into table in parallel.'
    log.info '    --needlestack                      FLAG                 Specify that calling was launched with needlestack, to use predifined features.'
    log.info ''
    exit 0
}
