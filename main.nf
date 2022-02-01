#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2021 IRB Barcelona

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
log.info "  SimulateLowCovSeq: nextflow pipeline to simulate low coverage seq data "
log.info "          and run a calling with strelka2, based on existing BAM files   "
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
    log.info '    --target_coverage              NUMBER         Target coverage for the downsampled BAM file.'
    log.info '    --ref                          FILE           Genome reference file.'
    log.info '    --avdb                         PATH           Path to annovar database.'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --tn_pairs                     FILE           Text file containing 2 columns: tumor=file name of tumor bams and normal=file name of normal bam (with colnames).'
    log.info '    --strelka2                     PATH           Strelka2 installation dir.'
    log.info '    --output_folder                FOLDER         Output folder (default: calling_lowcovWES).'
    log.info '    --cpu                          INTEGER        Number of cpu to use with strelka2 (default=2).'
    log.info '    --mem                          INTEGER        Memory to be used in GB (default=8).'
    log.info '    --genome                       STRING         Reference genome version (default=hg38).'
    log.info 'Flags:'
    log.info '    --no_calling                                  Do not perform strelka2 variant calling and annotation'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.bam_folder = null
params.target_coverage = null
params.ref = null
params.strelka2 = null
params.tn_pairs = null
params.avdb = null
params.output_folder = "calling_LowCovSeq"
params.cpu = 2
params.mem = 8
params.genome = "hg38"
params.no_calling = null

if(params.bam_folder == null | params.target_coverage == null |  params.ref == null |  (params.strelka2 == null & params.no_calling == null & params.tn_pairs == null & params.avdb == null) ){
  exit 1, "Please specify each of the following parameters: --bam_folder, --target_coverage, --ref, (--strelka2 and --tn_pairs and --avdb) or --no_calling"
}

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )

avdb = file(params.avdb)

if(params.no_calling == null){
  pairs = Channel.fromPath(params.tn_pairs).splitCsv(header: true, sep: '\t', strip: true)
    .map{ row -> [ file(params.bam_folder + "/" + row.tumor), file(params.bam_folder + "/" + row.tumor+'.bai'),
                   file(params.bam_folder + "/" + row.normal), file(params.bam_folder + "/" + row.normal+'.bai') ] }

  workflow = params.strelka2 + '/bin/configureStrelkaSomaticWorkflow.py'


  process samtoolsDownsampling {

    publishDir params.output_folder+"/BAM/", mode: 'copy', pattern: '*_DS.bam*'

    tag {bam_tag_t}

    input:
    file pair from pairs

    output:
    set val(bam_tag_t), file("tumor*.bam"), file("tumor*.bai"), file("normal*.bam"), file("normal*.bai") into ds_bambai

    shell:
    bam_tag_t = pair[0].baseName
    bam_tag_n = pair[2].baseName
    '''
    tumorbam=!{pair[0]}
    normalbam=!{pair[2]}
    for inputbam in ${tumorbam} ${normalbam}
    do
      echo "Starting processing bam file: " ${inputbam}
      if [[ "$inputbam" == "$tumorbam" ]]; then pref=tumor; else pref=normal; fi
      bamtag=$(basename "$inputbam" | cut -d. -f1)
      declare -i meanreadlength
      meanreadlength=`samtools view ${inputbam} | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | awk 'BEGIN {total=0} {total += $1} END { print int(total/NR) }'`
      echo "mean read length: " $meanreadlength

      declare -i numberreads
      numberreads=`samtools idxstats ${inputbam} | awk 'BEGIN {total=0} {total += $3} END {print total}'`
      echo "total number of reads: " $numberreads

      declare -i lengthsequence
      lengthsequence=`samtools idxstats ${inputbam} | awk 'BEGIN {total=0} {total += $2} END {print total}'`
      echo "total sequence length: " $lengthsequence

      meancov=$((meanreadlength * numberreads / lengthsequence))
      echo "estimated mean coverage: " $meancov
      targetcov=!{params.target_coverage}
      echo "requested target coverage: " $targetcov
      # samtools_ds=$((100 / (meancov / targetcov)))
      baseDir=!{baseDir}
      samtools_ds=`Rscript  ${baseDir}/bin/compute_downsampling_proportion.R $meancov $targetcov`
      echo "downsampling proportion in percent: " $samtools_ds
      
      len=`expr length "$samtools_ds"`
      if [ 2 -gt "$len" ]; then
        samtools_ds=0$samtools_ds
      fi

      echo "command samtools view -s 42.${samtools_ds} -b ${inputbam} -o ${bamtag}_DS.bam"
      samtools view -s 42.${samtools_ds} -b ${inputbam} -o ${pref}_${bamtag}_DS.bam
      samtools index ${pref}_${bamtag}_DS.bam
    done
    '''
  }
} else {

	  bams  = Channel.fromPath(params.bam_folder +'/*bam')

    process samtoolsDownsampling_calling {

    publishDir params.output_folder+"/BAM/", mode: 'copy', pattern: '*_DS.bam*'

    tag {bam_tag_t}

    input:
    file bam from bams

    output:
    file '*DS.bam*' into ds_bambai

    shell:
    bam_tag =bam.baseName
    '''
    declare -i meanreadlength
    meanreadlength=`samtools view file.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | awk 'BEGIN {total=0} {total += $1} END { print int(total/NR) }'`
    echo "mean read length: " $meanreadlength

    declare -i numberreads
    numberreads=`samtools idxstats file.bam | awk 'BEGIN {total=0} {total += $3} END {print total}'`
    echo "total number of reads: " $numberreads

    declare -i lengthsequence
    lengthsequence=`samtools idxstats file.bam | awk 'BEGIN {total=0} {total += $2} END {print total}'`
    echo "total sequence length: " $lengthsequence

    meancov=$((meanreadlength * numberreads / lengthsequence))
    targetcov=!{params.target_coverage}
    samtools_ds=$((meancov / targetcov))
    echo "final downsample proportion: " $samtools_ds
    
    samtools view -s 42.${samtools_ds} -b !{bam} -o !{bam_tag}_DS.bam
    samtools index !{bam_tag}_DS.bam
    '''
  }
}

if(params.no_calling == null){

  process strelka2Somatic {

       cpus params.cpu
       memory params.mem+'GB'

       publishDir params.output_folder+"/PASS/", mode: 'copy', pattern: '*PASS.vcf'

       tag {bam_tag_t}

       input:
       set val(tumor_id), file(tumor_bamds), file(tumor_baids), file(normal_bamds), file(normal_baids) from ds_bambai
       file fasta_ref
       file fasta_ref_fai

       output:
       file '*snvs.vcf.gz' into vcffiles
       file '*PASS.vcf' into passvcf

       shell:
       bam_tag_t = tumor_id
       '''
       !{workflow} --tumorBam !{tumor_bamds} --normalBam !{normal_bamds} --referenceFasta !{fasta_ref} --reportEVSFeatures --runDir strelkaAnalysis
       cd strelkaAnalysis
       ./runWorkflow.py -m local -j !{params.cpu} -g !{params.mem}
       cd ..
       mv strelkaAnalysis/results/variants/* .
       mv somatic.indels.vcf.gz !{tumor_bamds}_vs_!{normal_bamds}.somatic.indels.vcf.gz
       mv somatic.snvs.vcf.gz !{tumor_bamds}_vs_!{normal_bamds}.somatic.snvs.vcf.gz
       mv somatic.indels.vcf.gz.tbi !{tumor_bamds}_vs_!{normal_bamds}.somatic.indels.vcf.gz.tbi
       mv somatic.snvs.vcf.gz.tbi !{tumor_bamds}_vs_!{normal_bamds}.somatic.snvs.vcf.gz.tbi
       fixStrelkaOutput.sh *.vcf.gz

       vcf=`ls *snvs.vcf.gz`
       bcftools view -f PASS ${vcf} > "${vcf/.vcf.gz}_PASS.vcf"
       '''
  }

  process annotation {

      publishDir params.output_folder, mode: 'copy'

      tag {vcf_tag}

      input:
      file vcf from vcffiles
      file avdb

      output:
      file '*multianno.vcf.bgz' into bgzipvcfanno
      file '*multianno.vcf.bgz.tbi' into tabixvcfanno

      shell:
      vcf_tag=vcf.baseName
      '''
      table_annovar.pl !{vcf} !{avdb} -buildver !{params.genome} -out !{vcf_tag} -remove -protocol refGene -operation g -nastring . -vcfinput
      bgzip -c !{vcf_tag}.!{params.genome}_multianno.vcf > !{vcf_tag}.!{params.genome}_multianno.vcf.bgz
      tabix -p vcf !{vcf_tag}.!{params.genome}_multianno.vcf.bgz
      '''
  }

}
