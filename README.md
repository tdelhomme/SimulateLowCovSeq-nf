# SimulateLcWES

## Description
Nextflow pipeline to simulate low coverage WES from existing BAM files and call/annotate variants with Strelka2

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io).

2. External software:
- samtools
- strelka2
- annovar


## Input
  | Type      | Description     |
  |-----------|---------------|
  | bam_folder    | Input folder containing BAM files to downsample (should be indexed). |

  Specify the test files location


  log.info '    --bam_folder                   FOLDER         Input folder containing BAM files to downsample (should be indexed).'
  log.info '    --downsampling_prop            INTEGER        Proportion of reads that will be randomly selected from the BAM (between 1 and 100).'
  log.info '    --ref                          FILE           Genome reference file.'
  log.info '    --strelka2                     PATH           Strelka2 installation dir.'
  log.info '    --tn_pairs                     FILE           Text file containing 2 columns: tumor=file name of tumor bams and normal=file name of normal bam (with colnames).'
  log.info '    --avdb                         PATH           Path to annovar database.'
  log.info ''
  log.info 'Optional arguments:'
  log.info '    --output_folder                FOLDER         Output folder (default: calling_lowcovWES).'
  log.info '    --cpu                          INTEGER        Number of cpu to use with strelka2 (default=2).'
  log.info '    --mem                          INTEGER        Memory to be used in GB (default=8).'
  log.info '    --genome                       STRING         Reference genome version (default=hg38).'


## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --downsampling_prop    |            20 | Proportion of reads that will be randomly selected from the BAM (between 1 and 100). |
| --ref    |    /ref/hg38.fa | Genome reference file |
| --strelka2    |  /softs/strelka-2.9.10-0  | Strelka2 installation dir |
| --tn_pairs    |      pairs.txt | Text file containing 2 columns: tumor=file name of tumor bams and normal=file name of normal bam (with colnames) |
| --avdb    |      /db/annovar/humandb | Path to annovar database |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --output_folder   |      calling_lowcovWES | Output folder  |
| --cpu    |      2 | Number of cpu to use with strelka2 |
| --mem    |      8 | Memory to be used in GB  |
| --genome    |      hg38 | Reference genome version |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |


## Usage
  ```
  nextflow run tdelhomme/SimulateLcWES-nf --bam_folder BAM/ --downsampling_prop 20 --strelka2 /opt/anaconda3/envs/strelka/share/strelka-2.9.10-0 --ref /data/genome_references/GRCh38.d1.vd1.fa --tn_pairs pairs.txt --genome hg38
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | multianno.vcf.bgz | strelka2 VCF file annotated with annovar and compressed with bgzip |
  |  multianno.vcf.bgz.tbi   | tabix file |
