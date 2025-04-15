# Sequence data processing
This Markdown file describes the steps taken to process and align the raw reads obtained from whole genome sequencing island and mainland *P. maniculatus* genomes. 

## Sequencing details
83 whole genomes were sequenced to a target coverage of 30X using 2 NovaSeq X+ 10B flowcells plus a single lane. The 83 samples include 79 unique individuals plus 4 duplicated samples (FSB114-DUPLICATE, FSB42-DUPLICATE, FSB85-DUPLICATE, and FSB6-DUPLICATE).

Most of the samples are represented by 17 distinct lanes (2 full flowcells plus a single lane), but some are represented by 18 distinct lanes. For the purposes of data processing, I recoded the lanes (denoted by the string L00{X} in each FASTQ file name) to run from 1-17 (or 18) such that all FASTQ files posess unique names.

## Software information
All of the software used for this component are described in the [packages](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/packages) directory.

## Pipeline
This section details the exact commands and parameters used for read quality control and alignment. The following wildcards were used in constructing file names for the example commands:
- Sample name: `{sample}`
- Population name: `{pop}`
- Miscellaneous numbers: `{X}`
- Directory: `{dir}`
- Chromosome: `{chrom}`

#### Read quality control
First, we trimmed adapter sequences in all raw FASTQ files with BBduk. This was done using the paired-end setting of BBduk which takes in the forward and reverse read files and creates a trimmed output file for each. The [`adapters.fa`](https://github.com/PayseurLabUWMadison/gi_demography_inversions/blob/main/sequence_processing/adapters.fa) file contains possible adapter sequences.
```
bbduk.sh -Xmx1g in1={sample}_S{X}_L00{X}_R1_00{X}.fastq.gz in2={sample}_S{X}_L00{X}_R2_00{X}.fastq.gz out1={sample}_L00{X}_adaptrim_1P.fastq.gz out2={sample}_L00{X}_adaptrim_2P.fastq.gz  ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
```

After removing adapters, we trimmed back the low quality ends of reads with Trimmomatic. 
```
java -jar trimmomatic-0.39.jar PE {sample}_L00{X}_adaptrim_1P.fastq.gz {sample}_L00{X}_adaptrim_2P.fastq.gz  -baseout {sample}_L00{X}_adaptrim_qualtrim.fastq.gz SLIDINGWINDOW:4:15 MINLEN:36
```

#### Alignment
We used the BWA-MEM algorithm to map the processed, paired-end reads to the *P. maniculatus* reference genome [HU_Pman_2.1.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003704035.1). This step was performed separately for each pair of FASTQ files.
```
bash read_group.sh GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz {sample}_L00{X}_adaptrim_qualtrim_1P.fastq.gz {sample}_L00{X}_adaptrim_qualtrim_1P.fastq.gz {sample} {sample}_L00{X}_sorted.bam
```
>Note: The [`read_group.sh`](https://github.com/PayseurLabUWMadison/gi_demography_inversions/blob/main/sequence_processing/read_group.sh) script is required to properly format the RG tags in the output BAM file. The BWA-MEM algorithm is run from inside this bash script. The BWA command within this script has the following format: `bwa mem -t 8 -R $(echo "@RG\tID:${ID}\tPL:${PL}\tPU:${PU}\tLB:${LB}\tSM:${SM}") $REF $R1 $R2 | samtools sort -o $OUT`.

Then, for each sample, we merged together the lane-separated alignments such that each sample is subsequently represented by a single BAM file. 
```
samtools merge {sample}_sorted_merged.bam {sample}_L00{X}_sorted.bam {sample}_L00{Y}_sorted.bam {sample}_L00{Z}_sorted.bam
```

Finally, duplicate reads present in the alignments can be marked for removal with Picard MarkDuplicates.
```
gatk MarkDuplicates --TMP_DIR {dir} -I {sample}_sorted_merged.bam -O {sample}_sorted_merged_rmdups.bam
```
