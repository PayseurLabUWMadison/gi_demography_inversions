# Inversion identification, genotyping, and polarization
This Markdown file describes the steps taken to identify large-scale inversion polymorphisms from short-read whole genome sequences and to reconstruct the ancestral state of each rearrangement. 

## Methodological details
Following sequence data processing, variant calling, and variant filtering, we have three single-population VCFs that contain hard-filtered, biallelic SNPs found in each population sample as well as one multi-population VCF that contains the union of hard-filtered, biallelic SNPs found within each population sample. These steps are described [here](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/sequence_processing) and [here](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/variant_calling).

As detailed in the Materials and Methods section of the paper, we pruned each population sample for close relatives, giving us the following versions of the above callsets that exclude close relatives:
- multi-population: `unrelated_snpOnly_filteredPASS_gulf_islands_sorted_merged_rmdups.vcf.gz`
- single-population: `unrelated_snpOnly_filteredPASS_{pop}_sorted_merged_rmdups.vcf.gz` (where `{pop}` could be `saturna`, `pender`, or `maple_ridge`)

In the first step of this analysis, we demonstrate how polymorphic inversions can be detected based on the distortions they cause in allele frequency distributions and the decay of LD. We use PCA to then confirm their presence in the data. In the second step of this analysis, we leverage the reference assemblies of four additional *Peromyscus* species to polarize these segregating rearrangements in the *P. maniculatus* genome. 

## Software information
All of the software used for this component are described in the [packages](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/packages) directory.

## Code
This section details the commands and parameters used for inversion identification, genotyping, and polarization. The following wildcards were used in constructing file names for the example commands:
- Population name: `{pop}`
- Chromosome: `{chrom}`
- Species: `{species}`
- Miscellaneous numbers: `{X}`, `{Y}`, etc.

#### Inversion identification and genotyping
To assay the impact of inversions on allele frequencies and conduct the MAF-partitioned genotype analysis described in the paper to "screen" for segregating inversions, we used the VariantsToTable function in GATK to extract information about SNP positions, minor allele frequencies, and individual genotypes from each single-population VCF.
```
gatk VariantsToTable -V unrelated_snpOnly_filteredPASS_{pop}_sorted_merged_rmdups.vcf.gz -F CHROM -F POS -F REF -F ALT -F AN -F AC -GF GT -O {pop}_genotypes.table
```

To assay the impact of inversions of LD decay, we used PLINK to obtain information about the long-range decay of average pairwise r^2 across each chromosome in each population.
```
plink --vcf unrelated_snpOnly_filteredPASS_{pop}_sorted_merged_rmdups.vcf.gz --allow-extra-chr --chr {chrom} --geno 0 --mac 1 -r2 gz --ld-window 1000000000 --ld-window-kb 100000 --ld-window-r2 0 --out {pop}_{chrom}
```

To confirm that segregating inversions yield PCA clusters consistent with the expected homozygous and heterozygous genotypic classes, we conducted PCA with PLINK using SNPs within the bounds of previously-identified inversion breakpoints. This was performed using the combined, multi-population callset.
```
plink --vcf unrelated_snpOnly_filteredPASS_gulf_islands_sorted_merged_rmdups.vcf.gz --allow-extra-chr --chr {chrom} --geno 0 --from-bp {X} --to-bp {X} --pca 47 --out {chrom}_{X}-{Y}
```
Here, the interval `{X}`-`{Y}` refers to inversion breakpoints defined by [Harringmeyer and Hokestra (2022)](https://www.nature.com/articles/s41559-022-01890-0). 

#### Inversion polarization
To determine the state (i.e., ancestral vs derived) of segregating inversions in the *P. maniculatus* reference genome (HU_Pman_2.1.3), we inspected the status of these arrangements in the reference assemblies of four additional species:
1. *P. polionotus* ([HU_Ppol_1.3.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003704135.2/))
2. *P. leucopus* ([UCI_PerLeu_2.1](https://www.ncbi.nlm.nih.gov/assembly/8220561))
3. *P. eremicus* ([PerEre_H2_v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_949786415.1/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=refseq-release-220-20230911))
4. *P. californicus* ([ASM782708v3](https://www.ncbi.nlm.nih.gov/assembly/12613101))

Prior to conducting whole genome alignments between each of these species and *P. maniculatus*, we downloaded each FASTA file and used seqkit to recode NCBI RefSeq contigs to reflect their chromosome assignments.
```
cat {species}.fna | seqkit replace -p "\s.+" | seqkit replace -p '^(\S+)' -r '{kv}$2' -k {species}.alias | seqkit grep -r -p ^chr > {species}_chromID.fna
```
Here, the `{species}.alias` files are tab-separated lists organized as `[old contig name] [new contig name]`. These are present in this directory for all species with chromosome assignments.

Then we conducted the initial alignments with nucmer, filtered them, summarized them, and generated the corresponding dotplots.
```
nucmer --maxmatch -c 500 -b 500 -l 100 -p maniculatus_{species} maniculatus_chromID.fna {species}_chromID.fna
```
```
delta-filter -m -i 90 -l 100 maniculatus_{species}.delta > filtered_maniculatus_{species}.delta
```
```
show-coords -THrd filtered_maniculatus_{species}.delta > maniculatus_{species}.coords
```
```
dotplot -c maniculatus_{species}.coords -r maniculatus_chromID.fna -q {species}_chromID.fna -o maniculatus_{species}.pdf
```

As described in the Materials and Methods section of the paper, based on the resulting dotplots, we used SeqKit to correct reverse-strand errors and re-assign contig names in the *P. leucopus*, *P. eremicus*, and *P. californicus* references to reflect their homology with the *P. maniculatus* reference. The following sections describe these steps for each of these three species. 

###### *P. leucopus*
First, reassign contig names in the reference FASTA
```
cat leucopus_chromID.fna | seqkit replace -p '^(\S+)' -r '{kv}$2' -k leucopus_reAssigned.alias > leucopus_reAssigned.fna
```
Then, reverse complement chromosomes exhibiting reverse-strand errors.
```
cat leucopus_reAssigned.fna | seqkit grep -n -p "chr1" -p "chr5" -p "chr9" -p "chr11" -p "chr12" -p "chr17" -p "chr18" -p "chr20" -p "chrX" | seqkit seq -r -p > leucopus_reversedChroms.fna
```
Finally, remove these reverse complemented chromosomes from the full FASTA and concatenate with the corrected FASTA.
```
cat leucopus_reAssigned.fna | seqkit grep -n -p "chr1" -p "chr5" -p "chr9" -p "chr11" -p "chr12" -p "chr17" -p "chr18" -p "chr20" -p "chrX" -v > leucopus_correctChroms.fna
```
```
seqkit concat leucopus_correctChroms.fna leucopus_reversedChroms.fna --full > leucopus_alignReady.fna
```

###### *P. eremicus*
First, reassign contig names in the reference FASTA
```
cat eremicus_chromID.fna | seqkit replace -p '^(\S+)' -r '{kv}$2' -k eremicus_reAssigned.alias > eremicus_reAssigned.fna
```
Then, reverse complement chromosomes exhibiting reverse-strand errors.
```
cat eremicus_reAssigned.fna | seqkit grep -n -p "chr20" -p "chr1" -p "chr5" -p "chr9" -p "chr11" -p "chr12" -p "chr17" -p "chr18" -p "chrX" | seqkit seq -r -p > eremicus_reversedChroms.fna
```
Finally, remove these reverse complemented chromosomes from the full FASTA and concatenate with the corrected FASTA.
```
cat eremicus_reAssigned.fna | seqkit grep -n -p "chr20" -p "chr1" -p "chr5" -p "chr9" -p "chr11" -p "chr12" -p "chr17" -p "chr18" -p "chrX" -v > eremicus_correctChroms.fna
```
```
seqkit concat eremicus_correctChroms.fna eremicus_reversedChroms.fna --full > eremicus_alignReady.fna
```

###### *P. californicus*
First, reassign contig names in the reference FASTA
```
cat californicus.fna | seqkit replace -p "\s.+" | seqkit replace -p '^(\S+)' -r '{kv}$2' -k californicus.alias | seqkit grep -r -p ^chr > californicus_reAssigned.fna
```
Then, reverse complement chromosomes exhibiting reverse-strand errors.
```
cat californicus_reAssigned.fna | seqkit grep -n -p "chr1" -p "chr2" -p "chr5" -p "chr8" -p "chr9" -p "chr13" -p "chr16" -p "chr17" -p "chr22" -p "chrX" | seqkit seq -r -p > californicus_reversedChroms.fna
```
Finally, remove these reverse complemented chromosomes from the full FASTA and concatenate with the corrected FASTA.
```
cat californicus_chromID.fna | seqkit grep -n -p "chr1" -p "chr2" -p "chr5" -p "chr8" -p "chr9" -p "chr13" -p "chr16" -p "chr17" -p "chr22" -p "chrX" -v > californicus_correctChroms.fna
```
```
seqkit concat californicus_correctChroms.fna californicus_reversedChroms.fna --full > californicus_alignReady.fna
```

Following these corrections to the FASTA files, we re-ran the pairwise alignments using the same nucmer paramters described above. We then used these revised pairwise alignments as input to SyRI.
```
syri -c revised_maniculatus_{species}.coords -d filtered_revised_maniculatus_{species}.delta -r maniculatus_chromID.fna -q {species}_alignReady.fna --nc 8 --prefix maniculatus_{species}. --nosnp
```
And visualized the results of this synteny analysis with plotsr.


