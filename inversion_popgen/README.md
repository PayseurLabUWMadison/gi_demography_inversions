# Population genetic analysis of segregating inversions
This Markdown file describes the steps taken to compute population genetic summary statistics within and between arrangement classes. 

## Methodological details
In order to compare patterns of variation on inverted and standard haplotypes, we leveraged the subset of individuals that are either homozygous for the inverted haplotype or homozygous for the standard haplotype. Supplementary Table S3 tabulates these assignments at each locus in each individual. 

For each locus and each population, we used scikit-allel to compute nucleotide diversity within each arrangement class (i.e., inverted vs standard) and divergence between arrangements. These analyses were performed using the single-population callsets. We also computed divergence between like arrangements sampled from different populations using the multi-population callset.

## Software information
All of the software used for this component are described in the [packages](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/packages) directory.

## Code
This section details the commands and parameters used for population genetic analyses of segregating inversions. Relevant files are provided in this directory. The following wildcards were used in constructing file names for the example commands:
- Population name: `{pop}`
- Inversion locus: `{inv}`

Prior to summary statistic calculation, we constructed separate VCF files for each segregating inversion locus (3.0, 7.2, 14.0, 15.0, 21.0, and 22.0) using GATK's SelectVariants tool. The breakpoints for each locus are provided in the `inversion_breakpoints.txt` file, which is organized as `[chrom] [start] [stop] [locus]`. In addition to including only SNPs falling within the breakpoints of one inversion, the VCFs used as input in the commands below should exclude close relatives. 

To perform the within-population calculations of diversity and divergence, we used the `homozygotes_analysis.py` script which can be run with the following command.
```
python homozygotes_analysis.py --vcf {inv}_{pop}.vcf.gz --locus {inv} --pop {pop} --assignments genotype_assignments.txt --out {inv}_{pop}
```
This script takes in a single-population (`{pop}`), hard-filtered VCF file (see the [variant_calling](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/variant_calling) directory for details) which has been pruned for close relatives and is restricted to SNPs falling within the breakpoints of the defined inversion (`{inv}`). In addition, this script uses the genotype assignments present in the `genotype_assignments.txt` file. 

It outputs a `.pi` file formatted as `[subtype] [pi] [start] [stop] [n_bases] [counts] [sample_size] [pop] [locus]` and a `.dxy` file formatted as `[dxy] [start] [stop] [counts] [sample_size_inv] [sample_size_std] [pop] [locus]`.

To perform the between-population calculations of divergence between like arrangements, we used the `homozygotes_analysis_between_pops.py` script which can be run as follows.
```
python homozygotes_analysis_between_pops.py --vcf {inv}_gulf_islands.vcf.gz --locus {inv} --assignments genotype_assignments.txt --out {inv}_gulf_islands
```
This script takes in the combined multi-population VCF file (see the [variant_calling](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/variant_calling) directory for details) which has been pruned for close relatives and is restricted to SNPs falling within the breakpoints of the defined inversion (`{inv}`). This script also uses the genotype assignments present in the `genotype_assignments.txt` file. 

It outputs a `.dxy` file for each pair of populations that is present in the multi-population VCF. It is formatted as `[dxy] [start] [stop] [counts] [pop] [locus] [subtype]` where `[pop]` represents the population pair being compared. 



