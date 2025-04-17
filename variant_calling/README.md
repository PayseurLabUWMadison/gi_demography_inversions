# Variant calling and filtering
This Markdown file describes the steps taken to call and filter (based on site-level annotations) single nucleotide variants identified in whole genome sequences of island and mainland *P. maniculatus*. 

## Methodological details
Sampled mice span two islands (Saturna Island and Pender Island) and one mainland location (Maple Ridge). Since these sampling locations comprise distinct populations, joint calling of variants was conducted separately for each cohort (e.g., Saturna, Pender, and Maple Ridge). Unfortunately, this practice generates ambiguity about the state of variants that are present in one cohort but missing in others (i.e., we don't know whether they are monomorphic REF or whether there was insufficient information to make a confident call). To address this, we "rescue" genotypes at private variants by forcing GATK to output monomorphic REF records. 

## Software information
All of the software used for this component are described in the [packages](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/packages) directory.

## Code
This section details the exact commands and parameters used for variant calling and basic site-level filtering. The following wildcards were used in constructing file names for the example commands:
- Sample name: `{sample}`
- Population name: `{pop}`
- Directory: `{dir}`
- Chromosome: `{chrom}`

#### Variant calling
Before running GATK on the alignments, the merged and duplicate-marked BAM files need to be indexed with samtools. At this stage, each individual is represented by a single BAM file.
```
samtools index {sample}_sorted_merged_rmdups.bam
```

To compute genotype likelihoods at all sites and create GVCFs for each individual, we ran HaplotypeCaller on the GVCF Reference Confidence Model mode. This can be done on a per-individual basis, but we found it necessary to additionally split the input BAM for each individual into separate chromosomes, then make the GVCFs on a per-individual, per-chromosome basis.
```
gatk HaplotypeCaller --tmp-dir {dir} -R GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz -I {sample}_sorted_merged_rmdups.bam -L {chrom} -ERC GVCF -O {sample}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz
```
This produces a single GVCF for each individual, for each of their chromosomes.

In order to joint-call variants across individuals from a given population, the per-sample GVCFs need to be combined into a single GVCF that contains all individuals in the population. We did this on a per-chromosome basis for each population using CombineGVCFs.
```
gatk CombineGVCFs --tmp-dir {dir} -R GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz -V {sampleX}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -V {sampleY}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -V {sampleZ}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -L {chrom} -O {pop}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz
```
>Note: The number of input GVCFs specified with each `-V` flag will depend on the number of samples in that population. 
The output of this command is a single combined GVCFs for each chromosome, for each population.

We then performed joint-calling of variants within each population by running GenotypeGVCFs on a per-chromosome basis.
```
gatk GenotypeGVCFs --tmp-dir {dir} -R GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz -V {pop}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -L {chrom} -O {pop}_chr_{chrom}_sorted_merged_rmdups.vcf.gz
```

#### Variant filtering
To obtain a high-confidence variant callset, we applied a collection of site-level annotation filters using VariantFiltration on a per-population, per-chromosome basis.
```
gatk VariantFiltration --tmp-dir {dir} -V {pop}_chr_{chrom}_sorted_merged_rmdups.vcf.gz -L {chrom} -filter "QD < 5.0" --filter-name "QD5" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O filtered_{pop}_chr_{chrom}_sorted_merged_rmdups.vcf.gz
```

We then used SelectVariants to remove sites that did not pass the above filtering thresholds and to restrict our downstream callset to biallelic, SNP-type variants.
```
gatk SelectVariants --tmp-dir {dir} -V filtered_{pop}_chr_{chrom}_sorted_merged_rmdups.vcf.gz -L {chrom} --restrict-alleles-to BIALLELIC --select-type-to-include SNP --exclude-filtered true -O snpOnly_filteredPASS_{pop}_chr_{chrom}_sorted_merged_rmdups.vcf.gz
```

#### Rescuing genotypes at "private" variants
Each population-specific VCF will contain sites that are polymorphic only within the observed sample. This means that variants "private" to one population will not appear in the VCFs of other populations. This creates ambiguity as to whether such missing records are monomorphic versus uncalled. The steps below describe how we "rescued" the genotype calls at these "private" variants in order to properly merge the callsets across all three populations.

First we need to identify sites that are unique to each population. This can be accomplished with bcftools `isec`. For each focal population, we start by identifying records present in each of the *other* two populations, but *not* in our focal population.

This command prints records present in pop1 but *not* our focal population:
```
bcftools isec -C -Oz -o {focal}_unique_to_{pop1}_{chrom}.txt snpOnly_filteredPASS_{pop1}_chr_{chrom}_sorted_merged_rmdups.vcf.gz snpOnly_filteredPASS_{focal}_chr_{chrom}_sorted_merged_rmdups.vcf.gz
```

This command prints records present in pop2 but *not* our focal population:
```
bcftools isec -C -Oz -o {focal}_unique_to_{pop2}_{chrom}.txt snpOnly_filteredPASS_{pop2}_chr_{chrom}_sorted_merged_rmdups.vcf.gz snpOnly_filteredPASS_{focal}_chr_{chrom}_sorted_merged_rmdups.vcf.gz
```

We then merge these lists of private SNPs together:
```
cat {focal}_unique_to_{pop1}_{chrom}.txt {focal}_unique_to_{pop2}_{chrom}.txt | cut -f1,2 | sort -k1,1 -k2,2n | uniq > {focal}_to_force_call_{chrom}.txt
```
And use this combined list to create a [GATK-compatible interval file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists):
```
awk '{print $1":"$2"-"$2}' {focal}_to_force_call_{chrom}.txt > {focal}_to_force_call_{chrom}.list
```

We create this `{focal}_to_force_call_{chrom}.list` for *each* of the three populations and for every chromosome.

We can then pass these lists of variants that are "missing" from our focal population callset to GATK and use the `-all-sites` option to output records for these invariant sites:
```
gatk GenotypeGVCFs -R GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz -V {pop}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -L {pop}_to_force_call_{chrom}.list -ip 100 -all-sites -stand-call-conf 0.0 -O {chrom}_{pop}_force_calls.vcf.gz
```

We still want to filter these "rescued" records to ensure that we have confidence in the genotypes being called. However, in contrast to the filtering performed above, we omit any QUAL score-based filters since, by nature of these sites being invariant, they will yield poor QUAL scores:
```
gatk VariantFiltration -V {chrom}_{pop}_force_calls.vcf.gz -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O filtered_{chrom}_{pop}_force_calls.vcf.gz
```

Then we exclude those that failed the above filters and restrict our monomorphic callset to invariant SNPs using SelectVariants:
```
variant_calling/packages/gatk-4.2.0.0/gatk SelectVariants -V filtered_{chrom}_{pop}_force_calls.vcf.gz --select-type-to-include NO_VARIATION --exclude-filtered true -L {pop}_to_force_call_{chrom}.list -O monomorphicOnly_filteredPASS_{chrom}_{pop}_force_calls.vcf.gz
```

#### Merging callsets
Once these invariant calls have been "rescued", we can combine them across chromosomes and with the original variant callsets produced for each population.
```
gatk MergeVcfs I={pop}_polymorphic.list I={pop}_monomorphic.list O=allsites_filteredPASS_{pop}_sorted_merged_rmdups.vcf.gz
```
Here, the `{pop}_polymorphic.list` and `{pop}_monomorphic.list` represent lists of all the chromosome-split VCFs for the polymorphic and monomorphic calls.

After constructing these full genome, monomorphic+polymorphic callsets for each population, we can combine them to create a multi-population callset using bcftools `merge`. 
```
bcftools merge -m all -Oz -o allsites_filteredPASS_gulf_islands_sorted_merged_rmdups.vcf.gz allsites_filteredPASS_saturna_sorted_merged_rmdups.vcf.gz allsites_filteredPASS_pender_sorted_merged_rmdups.vcf.gz allsites_filteredPASS_maple_ridge_sorted_merged_rmdups.vcf.gz
``` 

Finally, to remove any "rescued" sites that are monomorphic in *all* populations, we can use GATK SelectVariants to restrict the combined callset to only include sites that are biallelic and polymorphic across the three populations.
```
gatk SelectVariants -V allsites_filteredPASS_gulf_islands_sorted_merged_rmdups.vcf.gz --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O snpOnly_filteredPASS_gulf_islands_sorted_merged_rmdups.vcf.gz
```
























