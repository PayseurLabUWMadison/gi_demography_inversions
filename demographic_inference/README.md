# Demographic inference for Gulf Islands populations
This Markdown file describes the steps taken to fit two- and three-population demographic models to the joint SFS of island and mainland populations.

## Methodological details
To obtain a high-confidence, putatively neutral SNP callset, we conducted individual-level filtering on coverage depth (DP) and genotype quality score (GQ). We also excluded sites in violation of Hardy-Weinberg expected genotype counts by filtering on GATK's ExcessHet annotation within each population cohort. Details about filtering thresholds are described in the Materials and Methods section of the paper. 

We additionally leveraged existin UCSC Genome Browser tracks to exclude SNPs falling within annotated repetitive elements, within 1 kbp of annotated genes, and within the boundaries of known inversion polymorphisms in *P. maniculatus*. Additional details about our masking approach can be found in the Materials and Methods section of the paper. 

Following these filtering and masking procedures, we excluded sites containing missing genotype calls for one or more individuals and we restricted downstream demographic inference to the unrelated subset of individuals identified in each population. We also restricted demographic inference to autosomal SNPs (i.e., we excluded X chromosome and mitochondrial SNPs). 

Altogether, this resulted in a single, multi-population autosomal callset (see [variant_calling](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/variant_calling) for details on multi-population callset construction) that is free from close relatives, missing data, low-confidence genotype calls, and non-neutral sites. In the following code, we call this inference-ready VCF file `gulf_islands_inference_ready.vcf.gz`. 

## Software information
All of the software used for this component are described in the [packages](https://github.com/PayseurLabUWMadison/gi_demography_inversions/tree/main/packages) directory.

## Code
This section details the commands and parameters used for demographic inference.

Before fitting demographic models, we constructed Moments frequency spectrum objects for both 2D (across all population pairs) and 2D jSFS. The script [`create_fs.py`](https://github.com/PayseurLabUWMadison/gi_demography_inversions/blob/main/demographic_inference/create_fs.py) can read in a multi-population VCF file and create either the marginal, joint (2D), or joint (3D) frequency spectrum for the samples/populations included in the specified population manifest. The script outputs the frequency spectrum object as a `.fs` file.
```
python create_fs.py --vcf gulf_islands_inference_ready.vcf.gz --prefix {prefix} --popfile {popfile} --dimension {1, 2, or 3}
```
Where the `--popfile` is the corresponding sample manifest in the format `[individual] [pop]`.

Tested 2D and 3D demographic models are specified in the [`fit_models.py`](https://github.com/PayseurLabUWMadison/gi_demography_inversions/blob/main/demographic_inference/fit_models.py) script. This python script `fit_models.py` takes in the `.fs` objects generated above, the name of the demographic model we wish to fit, the mutation rate and effective sequence length to use for parameter conversion, and a collection of optimization parameters. Using these optimization parameters, it then conducts a series of searches for the MLE of each parameter defined by the demographic model.
```
python fit_models.py --prefix {prefix} --fs {frequency spectrum}.fs --model {model} --opt_num {W} --fold_num {X} --outer_rep_num {Y} --inner_rep_num {Z} --mut {mutation rate} --L {effective sequence length}
```
For each of the optimizations defined by the input parameters, the script outputs the following:
```
[Model] [Param 1] [Param 2] ... [Theta] [Likelihood] [Popt]
```
Where [Popt] represents the un-converted parameter estimates to be more easily read in by downstream model evaluation scripts.

To assess model fit, we inspected residual plots, which compared differences between the jSFS expected under the given model to the observed jSFS for each bin of the frequency spectrum. The python script [`eval_models.py`](https://github.com/PayseurLabUWMadison/gi_demography_inversions/blob/main/demographic_inference/eval_models.py) uses the the best-fit, unscaled, parameter estimates obtained with the [`fit_models.py`](https://github.com/PayseurLabUWMadison/gi_demography_inversions/blob/main/demographic_inference/fit_models.py) script and performs this comparison, outputting a `.png` summarizing the results.
```
python eval_models.py --prefix {prefix} --fs {frequency spectrum}.fs --model {model} --popt {parameter estimates} --dimension {1, 2, or 3}
```

To estimate parameter uncertanties for our best-fitting three-population demographic model, we used the GIM method (Coffman et al. 2016) as an alternative to conventional bootstrapping. The python script [`est_uncert.py`](https://github.com/PayseurLabUWMadison/gi_demography_inversions/blob/main/demographic_inference/est_uncert.py) outputs standard deviations for each parameter estimate. 
```
python est_uncert.py --vcf gulf_islands_inference_ready.vcf.gz --popfile {popfile} --popt {parameter estimates}
```

The python script [`sim_data.py`](https://github.com/PayseurLabUWMadison/gi_demography_inversions/blob/main/demographic_inference/sim_data.py) uses msprime to simulate sequence data in VCF format under the best-fitting three-population demographic model.
```
python sim_data.py --model {model} --pop1_id {pop1 name} --pop1_nsamp {number of diploid individuals in pop1} --pop2_id {pop2 name} --pop2_nsamp {number of diploid individuals in pop2} --pop3_id {pop3 name} --pop3_nsamp {number of diploid individuals in pop3} --nsites {length of genomic element} --rep {index for given replicate} --recomb {per-bp, per-gen recombination rate} --mut {per-bp, per-gen recombination rate}
```

















