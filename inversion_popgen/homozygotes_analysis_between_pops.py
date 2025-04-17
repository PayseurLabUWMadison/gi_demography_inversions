#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: computes summary statistics within inv or std arrangement *between* populations
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python homozygotes_analysis_between_pops.py --vcf [in.vcf.gz] --locus [inv] --assignments [genotype assignments] --out [out]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import allel
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import itertools

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	print("--------------------------------------------------------------")
	print("Analyzing inputs...")
	print("--------------------------------------------------------------")

	print("Intializing variables...")

	# Parse command-line arguments
	args = parse_args()
	# VCF to analyze
	invcf = args.vcf
	# Output file prefix
	out = args.out
	# Locus name
	locus = args.locus
	# Pop name
	pop = args.pop
	# Assignments df
	assignments = pd.read_csv(args.assignments, delimiter=' ', header='infer')

	# Print out update message
	print('Analyzing data for locus ' + str(locus) + ' from population ' + str(pop) + '.')
	print('Subsetting assignments to reflect specified population and locus...')	

	# Collect a list of heterozygous individuals to exclude
	to_exclude = np.unique(assignments.loc[(assignments['Population'] == pop) & (assignments['Locus'] == locus) & (assignments['Genotype'] == 'inv/std'), 'Sample'].to_list())

	# Extract homozygous individuals pertaining to the specified population and inversion locus
	assignments_sub = assignments.loc[(assignments['Population'] == pop) & (assignments['Locus'] == locus) & (~assignments['Sample'].isin(to_exclude)), :]

	# Update message
	if len(to_exclude) > 0:
		print('Excluding heterozygous individuals:')
		print(to_exclude)
	else:
		print('There are no heterozygous individuals to exclude.')

	# Create a list of samples to include for the purposes of VCF processing
	to_include = np.unique(assignments_sub['Sample'])

	print('Keeping individuals: ')
	print(to_include)

	print("Reading in the VCF file...")

	# Read in the VCF file
	vcf = allel.read_vcf(invcf, samples = to_include, fields='*')
	# Get a list of samples in the VCF
	sample_vec = vcf['samples']
	# Extract the genotype field
	geno_vec = vcf['calldata/GT']
	# Get the vector of chromosomes (1st col in VCF)
	chrom_vec = vcf['variants/CHROM']
	# Get the vector of positions
	pos = vcf['variants/POS']
	# Get a list of the chromosomes represented in the VCF file
	chrom_ls = np.unique(chrom_vec)

	# Print summary message
	print("Found samples:")
	print(sample_vec)
	print("Found " + str(len(pos)) + " sites")
	print("Found chromosomes:")
	print(chrom_ls)

	#------------------------------------------------------------------------
	# Identify genomic subtypes and population subsets
	#------------------------------------------------------------------------

	print("--------------------------------------------------------------")
	print("Defining genomic subpopulations...")
	print("--------------------------------------------------------------")

	# Create the GenotypeArray and positions vector for this locus
	ga = allel.GenotypeArray(geno_vec)

	print('Identifying subpopulations...')

	# Store a list of populations for this locu
	pop_ls = np.unique(assignments_sub['Population'].to_list())

	print('Homozygous inversion individuals in each population:')

	# Initialize a dictionary where the population name is the key and the sample names are the values
	inv_genos = dict()
	inv_bool = dict()

	# Loop through and add the sample name array and the boolean array to the dict
	for pop in pop_ls:
		# Create a list of genotype names that were assigned as inversion homozygotes for this population
		inv_genos[pop] = assignments_sub.loc[(assignments_sub['Population'] == pop) & (assignments_sub['Genotype'] == 'inv/inv'), 'Sample'].to_list()
		print(inv_genos[pop])
		# Create a boolean list that tells us, for each element in the sample_vec, whether the genotype is found in this inversion suboppulation
		inv_bool[pop] = np.isin(sample_vec, inv_genos[pop])

	print('Homozygous standard individuals in each population:')

	# Initialize a dictionary where the population name is the key and the sample names are the values
	std_genos = dict()
	std_bool = dict()

	# Loop through and add the sample name array and the boolean array to the dict
	for pop in pop_ls:
		# Create a list of genotype names that were assigned as standard homozygotes for this population
		std_genos[pop] = assignments_sub.loc[(assignments_sub['Population'] == pop) & (assignments_sub['Genotype'] == 'std/std'), 'Sample'].to_list()
		print(std_genos[pop])
		# Create a boolean list that tells us, for each element in the sample_vec, whether the genotype is found in this standard suboppulation
		std_bool[pop] = np.isin(sample_vec, std_genos[pop])

	print('Subsetting genotypes by inversion status...')

	for pop1, pop2 in itertools.combinations(pop_ls, r=2):

		if (len(inv_genos[pop1]) > 0) and (len(std_genos[pop1]) > 0) and (len(inv_genos[pop2]) > 0) and (len(std_genos[pop2]) > 0):

			print('Comparing populations ' + pop1 + ' and ' + pop2 + '...')

			print('Measuring divergence between inverted subtypes...')
			
			# Subset the original GenotypeeArray to get the inversions
			pop1_ga_inv = ga.subset(sel0=None, sel1=inv_bool[pop1])
			pop2_ga_inv = ga.subset(sel0=None, sel1=inv_bool[pop2])

			print("Beginning Dxy analysis...")

			# Create df to store the results
			dxy_df = dxy(pop1_ga_inv, pop2_ga_inv, pos)
			# Add columns with additional identifying information
			dxy_df['pop'] = pop1 + "_" + pop2
			dxy_df['locus'] = locus
			dxy_df['subtype'] = 'inv'

			print('Writing output...')
			# Write the output to a tab-delimited CSV
			dxy_df.to_csv(("inv_" + pop1 + "_" + pop2 + "_" + out + '.dxy'), sep="\t", index=False)
			print('Done.')

			print(dxy_df)

			print('Measuring divergence between standard subtypes...')

			# Subset the original GenotypeArray to get the standards
			pop1_ga_std = ga.subset(sel0=None, sel1=std_bool[pop1])
			pop2_ga_std = ga.subset(sel0=None, sel1=std_bool[pop2])

			print("Beginning Dxy analysis...")

			# Create df to store the results
			dxy_df = dxy(pop1_ga_std, pop2_ga_std, pos)
			# Add columns with additional identifying information
			dxy_df['pop'] = pop1 + "_" + pop2
			dxy_df['locus'] = locus
			dxy_df['subtype'] = 'std'

			print('Writing output...')
			# Write the output to a tab-delimited CSV
			dxy_df.to_csv(("std_" + pop1 + "_" + pop2 + "_" + out + '.dxy'), sep="\t", index=False)
			print('Done.')

			print(dxy_df)
		
		else:
			print('Insufficient sample size for ' + pop1 + ' and ' + pop2 + '...')

#---------------------
# DEFINE FUNCTIONS
#---------------------

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# Passes VCF file to use (should be gzipped)
	parser.add_argument('--vcf', required=True)
	# Passes the locus name
	parser.add_argument('--locus', required=True, choices=['inv3.0', 'inv7.2', 'inv9.0', 'inv14.0', 'inv15.0', 'inv21.0', 'inv22.0'])
	# Passes the df containing genotype assignments
	parser.add_argument('--assignments', required=True)
	# Passes the prefix to use for output
	parser.add_argument('--out', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

def dxy(ga_inv, ga_std, pos):
	"""
	Compute Dxy in windows
	"""
	# Make an allele counts array for inv subtype
	ac_inv = ga_inv.count_alleles()
	# Make an allele counts array for std subtype
	ac_std = ga_std.count_alleles()

	# Compute dxy, creating relevant df columns
	dxy,windows,n_bases,counts = allel.windowed_divergence(pos, ac_inv, ac_std, size=5000, fill="NA")

	# Extract start and stop positions from windows output
	start = np.split(windows, 2, axis=1)[0].flatten()
	stop = np.split(windows, 2, axis=1)[1].flatten()

	# Create pandas df to hold the Fst results
	dxy_df = pd.DataFrame({'dxy': dxy, 'start': start, 'stop': stop, 'counts': counts})
	
	return(dxy_df)

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()












