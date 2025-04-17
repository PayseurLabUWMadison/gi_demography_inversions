#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: computes summary statistics within and between arrangement classes for a given population
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python homozygotes_analysis.py --vcf [in.vcf.gz] --locus [inv] --pop [saturna, pender, or maple_ridge] --assignments [genotype assignments] --out [out]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import allel
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	#-----------------
	# Set variables 
	#-----------------

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
	# Conduct PCA on haplotypes using only SNPs within "compatible" windows
	#------------------------------------------------------------------------

	print("--------------------------------------------------------------")
	print("Defining genomic subpopulations...")
	print("--------------------------------------------------------------")

	# Create the GenotypeArray and positions vector for this locus
	ga = allel.GenotypeArray(geno_vec)

	print('Identifying subpopulations...')

	print('Homozygous inversion subpopulation:')

	# Create a list of haplotype names that were assigned as inversions or 'left' clusters
	inv_genos = assignments_sub.loc[(assignments_sub['Genotype'] == 'inv/inv'), 'Sample'].to_list()
	print(inv_genos)

	# Create a boolean list that tells us, for each element in the sample_vec, whether the individual belongs to the inversion suboppulation
	inv_bool = np.isin(sample_vec, inv_genos)

	print('Homozygous standard subpopulation:')

	# Create a list of haplotype names that were assigned as standard or 'right' clusters
	std_genos = assignments_sub.loc[(assignments_sub['Genotype'] == 'std/std'), 'Sample'].to_list()
	print(inv_genos)

	# Create a boolean list that tells us, for each element in the sample_vec, whether the individual belongs to the standard suboppulation
	std_bool = np.isin(sample_vec, std_genos)

	print('Subsetting genotypes by inversion status...')

	# Subset the GenotypeArray to get the inversion homozygotes
	ga_inv = ga.subset(sel0=None, sel1=inv_bool)
	# Keep track of sample size
	n_inv = ga_inv.shape[1]
	print(n_inv)

	# Subset the GenotypeArray to get the standard homozygotes
	ga_std = ga.subset(sel0=None, sel1=std_bool)
	# Keep track of sample size
	n_std = ga_std.shape[1]
	print(n_std)

	print("--------------------------------------------------------------")
	print("Beginning nucleotide diversity analysis...")
	print("--------------------------------------------------------------")

	# Create dfs to store the nucleotide diversity results for each subtype
	inv_pi = pi('inv', ga_inv, pos, n_inv)
	std_pi = pi('std', ga_std, pos, n_std)
	# Combine subtype dfs into one
	pi_df = pd.concat([inv_pi, std_pi])
	# Add columns with additional identifying information
	pi_df['pop'] = pop
	pi_df['locus'] = locus

	print('Writing output...')
	# Write the output to a tab-delimited CSV
	pi_df.to_csv((out + '.pi'), sep="\t", index=False)
	print('Done.')

	print("--------------------------------------------------------------")
	print("Beginning Dxy analysis...")
	print("--------------------------------------------------------------")

	# Create df to store the results
	dxy_df = dxy(ga_inv, ga_std, pos, n_inv, n_std)
	# Add columns with additional identifying information
	dxy_df['pop'] = pop
	dxy_df['locus'] = locus

	print('Writing output...')
	# Write the output to a tab-delimited CSV
	dxy_df.to_csv((out + '.dxy'), sep="\t", index=False)
	print('Done.')

	print(dxy_df)

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
	parser.add_argument('--locus', required=True, choices=['inv3.0', 'inv7.2', 'inv14.0', 'inv15.0', 'inv21.0', 'inv22.0'])
	# Passes the pop name
	parser.add_argument('--pop', required=True, choices=['saturna', 'pender', 'maple_ridge'])
	# Passes the df containing genotype assignments
	parser.add_argument('--assignments', required=True)
	# Passes the prefix to use for output
	parser.add_argument('--out', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

def pi(subtype, ga, pos, n):
	""" 
	Compute windowed pi
	"""
	# Make an allele counts array from the GenotypeArray object
	ac = ga.count_alleles()
	
	# Make sure there is at least one variant site in this chromosome
	if(ga.n_variants > 1):

		# Compute windowed nucleotide diversity, creating relevant df columns
		pi,windows,n_bases,counts = allel.windowed_diversity(pos, ac, size=5000, fill="NA")

		# Extract start and stop positions from windows output
		start = np.split(windows, 2, axis=1)[0].flatten()
		stop = np.split(windows, 2, axis=1)[1].flatten()

		# Create pandas df to hold the pi results
		pi_df = pd.DataFrame({'subtype': subtype, 'pi': pi, 'start': start, 'stop': stop, 'n_bases': n_bases, 'counts': counts})
		pi_df['sample_size'] = n

		return(pi_df)

	else:
		print("Not enough variants")

def dxy(ga_inv, ga_std, pos, n_inv, n_std):
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
	dxy_df['sample_size_inv'] = n_inv
	dxy_df['sample_size_std'] = n_std
	
	return(dxy_df)

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()












