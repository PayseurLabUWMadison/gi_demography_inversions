#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: Simulate sequence data under the best-fit demographic model
# Approach: This script is meant to be run in parallel using CHTC. Towards this end, each run generates a single replicate (with chromosome number controlled by the --rep argument)
# Output: A formatted VCF file and a corresponding popfile
# Sampling: Number of diploid individuals can be controlled with the nsamp options to match observed sample sizes
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python sim_data.py --model [demographic model] --pop1_id [name of pop1] --pop1_nsamp [number of diploid individuals in pop1] --pop2_id [name of pop2] --pop2_nsamp [number of diploid individuals in pop2] --pop3_id [name of pop3] --pop3_nsamp [number of diploid individuals in pop3] --nsites [length of genomic element] --rep [index for given replicate] --recomb [per-bp, per-gen recombination rate] --mut [per-bp, per-gen recombination rate]
# Example run:
# python sim_data.py --model best_fit --pop1_id Saturna --pop1_nsamp 20 --pop2_id Pender --pop2_nsamp 10 --pop3_id MapleRidge --pop3_nsamp 17 --nsites 1000000 --rep 0 --recomb 5e-9 --mut 5e-9 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import msprime
import tskit
import pandas as pd
import numpy as np
import argparse

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	#-----------------
	# Set variables 
	#-----------------

	print("Reading in the data and intializing variables...")

	# Parse command-line arguments
	args = parse_args()

	# Set the demographic model for which we are simulating data
	model = args.model
	# Specify population names
	pop1_id = args.pop1_id
	pop2_id = args.pop2_id
	pop3_id = args.pop3_id
	# Specify diploid sample sizes for each population
	pop1_nsamp = int(args.pop1_nsamp)
	pop2_nsamp = int(args.pop2_nsamp)
	pop3_nsamp = int(args.pop3_nsamp)
	# Set the length of the independent genomic segments to simulate
	nsites = int(args.nsites)
	# Determine which chromosome/element number is to be generated
	rep = int(args.rep)
	# Set the recombination rate
	recomb_rate = float(args.recomb)
	# Set the mutation rate
	mut_rate = float(args.mut)
	# These are the names of the output files
	vcf = model + "_mut" + str(args.mut) + "_rec" + str(args.recomb) + "_rep" + str(args.rep) + ".vcf"
	popfile = model + "_mut" + str(args.mut) + "_rec" + str(args.recomb) + ".popfile"

	print("Done.")

	#---------------------------------
	# Construct the demographic model
	#---------------------------------

	print("Constructing demographic model...")

	# Construct the demographic model requested by user

	# Best-fit 3D model for Saturna, Pender, and Maple Ridge
	if (model == "best_fit"):
		demography = best_fit(pop1_id, pop2_id, pop3_id)
	else:
		raise NameError('Model not defined!')

	print("Done.")

	#------------------------------------------------------------------------
	# Simulate the tree sequences, add mutations, and convert to a dataframe
	#------------------------------------------------------------------------

	print("Simulating sequence data...")
	# Use msprime to generate simulated data under above demographic model; output as VCF-like pandas dataframe
	df = sim_data(demography, pop1_id, pop2_id, pop3_id, pop1_nsamp, pop2_nsamp, pop3_nsamp, nsites, rep, recomb_rate, mut_rate)
	print("Done.")

	#-------------------------------
	# Writing VCF-like output
	#-------------------------------

	# Need a properly formatted header line
	print("Creating metadata header...")
	header = get_header(rep, nsites)
	print("Done.")

	# Write full output to tab-delimited CSV file
	print("Writing VCF output...")
	with open(vcf, 'w') as outfile:
		outfile.write(header + "\n")
		df.to_csv(outfile, sep="\t", index=False)
	print("Done.")

	#-------------------------------
	# Creating popfile output
	#-------------------------------

	# Assuming simulated samples are ordered by population index, create a corresponding popfile
	print("Creating popfile...")
	# Get sample fields
	samples = df.columns.tolist()[9:]
	# Repeat the pop name appropriate number of times
	pops = np.repeat([pop1_id, pop2_id, pop3_id], [pop1_nsamp, pop2_nsamp, pop3_nsamp])
	# Create df for popfile
	pf = pd.DataFrame(list(zip(samples, pops)))
	pf.to_csv(popfile, sep="\t", index=False, header=False)
	print("Done.")

#---------------------
# DEFINE FUNCTIONS
#---------------------

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# Command-line variables

	# Determines which demographic model should be used
	parser.add_argument('--model', required=True)
	# Sets population names
	parser.add_argument('--pop1_id', required=True)
	parser.add_argument('--pop2_id', required=True)
	parser.add_argument('--pop3_id', required=True)
	# Determines how many diploid samples to simulate for each population
	parser.add_argument('--pop1_nsamp', required=True)
	parser.add_argument('--pop2_nsamp', required=True)
	parser.add_argument('--pop3_nsamp', required=True)	
	# Determines the length of each independent genomic segment
	parser.add_argument('--nsites', required=True)
	# Determines the ID of the chromosome/genomic element to simulate
	parser.add_argument('--rep', required=True)
	# Determines the per-bp per-gen recombination rate
	parser.add_argument('--recomb', required=True)
	# Determines the per-bp per-gen mutation rate
	parser.add_argument('--mut', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

#-----------------------------------------
# Define demographic models
#-----------------------------------------

def best_fit(pop1_id, pop2_id, pop3_id):
	"""
	Best-fit model obtained for Saturna (pop1), Pender (pop2), and Maple Ridge (pop3) 3D SFS.
	Specified model involves split of island-anc and mainland-anc from a common ancestor with symmetric migration between them. Extant island and mainland populations derive from these two ancestors and there is symmetric migration between the two islands following their split. 
	
	return: an msprime demography object
	"""
	demography = msprime.Demography()

	# Initializing population states

	# Initialize contmeporary Saturna population
	demography.add_population(name=pop1_id, initial_size=18607)
	# Initialize contmeporary Pender population
	demography.add_population(name=pop2_id, initial_size=14984)
	# Initialize contemporary Maple Ridge population
	demography.add_population(name=pop3_id, initial_size=69313)
	# Initialize epoch 1 populations
	demography.add_population(name="Mainland_Ancestor", initial_size=579644)
	demography.add_population(name="Island_Ancestor", initial_size=255958)
	# Initialize ancestral population
	demography.add_population(name="Ancestral", initial_size=129563)

	# Model contemporary migration between Pender and Saturna
	demography.set_symmetric_migration_rate(populations=[pop1_id, pop2_id], rate=1.63e-7)

	# Epoch 2 events

	# Model the split between Saturna and Pender from the island ancestor
	demography.add_population_split(time=5090, derived=[pop1_id, pop2_id], ancestral="Island_Ancestor")
	# Model the split of Maple Ridge from the mainland ancestor
	demography.add_population_split(time=5090, derived=[pop3_id], ancestral="Mainland_Ancestor")

	# Epoch 1 events

	# Model migration between the island and mainland ancestors
	demography.add_symmetric_migration_rate_change(time=5090, populations=["Mainland_Ancestor", "Island_Ancestor"], rate=6.75e-6)
	# Model the split between island and mainland ancestors from the ancestral population
	demography.add_population_split(time=2565592, derived=["Mainland_Ancestor", "Island_Ancestor"], ancestral="Ancestral")

	print(demography.debug())

	return demography

#-----------------------------------------
# Simulation functions
#-----------------------------------------	

def sim_data(demography, pop1_id, pop2_id, pop3_id, pop1_nsamp, pop2_nsamp, pop3_nsamp, nsites, rep, recomb_rate, mut_rate):
	"""
	This function uses msprime to simulate the tree sequences, add mutations, and convert to a dataframe

	demography: this is the msprime demography object that was created
	nsamples: this is the number of diploid samples to simulates data for
	nsites: this is the length of each genomic segment to simulate
	nreps: this is the number of independent genomic segments to simulate
	recomb_rate: this is the per-bp per-gen recombination rate experienced by each genomic segment
	mut_rate: this is the per-bp per-gen mutation rate experienced by each genomic segment

	return: a pandas df storing the error-free simulated data in a VCF-like format
	"""
	# Create a list of lists (where each element is a row of the VCF output)
	data = list()

	# For the given replicate, simulate the corresponding tree sequence objects and add mutations to each

	# Create the ts for the given replicate
	ts = msprime.sim_ancestry(samples={pop1_id: pop1_nsamp, pop2_id: pop2_nsamp, pop3_id: pop3_nsamp}, demography=demography, sequence_length=nsites, discrete_genome=True, recombination_rate=recomb_rate, ploidy=2)
	# Then add mutations on this ts object
	mts = msprime.sim_mutations(ts, rate=mut_rate)
	# Store the current segment in VCF format
	vcf = mts.as_vcf(contig_id=rep)
	# Turn the VCF output into a list of lines
	lines = vcf.split("\n")

	# Loop through the list of lines, breaking each into elements
	for line in lines:
		# Each element is a "cell" in the VCF
		elements = line.split("\t")

		# Only append rows that contain either the records or the column titles (we don't want the rest of the header info for now)
		# I also want to exclude rows/sites that have more than one derived allele (those will break the script)
		if (len(elements) > 1):

			# Count the number of ALT alleles (in case there are multi-allelic sites)
			if(len(elements[4].split(",")) < 2):
				data.append(elements)

	# Turn the data object that we generated for all replicates into a pandas df, using the first record as the column titles
	df = pd.DataFrame(data[1:], columns=data[0])

	return df

def get_header(rep, nsites):
	"""
	This function generates an appropriate VCF-like header based on input parameters
	"""

	# Create a list to store the bookended and contig-specific header lines
	header_info = list()
	# Add the first several header elements to this list
	header_info.append("##fileformat=VCFv4.2")
	header_info.append("##source=tskit 0.5.4")
	header_info.append("##FILTER=<ID=PASS,Description=\"All filters passed\">")

	# For the given chromosome, create a contig header line that specifies the length (given by nsites parameter)
	header_info.append("##contig=<ID=" + str(rep) + ",length=" + str(nsites) + ">")

	# Add the last header element(s) to this list
	header_info.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")

	# Join the header lines
	header = "\n".join(header_info)

	return header

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()