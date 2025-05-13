#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: estimate parameter uncertainties for best-fit model
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python est_uncert.py --vcf [filtered VCF used for original inference] --popfile [population manifest] --popt [param1 param2 ... paramN]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import pandas as pd
from random import sample
import moments
import dill as pickle 
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import pylab
import glob

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

	# Variable for holding the name of the input VCF file
	vcf = args.vcf
	# Variable for holding input popfile
	popfile = args.popfile
	# Create a variable to hold the parameter estimates
	popt_est = args.popt
	# Create a variable to hold the parameter estimates as a list of floats
	popt = [float(p) for p in popt_est]

	print("Done.")

	#---------------------------------------------
	# Constructing bootstrapped frequency spectra
	#---------------------------------------------

	# Run function to create bootstrapped FS from input VCF
	all_boot, original_fs = create_boots(vcf, popfile)

	#---------------------------------
	# Compute parameter uncertainties
	#---------------------------------

	# Determine which demographic model we're estimating uncertainties for
	model_func = tri_split_two_epoch_mig_minus2_D

	print("Original parameter estimates:")
	print(popt)

	print("Original theta:")
	# Determine the sample sizes from the original FS
	ns = original_fs.sample_sizes
	# Calculate the best-fit model FS
	model = model_func(popt, ns)
	# Calculate the corresponding theta value
	theta = moments.Inference.optimal_sfs_scaling(model, original_fs)
	print(theta)

	print("Estimating parameter uncertainties...")

	# Loop through different possible step sizes to ensure uncertainty estimates are stable
	for eps in [0.01, 0.001, 0.0001]:

		print("#------------------------------------------------------------")
		print("Using step size = " + str(eps))

		# Function returns standard deviation of parameter values (including theta– listed last) along with the full GIM to use in propogating uncertainties
		uncert, GIM = moments.Godambe.GIM_uncert(model_func, all_boot, popt, original_fs, log=False, multinom=True, eps=eps, return_GIM=True)

		print("Standard deviations of parameter values (plus theta):")
		print(uncert)

		print("Full GIM matrix:")
		print(GIM)

		print("Inverse of GIM matrix:")
		print(np.linalg.inv(GIM))
		print("#------------------------------------------------------------")


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

	# Determines which VCF to use
	parser.add_argument('--vcf', required=True)
	# Determines which popfile to use
	parser.add_argument('--popfile', required=True)
	# These are the best-fit parameter estimates (excluding theta and the ll) for the given model
	parser.add_argument("--popt", default=[], nargs="+")

	# Parse and return arguments
	args = parser.parse_args()
	return args

def create_boots(vcf, popfile):
	"""
	This function takes in a VCF dataset and creates bootstrap replicates to use in GIM functions
	"""
	# Read in the provided popfile as a pandas df
	popf_df = pd.read_csv(popfile, sep='\t', header=None)

	# Identify all unique populations in the popfile
	pop_ls = np.unique(popf_df.iloc[:,1].tolist())
	print("Populations detected:")
	print(pop_ls)

	# Construct a list for the samples that correspond to each unique population
	samp_ls = [popf_df[popf_df[1] == p][0].tolist() for p in pop_ls]
	print("Corresponding samples:")
	print(samp_ls)

	# Construct a list for the sample sizes to use in the projection argument
	proj_ls = [2*len(s) for s in samp_ls]
	print("Corresponding sample sizes:")
	print(proj_ls)

	print("Constructing data dictionary from input VCF...")

	# Read in the input vcf and construct moments dictionary object
	dd = moments.Misc.make_data_dict_vcf(vcf, popfile)

	print("Generating bootstraps...")

	# Create folded bootstrap replicats
	# And now write the output to a directory
	boot_dir = popfile + "_uncert_boots"
	moments.Misc.bootstrap(dd, pop_ls, proj_ls, polarized = False, num_boots=100, save_dir=boot_dir)

	# Load the saved bootstrap replicates
	boots_fids = glob.glob(boot_dir+'/*')
	all_boot = [moments.Spectrum.from_file(fid) for fid in boots_fids]

	print("Bootstraps generated.")

	print("Generating full frequency spectrum...")

	original_fs = moments.Spectrum.from_data_dict(dd, pop_ls, projections = proj_ls, polarized = False)

	print("Full frequency spectrum generated.")

	return all_boot, original_fs

# Define the demographic function for the best-fit model
def tri_split_two_epoch_mig_minus2_D(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, NULL, m23E2
	Excludes: m12E2, m13E2
	"""
	# Ancestral dynamics #

	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m23E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)

	# Epoch 1 dynamics #

	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations with symmetric migration rate m12E1
	fs.integrate([nu1E1, nu2E1], TE1, m=np.array([[0, m12E1], [m12E1, 0]]))

	# Epoch 2 dynamics #

	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates 0 between MR-PI, 0 between MR-SI, and m23E2 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, 0, 0], [0, 0, m23E2], [0, m23E2, 0]]))

	# Return the frequency spectrum obtained from the model
	return fs

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()



