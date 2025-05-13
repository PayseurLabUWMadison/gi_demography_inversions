#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: the purpose of this script is to inspect the fit of demographic parameter estimates to the real data
# Output: a comparison of the empirical to model SFS and the corresponding residuals
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: python eval_models.py --prefix [prefix to use for output/population(s) name(s)] --fs [fs file] --model [model name] --popt [param1 param2 param3] --dimension [1, 2, or 3]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------
# Import required packages
#--------------------------
import dadi
import moments
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import re
import pylab

#---------------------------------------------
# Loading the data and initializing variables
#---------------------------------------------

# Parse command line input
parser = argparse.ArgumentParser()

# This is the name of the population
parser.add_argument("--prefix")
# This is the fs file name
parser.add_argument("--fs")
# This is the demographic model we wish to fit 
parser.add_argument("--model")
# These are the best-fit parameter estimates (excluding theta and the ll) for the given model
parser.add_argument("--popt", default=[], nargs="+")
# This is the number of dimensions of the input FS
parser.add_argument("--dimensions")

# Parse args
args = parser.parse_args()

# Create variable to store the population name
pop_name = args.prefix
# Read in the input frequency spectrum
fs = moments.Spectrum.from_file(args.fs)
# Create a variable to hold the model ID
model_id = args.model
# Create a variable to hold the parameter estimates
popt_est = args.popt
# Create a variable to hold the parameter estimates as a list of floats
popt = [float(p) for p in popt_est]
# Create a variable to hold the number of dimensions
dim = int(args.dimensions)

#----------------------------------
# Specifying the demographic model
#----------------------------------

#------------------------
# Two-population models:
#------------------------

# Define the demographic function for a simple population split
def split(params, ns):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1 size to the ancestral population size
	2. nu2 = the ratio of population 2 size to the ancestral population size
	3. T = the time since the population split
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1, nu2, T = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# After the splitting T generations ago, the populations take constant sizes nu1 and nu2
	fs.integrate([nu1, nu2], T)
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split demographic model
def instantiate_split():
	# Set upper and lower bounds for the list of parameter choices (in order nu1, nu2, T)
	upper_bound = [50,50,10]
	lower_bound = [1e-4,1e-4,0]
	# Initial guess at parameter values (nu1, nu2, T)
	p0 = [1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split demographic model
def params_split(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1, nu2, and T parameter values 
	nu1 = float(popt[0]) * float(nuA)
	nu2 = float(popt[1]) * float(nuA)
	T = float(popt[2]) * float(nuA) * 2
	# Return these computed parameter values
	return nuA,nu1,nu2,T,theta


# Define the demographic function for a population split followed by exponential size change for both population 1 and 2
def split_exp_exp(params, ns):
	"""
	Parameter values:
	1. nu1B = the ratio of population 1's bottleneck size to the ancestral population size
	2. nu2B = the ratio of population 2's bottleneck size to the ancestral population size
	3. nu1C = the ratio of population 1's current size to the ancestral population size
	4. nu2C = the ratio of population 2's current size to the ancestral population size
	5. T = the time since the population split
	Overview: the ancestral population splits into population 1 with size nu1B and population 2 with size nu2B T generations from the present.
		After the split, the populations both change expoentially to sizes nu1C and nu2C, respectively.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1B, nu2B, nu1C, nu2C, T = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth from size nu1B to nu1C over T generations
	nu1_func = lambda t: nu1B * (nu1C/nu1B) ** (t/T)
	# Create the same func for other pop (describes exponential growth from size nu2B to nu2C over T generations)
	nu2_func = lambda t: nu2B * (nu2C/nu2B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1_func(t), nu2_func(t)]
	# After the splitting T generations ago, the populations both change exponentially according to their functions
	fs.integrate(nu_func, T)
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then exp/exp change demographic model
def instantiate_split_exp_exp():
	# Set upper and lower bounds for the list of parameter choices (in order nu1B, nu2B, nu1C, nu2C, T)
	upper_bound = [50,50,50,50,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,0]
	# Initial guess at parameter values (nu1B, nu2B, nu1C, nu2C, T)
	p0 = [1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then exp/exp demographic model
def params_split_exp_exp(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1B, nu2B, nu1C, nu2C, and T parameter values 
	nu1B = float(popt[0]) * float(nuA)
	nu2B = float(popt[1]) * float(nuA)
	nu1C = float(popt[2]) * float(nuA)
	nu2C = float(popt[3]) * float(nuA)
	T = float(popt[4]) * float(nuA) * 2
	# Return these computed parameter values
	return nuA,nu1B,nu2B,nu1C,nu2C,T,theta


# Define the demographic function for a population split followed by exponential change in population 1 and no change in population 2
def split_exp_const(params, ns):
	"""
	Parameter values:
	1. nu1B = the ratio of population 1's bottleneck size to the ancestral population size
	2. nu1C = the ratio of population 1's current size to the ancestral population size
	3. nu2C = the ratio of population 2's current size to the ancestral population size
	4. T = the time since the population split
	Overview: the ancestral population splits into population 1 with size nu1B and population 2 with size nu2C T generations from the present.
		After the split, the population 1 undergoes exponential change to nu1C and population stays constant at nu2C.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1B, nu1C, nu2C, T = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth in pop 1 from size nu1B to nu1C over T generations
	nu1_func = lambda t: nu1B * (nu1C/nu1B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1_func(t), nu2C]
	# After the splitting T generations ago, pop 1 changes exponentially according to function while pop 2 stays the same
	fs.integrate(nu_func, T)
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then exp/const change demographic model
def instantiate_split_exp_const():
	# Set upper and lower bounds for the list of parameter choices (in order nu1B, nu1C, nu2C, T)
	upper_bound = [50,50,50,10]
	lower_bound = [1e-4,1e-4,1e-4,0]
	# Initial guess at parameter values (nu1B, nu1C, nu2C, T)
	p0 = [1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then exp/const demographic model
def params_split_exp_const(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1B, nu1C, nu2C, and T parameter values 
	nu1B = float(popt[0]) * float(nuA)
	nu1C = float(popt[1]) * float(nuA)
	nu2C = float(popt[2]) * float(nuA)
	T = float(popt[3]) * float(nuA) * 2
	# Return these computed parameter values
	return nuA,nu1B,nu1C,nu2C,T,theta


# Define the demographic function for a population split followed by constant size in population 1 and exponential change in population 2
def split_const_exp(params, ns):
	"""
	Parameter values:
	1. nu1C = the ratio of population 1's current size to the ancestral population size
	2. nu2B = the ratio of population 2's bottleneck size to the ancestral populations size
	3. nu2C = the ratio of population 2's current size to the ancestral population size
	4. T = the time since the population split
	Overview: the ancestral population splits into population 1 with size nu1C and population 2 with size nu2B T generations from the present.
		After the split, the population 1 stays constant at nu1C while population 2 undergoes exponential change to nu2C.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1C, nu2B, nu2C, T = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth in pop 2 from size nu2B to nu2C over T generations
	nu2_func = lambda t: nu2B * (nu2C/nu2B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1C, nu2_func(t)]
	# After the splitting T generations ago, pop 2 changes exponentially according to function while pop 1 stays the same
	fs.integrate(nu_func, T)
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then const/exp change demographic model
def instantiate_split_const_exp():
	# Set upper and lower bounds for the list of parameter choices (in order nu1C, nu2B, nu2C, T)
	upper_bound = [50,50,50,10]
	lower_bound = [1e-4,1e-4,1e-4,0]
	# Initial guess at parameter values (nu1C, nu2B, nu2C, T)
	p0 = [1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then const/exp demographic model
def params_split_const_exp(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1C, nu2B, nu2C, and T parameter values 
	nu1C = float(popt[0]) * float(nuA)
	nu2B = float(popt[1]) * float(nuA)
	nu2C = float(popt[2]) * float(nuA)
	T = float(popt[3]) * float(nuA) * 2
	# Return these computed parameter values
	return nuA,nu1C,nu2B,nu2C,T,theta


# Define the demographic function for a population split followed by a post-split size change
def split_two_epoch(params, ns):
	"""
	Parameter values:
	1. nu1E1 = the ratio of population 1 size to the ancestral population size during first epoch
	2. nu2E1 = the ratio of population 2 size to the ancestral population size during the first epoch
	3. nu1E2 = the ratio of population 1 size to the ancestral population size during second epoch
	4. nu2E2 = the ratio of population 2 size to the ancestral population size during the second epoch
	5. TE1 = the duration of the first epoch
	6. TE2 = the duration of the second epoch
	Overview: the ancestral population splits into population 1 with size nu1E1 and population 2 with size nu2E1. After TE1 generations, the populations enter a second epoch with sizes nu1E1 and nu2E2 that lasts for TE2 generations.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# After the first epoch, the populations take constant sizes nu1E2 and nu2E2 for TE2 generations
	fs.integrate([nu1E2, nu2E2], TE2)
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split_two_epoch demographic model
def instantiate_split_two_epoch():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2)
	upper_bound = [50,50,50,50,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,0,0]
	# Initial guess at parameter values (nu1, nu2, T)
	p0 = [1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split_two_epoch demographic model
def params_split_two_epoch(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	TE1 = float(popt[4]) * float(nuA) * 2
	TE2 = float(popt[5]) * float(nuA) * 2
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,TE1,TE2,theta



#--------------------------------------
# Two population models with migration:
#--------------------------------------

# Define the demographic function for a population split followed by constant sizes for both population 1 and 2 with symmetric migration
def split_mig(params, ns):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1's size to the ancestral population size
	2. nu2 = the ratio of population 2's size to the ancestral population size
	3. T = the time since the population split
	4. m = the symmetric migration rate
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.
		After the split, there is symmetric migration at rate m.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1, nu2, T, m = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# After the splitting T generations ago, the populations have constant size and experience symmetric migration
	fs.integrate([nu1, nu2], T, m=np.array([[0, m], [m, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split with mig demographic model
def instantiate_split_mig():
	# Set upper and lower bounds for the list of parameter choices (nu1, nu2, T, m)
	upper_bound = [50,50,10,10]
	lower_bound = [1e-4,1e-4,0,0]
	# Initial guess at parameter values (nu1, nu2, T, m)
	p0 = [1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split with mig demographic model
def params_split_mig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1, nu2, T, m parameter values 
	nu1 = float(popt[0]) * float(nuA)
	nu2 = float(popt[1]) * float(nuA)
	T = float(popt[2]) * float(nuA) * 2
	m = float(popt[3])/(2 * float(nuA)) 
	# Return these computed parameter values
	return nuA,nu1,nu2,T,m,theta


# Define the demographic function for a population split followed by exponential size change for both population 1 and 2 with symmetric migration
def split_exp_exp_mig(params, ns):
	"""
	Parameter values:
	1. nu1B = the ratio of population 1's bottleneck size to the ancestral population size
	2. nu2B = the ratio of population 2's bottleneck size to the ancestral population size
	3. nu1C = the ratio of population 1's current size to the ancestral population size
	4. nu2C = the ratio of population 2's current size to the ancestral population size
	5. T = the time since the population split
	5. m = the symmetric migration rate
	Overview: the ancestral population splits into population 1 with size nu1B and population 2 with size nu2B T generations from the present.
		After the split, the populations both change expoentially to sizes nu1C and nu2C, respectively, and there is symmetric migration at rate m.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1B, nu2B, nu1C, nu2C, T, m = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth from size nu1B to nu1C over T generations
	nu1_func = lambda t: nu1B * (nu1C/nu1B) ** (t/T)
	# Create the same func for other pop (describes exponential growth from size nu2B to nu2C over T generations)
	nu2_func = lambda t: nu2B * (nu2C/nu2B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1_func(t), nu2_func(t)]
	# After the splitting T generations ago, the populations both change exponentially according to their functions with symmetric migration
	fs.integrate(nu_func, T, m=np.array([[0, m], [m, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then exp/exp change with mig demographic model
def instantiate_split_exp_exp_mig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1B, nu2B, nu1C, nu2C, T, m)
	upper_bound = [50,50,50,50,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,0,0]
	# Initial guess at parameter values (nu1B, nu2B, nu1C, nu2C, T)
	p0 = [1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then exp/exp with mig demographic model
def params_split_exp_exp_mig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1B, nu2B, nu1C, nu2C, T, and m parameter values 
	nu1B = float(popt[0]) * float(nuA)
	nu2B = float(popt[1]) * float(nuA)
	nu1C = float(popt[2]) * float(nuA)
	nu2C = float(popt[3]) * float(nuA)
	T = float(popt[4]) * float(nuA) * 2
	m = float(popt[5])/(2 * float(nuA)) 
	# Return these computed parameter values
	return nuA,nu1B,nu2B,nu1C,nu2C,T,m,theta


# Define the demographic function for a population split followed by exponential change in population 1 and no change in population 2
def split_exp_const_mig(params, ns):
	"""
	Parameter values:
	1. nu1B = the ratio of population 1's bottleneck size to the ancestral population size
	2. nu1C = the ratio of population 1's current size to the ancestral population size
	3. nu2C = the ratio of population 2's current size to the ancestral population size
	4. T = the time since the population split
	5. m = the symmetric migration rate
	Overview: the ancestral population splits into population 1 with size nu1B and population 2 with size nu2C T generations from the present.
		After the split, the population 1 undergoes exponential change to nu1C and population stays constant at nu2C with symmetric migration between them.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1B, nu1C, nu2C, T, m = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth in pop 1 from size nu1B to nu1C over T generations
	nu1_func = lambda t: nu1B * (nu1C/nu1B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1_func(t), nu2C]
	# After the splitting T generations ago, pop 1 changes exponentially according to function while pop 2 stays the same. And there is symmetric migration between them
	fs.integrate(nu_func, T, m=np.array([[0, m], [m, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then exp/const change with mig demographic model
def instantiate_split_exp_const_mig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1B, nu1C, nu2C, T, m)
	upper_bound = [50,50,50,10,10]
	lower_bound = [1e-4,1e-4,1e-4,0,0]
	# Initial guess at parameter values (nu1B, nu1C, nu2C, T, m)
	p0 = [1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then exp/const with mig demographic model
def params_split_exp_const_mig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1B, nu1C, nu2C, T, and m parameter values 
	nu1B = float(popt[0]) * float(nuA)
	nu1C = float(popt[1]) * float(nuA)
	nu2C = float(popt[2]) * float(nuA)
	T = float(popt[3]) * float(nuA) * 2
	m = float(popt[4])/(2 * float(nuA))  
	# Return these computed parameter values
	return nuA,nu1B,nu1C,nu2C,T,m,theta


# Define the demographic function for a population split followed by no change in population 1, exponential change in population 2, and symmetric migration
def split_const_exp_mig(params, ns):
	"""
	Parameter values:
	1. nu1C = the ratio of population 1's current size to the ancestral population size
	2. nu2B = the ratio of population 2's bottleneck size to the ancestral population size
	3. nu2C = the ratio of population 2's current size to the ancestral population size
	4. T = the time since the population split
	5. m = the symmetric migration rate
	Overview: the ancestral population splits into population 1 with size nu1C and population 2 with size nu2B T generations from the present.
		After the split, the population 1 stays constant at nu1C and population 2 undergoes exponential change to nu2C with symmetric migration between them.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1C, nu2B, nu2C, T, m = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth in pop 1 from size nu1B to nu1C over T generations
	nu2_func = lambda t: nu2B * (nu2C/nu2B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1C, nu2_func(t)]
	# After the splitting T generations ago, pop 1 stays the same while population 2 changes exponentially according to function. And there is symmetric migration between them
	fs.integrate(nu_func, T, m=np.array([[0, m], [m, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then const/exp change with mig demographic model
def instantiate_split_const_exp_mig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1C, nu2B, nu2C, T, m)
	upper_bound = [50,50,50,10,10]
	lower_bound = [1e-4,1e-4,1e-4,0,0]
	# Initial guess at parameter values (nu1C, nu2B, nu2C, T, m)
	p0 = [1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then const/exp with mig demographic model
def params_split_const_exp_mig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1B, nu1C, nu2C, T, and m parameter values 
	nu1C = float(popt[0]) * float(nuA)
	nu2B = float(popt[1]) * float(nuA)
	nu2C = float(popt[2]) * float(nuA)
	T = float(popt[3]) * float(nuA) * 2
	m = float(popt[4])/(2 * float(nuA))  
	# Return these computed parameter values
	return nuA,nu1C,nu2B,nu2C,T,m,theta


# Define the demographic function for a population split followed by a post-split size change with multiple mig parameters
def split_two_epoch_mig(params, ns):
	"""
	Parameter values:
	1. nu1E1 = the ratio of population 1 size to the ancestral population size during first epoch
	2. nu2E1 = the ratio of population 2 size to the ancestral population size during the first epoch
	3. nu1E2 = the ratio of population 1 size to the ancestral population size during second epoch
	4. nu2E2 = the ratio of population 2 size to the ancestral population size during the second epoch
	5. TE1 = the duration of the first epoch
	6. TE2 = the duration of the second epoch
	7. mE1 = the migration rate during the first epoch
	8. mE2 = the migration rate during the second epoch
	Overview: the ancestral population splits into population 1 with size nu1E1 and population 2 with size nu2E1 and migration rate mE1. After TE1 generations, the populations enter a second epoch with sizes nu1E1 and nu2E2 that lasts for TE2 generations with migration rate mE2.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations with migration rate mE1
	fs.integrate([nu1E1, nu2E1], TE1, m=np.array([[0, mE1], [mE1, 0]]))
	# After the first epoch, the populations take constant sizes nu1E2 and nu2E2 for TE2 generations with migration rate mE2
	fs.integrate([nu1E2, nu2E2], TE2, m=np.array([[0, mE2], [mE2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split_two_epoch demographic model with mig
def instantiate_split_two_epoch_mig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2)
	upper_bound = [50,50,50,50,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2)
	p0 = [1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split_two_epoch with mig demographic model
def params_split_two_epoch_mig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	TE1 = float(popt[4]) * float(nuA) * 2
	TE2 = float(popt[5]) * float(nuA) * 2
	mE1 = float(popt[6])/(2 * float(nuA))
	mE2 = float(popt[7])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,TE1,TE2,mE1,mE2,theta



#--------------------------------------------------
# Two population models with asymmetric migration:
#--------------------------------------------------

# Define the demographic function for a population split followed by constant sizes for both population 1 and 2 with asymmetric migration
def split_asymmig(params, ns):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1's size to the ancestral population size
	2. nu2 = the ratio of population 2's size to the ancestral population size
	3. T = the time since the population split
	4. m12 = the migration rate from population 2 into population 1 
	5. m21 = the migration rate from population 1 into population 2
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.
		After the split, there is asymmetric migration.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1, nu2, T, m12, m21 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# After the splitting T generations ago, the populations have constant size and experience asymmetric migration
	fs.integrate([nu1, nu2], T, m=np.array([[0, m12], [m21, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split with asymmig demographic model
def instantiate_split_asymmig():
	# Set upper and lower bounds for the list of parameter choices (nu1, nu2, T, m12, m21)
	upper_bound = [50,50,10,10,10]
	lower_bound = [1e-4,1e-4,0,0,0]
	# Initial guess at parameter values (nu1, nu2, T, m12, m21)
	p0 = [1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split with asymmig demographic model
def params_split_asymmig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1, nu2, T, m12, m21 parameter values 
	nu1 = float(popt[0]) * float(nuA)
	nu2 = float(popt[1]) * float(nuA)
	T = float(popt[2]) * float(nuA) * 2
	m12 = float(popt[3])/(2 * float(nuA)) 
	m21 = float(popt[4])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1,nu2,T,m12,m21,theta


# Define the demographic function for a population split followed by exponential size change for both population 1 and 2 with asymmetric migration
def split_exp_exp_asymmig(params, ns):
	"""
	Parameter values:
	1. nu1B = the ratio of population 1's bottleneck size to the ancestral population size
	2. nu2B = the ratio of population 2's bottleneck size to the ancestral population size
	3. nu1C = the ratio of population 1's current size to the ancestral population size
	4. nu2C = the ratio of population 2's current size to the ancestral population size
	5. T = the time since the population split
	6. m12 = the migration rate from population 2 into population 1
	7. m21 = the migration rate from population 1 into population 2
	Overview: the ancestral population splits into population 1 with size nu1B and population 2 with size nu2B T generations from the present.
		After the split, the populations both change expoentially to sizes nu1C and nu2C, respectively, and there is asymmetric migration at rate m.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1B, nu2B, nu1C, nu2C, T, m12, m21 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth from size nu1B to nu1C over T generations
	nu1_func = lambda t: nu1B * (nu1C/nu1B) ** (t/T)
	# Create the same func for other pop (describes exponential growth from size nu2B to nu2C over T generations)
	nu2_func = lambda t: nu2B * (nu2C/nu2B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1_func(t), nu2_func(t)]
	# After the splitting T generations ago, the populations both change exponentially according to their functions with asymmetric migration
	fs.integrate(nu_func, T, m=np.array([[0, m12], [m21, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then exp/exp change with asymmig demographic model
def instantiate_split_exp_exp_asymmig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1B, nu2B, nu1C, nu2C, T, m12, m21)
	upper_bound = [50,50,50,50,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,0,0,0]
	# Initial guess at parameter values (nu1B, nu2B, nu1C, nu2C, T, m12, m21)
	p0 = [1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then exp/exp with asymmig demographic model
def params_split_exp_exp_asymmig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1B, nu2B, nu1C, nu2C, T, and m12, m21 parameter values 
	nu1B = float(popt[0]) * float(nuA)
	nu2B = float(popt[1]) * float(nuA)
	nu1C = float(popt[2]) * float(nuA)
	nu2C = float(popt[3]) * float(nuA)
	T = float(popt[4]) * float(nuA) * 2
	m12 = float(popt[5])/(2 * float(nuA)) 
	m21 = float(popt[6])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1B,nu2B,nu1C,nu2C,T,m12,m21,theta


# Define the demographic function for a population split followed by exponential change in population 1 and no change in population 2
def split_exp_const_asymmig(params, ns):
	"""
	Parameter values:
	1. nu1B = the ratio of population 1's bottleneck size to the ancestral population size
	2. nu1C = the ratio of population 1's current size to the ancestral population size
	3. nu2C = the ratio of population 2's current size to the ancestral population size
	4. T = the time since the population split
	5. m12 = the migration rate from population 2 into population 1
	6. m21 = the migration rate from population 1 into population 2
	Overview: the ancestral population splits into population 1 with size nu1B and population 2 with size nu2C T generations from the present.
		After the split, the population 1 undergoes exponential change to nu1C and population stays constant at nu2C with symmetric migration between them.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1B, nu1C, nu2C, T, m12, m21 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth in pop 1 from size nu1B to nu1C over T generations
	nu1_func = lambda t: nu1B * (nu1C/nu1B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1_func(t), nu2C]
	# After the splitting T generations ago, pop 1 changes exponentially according to function while pop 2 stays the same. And there is symmetric migration between them
	fs.integrate(nu_func, T, m=np.array([[0, m12], [m21, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then exp/const change with asymmig demographic model
def instantiate_split_exp_const_asymmig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1B, nu1C, nu2C, T, m12, m21)
	upper_bound = [50,50,50,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,0,0,0]
	# Initial guess at parameter values (nu1B, nu1C, nu2C, T, m12, m21)
	p0 = [1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then exp/const with asymmig demographic model
def params_split_exp_const_asymmig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1B, nu1C, nu2C, T, and m12, m21 parameter values 
	nu1B = float(popt[0]) * float(nuA)
	nu1C = float(popt[1]) * float(nuA)
	nu2C = float(popt[2]) * float(nuA)
	T = float(popt[3]) * float(nuA) * 2
	m12 = float(popt[4])/(2 * float(nuA))  
	m21 = float(popt[5])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1B,nu1C,nu2C,T,m12,m21,theta


# Define the demographic function for a population split followed by no change in population 1, exponential change in population 2, and asymmetric migration
def split_const_exp_asymmig(params, ns):
	"""
	Parameter values:
	1. nu1C = the ratio of population 1's current size to the ancestral population size
	2. nu2B = the ratio of population 2's bottleneck size to the ancestral population size
	3. nu2C = the ratio of population 2's current size to the ancestral population size
	4. T = the time since the population split
	5. m12 = the migration rate from population 2 into population 1
	6. m21 = the migration rate from population 1 into population 2
	Overview: the ancestral population splits into population 1 with size nu1C and population 2 with size nu2B T generations from the present.
		After the split, the population 1 stays constant at nu1C and population 2 undergoes exponential change to nu2C with symmetric migration between them.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1C, nu2B, nu2C, T, m12, m21 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth in pop 1 from size nu1B to nu1C over T generations
	nu2_func = lambda t: nu2B * (nu2C/nu2B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1C, nu2_func(t)]
	# After the splitting T generations ago, pop 1 stays the same while population 2 changes exponentially according to function. And there is asymmetric migration between them
	fs.integrate(nu_func, T, m=np.array([[0, m12], [m21, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split then const/exp change with asymmig demographic model
def instantiate_split_const_exp_asymmig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1C, nu2B, nu2C, T, m12, m21)
	upper_bound = [50,50,50,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,0,0,0]
	# Initial guess at parameter values (nu1C, nu2B, nu2C, T, m12, m21)
	p0 = [1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split then const/exp with asymmig demographic model
def params_split_const_exp_asymmig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1B, nu1C, nu2C, T, and m12, m21 parameter values 
	nu1C = float(popt[0]) * float(nuA)
	nu2B = float(popt[1]) * float(nuA)
	nu2C = float(popt[2]) * float(nuA)
	T = float(popt[3]) * float(nuA) * 2
	m12 = float(popt[4])/(2 * float(nuA))  
	m21 = float(popt[5])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1C,nu2B,nu2C,T,m12,m21,theta


# Define the demographic function for a population split followed by a post-split size change with multiple asymmig parameters
def split_two_epoch_asymmig(params, ns):
	"""
	Parameter values:
	1. nu1E1 = the ratio of population 1 size to the ancestral population size during first epoch
	2. nu2E1 = the ratio of population 2 size to the ancestral population size during the first epoch
	3. nu1E2 = the ratio of population 1 size to the ancestral population size during second epoch
	4. nu2E2 = the ratio of population 2 size to the ancestral population size during the second epoch
	5. TE1 = the duration of the first epoch
	6. TE2 = the duration of the second epoch
	7. m12E1 = the migration rate during the first epoch from population 2 into population 1
	8. m21E1 = the migration rate during the first epoch from population 1 into population 2
	9. m12E2 = the migration rate during the second epoch from population 2 into population 1
	10. m21E2 = the migration rate during the second epoch from population 1 into population 2
	Overview: the ancestral population splits into population 1 with size nu1E1 and population 2 with size nu2E1 and migration rate mE1. After TE1 generations, the populations enter a second epoch with sizes nu1E1 and nu2E2 that lasts for TE2 generations with migration rate mE2.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, m12E1, m21E1, m12E2, m21E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations with migration rate mE1
	fs.integrate([nu1E1, nu2E1], TE1, m=np.array([[0, m12E1], [m21E1, 0]]))
	# After the first epoch, the populations take constant sizes nu1E2 and nu2E2 for TE2 generations with migration rate mE2
	fs.integrate([nu1E2, nu2E2], TE2, m=np.array([[0, m12E2], [m21E2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the split_two_epoch demographic model with asymmig
def instantiate_split_two_epoch_asymmig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, m12E1, m21E1, m12E2, m21E2)
	upper_bound = [50,50,50,50,10,10,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,0,0,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, m12E1, m21E1, m12E2, m21E2)
	p0 = [1,1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split_two_epoch with asymmig demographic model
def params_split_two_epoch_asymmig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, m12E1, m21E1, m12E2, m21E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	TE1 = float(popt[4]) * float(nuA) * 2
	TE2 = float(popt[5]) * float(nuA) * 2
	m12E1 = float(popt[6])/(2 * float(nuA))
	m21E1 = float(popt[7])/(2 * float(nuA))
	m12E2 = float(popt[8])/(2 * float(nuA))
	m21E2 = float(popt[9])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,TE1,TE2,m12E1,m21E1,m12E2,m21E2,theta



#--------------------------------------------------
# Three population model
#--------------------------------------------------

# The best-fit two population demographic models support the following features:
# 1) A more recent split between Saturna and Pender compared to Saturna–Maple Ridge and Pender–Maple Ridge
# 2) Variable post-split migration rates between Saturna–Pender and the two islands to the mainland

# Given this, define a three-population demographic model with two post-split epochs. 
# During the first epoch, allow distinct Ne for the island-ancestor and mainland-ancestor and symmetric migration between the two
# The second epoch is defined by an additional branching of the island-ancestor into the two extant island populations
# During the second epoch, allow distinct Ne for each extant island and mainland population and symmetric (but unique) migration rates among the three populations
def tri_split_two_epoch_mig(params, ns):
	"""
	Population order:
	1. Maple Ridge
	2. Pender
	3. Saturna

	Parameter values:
	1. nu1E1 = the ratio of the mainland ancestor's size to the ancestral population size during first epoch
	2. nu2E1 = the ratio of the island ancestor's size to the anestral population size during the first epoch
	3. nu1E2 = the ratio of the Maple Ridge size to the ancestral population size during the second epoch
	4. nu2E2 = the ratio of the Pender size to the ancestral population size during the second epoch
	5. nu3E2 = the ratio of the Saturna size to the ancestral population size during the second epoch
	6. TE1 = the duration of the first epoch during which there exists the mainland ancestor and island ancestor
	7. TE2 = the duration of the second epoch during which there exists all three contemporary populations
	8. m12E1 = the symmetric migration rate between the mainland ancestor and island ancestor during the first epoch
	9. m12E2 = the symmetric migration rate between Maple Ridge and Pender during the second epoch
	10. m13E2 = the symmetric migration rate between Maple Ridge and Saturna during the second epoch
	11. m23E2 = the symmetric migration rate between Saturna and Pender during the second epoch

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2, m23E2 = params
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
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates m12E2 between MR-PI, m13E2 between MR-SI, and m23E2 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, m12E2, m13E2], [m12E2, 0, m23E2], [m13E2, m23E2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2, m23E2)
	upper_bound = [50,50,50,50,50,10,10,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2, m23E2)
	p0 = [1,1,1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2, m23E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E1 = float(popt[7])/(2 * float(nuA))
	m12E2 = float(popt[8])/(2 * float(nuA))
	m13E2 = float(popt[9])/(2 * float(nuA))
	m23E2 = float(popt[10])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E1,m12E2,m13E2,m23E2,theta


#---------------------------------------------------------
# Three population models (minus 1 migration parameter)
#---------------------------------------------------------

def tri_split_two_epoch_mig_minus1_A(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, m13E2, m23E2
	Excludes: m12E1
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E2, m13E2, m23E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Epoch 1 dynamics #
	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# Epoch 2 dynamics #
	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates m12E2 between MR-PI, m13E2 between MR-SI, and m23E2 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, m12E2, m13E2], [m12E2, 0, m23E2], [m13E2, m23E2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus1_A():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, m13E2, m23E2)
	upper_bound = [50,50,50,50,50,10,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, m13E2, m23E2)
	p0 = [1,1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus1_A(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, m13E2, m23E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E2 = float(popt[7])/(2 * float(nuA))
	m13E2 = float(popt[8])/(2 * float(nuA))
	m23E2 = float(popt[9])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E2,m13E2,m23E2,theta

def tri_split_two_epoch_mig_minus1_B(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, m13E2, m23E2
	Excludes: m12E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m13E2, m23E2 = params
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
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates 0 between MR-PI, m13E2 between MR-SI, and m23E2 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, 0, m13E2], [0, 0, m23E2], [m13E2, m23E2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus1_B():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, m13E2, m23E2)
	upper_bound = [50,50,50,50,50,10,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, m13E2, m23E2)
	p0 = [1,1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus1_B(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, m13E2, m23E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E1 = float(popt[7])/(2 * float(nuA))
	m13E2 = float(popt[8])/(2 * float(nuA))
	m23E2 = float(popt[9])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E1,m13E2,m23E2,theta

def tri_split_two_epoch_mig_minus1_C(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, NULL, m23E2
	Excludes: m13E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m23E2 = params
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
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates m12E2 between MR-PI, 0 between MR-SI, and m23E2 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, m12E2, 0], [m12E2, 0, m23E2], [0, m23E2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus1_C():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, NULL, m23E2)
	upper_bound = [50,50,50,50,50,10,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, NULL, m23E2)
	p0 = [1,1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus1_C(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, NULL, m23E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E1 = float(popt[7])/(2 * float(nuA))
	m12E2 = float(popt[8])/(2 * float(nuA))
	m23E2 = float(popt[9])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E1,m12E2,m23E2,theta

def tri_split_two_epoch_mig_minus1_D(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2, NULL
	Excludes: m23E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2 = params
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
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates m12E2 between MR-PI, m13E2 between MR-SI, and 0 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, m12E2, m13E2], [m12E2, 0, 0], [m13E2, 0, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus1_D():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2, NULL)
	upper_bound = [50,50,50,50,50,10,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2, NULL)
	p0 = [1,1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus1_D(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, m13E2, NULL parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E1 = float(popt[7])/(2 * float(nuA))
	m12E2 = float(popt[8])/(2 * float(nuA))
	m13E2 = float(popt[9])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E1,m12E2,m13E2,theta

#---------------------------------------------------------
# Three population models (minus 2 migration parameters)
#---------------------------------------------------------

def tri_split_two_epoch_mig_minus2_A(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, m13E2, m23E2
	Excludes: m12E1, m12E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m13E2, m23E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Epoch 1 dynamics #
	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# Epoch 2 dynamics #
	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates 0 between MR-PI, m13E2 between MR-SI, and m23E2 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, 0, m13E2], [0, 0, m23E2], [m13E2, m23E2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus2_A():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, m13E2, m23E2)
	upper_bound = [50,50,50,50,50,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, m13E2, m23E2)
	p0 = [1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus2_A(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, m13E2, m23E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m13E2 = float(popt[7])/(2 * float(nuA))
	m23E2 = float(popt[8])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m13E2,m23E2,theta

def tri_split_two_epoch_mig_minus2_B(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, NULL, m23E2
	Excludes: m12E1, m13E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E2, m23E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Epoch 1 dynamics #
	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# Epoch 2 dynamics #
	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates m12E2 between MR-PI, 0 between MR-SI, and m23E2 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, m12E2, 0], [m12E2, 0, m23E2], [0, m23E2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus2_B():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, NULL, m23E2)
	upper_bound = [50,50,50,50,50,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, NULL, m23E2)
	p0 = [1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus2_B(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, NULL, m23E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E2 = float(popt[7])/(2 * float(nuA))
	m23E2 = float(popt[8])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E2,m23E2,theta

def tri_split_two_epoch_mig_minus2_C(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, m13E2, NULL
	Excludes: m12E1, m23E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E2, m13E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Epoch 1 dynamics #
	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# Epoch 2 dynamics #
	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates m12E2 between MR-PI, m13E2 between MR-SI, and 0 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, m12E2, m13E2], [m12E2, 0, 0], [m13E2, 0, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus2_C():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, m13E2, NULL)
	upper_bound = [50,50,50,50,50,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, m13E2, NULL)
	p0 = [1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus2_C(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, m13E2, NULL parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E2 = float(popt[7])/(2 * float(nuA))
	m13E2 = float(popt[8])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E2,m13E2,theta

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
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus2_D():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, NULL, m23E2)
	upper_bound = [50,50,50,50,50,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, NULL, m23E2)
	p0 = [1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus2_D(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, NULL, m23E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E1 = float(popt[7])/(2 * float(nuA))
	m23E2 = float(popt[8])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E1,m23E2,theta

def tri_split_two_epoch_mig_minus2_E(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, m13E2, NULL
	Excludes: m12E2, m23E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m13E2 = params
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
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates 0 between MR-PI, m13E2 between MR-SI, and 0 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, 0, m13E2], [0, 0, 0], [m13E2, 0, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus2_E():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, m13E2, NULL)
	upper_bound = [50,50,50,50,50,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, m13E2, NULL)
	p0 = [1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus2_E(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, m13E2, NULL parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E1 = float(popt[7])/(2 * float(nuA))
	m13E2 = float(popt[8])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E1,m13E2,theta

def tri_split_two_epoch_mig_minus2_F(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, NULL, NULL
	Excludes: m13E2, m23E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2 = params
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
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates m12E2 between MR-PI, 0 between MR-SI, and 0 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, m12E2, 0], [m12E2, 0, 0], [0, 0, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus2_F():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, NULL, NULL)
	upper_bound = [50,50,50,50,50,10,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, NULL, NULL)
	p0 = [1,1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus2_F(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, m12E2, NULL, NULL parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E1 = float(popt[7])/(2 * float(nuA))
	m12E2 = float(popt[8])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E1,m12E2,theta

#---------------------------------------------------------
# Three population models (minus 3 migration parameters)
#---------------------------------------------------------

def tri_split_two_epoch_mig_minus3_A(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, NULL, m23E2
	Excludes: m12E1, m12E2, m13E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m23E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Epoch 1 dynamics #
	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# Epoch 2 dynamics #
	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates 0 between MR-PI, 0 between MR-SI, and m23E2 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, 0, 0], [0, 0, m23E2], [0, m23E2, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus3_A():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, NULL, m23E2)
	upper_bound = [50,50,50,50,50,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, NULL, m23E2)
	p0 = [1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus3_A(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, NULL, m23E2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m23E2 = float(popt[7])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m23E2,theta

def tri_split_two_epoch_mig_minus3_B(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, m13E2, NULL
	Excludes: m12E1, m12E2, m23E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m13E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Epoch 1 dynamics #
	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# Epoch 2 dynamics #
	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates 0 between MR-PI, m13E2 between MR-SI, and 0 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, 0, m13E2], [0, 0, 0], [m13E2, 0, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus3_B():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, m13E2, NULL)
	upper_bound = [50,50,50,50,50,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, m13E2, NULL)
	p0 = [1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus3_B(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, m13E2, NULL parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m13E2 = float(popt[7])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m13E2,theta

def tri_split_two_epoch_mig_minus3_C(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, NULL, NULL
	Excludes: m12E1, m13E2, m23E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Epoch 1 dynamics #
	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# Epoch 2 dynamics #
	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates m12E2 between MR-PI, 0 between MR-SI, and 0 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, m12E2, 0], [m12E2, 0, 0], [0, 0, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus3_C():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, NULL, NULL)
	upper_bound = [50,50,50,50,50,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, NULL, NULL)
	p0 = [1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus3_C(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, m12E2, NULL, NULL parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E2 = float(popt[7])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E2,theta

def tri_split_two_epoch_mig_minus3_D(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, NULL, NULL
	Excludes: m13E2, m23E2, m12E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1 = params
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
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations with migration rates 0 between MR-PI, 0 between MR-SI, and 0 between PI-SI
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2, m=np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus3_D():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, NULL, NULL)
	upper_bound = [50,50,50,50,50,10,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, NULL, NULL)
	p0 = [1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus3_D(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, m12E1, NULL, NULL, NULL parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	m12E1 = float(popt[7])/(2 * float(nuA))
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,m12E1,theta

#------------------------------------------------------------
# Three population models (minus all 4 migration parameters)
#------------------------------------------------------------

def tri_split_two_epoch_mig_minus4_A(params, ns):
	"""
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, NULL, NULL
	Excludes: m12E1, m13E2, m23E2, m12E2
	"""
	# Ancestral dynamics #
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2 = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Epoch 1 dynamics #
	# Create a 2D FS that represents the mainland ancestor and the island ancestor
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], (ns[1] + ns[2]))
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations
	fs.integrate([nu1E1, nu2E1], TE1)
	# Epoch 2 dynamics #
	# Create a 3D FS that represents the extant island and mainland populations
	# Do this by splitting the island ancestor into the contemporary Pender and Saturn populations
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	# After the second epoch, the populations take constant sizes nu1E2, nu2E2, and nu3E2 (for MR, PI, and SI) for TE2 generations 
	fs.integrate([nu1E2, nu2E2, nu3E2], TE2)
	# Return the frequency spectrum obtained from the model
	return fs
# Function to gather optimization starting information for the tri_split_two_epoch demographic model with mig
def instantiate_tri_split_two_epoch_mig_minus4_A():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, NULL, NULL)
	upper_bound = [50,50,50,50,50,10,10]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, NULL, NULL)
	p0 = [1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the tri_split_two_epoch with mig demographic model
def params_tri_split_two_epoch_mig_minus4_A(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, nu3E2, TE1, TE2, NULL, NULL, NULL, NULL parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	nu3E2 = float(popt[4]) * float(nuA)
	TE1 = float(popt[5]) * float(nuA) * 2
	TE2 = float(popt[6]) * float(nuA) * 2
	# Return these computed parameter values
	return nuA,nu1E1,nu2E1,nu1E2,nu2E2,nu3E2,TE1,TE2,theta


#--------------------------------------------
# Edit some internal dadi plotting functions
#--------------------------------------------
def plot_1d_comp_Poisson(model, data, pop_name, model_id, ns, fig_num=None, residual='Anscombe',
						 plot_masked=False, show=True):
	"""
	Poisson comparison between 1d model and data.


	model: 1-dimensional model SFS
	data: 1-dimensional data SFS
	fig_num: Clear and use figure fig_num for display. If None, an new figure
			 window is created.
	residual: 'Anscombe' for Anscombe residuals, which are more normally
			  distributed for Poisson sampling. 'linear' for the linear
			  residuals, which can be less biased.
	plot_masked: Additionally plots (in open circles) results for points in the 
				 model or data that were masked.
	show: If True, execute pylab.show command to make sure plot displays.
	"""
	if fig_num is None:
		f = pylab.gcf()
	else:
		f = pylab.figure(fig_num, figsize=(7,7))
	pylab.clf()

	if data.folded and not model.folded:
		model = model.fold()

	masked_model, masked_data = dadi.Numerics.intersect_masks(model, data)

	ax = pylab.subplot(2,1,1)
	pylab.semilogy(masked_data, '-ob', label='data')
	pylab.semilogy(masked_model, '-or', label='model')

	if plot_masked:
		pylab.semilogy(masked_data.data, '--ob', mfc='w', zorder=-100)
		pylab.semilogy(masked_model.data, '--or', mfc='w', zorder=-100)

	ax.legend(loc='lower left')

	res_ax = pylab.subplot(2,1,2, sharex = ax)
	if residual == 'Anscombe':
		resid = dadi.Inference.Anscombe_Poisson_residual(masked_model, masked_data)
	elif residual == 'linear':
		resid = dadi.Inference.linear_Poisson_residual(masked_model, masked_data)
	else:
		raise ValueError("Unknown class of residual '%s'." % residual)
	pylab.plot(resid, '-og')
	if plot_masked:
		pylab.plot(resid.data, '--og', mfc='w', zorder=-100)

	# Edited portion (4/19/23)-------------------------------------------------------------------
	ax.set_xlim(0, data.shape[0]-(ns/2)) 

	f.suptitle((str(pop_name) + "; " + str(model_id) + " model fit"))
	ax.set_xlabel('Folded SFS bin')
	ax.set_ylabel('Number of SNPs')
	res_ax.set_xlabel('Folded SFS bin')
	res_ax.set_ylabel('Residuals')

	pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
	#----------------------------------------------------------------------------------------------

	if show:
		pylab.show()

def plot_1d_comp_multinom(model, data, pop_name, model_id, ns, fig_num=None, residual='Anscombe',
						  plot_masked=False, show=True):
	"""
	Mulitnomial comparison between 1d model and data.


	model: 1-dimensional model SFS
	data: 1-dimensional data SFS
	fig_num: Clear and use figure fig_num for display. If None, an new figure
			 window is created.
	residual: 'Anscombe' for Anscombe residuals, which are more normally
			  distributed for Poisson sampling. 'linear' for the linear
			  residuals, which can be less biased.
	plot_masked: Additionally plots (in open circles) results for points in the 
				 model or data that were masked.
	show: If True, execute pylab.show command to make sure plot displays.

	This comparison is multinomial in that it rescales the model to optimally
	fit the data.
	"""
	model = dadi.Inference.optimally_scaled_sfs(model, data)

	plot_1d_comp_Poisson(model, data, pop_name, model_id, ns, fig_num, residual,
						 plot_masked, show)

#------------------------
# Model comparison
#------------------------

# Print out the best-fit model parameters
print("Best-fit " + str(pop_name) + " parameters for the " + str(model_id) + " model: " + ", ".join(map(str, popt_est)))

# Determine what the sample size is based on the input fs
ns = fs.sample_sizes

# Need to determine pop_ids based on input fs
pop_ids = fs.pop_ids

# Wrap our demographic function based on user specified model
if (model_id=="split"):
	model_func = split
elif (model_id=="split_exp_exp"):
	model_func = split_exp_exp
elif (model_id=="split_exp_const"):
	model_func = split_exp_const
elif (model_id=="split_const_exp"):
	model_func = split_const_exp
elif (model_id=="split_two_epoch"):
	model_func = split_two_epoch
#
elif (model_id=="split_mig"):
	model_func = split_mig
elif (model_id=="split_exp_exp_mig"):
	model_func = split_exp_exp_mig
elif (model_id=="split_exp_const_mig"):
	model_func = split_exp_const_mig
elif (model_id=="split_const_exp_mig"):
	model_func = split_const_exp_mig
elif (model_id=="split_two_epoch_mig"):
	model_func = split_two_epoch_mig
#
elif (model_id=="split_asymmig"):
	model_func = split_asymmig
elif (model_id=="split_exp_exp_asymmig"):
	model_func = split_exp_exp_asymmig
elif (model_id=="split_exp_const_asymmig"):
	model_func = split_exp_const_asymmig
elif (model_id=="split_const_exp_asymmig"):
	model_func = split_const_exp_asymmig
elif (model_id=="split_two_epoch_asymmig"):
	model_func = split_two_epoch_asymmig
#
elif (model_id=="tri_split_two_epoch_mig"):
	model_func = tri_split_two_epoch_mig
#
elif (model_id=="tri_split_two_epoch_mig_minus1_A"):
	model_func = tri_split_two_epoch_mig_minus1_A
elif (model_id=="tri_split_two_epoch_mig_minus1_B"):
	model_func = tri_split_two_epoch_mig_minus1_B
elif (model_id=="tri_split_two_epoch_mig_minus1_C"):
	model_func = tri_split_two_epoch_mig_minus1_C
elif (model_id=="tri_split_two_epoch_mig_minus1_D"):
	model_func = tri_split_two_epoch_mig_minus1_D
#
elif (model_id=="tri_split_two_epoch_mig_minus2_A"):
	model_func = tri_split_two_epoch_mig_minus2_A
elif (model_id=="tri_split_two_epoch_mig_minus2_B"):
	model_func = tri_split_two_epoch_mig_minus2_B
elif (model_id=="tri_split_two_epoch_mig_minus2_C"):
	model_func = tri_split_two_epoch_mig_minus2_C
elif (model_id=="tri_split_two_epoch_mig_minus2_D"):
	model_func = tri_split_two_epoch_mig_minus2_D
elif (model_id=="tri_split_two_epoch_mig_minus2_E"):
	model_func = tri_split_two_epoch_mig_minus2_E
elif (model_id=="tri_split_two_epoch_mig_minus2_F"):
	model_func = tri_split_two_epoch_mig_minus2_F
#
elif (model_id=="tri_split_two_epoch_mig_minus3_A"):
	model_func = tri_split_two_epoch_mig_minus3_A
elif (model_id=="tri_split_two_epoch_mig_minus3_B"):
	model_func = tri_split_two_epoch_mig_minus3_B
elif (model_id=="tri_split_two_epoch_mig_minus3_C"):
	model_func = tri_split_two_epoch_mig_minus3_C
elif (model_id=="tri_split_two_epoch_mig_minus3_D"):
	model_func = tri_split_two_epoch_mig_minus3_D
#
elif (model_id=="tri_split_two_epoch_mig_minus4_A"):
	model_func = tri_split_two_epoch_mig_minus4_A

# Calculate the best-fit model FS
model = model_func(popt, ns)

# Calculate the likelihood of the data given the model FS
ll_model = moments.Inference.ll_multinom(model, fs)

# Calculate the corresponding theta value
theta = moments.Inference.optimal_sfs_scaling(model, fs)

# Print out the computed likelihood
print("Maximum composite likelihood: " + str(ll_model))

# Print out the optimal value of theta
print("Optimal value of theta: " + str(theta))

# Plot a comparison of the resulting model FS with the data

# For 1 dimensional analysis
if (dim==1):
	# Open the figure
	pylab.figure(figsize=(7,6))
	# Plot the 1D comparison using the edited function
	plot_1d_comp_multinom(model, fs, pop_name, model_id, ns, show=False)
	# Save the figure
	pylab.savefig((str(pop_name) + "_" + str(model_id) + "_residuals.png"), dpi=250)

# For 2 dimensional analysis
elif (dim==2):
	# Open the figure
	pylab.figure(figsize=(7,6))
	# Plot the 2D comparison
	dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=None, vmax=None, show=False)
	# Save the figure
	pylab.savefig((str(pop_name) + "_" + str(model_id) + "_residuals.png"), dpi=250)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the second population
	data_fs1 = fs.marginalize([1])
	model_fs1 = model.marginalize([1])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs1, data_fs1, "Population 1", model_id, data_fs1.sample_sizes, show=False)
	# Save the figure
	pylab.savefig((str(pop_name) + "_" + str(model_id) + "_pop1_1Dresiduals.png"), dpi=250)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the first population
	data_fs2 = fs.marginalize([0])
	model_fs2 = model.marginalize([0])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs2, data_fs2, "Population 2", model_id, data_fs2.sample_sizes, show=False)
	# Save the figure
	pylab.savefig((str(pop_name) + "_" + str(model_id) + "_pop2_1Dresiduals.png"), dpi=250)

# For 3 dimensional analysis
elif (dim==3):
	# Open the figure
	pylab.figure(figsize=(15,15))
	# Plot the 2D comparison
	dadi.Plotting.plot_3d_comp_multinom(model, fs, vmin=1, vmax=None)
	# Save the figure
	pylab.savefig((str(pop_name) + "_" + str(model_id) + "_residuals.png"), dpi=250)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the second and third populations
	data_fs1 = fs.marginalize([1,2])
	model_fs1 = model.marginalize([1,2])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs1, data_fs1, "Population 1", model_id, data_fs1.sample_sizes, show=False)
	# Save the figure
	pylab.savefig((str(pop_name) + "_" + str(model_id) + "_pop1_1Dresiduals.png"), dpi=250)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the first and third populations
	data_fs2 = fs.marginalize([0,2])
	model_fs2 = model.marginalize([0,2])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs2, data_fs2, "Population 2", model_id, data_fs2.sample_sizes, show=False)
	# Save the figure
	pylab.savefig((str(pop_name) + "_" + str(model_id) + "_pop2_1Dresiduals.png"), dpi=250)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the first and second populations
	data_fs3 = fs.marginalize([0,1])
	model_fs3 = model.marginalize([0,1])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs3, data_fs3, "Population 3", model_id, data_fs3.sample_sizes, show=False)
	# Save the figure
	pylab.savefig((str(pop_name) + "_" + str(model_id) + "_pop3_1Dresiduals.png"), dpi=250)










	