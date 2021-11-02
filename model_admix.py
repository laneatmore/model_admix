#!/usr/bin/env python3

##########################################################################################
# Model for simulating admixture through msprime. Generates plink analysis and ADMIXTURE #
# results for simulated admixed population. Input is 3 populations - 2 parent, 1 admixed #
# 																						 #
# To begin simulation: run model_admix.py in shell with the following parameters:        #
# Pop1 Current size, Pop2 Current size, Pop3(admixed) Current size, time admixture(ybp), #
# prop of admixed population from pop1 (decimal), prop of admixed population from pop2,  #
# length of sequence sampled	                                                         #
#														                                 #
#            Before using this model it is necessary to install msprime! 				 #
#            https://msprime.readthedocs.io/en/stable/introduction.html					 #
#								Written by Lane M Atmore 2019 							 #	
#									Updated Oct 2021				             		 #
##########################################################################################


#NOTE: This program was written with msprime0.7
#msprime 1.0 is now out and is much more user-friendly
#however, the developers offer indefinite backwards compatibility with old APIs, therefore
#this has not been updated. The msprime.simulate() option is analogous to the new
#.demography and .sim_ancestry options

#In 0.7 msprime also produced some VCF quirks which were not compatible with downstream
#packages, therefore some scripts here were developed to modify the VCF output from msprime
#To my knowledge this is no longer necessary with the new version

##required packages
import numpy as np
import subprocess
import sys
import os
import pandas as pd
import msprime

###Updated for 2021 manuscript Oct 2021####

##define arguments and argument order
##arg1 - population 1 initial size
##arg2 - population 2 initial size
##arg3 - admixed population (pop3) initial size
##arg4 - time of admixture in years ago
##arg5 - proportion of admixed population from population 1 as decimal
##arg6 - proportion of admixed population from population 2 as decimal 
##arg7 - length of genome you want to simulate

#for the CAHG we want to model a few different things
#two ancestral populations -- ancient Pygmies and Bantu, plus current slightly admixed CAHG
#edit the initial population sizes where necessary
#currently optimized for modeling admixture between CAHG and Bantu populations

#chromosome links to the chromosome lengths and recombination maps

#dem_option will tell the program to model constant size, recent collapse, recent expansion


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

pop1 = int(sys.argv[1]) #pop1 initial size
pop2 = int(sys.argv[2]) #pop2 initial size
pop3 = int(sys.argv[3]) #adm initial size
time_admix = int(sys.argv[4])
prop_pop1 = float(sys.argv[5]) #pop1 admixture proportion
prop_pop2 = float(sys.argv[6]) #pop2 admixture proportion
chrom = int(sys.argv[7]) #chromosome
dem_option = sys.argv[8] #which model?
sample_pop1 = int(sys.argv[9]) #specify sample size for ancestor 1
sample_pop2 = int(sys.argv[10]) #specify sample size for ancestor 2
sample_pop3 = int(sys.argv[11]) #specify sample size for admixed population

if not os.path.exists("Admixture/"):
	os.mkdir("Admixture/")
	
os.chdir("Admixture/")

#Simulation function

def model_admix_constant(pop1, pop2, pop3, time_admix, prop_pop1, prop_pop2, chrom, sample_pop1, sample_pop2, sample_pop3):

	
	#defining the variables for the simulation by scaling to args			
	Pop1 = pop1
	Pop2 = pop2
	Pop3 = pop3
	generation_time = 30
	T_Pop1 = 10000000/generation_time  #pygmies
	T_Pop2 = 5000000/generation_time  #bantu
	T_Admix = time_admix/generation_time
	r_Pop3 = (pop3/100)**(1/T_Admix) - 1
	r_Pop2_1 = (pop2/1000)**(1/T_Admix) - 1
	r_Pop2_2 = (0.5*pop2 / 20)**(1 / (T_Pop2 - T_Admix)) - 1
	r_Pop1 = (pop3/100)**(1/T_Pop1)-1
	s_Pop2_0 = 500
	s_Pop3_0 = 100
	m_Pop1 = 0.1*prop_pop1
	m_Pop2 = 0.1*prop_pop2
	
	chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','18',
	'19','20','21','22']
	
	length = ['248956422','242193529','198295559','190214555','181538259','170805979',
	'159345973','145138636','138394717','133797422','135086622','133275309','114364328',
	'107043718','101991189','90338345','83257441','80373285','58617616','64444167',
	'46709983','50818468']
	
	recomb_rate = ['1.14856e-08','1.10543e-08', '1.12796e-08','1.12312e-08','1.12809e-08',
	'1.12229e-08','1.17646e-08','1.14785e-08','1.17807e-08','1.33651e-08','1.17193e-08',
	'1.30502e-08','1.09149e-08','1.11973e-08','1.38358e-08','1.48346e-08','1.58249e-08',
	'1.5076e-08','1.82201e-08','1.71783e-08','1.30452e-08','1.4445e-08']
	
	chrom_tuple = list(zip(chroms, length, recomb_rate))
	chrom_map = pd.DataFrame(chrom_tuple, columns = ['chroms', 'length', 'recomb_rate'])
	
	chrom_map_pos = (1-int(chrom))
	
	chrom_recomb_rate = chrom_map.iloc[chrom_map_pos]['recomb_rate']
	chrom_length = chrom_map.iloc[chrom_map_pos]['length']
	
	print("BEGINNING SIMULATION:", str(sys.argv), flush = True)

	#begin simulation
	sim = msprime.simulate(
		length = int(chrom_length), 
		recombination_rate = float(chrom_recomb_rate),
		mutation_rate = 1.29e-8, #human mutation rate
		population_configurations = [
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop1)), 
			initial_size = Pop1,
			growth_rate = r_Pop1
			), 
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop2)), 
			initial_size = Pop2, 
			growth_rate = r_Pop2_1
			), 
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop3)), 
			initial_size = Pop3
			)], #no growth rate for Pop3, constant population size
		migration_matrix = [
			[0, 0, 0],
			[0, 0, 0],
			[0, 0, 0]
			],
		demographic_events = [
			msprime.MigrationRateChange(
				time = ((time_admix - 60)/generation_time),
				rate = 0
				),
			msprime.MigrationRateChange(
				time = T_Admix, 
				rate = 1.5*prop_pop1,
				matrix_index = [2,0]
				),
			msprime.MigrationRateChange(
				time = T_Admix,
				rate = 1.5*prop_pop2,
				matrix_index = [2,1]
				),
			msprime.PopulationParametersChange(
				time = T_Admix, 
				initial_size = s_Pop3_0, 
				population_id = 2
				),
			msprime.MassMigration(
				time=T_Admix, 
				source = 2,
				dest = 1, 
				proportion = prop_pop2
				),
			msprime.MassMigration(
				time=T_Admix, 
				source = 2, 
				dest = 0, 
				proportion = prop_pop1
				),
	 		msprime.PopulationParametersChange(
	 			time = T_Admix, 
	 			growth_rate = r_Pop2_2,
	 			population_id = 1
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Admix, 
	 			growth_rate = 0, 
	 			population_id = 2
	 			),
	 		msprime.MigrationRateChange(
	 			time=T_Admix, 
	 			rate = 0
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop2, 
	 			initial_size = s_Pop2_0,
	 			population_id = 1
	 			),
	 		msprime.MassMigration(
	 			time=T_Pop2, 
	 			source = 1, 
	 			dest = 0, 
	 			proportion = 1.0
	 			),
			msprime.MigrationRateChange(
				time = T_Pop2, 
				rate = 1
				),
			msprime.PopulationParametersChange(
				time = T_Pop2, 
				growth_rate = 0, 
				population_id = 2
				),
			msprime.PopulationParametersChange(
				time = T_Pop2, 
				growth_rate = 0, 
				population_id = 1
				),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop1,
	 			population_id = 0
	 			),
	 		msprime.MigrationRateChange(
	 			time = T_Pop1, 
	 			rate = 1
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop3, 
	 			population_id = 2
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop2_1, 
	 			population_id = 1
	 			)
	 			]
	 			)
	 	 		
	print("SIMULATION COMPLETE", flush = True)

	#Make VCF of simulated data
	with open("snps_" + str(chrom) + ".vcf", "w") as vcf_file: 
		sim.write_vcf(vcf_file, 2)

	#msprime makes files that aren't quite compatible with what we need to do
	#Therefore we'll need to clean these files up before proceeding
	with open("snps_" + str(chrom) + ".vcf", "r") as file: 
		filedata = file.read()
		
	filedata = filedata.replace("tsk_0", "tsk_00")
	
	with open ("snps_" + str(chrom) + ".vcf", "w") as file: 
		file.write(filedata)
	
	print("VCF FIXED", flush = True)
	
	#Make files compatible for plink	
	os.system("../Dependencies/pop_info_generator.py " + str(sample_pop1) + ' ' + 
	str(sample_pop2) + ' ' + str(sample_pop3))
		
	make_files = subprocess.Popen(
		"plink --vcf snps_" + str(chrom) + ".vcf --make-bed --out model_" + str(chrom), 
		shell=True
		)

	#Create population information for the simulated data
	print("POP INFO CREATED", flush = True)
	
	make_files.communicate()
	print("FILES CREATED", flush = True)


def model_admix_expansion(pop1, pop2, pop3, time_admix, prop_pop1, prop_pop2, chrom, sample_pop1, sample_pop2, sample_pop3):

	
	#defining the variables for the simulation by scaling to args			
	Pop1 = pop1
	Pop2 = pop2
	Pop3 = pop3
	generation_time = 30
	T_Pop1 = 10000000/generation_time  #pygmies
	T_Pop2 = 5000000/generation_time  #bantu
	T_Admix = time_admix/generation_time
	r_Pop2_1 = (pop2/1000)**(1/T_Admix) - 1
	r_Pop3 = (pop3/100)**(1/T_Admix) - 1
	r_Pop2_2 = (0.5*pop2 / 20)**(1 / (T_Pop2 - T_Admix)) - 1
	r_Pop1 = (pop3/100)**(1/T_Pop1)-1
	s_Pop2_0 = 500
	s_Pop3_0 = 100
	m_Pop1 = 0.1*prop_pop1
	m_Pop2 = 0.1*prop_pop2
	
	chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','18',
	'19','20','21','22']
	
	length = ['248956422','242193529','198295559','190214555','181538259','170805979',
	'159345973','145138636','138394717','133797422','135086622','133275309','114364328',
	'107043718','101991189','90338345','83257441','80373285','58617616','64444167',
	'46709983','50818468']
	
	recomb_rate = ['1.14856e-08','1.10543e-08', '1.12796e-08','1.12312e-08','1.12809e-08',
	'1.12229e-08','1.17646e-08','1.14785e-08','1.17807e-08','1.33651e-08','1.17193e-08',
	'1.30502e-08','1.09149e-08','1.11973e-08','1.38358e-08','1.48346e-08','1.58249e-08',
	'1.5076e-08','1.82201e-08','1.71783e-08','1.30452e-08','1.4445e-08']
	
	chrom_tuple = list(zip(chroms, length, recomb_rate))
	chrom_map = pd.DataFrame(chrom_tuple, columns = ['chroms', 'length', 'recomb_rate'])
	
	chrom_map_pos = (1-int(chrom))
	
	chrom_recomb_rate = chrom_map.iloc[chrom_map_pos]['recomb_rate']
	chrom_length = chrom_map.iloc[chrom_map_pos]['length']
	
	print("BEGINNING SIMULATION:", str(sys.argv), flush = True)

	#begin simulation
	
	sim = msprime.simulate(
		length = int(chrom_length), 
		recombination_rate = float(chrom_recomb_rate),
		mutation_rate = 1.29e-8, #human mutation rate
		population_configurations = [
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop1)), 
			initial_size = Pop1,
			growth_rate = r_Pop1
			), 
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop2)), 
			initial_size = Pop2, 
			growth_rate = r_Pop2_1
			), 
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop3)), 
			initial_size = Pop3, 
			growth_rate = r_Pop3 #pop3 is growing in the recent past
			)],
		migration_matrix = [
			[0, 0, 0],
			[0, 0, 0],
			[0, 0, 0]
			],
		demographic_events = [
			msprime.MigrationRateChange(
				time = ((time_admix - 60)/generation_time),
				rate = 0
				),
			msprime.MigrationRateChange(
				time = T_Admix, 
				rate = 1.5*prop_pop1,
				matrix_index = [2,0]
				),
			msprime.MigrationRateChange(
				time = T_Admix,
				rate = 1.5*prop_pop2,
				matrix_index = [2,1]
				),
			msprime.PopulationParametersChange(
				time = T_Admix, 
				initial_size = s_Pop3_0, 
				population_id = 2
				),
			msprime.MassMigration(
				time=T_Admix, 
				source = 2,
				dest = 1, 
				proportion = prop_pop2
				),
			msprime.MassMigration(
				time=T_Admix, 
				source = 2, 
				dest = 0, 
				proportion = prop_pop1
				),
	 		msprime.PopulationParametersChange(
	 			time = T_Admix, 
	 			growth_rate = r_Pop2_2,
	 			population_id = 1
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Admix, 
	 			growth_rate = 0, 
	 			population_id = 2
	 			),
	 		msprime.MigrationRateChange(
	 			time=T_Admix, 
	 			rate = 0
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop2, 
	 			initial_size = s_Pop2_0,
	 			population_id = 1
	 			),
	 		msprime.MassMigration(
	 			time=T_Pop2, 
	 			source = 1, 
	 			dest = 0, 
	 			proportion = 1.0
	 			),
			msprime.MigrationRateChange(
				time = T_Pop2, 
				rate = 1
				),
			msprime.PopulationParametersChange(
				time = T_Pop2, 
				growth_rate = 0, 
				population_id = 2
				),
			msprime.PopulationParametersChange(
				time = T_Pop2, 
				growth_rate = 0, 
				population_id = 1
				),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop1,
	 			population_id = 0
	 			),
	 		msprime.MigrationRateChange(
	 			time = T_Pop1, 
	 			rate = 1
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop3, 
	 			population_id = 2
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop2_1, 
	 			population_id = 1
	 			)
	 			]
	 			)
	 	 		
	print("SIMULATION COMPLETE", flush = True)

	#Make VCF of simulated data
	with open("snps_" + str(chrom) + ".vcf", "w") as vcf_file: 
		sim.write_vcf(vcf_file, 2)

	#msprime makes files that aren't quite compatible with what we need to do
	#Therefore we'll need to clean these files up before proceeding
	with open("snps_" + str(chrom) + ".vcf", "r") as file: 
		filedata = file.read()
		
	filedata = filedata.replace("tsk_0", "tsk_00")
	
	with open ("snps_" + str(chrom) + ".vcf", "w") as file: 
		file.write(filedata)
	
	print("VCF FIXED", flush = True)
	
	#Make files compatible for plink	
	os.system("../Dependencies/pop_info_generator.py " + str(sample_pop1) + ' ' + 
	str(sample_pop2) + ' ' + str(sample_pop3))
		
	make_files = subprocess.Popen(
		"plink --vcf snps_" + str(chrom) + ".vcf --make-bed --out model_" + str(chrom), 
		shell=True
		)

	#Create population information for the simulated data
	print("POP INFO CREATED", flush = True)
	
	make_files.communicate()
	print("FILES CREATED", flush = True)



def model_admix_collapse(pop1, pop2, pop3, time_admix, prop_pop1, prop_pop2, chrom, sample_pop1, sample_pop2, sample_pop3):

	
	#defining the variables for the simulation by scaling to args			
	Pop1 = pop1
	Pop2 = pop2
	Pop3 = pop3
	generation_time = 30
	T_Pop1 = 10000000/generation_time  #pygmies
	T_Pop2 = 5000000/generation_time  #bantu
	T_Admix = time_admix/generation_time
	r_Pop2_1 = (pop2/1000)**(1/T_Admix) - 1
	r_Pop3 = (pop3/100)**(1/T_Admix) - 1
	r_Pop2_2 = (0.5*pop2 / 20)**(1 / (T_Pop2 - T_Admix)) - 1
	r_Pop1 = (pop3/100)**(1/T_Pop1)-1
	s_Pop2_0 = 500
	s_Pop3_0 = 100
	m_Pop1 = 0.1*prop_pop1
	m_Pop2 = 0.1*prop_pop2
	
	chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','18',
	'19','20','21','22']
	
	length = ['248956422','242193529','198295559','190214555','181538259','170805979',
	'159345973','145138636','138394717','133797422','135086622','133275309','114364328',
	'107043718','101991189','90338345','83257441','80373285','58617616','64444167',
	'46709983','50818468']
	
	recomb_rate = ['1.14856e-08','1.10543e-08', '1.12796e-08','1.12312e-08','1.12809e-08',
	'1.12229e-08','1.17646e-08','1.14785e-08','1.17807e-08','1.33651e-08','1.17193e-08',
	'1.30502e-08','1.09149e-08','1.11973e-08','1.38358e-08','1.48346e-08','1.58249e-08',
	'1.5076e-08','1.82201e-08','1.71783e-08','1.30452e-08','1.4445e-08']
	#chromosome lengths and recombination rates from stdpopsim catalogue
	
	chrom_tuple = list(zip(chroms, length, recomb_rate))
	chrom_map = pd.DataFrame(chrom_tuple, columns = ['chroms', 'length', 'recomb_rate'])
	
	chrom_map_pos = (1-int(chrom))
	
	chrom_recomb_rate = chrom_map.iloc[chrom_map_pos]['recomb_rate']
	chrom_length = chrom_map.iloc[chrom_map_pos]['length']
	
	print("BEGINNING SIMULATION:", str(sys.argv), flush = True)

	#begin simulation
	sim = msprime.simulate(
		length = int(chrom_length), 
		recombination_rate = float(chrom_recomb_rate),
		mutation_rate = 1.29e-8, #human mutation rate
		population_configurations = [
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop1)), 
			initial_size = Pop1,
			growth_rate = r_Pop1
			), 
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop2)), 
			initial_size = Pop2, 
			growth_rate = r_Pop2_1
			), 
		msprime.PopulationConfiguration(
			sample_size = (2*int(sample_pop3)), 
			initial_size = Pop3, 
			growth_rate = r_Pop3
			)],
		migration_matrix = [
			[0, 0, 0],
			[0, 0, 0],
			[0, 0, 0]
			],
		demographic_events = [
			msprime.PopulationParametersChange(
				time = (time_admix - 5820)/generation_time,
				population_id = 2,
				growth_rate = r_Pop3), #pop3 grows back from the bottleneck
			msprime.PopulationParametersChange(
				time = (time_admix - 5760)/generation_time,
				population_id = 2,
				initial_size = 100), #pop3 experienced a bottleneck 8 generations ago that lasted for 2 generations
			msprime.MigrationRateChange(
				time = ((time_admix - 60)/generation_time),
				rate = 0
				),
			msprime.MigrationRateChange(
				time = T_Admix, 
				rate = 1.5*prop_pop1,
				matrix_index = [2,0]
				),
			msprime.MigrationRateChange(
				time = T_Admix,
				rate = 1.5*prop_pop2,
				matrix_index = [2,1]
				),
			msprime.PopulationParametersChange(
				time = T_Admix, 
				initial_size = s_Pop3_0, 
				population_id = 2
				),
			msprime.MassMigration(
				time=T_Admix, 
				source = 2,
				dest = 1, 
				proportion = prop_pop2
				),
			msprime.MassMigration(
				time=T_Admix, 
				source = 2, 
				dest = 0, 
				proportion = prop_pop1
				),
	 		msprime.PopulationParametersChange(
	 			time = T_Admix, 
	 			growth_rate = r_Pop2_2,
	 			population_id = 1
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Admix, 
	 			growth_rate = 0, 
	 			population_id = 2
	 			),
	 		msprime.MigrationRateChange(
	 			time=T_Admix, 
	 			rate = 0
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop2, 
	 			initial_size = s_Pop2_0,
	 			population_id = 1
	 			),
	 		msprime.MassMigration(
	 			time=T_Pop2, 
	 			source = 1, 
	 			dest = 0, 
	 			proportion = 1.0
	 			),
			msprime.MigrationRateChange(
				time = T_Pop2, 
				rate = 1
				),
			msprime.PopulationParametersChange(
				time = T_Pop2, 
				growth_rate = 0, 
				population_id = 2
				),
			msprime.PopulationParametersChange(
				time = T_Pop2, 
				growth_rate = 0, 
				population_id = 1
				),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop1,
	 			population_id = 0
	 			),
	 		msprime.MigrationRateChange(
	 			time = T_Pop1, 
	 			rate = 1
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop3, 
	 			population_id = 2
	 			),
	 		msprime.PopulationParametersChange(
	 			time = T_Pop1, 
	 			growth_rate = r_Pop2_1, 
	 			population_id = 1
	 			)
	 			]
	 			)
	 	 		
	print("SIMULATION COMPLETE", flush = True)

	#Make VCF of simulated data
	with open("snps_" + str(chrom) + ".vcf", "w") as vcf_file: 
		sim.write_vcf(vcf_file, 2)

	#msprime makes files that aren't quite compatible with what we need to do
	#Therefore we'll need to clean these files up before proceeding
	with open("snps_" + str(chrom) + ".vcf", "r") as file: 
		filedata = file.read()
		
	filedata = filedata.replace("tsk_0", "tsk_00")
	
	with open ("snps_" + str(chrom) + ".vcf", "w") as file: 
		file.write(filedata)
	
	print("VCF FIXED", flush = True)
	
	#Make files compatible for plink	
	os.system("../Dependencies/pop_info_generator.py " + str(sample_pop1) + ' ' + 
	str(sample_pop2) + ' ' + str(sample_pop3))
		
	make_files = subprocess.Popen(
		"plink --vcf snps_" + str(chrom) + ".vcf --make-bed --double-id --out model_" + str(chrom), 
		shell=True
		)

	#Create population information for the simulated data
	print("POP INFO CREATED", flush = True)
	
	make_files.communicate()
	print("FILES CREATED", flush = True)
	
#We'll need to do a bit more file prep before we're ready to get moving on analysis
def fam_fix():
	fam_fix = subprocess.Popen(
		"../Dependencies/fam_fix.pl ../Admixture/model_" + str(chrom) + ".fam ../Admixture/model_" + str(chrom) + "_fixed.fam", 
		shell=True
		)	
	fam_fix.communicate()
	print("FAM FIXED", flush = True) 

	os.rename("model_" + str(chrom) + "_fixed.fam", "model_" + str(chrom) + ".fam")
	print("FAM REPLACED", flush = True)


#Add SNP IDs to the VCF
def bim_fix():
	bim_fix = subprocess.Popen(
			"../Dependencies/bim_fix.py ../Admixture/model_" + str(chrom) + " ../Admixture/model_" + str(chrom) + "_fixed.bim",
			shell=True
			)
	bim_fix.communicate()
	os.system('mv ../Admixture/model_' + str(chrom) + '_fixed.bim ../Admixture/model_' + str(chrom) + '.bim')
	print("BIM FIXED", flush = True)

def new_vcf():
	new_vcf = subprocess.Popen(
		"plink --bfile model_" + str(chrom) + " --recode vcf-iid --double-id --out snps_" + str(chrom), 
		shell = True
		)
	new_vcf.communicate()
	print("NEW VCF CREATED", flush = True) 

#def snp_id():
#	snp_id = subprocess.Popen(
#		"../Dependencies/snp_id.py",
#		shell=True
#		)

#Finally ready for plink	
#Perform pca and plot
def pca_test():
	pca = subprocess.Popen(
		"plink --bfile model_" + str(chrom) + " --pca --out model_" + str(chrom) + "_pca", 
		shell=True
		)
	pca.communicate()
	
	#plot = subprocess.Popen(
	#	"Rscript ../Dependencies/pca_plot.R", 
	#	shell=True
	#	)
	#plot.communicate()
	print("PCA FINISHED", flush = True)

#Prune the dataset
def prune():	
	prune = subprocess.Popen(
		"plink --bfile model_" + str(chrom) + " --indep-pairwise 50 2 0.8 --double-id --out " + str(chrom),
		shell=True
		)
	prune.communicate()

#New dataset for pruned beds
def make_beds():
	new_bed = subprocess.Popen(
		"plink --bfile model_" + str(chrom) + " --extract " + str(chrom) + ".prune.in --make-bed --double-id --out pruned_model_" + str(chrom), 
		shell=True
		)
	new_bed.communicate()
	
	new_bedout = subprocess.Popen(
		"plink --bfile model_" + str(chrom) + " --extract " + str(chrom) + ".prune.out --make-bed --double-id --out outpruned_model_" + str(chrom), 
		shell=True
		)
	new_bedout.communicate()

def prune_mp():
	m_p = subprocess.Popen(
		"plink --bfile pruned_model_" + str(chrom) + " --recode --double-id --out pruned_model_" + str(chrom), 
		shell=True
		)
	m_p.communicate()
	
#run ADMIXTURE
def admixture_test():
	admix = subprocess.Popen(
		"for K in 1 2 3; \
		do admixture --cv pruned_model_" + str(chrom) + ".bed $K | tee log${K}_" + str(chrom) + ".out; done", 
		shell=True
		)
	admix.communicate()
	
	cv_grab = subprocess.Popen(
		"grep -h CV log*_" + str(chrom) + ".out > cv_error_" + str(chrom) + ".txt",
		shell = True
		)
	cv_grab.communicate()
	
	#cv_error = subprocess.Popen(
	#	"Rscript ../Dependencies/cv_error_plot.R",
	#	shell=True
	#	)
	#cv_error.communicate()
	
	#admix_plot = subprocess.Popen(
	#	"Rscript ../Dependencies/admix_plot.R", 
	#	shell = True
	#	)
	#admix_plot.communicate()
	
	print("ADMIXTURE FINISHED", flush = True)

#MAF removal
def freq():
	freq = subprocess.Popen(
		"plink --bfile model_" + str(chrom) + " --maf 0.05 --double-id --make-bed --out sims_" + str(chrom),
		shell=True
		)
	freq.communicate()

def main(pop1, pop2, pop3, time_admix, prop_pop1, prop_pop2, chrom, dem_option, sample_pop1, sample_pop2, sample_pop3):
	if (len(sys.argv) >= 13):
		print("Too many arguments specified", flush = True)
		exit()
	elif (len(sys.argv) < 12):
		print("Too few arguments specified", flush = True)
		exit()
	else:
		pass
		
	#for file in os.listdir("."):
	#	if file.endswith(".pdf"):
	#		print("Clear working directory before simulating new data", flush = True) 
	#		exit()
	#	else:
	#		pass 

	if sys.argv[8] == 'constant':
		model_admix_constant(pop1, pop2, pop3, time_admix, prop_pop1, prop_pop2, chrom, 
		sample_pop1, sample_pop2, sample_pop3)
	elif sys.argv[8] == 'collapse':
		model_admix_collapse(pop1, pop2, pop3, time_admix, prop_pop1, prop_pop2, chrom, 
		sample_pop1, sample_pop2, sample_pop3)
	elif sys.argv[8] == 'expansion':
		model_admix_expansion(pop1, pop2, pop3, time_admix, prop_pop1, prop_pop2, chrom, 
		sample_pop1, sample_pop2, sample_pop3)
	else:
		sys.exit('Did you specify constant, collapse, or expansion models?', flush = True)
		
	fam_fix()
	bim_fix()
	new_vcf()
	#snp_id()
	pca_test()
	prune()
	make_beds()
	prune_mp()
	admixture_test()
	freq()


if __name__ == '__main__':
	main(pop1, pop2, pop3, time_admix, prop_pop1, prop_pop2, chrom, dem_option, sample_pop1, sample_pop2, sample_pop3)
