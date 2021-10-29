# model_admix
Modeling admixture

This is a python program I put together during my MPhil that uses msprime to simulate whole-genome data.
It's a bit clunky and uses msprime0.7, which is deprecated, but the developers of msprime are allowing
backwards-compatibility indefinitely. 

You can use the program as follows:

python model_admix.py $pop1_initial_size $pop2_initial_size $pop3_initial_size $time_of_admixture_in_years $proportion_of_pop1_ancestry $proportion_of_pop2_ancestry $chromosome $demographic_model_option $sample_pop1 $sample_pop2 $sample_pop3
  
where: \
  pop1 initial size: (integer) size of ancestor 1 at the beginning of its existence (time designated in the program) \
  pop2 initial size: (integer) size of ancestor 2 at the beginning of its existence \
  pop3 initial size: (integer) size of admixed population at beginning of its existence \
  time of admixture in years: (integer) time at which the admixture event occurred \
  proportion of pop1 ancestry: (float) proportion of population 1 ancestry found in admixed population \
  proportion of pop2 ancestry: (float) proportion of population 2 ancestry found in admixed population \
  chromosome: (integer) specify the chromosome to set recombination rate and chromosome length \
  demographic model option: (string) collapse, expansion, or constant \
  sample_pop1: sample size desired pop1 \
  sample_pop2: sample size desired pop2 \
  sample_pop3: sample size desired pop3 

And all arguments are required.

Each demographic model accounts for one admixture event and two generations of ongoing migration (0.1*proportion of ancestor) before removing migration between the three populations. The constant population size model does not specify a population growth rate for the admixed population. The collapse model specifies a bottleneck in the admixed population for two generations (currently set at 8 generations ago provided the time to admixture is 6000 years ago -- the user may change this where desired (edit the time parameters in lines 540-545). The population growth model allows the admixed population to grow after the time of admixture according to the following formula: (pop3/100)**(1/T_Admix) - 1

The script contains dataframes for each human chromosome, but these can easily be edited to be for a different species. Just
change the pandas dataframes to whatever species' info you need. 

Be sure to also check things like initial pop size and the timing of the population collapse parameter to make sure it works for your desired model.

Dependencies: \
pandas
pandas_plink
msprime

*remember to download the dependencies folder -- msprime 0.x wasn't the greatest for VCFs that were compatible with downstream programs, so these scripts will help clean those up and run ADMIXTURE and do a quick PCA on your simulation*

*run the program in the same directory as your data!*
