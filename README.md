# phyOverlap
Method for testing for phylogenetic overlap of mutations in a gene and a phenotype of interest.

phyOverlap is written for use in python 2.7.X

Required outside python packages:
	
	dendropy v3.12
	Biopython
	Scipy
	sys
	multiprocessing
	
Python packages included:
	treeFunctions

If you install everything correctly, you should be able to navigate into the Example_Input_Files folder, UNZIP THE .FAS FILE, run the following command, and generate the p-values used in Figure 1C

python ../phyOverlap.py -s INH_sens.txt -r INH_res.txt -f 549str_all_SNPs.filter95.Q.fas -t 549strain_Tree.nwk -i Sites_Steps_Ancestors.txt -g Rv_000962.3.features.txt -S Synon_Sites.txt -o Test.txt -p 4
