#! /usr/bin/env python

import treeFunctions
import dendropy
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import scipy.stats
import sys
import multiprocessing as mp
import optparse
	
if __name__ == "__main__":

	usage="usage: %prog [options]"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option("-s", "--sensitive", dest="sensitive",
			help="List of taxa without attribute",
			type="string")
        parser.add_option("-r", "--resistant", dest="resistant",
                        help="List of taxa with attribute",
                        type="string")
        parser.add_option("-f", "--fasta", dest="fasta",
                        help="SNP file in FASTA format. Names must match tree. Unknown positions M or ?",
                        type="string")
        parser.add_option("-t", "--tree", dest="inTree",
                        help="Phylogenetic tree in nwk format. Rooted.",
                        type="string")
        parser.add_option("-i", "--siteInfo", dest="site_info",
                        help="File containing one row for each SNP position detailing genomic site, number of steps on tree, and ancestral state",
                        type="string")
        parser.add_option("-g", "--genes", dest="geneFile",
                        help="List of genomic intervals to test, see README for description",
                        type="string")
        parser.add_option("-S", "--synonymous", dest="synon_sites",
                        help="List of sites where mutations are synonymous",
                        type="string")
        parser.add_option("-o", "--output", dest="outfile",
                        help="The name of an output file locations",
                        type="string")
        parser.add_option("-p", "--processors", dest="nProc",
                        help="The name of an output file locations",
                        type="int", default=1)
	options, args = parser.parse_args()

	mandatories = ["sensitive", "resistant", "fasta", "inTree", "site_info", "geneFile", "synon_sites", "outfile"]
	for m in mandatories:
		if not getattr(options, m, None):
			print "\nMust provide %s.\n" %m
			parser.print_help()
			exit(-1)
	
	mySensTaxa = treeFunctions.getSampleNames(options.sensitive)
	myResTaxa = treeFunctions.getSampleNames(options.resistant)
	myTaxa = mySensTaxa + myResTaxa
	
	#Getting from file now of ancestral reconstructions
	#refStates = treeFunctions.getMajorAlleles(sys.argv[3])
	
	print 'gotRefs'
	

	
	hashFasta = SeqIO.to_dict(SeqIO.parse(open(options.fasta),"fasta"))

	iSites = []
	iSteps = []
	refStates = []
	
	pNames, pAlts, pCoords = treeFunctions.getProteinSimple(options.geneFile)
	
	for strLine in open(options.site_info):
		astrLine = strLine.strip().split('\t')
		
		iSites.append(int(astrLine[0]))
		iSteps.append(int(astrLine[1]))
		refStates.append(astrLine[2])
	#	MinorAllele.append(int(astrLine[2]))


	synonSNPs = []
	
	for inline in open(options.synon_sites):
		synonSNPs.append(int(inline.strip()))

	prot_dict = {}

	for i, prot in enumerate(pNames):
		coords = str(pCoords[i])
		coordList = coords.split('.')
		start = int(coordList[0])
		stop = int(coordList[2])
		
		snpNumList = []
		siteList = []
		stepsList = []
		refsList = []
		numAmbiguous = 0
		for iNum, iSite in enumerate(iSites):
				if int(iSite) >= start and int(iSite) <= stop:
					if iSites[iNum] in synonSNPs:
						continue
					if len(refStates[iNum]) > 1:
						numAmbiguous += 1
						continue	
					if iSteps[iNum] == 0:
						continue		
					siteList.append(iSites[iNum])
					snpNumList.append(iNum)
					stepsList.append(iSteps[iNum])
					refsList.append(refStates[iNum])
		prot_dict[prot] = [snpNumList, stepsList, refsList, siteList, numAmbiguous]
		
	print prot_dict[pNames[0]]

	actualTestVals = []
	proteinSiteNums = []
	proteinStepNums = []
	for prot in pNames:	
		pSNPnums, pSteps, pRefs = prot_dict[prot][:3]
		
		if sum(pSteps) == 0:
			actualTestVals.append('blah')

			prot_dict[prot] = False
			continue
			
		pTestStats = []
		print prot	
		numTested = 0
		for pNum, pSNPnum in enumerate(pSNPnums):
	
			numSens_major = 0
			numRes_major = 0
			numSens_minor = 0
			numRes_minor = 0 
			numSteps = pSteps[pNum] 
			
#			if numSteps == 0: 
#				continue
			numTested += 1
			for strTaxaSens in mySensTaxa:
				myCall = hashFasta[strTaxaSens].seq[pSNPnum]
				if myCall == pRefs[pNum] or myCall == '?' or myCall == 'M':
					numSens_major += 1
				else:
					numSens_minor += 1
				
			for strTaxaRes in myResTaxa:
				myCall = hashFasta[strTaxaRes].seq[pSNPnum]
				if myCall == pRefs[pNum] or myCall == '?' or myCall == 'M':
					numRes_major += 1
				else:
					numRes_minor += 1	

			testStat = float(numRes_minor) / (numRes_minor + numSens_minor)
			testStatList = [testStat] * numSteps
			pTestStats += testStatList
		
		proteinSiteNums.append(numTested)
		proteinStepNums.append(len(pTestStats))
		pTestComposite = sum(pTestStats) / float(len(pTestStats))
		actualTestVals.append(pTestComposite)


	ntests = 0

	outReal = open(options.outfile, 'w')

	tree = dendropy.Tree.get_from_path(options.inTree, 'newick')
	cleanTree = treeFunctions.trimTree(tree, myTaxa)
	treeFunctions.labelTree(cleanTree)
	nodes, weights = treeFunctions.buildWeightedArray(cleanTree)	
	sensNodes, resNodes, nodeDict = treeFunctions.getNodePhenoLists(cleanTree, myResTaxa)
	print len(resNodes)
	print len(actualTestVals)
	print len(pNames)
	
	myPool = mp.Pool(processes=options.nProc)
	
	jobInputs = []
	
	proteinsToWrite = []
	
	for i,prot in enumerate(pNames):
		if prot_dict[prot]:
			actualOR = actualTestVals[i]
			stepsVector = prot_dict[prot][1]
			numSkipped = prot_dict[prot][4]
			proteinsToWrite.append([prot, actualOR])
			jobInputs.append([nodes, weights, stepsVector, 50000, cleanTree, nodeDict, sensNodes, resNodes, actualOR])
		
	myOutputs = myPool.map(treeFunctions.newScoringAlgorithm_vFourPntOne_mp, jobInputs)
	myPool.close()
	myPool.join()

	for i in range(len(myOutputs)):
		pInfo = proteinsToWrite[i]
		statInfo = myOutputs[i]
		sites = [proteinSiteNums[i]]
		steps = [proteinStepNums[i]]
		myLine = pInfo + sites + steps + statInfo
		outReal.write('\t'.join(map(str, myLine)) + '\n')
	outReal.close()
