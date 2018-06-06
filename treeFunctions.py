#! /usr/bin/env python

"""
Tree functions for doing an in house phyC set of calculations
"""

import random
import dendropy
import scipy.stats
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import numpy



def getSampleNames(textFile):
        outNames = []
        for strLine in open(textFile):
                astrLine = strLine.strip()
                outNames.append(astrLine)
        return outNames


#Give a dendropy tree and a list of taxa on the nodes.  Trims all other taxa off the tree
#Everything in the taxon list should be spaces because dendropy annoying converts _ to space
def trimTree(tree, taxonList):
    taxa = tree.taxon_set
    prunedTree = dendropy.Tree(tree, taxon_set=taxa)
#    prunedTree.retain_taxa_with_labels(taxonList, update_splits=True)
    prunedTree.retain_taxa_with_labels(taxonList)
    return prunedTree

# Takes a tree object and labels all of the nodes
def labelTree(tree):
    nodeNum = 1
    for myNode in tree.preorder_node_iter():
        myNode.label = str(nodeNum)
        nodeNum +=1

# Builds a list of weights for each node in a tree.  Discards any with NO length
# (This occurs at the outgroup, only one branch has a length)
def buildWeightedArray(treeObject):
    nodeList = []
    nodeWeights = []
    for myNode in treeObject.preorder_node_iter():
        if myNode.edge_length:
            nodeList.append(myNode.label)
            nodeWeights.append(myNode.edge_length)
    return [nodeList, nodeWeights]

# Takes the lists from buildWeightedArray, a number of steps to put on tree and
# a number for iterating, and then makes a list of lists of node draws (with replacement)
def randomDrawMachine(nodes, weights, iSteps, numIters):
    aaDraws = []
    for iPermutation in range(numIters):
        aPermDraws = []
        for iStep in range(iSteps):
            rnd = random.random() * sum(weights)
            for i, w in enumerate(weights):
                rnd -= w
                if rnd < 0:
                    aPermDraws.append(nodes[i])
                    break        
        aaDraws.append(aPermDraws)
    return aaDraws

def getMajorAlleles_plus(multiFasta):
    hashFasta = SeqIO.to_dict(SeqIO.parse(open(multiFasta),"fasta"))
    
    majorAlleles = []
    minorAlleleFreq = []
    
    myKeys = hashFasta.keys()
    
    iFastaLength = len(hashFasta[myKeys[0]].seq)

    for iPos in range(iFastaLength):
        SNPcalls = []

        for myKey in myKeys:
            SNPcalls.append(hashFasta[myKey].seq[iPos])
        
        SNP_callsClean = [elem for elem in SNPcalls if elem != "?" and elem != "M"]
        callSet = set(SNP_callsClean)
        callDict = {}
        for baseCall in callSet:
            callDict[baseCall] = SNP_callsClean.count(baseCall)

        majorAllele = max(callDict.iteritems(), key=operator.itemgetter(1))[0]
        majorAlleles.append(majorAllele)
        numMajor = SNPcalls.count(majorAllele) + SNPcalls.count('?')
        numMinor = len(SNPcalls) - numMajor
        minorAlleleFreq.append(numMinor)
    
    return [majorAlleles, minorAlleleFreq]

def getProteinSimple(myGenesFile):
	with open(myGenesFile) as f:
		protein_Names = []
                protein_Alts = []
                protein_Coords = []
		for inline in f:
			info = inline.strip().split('\t')
			protein_Names.append(info[0])
			protein_Alts.append(info[1])
			protein_Coords.append(info[2])
		return (protein_Names, protein_Alts, protein_Coords) 
def getProteinInfo_fromJian(myGenesFile):
        with open(myGenesFile) as f:
                protein_Names = []
                protein_Alts = []
                protein_Coords = []

                for inline in f:
                        lineList = inline.strip().split('\t') 
      			protein_Name = lineList[3]
			if protein_Name in protein_Names:
				continue
			else:
				protein_Names.append(protein_Name)
	                        protein_Alt = lineList[6]
				if protein_Alt:
					protein_Alts.append(protein_Alt)
				else:
					protein_Alts.append(protein_Name)
        	                protein_Coords.append(lineList[4]) 
                return (protein_Names, protein_Alts, protein_Coords)

def getNodePhenoLists(cleanTree, resTaxa):
    resNodes = []
    sensNodes = []
    nodeDict = {}
    for myNode in cleanTree.leaf_iter():
        if str(myNode.taxon) in resTaxa:
            resNodes.append(myNode)
        else:
            sensNodes.append(myNode)
    for myNode in cleanTree.preorder_node_iter():
        nodeDict[myNode.label] = myNode
        
    return [sensNodes, resNodes, nodeDict]

def newScoringAlgorithm_vFourPntOne_mp(myData):
	nodes, weights, stepsVector, numIters, cleanTree, nodeDict, sensNodes, resNodes, actualOR = myData
	numSens = len(sensNodes)
	numRes = len(resNodes)
	halfSize = float(numSens + numRes)/2
	numTests = 0
	numExtreme = 0
    
	for iPermutation in range(numIters):
		numTests += 1
		aPermDraws = []
		for iSiteSteps in stepsVector:
			aSiteDraws = []
			drawCount = 0
#            for iStep in range(iSiteSteps)
			while drawCount < iSiteSteps and drawCount < 1000:
				rnd = random.random() * sum(weights)
				for i, w in enumerate(weights):
					rnd -= w
					if rnd < 0:
#						if nodes[i] not in aSiteDraws:
						aSiteDraws.append(nodes[i])
						drawCount += 1
						break
#						else:
#							break        
			aPermDraws.append(aSiteDraws)
            
		pTestVals = []
        
		totalSNPs = 0
		for SNPdraw in aPermDraws:

			snpPheno1Minor = []
			snpPheno2Minor = []
            
			if not SNPdraw:
				continue
			totalSNPs +=1
			for nodeDraw in SNPdraw:
				myNode = nodeDict[nodeDraw]
				for myTip in myNode.leaf_iter():
					if myTip in sensNodes:
						snpPheno1Minor.append(myTip)
					elif myTip in resNodes:
						snpPheno2Minor.append(myTip)
					else:
						raise Exception("shitty lookup")
            
			iPheno1Minor = len(set(snpPheno1Minor))
			iPheno2Minor = len(set(snpPheno2Minor))	            
			numMinor = iPheno1Minor + iPheno2Minor
			
			testStat = float(iPheno2Minor)/ numMinor
			testStatList = [testStat] * len(SNPdraw)
			pTestVals += testStatList
		overallStat = sum(pTestVals) / float(len(pTestVals))
		if overallStat >= actualOR:
			numExtreme += 1
		if numTests > 500:
			if float(numExtreme)/numTests > 0.10:
				break
				
	return [float(numExtreme)/numTests, numTests]
	
