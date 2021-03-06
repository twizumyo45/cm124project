import numpy as np
import scipy.stats as stat 
import math
from copy import copy, deepcopy

#variable parameters
numRows = 50
numCols = 7
numParents = 3 #start counting from 0; set to be about %% of number of individuals

#function to count the number of genotypes this phase pair explains
def explainsGenotypes(phase1, phase2, genotypes, alrdyExplained):
    maxExplainedCurr = 0
    #get the genotype this haplotype pair is 
    pairsGenotype = np.random.random_integers(0,1,numCols)
    for l in range(0, numCols):
        if (phase1[l] == phase2[l]):
            pairsGenotype[l] = phase1[l]
        else:
            pairsGenotype[l] = 2
            
    #loop through and count the number of genotypes it explains
    for i in range(0, numRows/2):
        if (alrdyExplained[i] == False and (pairsGenotype == genotypes[i]).all()):
            maxExplainedCurr += 1
    
    return maxExplainedCurr            

#function to translate number into binary (with array of 0/1's)
#m is number of digits
#note output is displayed in reverse digit place
def toBinary(num, m):
    binary = np.random.random_integers(0,1,numCols)
    for i in range (m-1, -1, -1):
        toSubtract = 2 ** i
        if (toSubtract <= num):
            num = num - toSubtract
            binary[i] = 1
        else:
            binary[i] = 0
    return binary

#helper function for hapCount
def hapContains (dataInRow, hapList):
    for i in range (0, numRows):
        if (dataInRow == hapList[i]).all():
            return i
    return -1

#set of haplotypes + counts (first slot is the count, rest is the hap rep)
def hapCount (dataIn):
    countList = [0]*numRows
    hapList = [[5]*numCols]*numRows
    currSize = 0
    #go thru and add to counts if alrdy exists, or else create
    for i in range (0, numRows):
        ind = hapContains(dataIn[i], hapList)
        if (ind != -1):
            countList[ind] += 1
        else:
            hapList[currSize] = dataIn[i]
            countList[currSize] += 1
            currSize += 1
    print countList
    print hapList  
    print currSize  
    

#build genotypes from true haplotypes function
def getGenotypes (dataIn):
    genotypes = np.random.random_integers(0,1,(numRows/2, numCols))
    for r in range(0, numRows/2):
        for l in range(0, numCols):
            if dataIn[2*r][l] == dataIn[2*r+1][l]:
                genotypes[r][l] = dataIn[2*r][l]
            else:
                genotypes[r][l] = 2
    return genotypes

#generator function 
#[0, i] is range of indexes to draw from
def buildHaplotypes (dataIn, i):
    for k in range(i/2+1, numRows/2):
        toCopy = np.random.random_integers(0,i,2)
        dataIn[2*k] = dataIn[toCopy[0]]
        dataIn[2*k+1] = dataIn[toCopy[1]]

#each row is a person, each col is a SNP
#so 1000 people and 15 SNPs
dataIn = np.random.random_integers(0,1,(numRows, numCols))


buildHaplotypes(dataIn, numParents)

#create genotypes from true haplotypes
genotypes = getGenotypes(dataIn)

#ALGORITHM START
print "----initializing----"
currIndex = 0
alrdyExplained = [False]*(numRows/2)
hapSetSol = [[5]*numCols]*numRows

#initialize like clarcks - for all unambiguous genotypes, add to the haplotype list, and rmv from the genotype list
#remember also rmv genotypes it might explain
for i in range(0, numRows/2):
    if (not 2 in genotypes[i]):
        alrdyExplained[i] = True
        hapSetSol[2*i] = genotypes[i]
        hapSetSol[2*i+1] = genotypes[i]
        currIndex += 1

print "-----after init-----"
print alrdyExplained

#start greedy portion
numIterations = 0
#do while there are still genotypes to be explained
while (False in alrdyExplained):
    #get the solution with the max number of explanations
    maxExplainedSoFar = -1
    maxPhase1 = np.random.random_integers(0,1,numCols)
    maxPhase2 = np.random.random_integers(0,1,numCols)
    #loop through all possible 2^m haplotype phases
    for i in range(0, (2**numCols) - 1):
        currPhase1 = toBinary(i, numCols)
        for k in range(i, (2**numCols) - 1):
            currPhase2 = toBinary(k, numCols)
            #for each possible solution phase pair, count all the number of genotypes it might explain
            maxExplainedCurr = explainsGenotypes(currPhase1, currPhase2, genotypes, alrdyExplained)
            if (maxExplainedSoFar < maxExplainedCurr):
                maxPhase1 = currPhase1
                maxPhase2 = currPhase2
                maxExplainedSoFar = maxExplainedCurr
                
                
    #add to the haplotype list
    #remove explained genotypes from the genotype list
    #get the genotype this haplotype pair is 
    pairsGenotype = np.random.random_integers(0,1,numCols)
    for l in range(0, numCols):
        if (maxPhase1[l] == maxPhase2[l]):
            pairsGenotype[l] = maxPhase1[l]
        else:
            pairsGenotype[l] = 2
            
    #remove explained genotypes from the genotype list
    for i in range(0, numRows/2):
        if (alrdyExplained[i] == False and (pairsGenotype == genotypes[i]).all()):
            alrdyExplained[i] = True
            hapSetSol[2*i] = maxPhase1
            hapSetSol[2*i+1] = maxPhase2
    
    print "-----new iteration-----"
    print alrdyExplained
    print maxPhase1
    print maxPhase2
    print pairsGenotype
    print maxExplainedSoFar
    numIterations += 1

print "-----DONE-----"
print numIterations
print hapCount(hapSetSol)