import numpy as np
import scipy.stats as stat 
import math
from copy import copy, deepcopy

numRows = 10
numCols =10

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
    

#breeding function
def breed (dataIn,i,j): 
    ''' i mates j, make person k '''
    personI = dataIn[i,:]
    personJ = dataIn[j,:]
    ''' cross over happens '''
    flip = personI == personJ # same allele will be passed on. 
#    flip = np.where(flip==False) # any position with different allele, the allele will be chosen randomly 
    personK = deepcopy(personI) ## WILL NOT WORK . if you do personK = personI without deepcopy 
#!! this is a bit incorrect, the true process of 'cross-over' is not simple fliping of bits, usually, "chunks" are cross-over 
#!! you can modify this process. 
#  ''' much more complicated if alleles are not independent '''
# personK[flip] = np.random.binomial( 1,.5,np.shape(flip)[1] )
    for m in range(0, numCols/5):
        for l in range(0, 5):
            IOrJ = np.random.random_integers(0,1,1)
            person = personJ
            if (IOrJ == 0):
                person = personI
            if (flip[m*5 + l] == False):
                personK[m*5 + l] = person[m*5 + l]
    return personK

#build genotype function
def getGenotypes (dataIn):
    genotypes = np.random.random_integers(0,1,(numRows/2, numCols))
    for r in range(0, numRows/2):
        for l in range(0, numCols):
            if dataIn[2*r][l] == dataIn[2*r+1][l]:
                genotypes[r][l] = dataIn[2*r][l]
            else:
                genotypes[r][l] = 2
    return genotypes

#each row is a person, each col is a SNP
#so 1000 people and 15 SNPs
dataIn = np.random.random_integers(0,1,(numRows, numCols))

#breed
for i in range(2, numRows):
    dataIn[i] = breed(dataIn, 0, 1)
    
#hapCount(dataIn)
genotypes = getGenotypes(dataIn)