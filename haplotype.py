import numpy as np
import scipy.stats as stat 
import math
from copy import copy, deepcopy

#breeding function
def breed (dataIn,i,j): 
    ''' i mates j, make person k '''
    personI = dataIn[i,:]
    personJ = dataIn[j,:]
    ''' cross over happens '''
    flip = personI == personJ # same allele will be passed on. 
    flip = np.where(flip==False) # any position with different allele, the allele will be chosen randomly 
    personK = deepcopy(personI) ## WILL NOT WORK . if you do personK = personI without deepcopy 
#!! this is a bit incorrect, the true process of 'cross-over' is not simple fliping of bits, usually, "chunks" are cross-over 
#!! you can modify this process. 
    ''' much more complicated if alleles are not independent '''
    personK[flip] = np.random.binomial( 1,.5,np.shape(flip)[1] )
    return personK

#each row is a person, each col is a SNP
#so 1000 people and 15 SNPs
dataIn = np.random.random_integers(0,1,(1000,15))
