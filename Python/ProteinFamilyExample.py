# Example python code for loading and matching sequence and contact map

import numpy
import pickle

# Data = pickle.load(open('MSATestProteins.pkl','rb')) # large array  (takes time to load) Data['PF00034.20']


# load small array of distance matrices 
Data1 = pickle.load(open('MSATestProteins1.pkl','rb'))
PF = Data1['PF03061.21'] # take one protein family

n_family = len(PF) # number of proteins in family (alignment)

P = PF[0] # take first protein in family
DistanceMat = P[0] # distance matrix (2d array)
NameChain = P[1]
UniprotName = P[2]
AlignmentInds = P[3]  # from 50 to 127 
AlignmentInds = AlignmentInds.split('-')  
AlignmentInds = [int(AlignmentInds[0]), int(AlignmentInds[1])] # convert to integers 

# load alignment
Align = numpy.load('Pfam/alignments/PF03061.21')
Align.Sequence[0] # view first protein in alignment 
num_proteins = len(Align.Sequence) # get # proteins 

Align.Sequence[NameChain] # entire sequence for first protein: 

# extract sequence and distance matrix for a specific protein: 'ENTH_ECOLI/50-127' in family: 'PF03061.21'
S = Align.Sequence[NameChain].replace('-', '') # entire sequence without gaps. Of length 78=127-50+1
D = DistanceMat[AlignmentInds[0]:AlignmentInds[1]+1,AlignmentInds[0]:AlignmentInds[1]+1]



