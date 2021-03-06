import pandas as pd
import numpy as np
import pickle

import Bio.PDB
import numpy as np


amino = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX',
         'CYS', 'GLN', 'GLU', 'GLX', 'GLY',
         'HIS', 'ILE', 'LEU', 'LYS', 'MET',
         'PHE', 'PRO', 'SER', 'THR', 'TRP',
         'TYR', 'UNK', 'VAL']




def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""
    try:
        onecoord = residue_one['CA'].coord
    except KeyError:
        onecoord = np.random.choice(residue_one.child_list).coord


    try:
        twocoord = residue_two['CA'].coord
    except KeyError:
        twocoord = np.random.choice(residue_two.child_list).coord


    diff_vector = onecoord - twocoord

    return np.sqrt(np.sum(diff_vector * diff_vector))


def model(pdb_filename):
    """returns the model of a pdb file"""
    pdb_code = pdb_filename[0:4]
    structure = Bio.PDB.PDBParser().get_structure(pdb_code,pdb_filename)
    model = structure[0]
    return(model)


def calc_dist_matrix(chain):
    """Returns a matrix of C-alpha distances between two chains"""
    ids1 = list()
    for res in chain:
        if res.resname not in amino:
            ids1.append(res.id)
    for id in ids1:
        chain.detach_child(id)


    answer = np.zeros((len(chain), len(chain)), np.float)
    for i in range(len(chain.child_list)):
        for j in range(i,len(chain.child_list)):
            answer[i, j] = calc_residue_dist(chain.child_list[i], chain.child_list[j])
            answer[j, i] = answer[i, j]
    return answer



def contactmap(model,chain,MapThreshold = 12.0):
    """Generates a chain's contact map and distance matrix"""
    Chain = model.child_dict[chain]
    dist_matrix = calc_dist_matrix(Chain)
    contact_map = dist_matrix < MapThreshold
    return dist_matrix



def MathcingMapAlign(mapind,Alignmentname,Data1):
    '''this fucntion take the Alignment name and the maps map index (in the MSATestProteins.pkl dict),
     returns the sub-map and mathcing part in the alignment '''
    Align = np.load('/cs/cbio/orzuk/projects/ContactMaps/data/Pfam/alignments/%s'%Alignmentname)
    Data1list = Data1[Alignmentname][mapind]
    rowname = Data1list[1]
    mapindices = Data1list[3][:-1].split('-')
    SubMap = Data1list[0][int(mapindices[0]):int(mapindices[1]),int(mapindices[0]):int(mapindices[1])]
    ind = [a ==rowname for a in Align.index]
    RowStringlist = list(list(Align[ind].iloc[0])[0])
    MapAlignIndices = [x for x in range(len(RowStringlist)) if RowStringlist[x] != '-']
    d = [Align.Sequence.str[k] for k in MapAlignIndices]
    AlignMap = pd.concat(d,axis=1)
    AlignMap.columns = MapAlignIndices
    return [AlignMap,SubMap]


'''The following part generates dictionairy to match distance matrices to Pfam alignments'''


PfamFamTest = pickle.load(open('/cs/cbio/orzuk/projects/ContactMaps/data/Pfam/PfamFamsTest.p', 'rb'))
PfamFamTest = PfamFamTest[1:]
for i in range(len(PfamFamTest)):
    PfamFamTest[i] = str(PfamFamTest[i])

file = open('/cs/cbio/orzuk/projects/ContactMaps/data/Pfam/PfamPDB')
MSATestProteins = pickle.load(open('/cs/cbio/orzuk/projects/ContactMaps/data/MSATestProteins.pkl','rb'),encoding = 'latin1')
Ks = MSATestProteins.keys()
lines = file.readlines()

count = 0
for line in lines:
    if line.startswith('#=GF AC'):
        key = line.split(' ')
        key = list(filter(None, key))[2]
        key = key[:-1]
        if key not in Ks:
            MSATestProteins[key] = []
    if key not in Ks:

        if key.split('.')[0] in PfamFamTest:

            if line.startswith('#=GS'):
                a = line.replace(';', '').split(' ')
                a = list(filter(None, a))
                if len(a)==7:
                    Uniprot = a[1]
                    protein = a[4].lower()
                    protfile = protein[1:3]
                    chain = a[5]
                    indices = a[6]
                    try:
                        DistMat = contactmap(model(protein+'.pdb'), chain)
                        MSATestProteins[key].append([DistMat,Uniprot,protein+chain,indices])
                    except IOError:
                        print('pick')
                    except KeyError:
                        print('pick2')
    print(count)
    if count==10000:

        pickle.dump(MSATestProteins,open('/cs/cbio/orzuk/projects/ContactMaps/data/MSATestProteins.pkl','wb'))
        count=0
    count +=1



def PrtoteinFeatures(ProteinDataDict):
    '''this function generates the features per protein our of jinbo's data'''
    ACC1 = ['ACC1' + '.%s' % i for i in range(3)]
    ACC2 = ['ACC2' + '.%s' % i for i in range(3)]
    PSFM1 = ['PSFM1' + '.%s' % i for i in range(20)]
    PSFM2 = ['PSFM2' + '.%s' % i for i in range(20)]
    DISO1 = ['DISO1' + '.%s' % i for i in range(1)]
    DISO2 = ['DISO2' + '.%s' % i for i in range(1)]
    PSSM1 = ['PSSM1' + '.%s' % i for i in range(20)]
    PSSM2 = ['PSSM2' + '.%s' % i for i in range(20)]
    SS31 = ['SS31' + '.%s' % i for i in range(3)]
    SS32 = ['SS32' + '.%s' % i for i in range(3)]
    SS81 = ['SS81' + '.%s' % i for i in range(8)]
    SS82 = ['SS82' + '.%s' % i for i in range(8)]

    amino = 'ACDEFGHIKLMNPQRSTVWY'
    Amino = list(amino)
    Amino2 = [Amino[i]+str(2) for i in range(len(Amino))]
    sequence = ProteinDataDict['sequence']
    sequence = list(sequence)
    L = len(sequence)
    Seqindices1D = np.arange(0, L)
    '''Pairs of indices for the 1D features'''
    Seqindices1D = list(itertools.combinations(Seqindices1D, 2))
    '''All pairs of aminos (on the upper triangular)'''
    aminos = list(map(list, zip(*list(itertools.combinations(sequence, 2)))))
    aminos = pd.concat([pd.get_dummies(pd.Categorical(aminos[0], categories=Amino), drop_first = True)
                           ,pd.get_dummies(pd.Categorical(aminos[1], categories=Amino), drop_first = True)],axis=1)
    '''upper triangular indices for pairwise features'''
    index = np.triu_indices(L, k=1)
    n=len(index[1])
    SequneceDistance = (index[1]-index[0]).reshape((n,1))
    Contact = ProteinDataDict['contactMatrix'][index].reshape((n,1))
    ProteinName = np.array([ProteinDataDict['name']]*n).reshape((n,1))
    Length = np.array([L]*n).reshape((n,1))
    ccmpred = ProteinDataDict['ccmpredZ'][index].reshape((n,1))

    psi = ProteinDataDict['psicovZ'][index].reshape((n,1))

    OtherPairs = ProteinDataDict['OtherPairs']
    information = OtherPairs[:, :, 1][index].reshape((n,1))

    potential = OtherPairs[:, :, 2][index].reshape((n,1))

    ACC=ProteinDataDict['ACC'][Seqindices1D,]
    TwoDshape = [ACC.shape[0],ACC.shape[2]*2]
    ACC = ACC.reshape(TwoDshape)
    PSFM = ProteinDataDict['PSFM'][Seqindices1D,]
    TwoDshape = [PSFM.shape[0],PSFM.shape[2]*2]
    PSFM = PSFM.reshape(TwoDshape)



    DISO = ProteinDataDict['DISO'][Seqindices1D,]
    TwoDshape = [DISO.shape[0],DISO.shape[2]*2]
    DISO = DISO.reshape(TwoDshape)

    PSSM = ProteinDataDict['PSSM'][Seqindices1D,]
    TwoDshape = [PSSM.shape[0],PSSM.shape[2]*2]
    PSSM = PSSM.reshape(TwoDshape)

    SS3 = ProteinDataDict['SS3'][Seqindices1D,]
    TwoDshape = [SS3.shape[0],SS3.shape[2]*2]
    SS3 = SS3.reshape(TwoDshape)

    SS8 = ProteinDataDict['SS8'][Seqindices1D,]
    TwoDshape = [SS8.shape[0],SS8.shape[2]*2]
    SS8 = SS8.reshape(TwoDshape)

    Features = [aminos,ccmpred,psi,information,potential,Length,SequneceDistance,Contact,ACC,PSFM,DISO,PSSM,SS3,SS8]
    FeatureNames = Amino[1:]+Amino2[1:]+ ['ccmpredZ', 'psicovZ', 'mutual information', 'pairwise contact potential',
             'Length', 'SequneceDistance',
            'Contact'] + ACC1 + PSFM1 + DISO1 + PSSM1 + SS31 + SS81 + ACC2 + PSFM2 + DISO2 + PSSM2 + SS32 + SS82

    Data = pd.DataFrame(np.column_stack(Features),columns = FeatureNames).apply(pd.to_numeric)
    ind = (Data['Contact'] != -1)&(Data['SequneceDistance']>5)
    Data  = Data[ind]
    return (Data)



def Sequence(pdb,chain,path):
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.Polypeptide import three_to_one
    from Bio.PDB.Polypeptide import is_aa
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    #path = '/cs/cbio/orzuk/projects/ContactMaps/data/pdb/'

    ## First, open and parse the protein file
    p = PDBParser(PERMISSIVE=1)
    structure = p.get_structure(path+pdb+'.pdb',path+ pdb+'.pdb')
    model = structure[0]
    chain=model[chain]

    ## Now go through the hierarchy of the PDB file
    ##
    ## 1- Structure
    ##      2- Model
    ##          3- Chains
    ##              4- Residues
    ##



    seq = list()
    chainID = chain.get_id()

    for residue in chain:
                ## The test below checks if the amino acid
                ## is one of the 20 standard amino acids
                ## Some proteins have "UNK" or "XXX", or other symbols
                ## for missing or unknown residues
        if is_aa(residue.get_resname(), standard=True):
            seq.append(three_to_one(residue.get_resname()))
        else:
            seq.append("X")


    myProt = Seq(str(''.join(seq)), IUPAC.protein)
    #seqObj = SeqRecord(myProt, id=chainID, name="", description="")
    return str(myProt)


from os import listdir

files = listdir('/cs/cbio/orzuk/projects/ContactMaps/data/pdb')
path ='/cs/cbio/orzuk/projects/ContactMaps/data/pdb/'
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
p = PDBParser(PERMISSIVE=1)
SequnecesDict={}
count = 1.0
k=1
for pdb in files:
    print(count)
    try:
        structure = p.get_structure(path+pdb,path+ pdb)
        model = structure[0]
        chains =[]

        for i in model.child_list:
            if str(i)[-2]!=' ':
                chains.append(str(i)[-2])
        for chain in chains:
            SequnecesDict[pdb[:4]+chain] = Sequence(pdb[:4],chain,path).strip('X')
        count+=1
    except UnicodeDecodeError:
        'pcik'


    if count==k*100:
        print(count/len(files))
        k+=1

##Generate phylogenetic tree
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from Bio.Phylo.TreeConstruction import DistanceCalculator

import numpy as np
import pickle

Data1 =  pickle.load(open('/cs/cbio/orzuk/projects/ContactMaps/data/MSATestProteins1.pkl','rb'),encoding = 'latin1')
Y='PF03061.21'
Align = np.load('/cs/cbio/orzuk/projects/ContactMaps/data/Pfam/alignments/%s'%Y)
ind = [a =='ENTH_ECOLI/50-127' for a in Align.index]
Align[ind]
'''Take the alignment and tuen it to MultipleSeqAlignment object'''
def MSAOBJ(Align):
    calculator = DistanceCalculator('identity')

    MSAlst = []
    for indx  in Align.index:
        ind = [a == indx for a in Align.index]
        seq = Seq(list(Align[ind].iloc[0])[0])
        MSAlst.append(SeqRecord(seq,id=indx))
        align = MultipleSeqAlignment(MSAlst)
        dm = calculator.get_distance(align)
        return(align,dm)


'''now for the tree'''
from Bio.Phylo.TreeConstruction import DistanceCalculator

path ='/cs/cbio/orzuk/projects/ContactMaps/data/pdb/'
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

from Bio.SeqRecord import SeqRecord


import numpy as np
import pickle


Data1 =  pickle.load(open('/cs/cbio/orzuk/projects/ContactMaps/data/MSATestProteins.pkl','rb'),encoding = 'latin1')
Y='PF00036.31'
Align = np.load('/cs/cbio/orzuk/projects/ContactMaps/data/Pfam/alignments/%s'%Y)
ind = [a =='ENTH_ECOLI/50-127' for a in Align.index]
Align[ind]
Yprots =[]
Uniprots = []
for i in range(len(Data1[Y])):
    Yprots.append([Data1[Y][i][2][:4]+'.pdb',Data1[Y][i][2][-1]])
    Uniprots.append(Data1[Y][i][1])




'''Starting with an alignment'''
JAlign = open('/cs/cbio/orzuk/projects/ContactMaps/data/272871.a2m')
#JAlign = open('272871.all_in_one/272871.a2m')

MSA0 = JAlign.readlines()
'''Next turn the line in to a sequence object'''
MSA = []
MSAList =[]
count=1
for seq in MSA0:
    Seqlist = list(seq[:-1])
    Seq1 = Seq(''.join(Seqlist))
    Seqlist = [Seqlist[x] for x in range(len(Seqlist)) if (Seqlist[x] != '-' and Seqlist[x].upper()==Seqlist[x])]
    SEQ = Seq(''.join(Seqlist))
    MSA.append(SeqRecord(Seq1,id=str(count)))
    MSAList.append(SEQ)
    count+=1

MSA = MultipleSeqAlignment(MSA)


'''indices for the maps'''
def AlignIndices(MSAline):
    MSAline = list(MSAline)
    Indices = [x for x in range(len(MSAline)) if (MSAline[x] != '-' and MSAline[x].upper()==MSAline[x])]
    return(Indices)


def TreeDistMat(AlignObject):
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(AlignObject)

    return(dm)

def FindMaps(AlignObject,protein,pdblist):
    dm =TreeDistMat(AlignObject)
    mat = dm.matrix[protein]
    MapsDict = {}
    Distances =[]
    for Protein in pdblist:
        MapsDict[Protein] = #The protein's map
        Distances.append(mat[protein,Protein])
    PDBsDistance = pd.DataFrame({'Distance':Distances,"PDB'S":pdblist})
    PDBsDistance = PDBsDistance.sort_values('Distance',1)
    '''the '''
    Outcontact.iloc[AlignIndices,AlignIndices]





SequnecesDict=pickle.load(open('/cs/cbio/orzuk/projects/ContactMaps/data/SequnecesDict','rb'))

for i in range(len(MSAList)):
    if MSAList[i] in SequnecesDict.values():
        print('pick')







'''We need the sequences dict to match sequnece to proteins'''


from Bio.Phylo.TreeConstruction import DistanceCalculator

path ='/cs/cbio/orzuk/projects/ContactMaps/data/pdb/'
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq



import numpy as np
import pickle

from Bio import SeqIO

from os import listdir

files = listdir('/cs/cbio/orzuk/projects/ContactMaps/data/pdb')
files = [file for file in files if file.endswith('.pdb') ]
p = PDBParser(PERMISSIVE=1)
SequnecesDict={}
count = 1.0
k=1

'''Sequences for pdbs'''
for pdb in files:
    try:

        handle = open(path+pdb, "r")

        for record in SeqIO.parse(handle, "pdb-seqres"):
            SequnecesDict[record.id]=record.seq

        count+=1
    except UnicodeDecodeError:
        'pick'


    if count==k*round(len(files)/100.0):
        print('%i precent completed'%k)
        k+=1

pickle.dump(SequnecesDict,open('/cs/cbio/orzuk/projects/ContactMaps/data/SequnecesDict','wb'))













for indx  in Uniprots:
    ind = [a == indx for a in Align.index]
    lst = list(list(Align[ind].iloc[0])[0])
    c = [lst[x] for x in range(len(lst)) if lst[x] != '-' and lst[x].upper()==lst[x]]
    seq=''
    for j in c:
        seq+=j
    seq = Seq(seq)


'''parsing the pdb-uniprot dict'''
import pandas as pd
from Bio import SeqIO

Lines=open('pdbtosp.txt').readlines()
file =''
def Prse(line):
    Linelst = line.split(' ')
    Linelst =list(filter(None, Linelst))
    try:
        Linelst.remove(',')
    except ValueError:
        'pick'
    try:
        Linelst[2] = Linelst[2] + Linelst[3]
    except IndexError:
        'pick'
    String = ','.join(Linelst)
    return (String)
for line in Lines:
    file=file+'\n'+Prse(line)
PDBUnip =pd.read_csv("Output.txt",error_bad_lines=False)
def matchpdbs(alignmentfile):
    align=SeqIO.parse(alignmentfile,'stockholm')
    Align = list(align)
    names = []
    for i in range(len(Align)):
        names.append(Align[i].id.split('|'))
    dat = pd.DataFrame(names)
