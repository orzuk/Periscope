import pandas as pd
import numpy as np
import pickle
import itertools
from os import listdir
import tensorflow as tf


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
ColNames = ['Amino1','Amino2','ccmpredZ', 'psicovZ', 'mutual_information', 'pairwise_contact_potential',
 'SequneceDistance'] + ACC1 + PSFM1 + DISO1 + PSSM1 + SS31 + SS81 + ACC2 + PSFM2 + DISO2 + PSSM2 + SS32 + SS82

FeatCols = []
for i in ColNames:
    if i.startswith('Amino'):
        FeatCols.append(tf.feature_column.categorical_column_with_vocabulary_list(i,Amino))
    else:
        FeatCols.append(tf.feature_column.numeric_column(i))

def input_fn_fen(ProteinDataDict):
    Feat = PrtoteinFeatures(ProteinDataDict)


    y=pd.DataFrame(Feat['Contact'])
    temp = Feat.pop('Contact')
    x = pd.DataFrame(Feat)

    return(tf.estimator.inputs.pandas_input_fn(x,y,shuffle=False))

path = input('Please enter path to the train data (pdb25-6767-train.release.contactFeatures.pkl) :\nIf you dont have it just use:\n'
             'wget "http://raptorx.uchicago.edu/download/b3IuenVrQG1haWwuaHVqaS5hYy5pbA==/26/" -O pdb25-6767-train.release.contactFeatures.pkl\npath(including the file name): ')
data = pickle.load(open(path,'rb'),encoding='latin1')

def PrtoteinFeatures(ProteinDataDict):
    print('Generating Feature for protein %s chain %s\n'%(protein['name'][:-1],protein['name'][-1]))

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

    aminos = pd.DataFrame([pd.Categorical(aminos[0], categories=Amino), pd.Categorical(aminos[1], categories=Amino)])
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

    FeatureNames = ['Amino1','Amino2','ccmpredZ', 'psicovZ', 'mutual_information', 'pairwise_contact_potential',
             'Length', 'SequneceDistance',
            'Contact'] + ACC1 + PSFM1 + DISO1 + PSSM1 + SS31 + SS81 + ACC2 + PSFM2 + DISO2 + PSSM2 + SS32 + SS82
    Data = pd.DataFrame(np.column_stack(Features),columns = FeatureNames)
    ind = (Data['Contact'] != -1)&(Data['SequneceDistance']>5)
    Data  = Data[ind]
    #return (Data)
    return (Data.to_dict(orient='list'))


files = listdir('/cs/cbio/orzuk/projects/ContactMaps/data/ProteinTrainData')



'''TF-Train'''
run_config = tf.estimator.RunConfig().replace(session_config=tf.ConfigProto(device_count={'GPU': 0}))
path2 = input('Please enter path for saving the model')
model = tf.estimator.LinearClassifier(feature_columns=FeatCols,model_dir=path2, config=run_config)
for protein in data:
    print('Training protein %s chain %s\n'%(protein['name'][:-1],protein['name'][-1]))

    model.train(input_fn=input_fn_fen(protein))
