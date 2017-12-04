from sklearn.linear_model import LogisticRegression
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


data = pickle.load(open('/cs/cbio/orzuk/projects/ContactMaps/data/Raptor/pdb25-6767-train.release.contactFeatures.pkl','rb'),encoding='latin1')

def PrtoteinFeatures(ProteinDataDict):

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
    #aminos = pd.concat([pd.get_dummies(pd.Categorical(aminos[0], categories=Amino), drop_first = True)
    #                       ,pd.get_dummies(pd.Categorical(aminos[1], categories=Amino), drop_first = True)],axis=1)
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
    #FeatureNames = Amino[1:]+Amino2[1:]+ ['ccmpredZ', 'psicovZ', 'mutual information', 'pairwise contact potential',
    #         'Length', 'SequneceDistance',
    #        'Contact'] + ACC1 + PSFM1 + DISO1 + PSSM1 + SS31 + SS81 + ACC2 + PSFM2 + DISO2 + PSSM2 + SS32 + SS82


    FeatureNames = ['Amino1','Amino2','ccmpredZ', 'psicovZ', 'mutual_information', 'pairwise_contact_potential',
             'Length', 'SequneceDistance',
            'Contact'] + ACC1 + PSFM1 + DISO1 + PSSM1 + SS31 + SS81 + ACC2 + PSFM2 + DISO2 + PSSM2 + SS32 + SS82
    Data = pd.DataFrame(np.column_stack(Features),columns = FeatureNames)
    ind = (Data['Contact'] != -1)&(Data['SequneceDistance']>5)
    Data  = Data[ind]
    return (Data.to_dict(orient='list'))
count = 1


files = listdir('/cs/cbio/orzuk/projects/ContactMaps/data/ProteinTrainData')



'''TF-Train'''
run_config = tf.estimator.RunConfig().replace(session_config=tf.ConfigProto(device_count={'GPU': 0}))

model = tf.estimator.LinearClassifier(feature_columns=FeatCols,model_dir='/cs/cbio/orzuk/projects/ContactMaps/data/Jinbo', config=run_config)
count =1
for protein in data:
    print(count)

    model.train(input_fn=input_fn_fen(protein))
    count+=1

'''Testing'''

with open('/cs/cbio/orzuk/projects/ContactMaps/data/PfamTestData.p','rb') as f:
    PfamTestData = pickle._Unpickler(f)
    PfamTestData.encoding = 'latin1'
    PfamTestData = PfamTestData.load()
with open('/cs/cbio/orzuk/projects/ContactMaps/data/Jinbo/LogRegModelNew2','rb') as f:
    LogRegModel50 = pickle._Unpickler(f)
    LogRegModel50.encoding = 'latin1'
    LogRegModel50 = LogRegModel50.load()

with open('LogRegModelNew2','rb') as f:
    LogRegModel50 = pickle._Unpickler(f)
    LogRegModel50.encoding = 'latin1'
    LogRegModel50 = LogRegModel50.load()

with open('PfamTestData.p','rb') as f:
    PfamTestData = pickle._Unpickler(f)
    PfamTestData.encoding = 'latin1'
    PfamTestData = PfamTestData.load()



PfamTestData = pickle.load(open('PfamTestData.p','rb'),encoding='latin1')

PfamTestData = pickle.load(open('/cs/cbio/orzuk/projects/ContactMaps/data/PfamTestData.p','rb'))
LogRegModel50 = pickle.load(open('/cs/cbio/orzuk/projects/ContactMaps/data/Jinbo/LogRegModelNew2','rb'))
PPVS = pd.DataFrame(columns = ['L','L/2','L/5','L/10'])
PPVM = pd.DataFrame(columns = ['L','L/2','L/5','L/10'])
PPVL = pd.DataFrame(columns = ['L','L/2','L/5','L/10'])
FinalPPV  = [PPVS,PPVM,PPVL]
i=0
Short = []
Medium = []
Long = []
ROCData = pd.DataFrame(columns = ['Y','Yhat'])

for testprot in PfamTestData:
    sequence = testprot['sequence']
    sequence = list(sequence)
    L0 = len(sequence)
    Data0 = PrtoteinFeatures(testprot)
    Y = np.array(Data0['Contact'])
    Yhat = LogRegModel50.predict_proba(Data0.drop(['Contact'], axis=1))[:,1]

    shortind = Data0['SequneceDistance']<=11
    mediumind = np.logical_and(Data0['SequneceDistance']>11,Data0['SequneceDistance']<=23)
    longind = Data0['SequneceDistance']>=24
    Distind = [shortind,mediumind,longind]
    Ys = pd.DataFrame([Y,Yhat,shortind,mediumind,longind],
                      index = ['Y','Yhat','Short','Medium','Long']).transpose()
    YsSorted= Ys.sort_values('Yhat',ascending=False)



    ks = [1.0, 2.0, 5.0, 10.0]

    lengths = [L0,round(L0/2.0),
               round(L0/5.0),round(L0/10.0)]
    PPV0 = list()
    PPV1 = list()
    PPV2 = list()
    for length in lengths:
        S=YsSorted[:int(length)][YsSorted['Short']]['Y'].mean()
        if S==0:
            PPV0.append(1/length)
        else:
            PPV0.append(S)

        M=YsSorted[:int(length)][YsSorted['Medium']]['Y'].mean()
        if  M== 0:
            PPV1.append(1/length)
        else:
            PPV1.append(M)


        L= YsSorted[:int(length)][YsSorted['Long']]['Y'].mean()
        if L== 0:
            PPV2.append(1/length)
        else:
            PPV2.append(L)

    FinalPPV[0] = pd.concat((FinalPPV[0],pd.DataFrame(PPV0,index = ['L','L/2','L/5','L/10']).transpose()),axis=0)
    FinalPPV[1] = pd.concat((FinalPPV[1],pd.DataFrame(PPV1,index = ['L','L/2','L/5','L/10']).transpose()),axis=0)
    FinalPPV[2] = pd.concat((FinalPPV[2],pd.DataFrame(PPV2,index = ['L','L/2','L/5','L/10']).transpose()),axis=0)
    ROCData = pd.concat((ROCData,YsSorted))

Lengths = ['Short','Medium','Long']
i=0
for results in FinalPPV:
    print(Lengths[i])
    print(np.mean(results))
    i+=1



'''Roc curve'''
from sklearn import metrics


fprS, tprS, thresholdsS = metrics.roc_curve(list(ROCData['Y'][ROCData['Short']]), list(ROCData['Yhat'][ROCData['Short']]), pos_label=1)
fprM, tprM, thresholdsM = metrics.roc_curve(list(ROCData['Y'][ROCData['Medium']]), list(ROCData['Yhat'][ROCData['Medium']]), pos_label=1)
fprL, tprL, thresholdsL = metrics.roc_curve(list(ROCData['Y'][ROCData['Long']]), list(ROCData['Yhat'][ROCData['Long']]), pos_label=1)



FinalRocData = [[fprS, tprS, thresholdsS],[fprM, tprM, thresholdsM],[fprL, tprL, thresholdsL]]


import matplotlib.pyplot as plt

plt.figure()
lw = 2
plt.plot(fprS, tprS, color='darkorange',
         lw=lw,label = 'Short' )
plt.plot(fprM, tprM, color='green',
         lw=lw,label = 'Medium' )
plt.plot(fprL, tprL, color='red',
         lw=lw,label = 'Long' )
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Roc for the logisitic regression model')
plt.legend(loc="lower right")
plt.show()

'''Old Training'''
for protein in data:
    print(count)
    if protein['name'] not in files:
        Feat = PrtoteinFeatures(protein)
        pickle.dump(Feat, open('/cs/cbio/orzuk/projects/ContactMaps/data/ProteinTrainData/%s'%protein['name'], 'wb'))
    count+=1

    if len(list(set(Feat['Contact'])))<=1:
        print('No contact')
    else:
        Model.fit(Feat.drop(['Contact'], axis=1), Feat['Contact'])
    pickle.dump(Model, open('Jinbo/LogRegModelNew2', 'wb'))

    count+=1

mshelffile='/cs/cbio/orzuk/projects/ContactMaps/data/Raptor/pdb25-6767-train.release.contactFeatures.pkl'
with open(mshelffile, 'rb') as f:
    d = pickle.load(f, encoding='latin1')
Model = LogisticRegression(solver = 'sag', warm_start = True)
