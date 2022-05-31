import random
import numpy
import pandas as pd
import sklearn
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import *
from sklearn.model_selection import StratifiedKFold
from Dataset_Utils.utils import *
from Evaluation_Measures.Evaluation_Metrics import Evalutation_Metrics
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import os
def random_seed(seed_value, use_cuda):
    numpy.random.seed(seed_value) # cpu vars
    random.seed(seed_value) # Python
random_seed(42, False)
ev=Evalutation_Metrics()
# set the folder directory...
dir = "./LGCA-PPI-A-Local-Global-Context-Aware-Virus-Host-Protein-Protein-Interaction-Predictor"
# set the datasets path....
data_path="/home/rehab/Videos/DeNovo_Dataset/DeNovo_full_train.csv"
save_path=os.path.join(dir,"encodings")
# set th eencoding type..
encoding_type=["Global"]


# define the number of folds...
kfold=5

chunk=4
cmb = [100,0,100,0]
columns = ["seq_x","seq_y","class"]
dataset="Berman"



for enc_type in  encoding_type:

    print("processing train data")
    data_train,df=process_Data(data_path,cmb,columns,"tr", chunk, enc_type)
    df.to_csv(save_path+"/"+str("-".join([str(i) for i in cmb]))+dataset+"_"+enc_type+"-features.csv", index=False, header=False)
    averageScores = {"acc": 0, "pre": 0, "rec": 0, "f1": 0, "mcc":0, "spec":0}

    kf = StratifiedKFold(n_splits=kfold, random_state=42, shuffle=True)
    fold=1
    actual=[]
    predict_proba=[]
    for train_index, test_index in kf.split(data_train[:, 1:], data_train[:, 0]):
        print("Fold No.",fold)

        X_train, X_test = data_train[:, 1:][train_index],data_train[:, 1:][test_index]
        y_train, y_test = data_train[:, 0][train_index],data_train[:, 0][test_index]
        clf = RandomForestClassifier(max_depth=702, n_estimators=291, random_state=42)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        y_pred_proba = clf.predict_proba(X_test)

        actual.extend(y_test)
        predict_proba.extend(y_pred_proba)

        acc, pre, rec, f1, mcc, spec, sensitivity = ev.fundamentalEvalMetrics(y_pred, y_test, task="multi-class")

        print("Accuracy", acc)
        print("Precision", pre)
        print("Recall", rec)
        print("F1", f1)
        print("MCC", mcc)
        print("SPEC", spec)

        averageScores["acc"] = averageScores["acc"] + acc
        averageScores["pre"] = averageScores["pre"] + pre
        averageScores["rec"] = averageScores["rec"] + rec
        averageScores["f1"] = averageScores['f1'] + f1
        averageScores["mcc"] = averageScores['mcc'] + mcc
        averageScores["spec"] = averageScores['spec'] + spec

        fold += 1
    averageScores["acc"] = averageScores["acc"] / kfold
    averageScores["pre"] = averageScores["pre"] / kfold
    averageScores["rec"] = averageScores["rec"] / kfold
    averageScores["f1"] = averageScores['f1'] / kfold
    averageScores["mcc"] = averageScores['mcc'] / kfold
    averageScores["spec"] = averageScores['spec'] / kfold
    numpy.savez(save_path+"/"+enc_type+"-Encod.npz", actual, numpy.array(predict_proba))
    plot_roc(actual, numpy.array(predict_proba),dir+"/"+str("-".join([str(i) for i in cmb]))+"-"+"roc.png" )
    plot_prc(actual, numpy.array(predict_proba),dir+"/"+str("-".join([str(i) for i in cmb]))+"-"+"prc.png")
    df=pd.DataFrame()
    df["Acc"]=[str(round(averageScores["acc"], 4))]
    df["pre"]=[str(round(averageScores["pre"], 4))]
    df["f1"]=[str(round(averageScores["f1"], 4))]
    df["Rec"]=[str(round(averageScores["rec"], 4))]
    df["mcc"]=[str(round(averageScores["mcc"], 4))]
    df["spec"]=[str(round(averageScores["spec"], 4))]
    df["Encoding_type"]=[enc_type]
    df["cmb"]=[str("-".join([str(i) for i in cmb]))]
    df.to_csv(os.path.join(dir,"Results_Berman.csv"),index=False)
