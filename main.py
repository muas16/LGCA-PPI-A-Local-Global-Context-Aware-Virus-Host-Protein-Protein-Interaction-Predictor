import os
import pickle
from pathlib import Path
import numpy
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import *
from  Dataset_Utils.utils import *
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
# set the folder directory...
dir = "./LGCA-PPI-A-Local-Global-Context-Aware-Virus-Host-Protein-Protein-Interaction-Predictor"
# set the datasets path....
train_data = "PVI_Datasets/DeNovo_Dataset/DeNovo_full_train.csv"
test_data = "PVI_Datasets/DeNovo_Dataset/DeNovo_full_test.csv"

save_path = os.path.join(dir, "encodings")

dataset="denovo"

if not os.path.exists(save_path):

    os.mkdir(save_path)
# define encoding type
encoding_type=["Local"]
cmb = [100, 0, 100, 0]
columns = ["seq_x", "seq_y", "class"]

chunk=4
for enc_type in encoding_type:

    print("processing train data")
    data_train, df_train=process_Data(train_data,cmb,columns,"tr", chunk=chunk, enc_type=enc_type)
    print("processing test data")
    data_test, df_test=process_Data(test_data,cmb,columns,"tr", chunk=chunk, enc_type=enc_type)
    df=pd.concat([df_train,df_test],ignore_index=False)
    df.to_csv(save_path+"/"+"_"+str("-".join([str(i) for i in cmb]))+dataset+"_"+enc_type+"-features.csv", index=False, header=False)
    train_labels = df_train.iloc[:, -1:]
    df_train.drop(df_train.columns[-1], axis=1, inplace=True)
    test_labels = df_test.iloc[:, -1:].iloc[:, 0].tolist()
    df_test.drop(df_test.columns[-1], axis=1, inplace=True)
    train_features = np.array(df_train)
    test_features = np.array(df_test)
    clf = RandomForestClassifier(random_state=42)
    clf.fit(train_features, train_labels.iloc[:, 0].tolist())
    y_pred = clf.predict(test_features)
    y_pred_proba=clf.predict_proba(test_features)
    print("Accuracy", str(round(accuracy_score(y_pred, test_labels),4)))
    print("Precision", str(round(precision_score(test_labels, y_pred, average="weighted"),4)))
    print("Recall", str(round(recall_score(test_labels, y_pred, average="weighted"),4)))
    print("F1", str(round(f1_score(test_labels, y_pred, average="weighted"),4)))
    print("MCC", str(round(matthews_corrcoef(test_labels, y_pred),4)))
    tn, fp, fn, tp = confusion_matrix(test_labels, y_pred).ravel()
    print("Specificity", str(round(tn / (tn + fp),4)))
    numpy.savez(save_path+"/"+enc_type+"-Encod.npz", test_labels, y_pred_proba)

    plot_roc(testY, y_pred_proba,"Denovo/"+str("-".join([str(i) for i in cmb]))+"roc.png" )
    plot_prc(testY, y_pred_proba,"Denovo/"+str("-".join([str(i) for i in cmb]))+"prc.png" )

    df = pd.DataFrame()
    df["Acc"] = [str(round(accuracy_score(y_pred, test_labels), 4))]
    df["pre"] = [str(round(precision_score(test_labels, y_pred, average="weighted"), 4))]
    df["f1"] = [str(round(f1_score(test_labels, y_pred, average="weighted"), 4))]
    df["rec"] = [str(round(recall_score(test_labels, y_pred, average="weighted"), 4))]
    df["mcc"] = [ str(round(matthews_corrcoef(test_labels, y_pred), 4))]
    df["spec"] = [str(round(tn / (tn + fp),4))]
    df["Encoding_type"] = [enc_type]
    df["cmb"] = [str("-".join([str(i) for i in cmb]))]
    df.to_csv(os.path.join(dir, "Results.csv"), index=False)



