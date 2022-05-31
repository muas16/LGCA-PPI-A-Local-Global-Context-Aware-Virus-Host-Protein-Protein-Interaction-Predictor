from math import exp
from numpy.linalg import norm

import numpy
import pandas as pd
from numpy.linalg import norm
from sklearn.decomposition import PCA
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.preprocessing import StandardScaler

from src.Encoding import Encoding


def generatesequences(sepearteddataset, columns, protAstart, protAend, protBstart, protBend):

    final_df = pd.read_csv(sepearteddataset, names=columns, sep=',', skiprows=1)
    x=final_df["seq_x"]
    y=final_df["seq_y"]

    final_df["seq_x"]=x
    final_df["seq_y"]=y

    for iterator in range(len(final_df)):
        sequence_1 = ""
        sequence_2 = ""

        if protAstart==100 and protAend==0:
            sequence_1=final_df[columns[0]].iloc[iterator][:]

        if protAstart != 100:
            sequence_1 += final_df[columns[0]].iloc[iterator][:protAstart]
        if protAend != 0:
            sequence_1 += final_df[columns[0]].iloc[iterator][-protAend:]

        if protBstart==100 and protBend ==0:
            sequence_2=final_df[columns[1]].iloc[iterator][:]

        if protBstart != 100:
            sequence_2 += final_df[columns[1]].iloc[iterator][:protBstart]
        if protBend != 0:
            sequence_2 += final_df[columns[1]].iloc[iterator][-protBend:]


        final_df["seq_x"].iloc[iterator] = sequence_1
        final_df["seq_y"].iloc[iterator] = sequence_2
    return final_df

def extract_feature(P_protein_a, P_protein_b, N_protein_a, N_protein_b, chunk, enc_type):

    N_A_feature = []
    N_B_feature = []
    P_A_feature = []
    P_B_feature = []



    m = len(N_protein_a)
    n = len(P_protein_a)

    for i in range(m):


        SEQ = N_protein_a[i]

        enc=Encoding()

        if enc_type=="Global":
            C = enc.GlobalEncoding(seq=SEQ)
            C = list(C.transpose())
            N_A_feature.append(list(C))
        if enc_type == "Local":
            D = enc.Sline(seq=SEQ)
            N_A_feature.append(list(D))
        if enc_type=="Local_Global":
            C = enc.GlobalEncoding(seq=SEQ)
            C = list(C.transpose())
            D = enc.Sline(seq=SEQ)
            N_A_feature.append(list(list(C) + D))

        # D = enc.Sline(seq=SEQ)
        # N_A_feature.append(list(D))

    N_A_feature=numpy.array(N_A_feature)
    for i in range(m):
        SEQ = N_protein_b[i]

        enc = Encoding()

        if enc_type == "Global":
            C = enc.GlobalEncoding(seq=SEQ)
            C = list(C.transpose())
            N_B_feature.append(list(C))
        if enc_type == "Local":
            D = enc.Sline(seq=SEQ)
            N_B_feature.append(list(D))
        if enc_type == "Local_Global":
            C = enc.GlobalEncoding(seq=SEQ)
            C = list(C.transpose())
            D = enc.Sline(seq=SEQ)
            N_B_feature.append(list(list(C) + D))
    N_B_feature=numpy.array(N_B_feature)

    for i in range(n):


        SEQ = P_protein_a[i]

        enc = Encoding()

        if enc_type == "Global":
            C = enc.GlobalEncoding(seq=SEQ)
            C = list(C.transpose())
            P_A_feature.append(list(C))
        if enc_type == "Local":
            D = enc.Sline(seq=SEQ)
            P_A_feature.append(list(D))
        if enc_type == "Local_Global":
            C = enc.GlobalEncoding(seq=SEQ)
            C = list(C.transpose())
            D = enc.Sline(seq=SEQ)
            P_A_feature.append(list(list(C) + D))
    P_A_feature = numpy.array(P_A_feature)

    for i in range(n):


        SEQ = P_protein_b[i]

        enc = Encoding()

        if enc_type == "Global":
            C = enc.GlobalEncoding(seq=SEQ)
            C = list(C.transpose())
            P_B_feature.append(list(C))
        if enc_type == "Local":
            D = enc.Sline(seq=SEQ)
            P_B_feature.append(list(D))
        if enc_type == "Local_Global":
            C = enc.GlobalEncoding(seq=SEQ)
            C = list(C.transpose())
            D = enc.Sline(seq=SEQ)
            P_B_feature.append(list(list(C) + D))


    P_B_feature = numpy.array(P_B_feature)


    return N_A_feature,N_B_feature,P_A_feature,P_B_feature

def write_features(N_A_feature,N_B_feature,P_A_feature,P_B_feature,type,):

    pd.DataFrame(N_A_feature).to_csv(type+"_N_A_feature.csv",index=False, header=False)
    pd.DataFrame(N_B_feature).to_csv(type+"_N_B_feature.csv",index=False, header=False)
    pd.DataFrame(P_A_feature).to_csv(type+"_P_A_feature.csv",index=False, header=False)
    pd.DataFrame(P_B_feature).to_csv(type+"_P_B_feature.csv",index=False, header=False)

def process_Data(dataset_path, cmb, columns,dt,chunk,enc_type):

    df = generatesequences(dataset_path, columns, cmb[0], cmb[1], cmb[2], cmb[3])

    classes = set(df["class"].values.tolist())
    Encoding = []

    P_protein = df[df["class"] == 1]
    # N_protein=df[df["class"] == -1]
    N_protein = df[df["class"] == 0]

    P_protein_A = P_protein["seq_x"].values.tolist()
    P_protein_B = P_protein["seq_y"].values.tolist()

    N_protein_A = N_protein["seq_x"].values.tolist()
    N_protein_B = N_protein["seq_y"].values.tolist()

    N_A_feature, N_B_feature, P_A_feature, P_B_feature = extract_feature(P_protein_A, P_protein_B, N_protein_A,
                                                                         N_protein_B, chunk, enc_type)

    PP = numpy.array([m + n for m, n in zip(P_A_feature.tolist(), P_B_feature.tolist())])
    NN = numpy.array([m + n for m, n in zip(N_A_feature.tolist(), N_B_feature.tolist())])
    Elegan_data_GE=numpy.concatenate((PP, NN))
    m = len(N_A_feature)
    n = len(P_A_feature)
    p_label= numpy.array([1]*n)
    n_label =numpy.array([-1]*m)
    Elegan_label=numpy.concatenate((p_label, n_label))

    df1 = pd.DataFrame(Elegan_data_GE)
    df2 = pd.DataFrame({"class": Elegan_label})
    df = pd.concat([df1, df2], axis=1)
    data = numpy.array([[m] + n for m, n in zip(Elegan_label.tolist(), Elegan_data_GE.tolist())])

    return data, df

def plot_roc(testY, y_pred_proba, path):
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc

    plt.figure()

    fpr, tpr, thresholds = roc_curve(testY, y_pred_proba[:, 1])
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, lw=1,
             label='  (area = %0.2f)' % roc_auc)


    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right", prop={'size': 20})
    plt.savefig(path)

    plt.show()
    print()
    return roc_auc

def plot_prc(testY, y_pred_proba,path):
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc

    plt.figure()

    pre, rec, thresholds = precision_recall_curve(testY, y_pred_proba[:, 1])
    ap = average_precision_score(testY, y_pred_proba[:, 1])

    plt.plot(pre, rec, lw=1,
             label='  (area = %0.2f)' % ap)

    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    # plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right", prop={'size': 20})
    plt.savefig(path)

    plt.show()
    return ap

