from sklearn import metrics
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef, confusion_matrix


class Evalutation_Metrics:
    def __init__(self):
        print("eval metrics")

    def fundamentalEvalMetrics(self, y_pred, y_test, task='multi-class'):
        if task=="multi-class":
            acc=accuracy_score(y_pred, y_test)
            pre=precision_score(y_test, y_pred, average="weighted")
            rec=recall_score(y_test, y_pred, average="weighted")
            f1=f1_score(y_test, y_pred, average="weighted")
            mcc= round(matthews_corrcoef(y_test, y_pred), 4)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
            spec=round(tn / (tn + fp), 4)
            sensitivity = round( tp / (tp + fn),4)
            return acc, pre, rec, f1, mcc, spec, sensitivity
