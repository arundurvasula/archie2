"""
|-----------------------------------+---------------------------------|
| Name                              | Definition                      |
|-----------------------------------+---------------------------------|
| Accuracy                          | (TP + TN) / (FP + TP + FN + TN) |
| True positive rate (TPR)          | TP / (TP + FN)                  |
| Sensitivity                       | see TPR                         |
| Recall                            | see TPR                         |
| Hit rate                          | see TPR                         |
| Probability of detection          | see TPR                         |
| True negative rate (TNR)          | TN / (TN + FP)                  |
| Specificity                       | see TNR                         |
| Selectivity                       | see TNR                         |
| Positive predictive value (PPV)   | TP / (TP + FP)                  |
| Precision                         | see PPV                         |
| Negative predictive value (NPV)   | TN / (TN + FN)                  |
| False negative rate (FNR)         | FN / (FN + TP)                  |
| Miss rate                         | see FNR                         |
| False positive rate (FPR)         | FP / (FP + TN)                  |
| Fall out                          | see FPR                         |
| False discovery rate (FDR)        | FP / (FP + TP)                  |
| False omission rate (FOR)         | FN / (FN + TN)                  |
| F1 score                          | 2TP / (2TP + FP + FN)           |
| F-score                           | see F1 score                    |
| F-measure                         | see F1 score                    |
| Sorensen-Dice coefficient         | see F1 score                    |
| Dice similarity coefficient (DSC) | see F1 score                    |
|-----------------------------------+---------------------------------|

Table from https://dzone.com/articles/accuracy-precision-and-recall
"""

import pandas
import numpy
import sys

mus = ["0.000000000125",
       "0.000000000625",
       "0.00000000125",
       "0.00000000625",
       "0.0000000125",
       "0.0000000625",
       "0.000000125"]

rs = ["0.0000000001",
      "0.0000000005",
      "0.000000001",
      "0.000000005",
      "0.00000001",
      "0.00000005",
      "0.0000001"]

thresholds = {}
thresholds["T_0_10"] = 0.10
thresholds["T_0_15"] = 0.15
thresholds["T_0_20"] = 0.20
thresholds["T_0_25"] = 0.25
thresholds["T_0_30"] = 0.30
thresholds["T_0_35"] = 0.35
thresholds["T_0_40"] = 0.40
thresholds["T_0_45"] = 0.45
thresholds["T_0_50"] = 0.50
thresholds["T_0_55"] = 0.55
thresholds["T_0_60"] = 0.60
thresholds["T_0_65"] = 0.65
thresholds["T_0_70"] = 0.70
thresholds["T_0_75"] = 0.75
thresholds["T_0_80"] = 0.80
thresholds["T_0_85"] = 0.85
thresholds["T_0_90"] = 0.90
thresholds["T_0_95"] = 0.95


df = pandas.read_csv("archie_statistics_per_snp.csv")

cols = ["SET", "MU", "MU_TEST2TRAIN", "R", "R_TEST2TRAIN",
        "THRESHOLD", "LOO_REPL",
        "TP", "FN", "FP", "TN",
        "TRUE_POS", "TRUE_NEG", "LABL_POS", "LABL_NEG", "ALL",
        "ACCURACY", "TNR", "FPR", "FNR",
        "RECALL", "PRECISION", "FDR", "FOR"]
rows = []

for i in range(1, 50):
    print("Processing %d" % i)
    for t in thresholds.keys():
        for repl in range(1, 101):
            # print("LOO replicate %d" % repl)

            xdf = df[(df["SET"] == i) & 
                     (df["THRESHOLD"] == t) &
                     (df["REPL"] != repl)]

            if xdf.shape[0] != 99:
                print(xdf)
                sys.exit("Wrong number of replicates!")

            tp = numpy.sum(xdf.TP.values)
            fn = numpy.sum(xdf.FN.values)
            fp = numpy.sum(xdf.FP.values)
            tn = numpy.sum(xdf.TN.values)

            true_pos = tp + fn
            true_neg = fp + tn
            labl_pos = tp + fp
            labl_neg = tn + fn
            tot = true_pos + true_neg

            row = {}
            row["SET"] = i
            row["MU"] = xdf.MU.values[0]
            row["R"] = xdf.R.values[0]
            row["MU_TEST2TRAIN"] = row["MU"] / 0.0000000125
            row["R_TEST2TRAIN"] = row["R"] / 0.00000001
            row["THRESHOLD"] = thresholds[t]
            row["LOO_REPL"] = repl

            row["TP"] = tp
            row["FN"] = fn
            row["FP"] = fp
            row["TN"] = tn

            row["TRUE_POS"] = true_pos
            row["TRUE_NEG"] = true_neg
            row["LABL_POS"] = labl_pos  # sometimes 0
            row["LABL_NEG"] = labl_neg
            row["ALL"] = tot

            if true_pos == 0:
              print("  %d / %d" % (true_pos, tot))

            row["ACCURACY"] = (tp + tn) / tot

            if true_neg > 0:
                row["TNR"] = tn / true_neg
                row["FPR"] = fp / true_neg
            if true_pos > 0:
                row["FNR"] = fn / true_pos
                row["RECALL"] = tp / true_pos  # TPR
            if labl_pos > 0:
                row["PRECISION"] = tp / labl_pos  # sometimes does not exist
                row["FDR"] = fp / labl_pos        # sometimes does not exist
            if labl_neg > 0:
                row["FOR"] = fn / labl_neg
            rows.append(row)

df = pandas.DataFrame(rows, columns = cols)
df.to_csv("archie_statistics_per_snp_loo.csv", 
          sep="," , index=None, na_rep="NA")
