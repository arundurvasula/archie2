import pandas
import numpy
import sys

mus = [0.000000000125,
       0.000000000625,
       0.00000000125,
       0.00000000625,
       0.0000000125,
       0.0000000625,
       0.000000125]

rs = [0.0000000001,
      0.0000000005,
      0.000000001,
      0.000000005,
      0.00000001,
      0.00000005,
      0.0000001]

thresholds = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 
              0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
              0.70, 0.75, 0.80, 0.85, 0.90, 0.95]


df = pandas.read_csv("archie_statistics_per_snp_loo.csv")

mets = ["TP", "FN", "FP", "TN",
        "TRUE_POS", "TRUE_NEG", "LABL_POS", "LABL_NEG", "ALL",
        "ACCURACY", "TNR", "FPR", "FNR",
        "RECALL", "PRECISION", "FDR", "FOR"]

cols = ["SET", "MU", "MU_TEST2TRAIN", "R", "R_TEST2TRAIN", "THRESHOLD"]
for met in mets:
    cols.append("AV_" + met)
    cols.append("SD_" + met)
    cols.append("SE_" + met)
    cols.append("NR_" + met)
rows = []

for i in range(1, 50):
    print("Processing %d" % i)
    for t in thresholds:
        xdf = df[(df["SET"] == i) & 
                 (df["THRESHOLD"] == t)]

        if xdf.shape[0] != 100:
            print(xdf)
            sys.exit("Wrong number of replicates!")

        row = {}
        row["SET"] = i
        row["MU"] = xdf.MU.values[0]
        row["R"] = xdf.R.values[0]
        row["MU_TEST2TRAIN"] = row["MU"] / 0.0000000125
        row["R_TEST2TRAIN"] = row["R"] / 0.00000001
        row["THRESHOLD"] = t

        for met in mets:
            vals = xdf[met].values
            vals = vals[~(numpy.isnan(vals))]
            if vals.size > 0:
                row["AV_" + met] = numpy.mean(vals)
                row["SD_" + met] = numpy.std(vals)
                row["SE_" + met] = row["SD_" + met] / numpy.sqrt(vals.size)
            row["NR_" + met] = vals.size

        rows.append(row)
            

df = pandas.DataFrame(rows, columns = cols)
df.to_csv("archie_statistics_per_snp_loo_summary.csv", 
          sep="," , index=None, na_rep="NA")
