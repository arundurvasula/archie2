import pandas
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

path = "../data/testing-data"

dfs = []

i = 1
for r in rs:
    for mu in mus:
        modl = "msmodified.set-" + str(i) + ".mu-" + mu + ".r-" + r
        for repl in [str(x) for x in range(1, 101)]:
            rdir = repl.zfill(4)
            fnam = path + '/' + modl + '/' + rdir + "/archie_statistics_per_snp.txt"
            df = pandas.read_csv(fnam, sep=" ")
            df.insert(loc=0, column="REPL", value=[repl] * 18)
            df.insert(loc=0, column="R", value=[r] * 18)
            df.insert(loc=0, column="MU", value=[mu] * 18)
            df.insert(loc=0, column="SET", value=[i] * 18)
            dfs.append(df)

        i += 1


df = pandas.concat(dfs, sort=False)
df.to_csv(path + "/archie_statistics_per_snp.csv", 
          sep=",", index=None, na_rep="NA")
