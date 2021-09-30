import argparse
import numpy
import pandas
import sys


def loo_test(x1, x2):
    n1 = len(x1)
    n2 = len(x2)
    if n1 != n2:
        sys.exit("Error: Input arrays have different sizes!\n")

    fracs = []
    for i in range(n1):
        tmp1 = x1[i]
        tmp2 = x2[i]

        x1[i] = 0
        x2[i] = 0
        sum1 = numpy.sum(x1)
        sum2 = numpy.sum(x2)

        if sum1 + sum2 > 0:
            frac = sum1 / (sum1 + sum2)
            fracs.append(frac)

        x1[i] = tmp1
        x2[i] = tmp2

    fracs = numpy.array(fracs)
    n = len(fracs)

    avg_frac = numpy.mean(fracs)
    mysum = 0
    for frac in fracs:
        mysum = mysum + ((frac - avg_frac) ** 2)
    var_frac = numpy.sqrt( ((n - 1.0) / n) * mysum)

    return [avg_frac, var_frac, n]


def main(args):
    # Read inputs
    nins = len(args.inputs)
    dfs = []
    for i in range(nins):
        df = pandas.read_csv(args.inputs[i], sep=" ")
        [nrow, ncol] = df.shape
        df.insert(0, "REPL", [i] * nrow)
        dfs.append(df)

    df = pandas.concat(dfs, ignore_index=True)

    thresholds = list(set(df.THRESHOLD.values))

    cols = ["THRESHOLD", "NREPS", "PRECISION", "RECALL",
            "LOO_AVG_PRECISION", "LOO_VAR_PRECISION", "LOO_N_PRECISION",
            "LOO_AVG_RECALL", "LOO_VAR_RECALL", "LOO_N_RECALL"]
    rows = []

    for threshold in thresholds:
        xdf = df[df["THRESHOLD"] == threshold]

        tps = []
        fns = []
        fps = []
        tns = []

        for i in range(nins):
            ydf = xdf[xdf["REPL"] == i]

            tps.append(ydf.TP.values[0])
            fns.append(ydf.FN.values[0])
            fps.append(ydf.FP.values[0])
            tns.append(ydf.TN.values[0])

        tps = numpy.array(tps)
        fns = numpy.array(fns)
        fps = numpy.array(fps)
        tns = numpy.array(tns)
        n = len(tps)

        # Overall all replicates
        tp = numpy.sum(tps)
        fn = numpy.sum(fns)
        fp = numpy.sum(fps)
        tn = numpy.sum(tns)

        pos = tp + fn
        if pos > 0:
            recall = tp / pos
            [loo_avg_r, loo_var_r, loo_n_r] = loo_test(tps, fns)
        else:
            recall = "NA"
            loo_avg_r = "NA"
            loo_var_r = "NA"
            loo_n_r = 0

        pos = tp + fp
        if pos > 0:
            precision = tp / (tp + fp)
            [loo_avg_p, loo_var_p, loo_n_p] = loo_test(tps, fps)
        else:
            precision = "NA"
            loo_avg_p = "NA"
            loo_var_p = "NA"
            loo_n_p = 0

        row = {}
        row["THRESHOLD"] = threshold

        row["NREPS"] = n
        row["PRECISION"] = precision
        row["RECALL"] = recall

        row["LOO_AVG_PRECISION"] = loo_avg_p
        row["LOO_VAR_PRECISION"] = loo_var_p
        row["LOO_N_PRECISION"] = loo_n_p

        row["LOO_AVG_RECALL"] = loo_avg_r
        row["LOO_VAR_RECALL"] = loo_var_r
        row["LOO_N_RECALL"] = loo_n_r

        rows.append(row)

    df = pandas.DataFrame(rows, columns = cols)
    df.to_csv(args.output, sep=" ", index=None, na_rep="NA")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", type=str, nargs='+', required=True,
                        help="Input files (one per replicate)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output file")
    main(parser.parse_args()) 

