import argparse
import numpy
import pandas
import sys


def main(args):
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

    threshold_keys = thresholds.keys()
    category_keys = ["TP", "FN", "FP", "TN"]

    replicate = {}
    for t in threshold_keys:
        replicate[t] = {}
        for c in category_keys:
           replicate[t][c] = 0

    with open(args.input, 'r') as fi:
        line = fi.readline()  # header

        for line in fi:
            data = line.split()

            pp = float(data[2])
            true_anc = int(data[3])

            if true_anc == 1:
                for t in threshold_keys:
                    if pp > thresholds[t]:
                       replicate[t]["TP"] = replicate[t]["TP"] + 1
                    else:
                       replicate[t]["FN"] = replicate[t]["FN"] + 1
            else:
                for t in threshold_keys:
                    if pp > thresholds[t]:
                        replicate[t]["FP"] = replicate[t]["FP"] + 1
                    else:
                        replicate[t]["TN"] = replicate[t]["TN"] + 1

        with open(args.output, 'w') as fo:
            fo.write("THRESHOLD %s\n" % " ".join(category_keys))
            for t in threshold_keys:
                fo.write("%s" % t)
                for c in category_keys: 
                    fo.write(" %d" % replicate[t][c])
                fo.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input file")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output file")
    main(parser.parse_args()) 

