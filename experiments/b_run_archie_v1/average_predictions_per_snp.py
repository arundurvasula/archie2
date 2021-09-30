import argparse
import numpy
import pandas
import sys


def main(args):
    # Read predictions per window
    df = pandas.read_csv(args.input, sep=" ")

    # Read SNP positions
    rows = []
    with open(args.snp, 'r') as fsnp, open(args.anc, 'r') as fanc, open(args.output, 'w') as fout:
        fout.write("SNP HAP MEAN_PRED TRUE_ANC\n")

        for lsnp in fsnp:
            lanc = fanc.readline().rstrip()

            snp = int(lsnp.split()[3])
            anc = list(lanc)
 
            xdf = df[(df["WSTART"] <= snp) & (df["WEND"] >= snp)]
            wstarts = numpy.unique(xdf.WSTART.values)
            nwindows = len(wstarts)
            nhaps = int(xdf.shape[0] / nwindows)

            preds = numpy.zeros(nhaps)
            for wstart in wstarts:
                preds += numpy.array(xdf[xdf["WSTART"] == wstart].PREDICTION.values)
            preds = preds / float(nwindows)

            haps = xdf[xdf["WSTART"] == wstarts[0]].HAP.values

            for i in range(nhaps):
                fout.write("%d %s %1.12f %s\n" % (snp, haps[i], preds[i], anc[i]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input file with prediction per window")
    parser.add_argument("-s", "--snp", type=str, required=True,
                        help="File with SNP information (.snp)")
    parser.add_argument("-a", "--anc", type=str, required=False,
                        help="File with ancestry (.anc)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output file")
    main(parser.parse_args()) 

