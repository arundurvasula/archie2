import argparse
import numpy
import random


def main(args):
    seeds = numpy.random.randint(1, 2147483647, size=args.size)
    with open(args.output, 'w') as f:
        f.write('\n'.join([str(x) for x in seeds]))
        f.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--size", type=int,
                        default=100, help="Number of seeds")
    parser.add_argument("-o", "--output", type=str,
                        default="seeds.txt", help="Output file")
    main(parser.parse_args()) 

