#!/usr/bin/env python

"""Plot the output of signal_at_orf.py, after grouping by position and calculating the mean"""

import pandas as pd
import pylab as pl

def main():
    print()
    print("---------------------------------------------------------")
    print("Computing average signal by position for data in file(s):")
    print(args.input_data_a)
    if args.input_data_b:
        print(args.input_data_b)
    print()

    data_a = pd.read_csv(args.input_data_a, sep='\t')
    grouped_a = data_a.groupby(['position'], as_index=False).mean()
    xa = grouped_a['position']
    ya = grouped_a['signal']

    if args.input_data_b:
        data_b = pd.read_csv(args.input_data_b, sep='\t')
        grouped_b = data_b.groupby(['position'], as_index=False).mean()
        xb = grouped_b['position']
        yb = grouped_b['signal']

    print('Plotting...')
    print("---------------------------------------------------------")

    pl.clf()
    pl.plot(xa, ya, label=args.input_data_a)
    if args.input_data_b:
        pl.plot(xb, yb, label=args.input_data_b)
    pl.legend()
    pl.show()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Plot output of signal_at_orf.py")
    parser.add_argument('-a', '--input_data_a', help=("path to a TSV file"), required=True)
    parser.add_argument('-b', '--input_data_b', help=("path to a TSV file"))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
    
    # Also print help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    main()
