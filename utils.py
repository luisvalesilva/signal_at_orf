
"""
    signal_at_orf utils
    ~~~~~~~~~~~~~~~~~~~~~
    Helper functions for signal_at_orf.

    :copyright: (c) 2016 by Luis Vale Silva.
    :license: MIT, see LICENSE for more details.
"""

__author__ = "Luis Vale Silva"
__status__ = "Development"

import sys
import os
import glob
import re
import time
import tqdm

import pandas as pd


def print_elapsed_time(t0):
    """Compute and print time elapsed from a start time (time.time() object).

    Parameters
    ----------
    t0: output of time.time()
        Start time to compute elapsed time from (no default)
    Returns
    -------
    Pandas data frame concatenating all read files"""
    elapsed_time = time.time() - t0

    if elapsed_time < 60:
        print("{:2.1f} sec.".format(elapsed_time))
    elif 60 < elapsed_time < 3600:
        print("{:2.1f} min.".format(elapsed_time / 60))
    else:
        print("{:2.1f} hr.".format(elapsed_time / 3600))


def read_wiggle(path, use_pbar=True):
    """Load wiggle files as dictionary of pandas data frames.

    Parameters
    ----------
    path: string
        Path to folder containing the wiggle files generated using the lab's
        ChIP-seq analysis pipelines
    use_pbar: bool, default True
        Whether to use a progress bar (requires tqdm package)
    Returns
    -------
    Dictionary of chromosome names and wiggle data (as pandas DataFrame)"""
    print('Reading wiggle data from:\n{}'.format(path))

    # Get all file names
    if os.path.isdir(path):
        all_files = glob.glob(path + '/*')
    else:
        sys.exit('Error: Incorrect path.')

    # Start dict to collect data
    out_dict = dict()
    cols = ['position', 'signal']

    if use_pbar:
        pbar = tqdm.tqdm(total=100)

    for file in all_files:
        # Skip the file containing all chromosomes, if present
        if re.search('all', file):
            continue

        # Get chr name string from file name
        chr_name = re.search('(?<=chr)[IXV]+', file)
        chr_name = 'chr' + chr_name.group(0)

        # Load data as pandas df and add to dict
        df = pd.read_table(file, header=None, names=cols, sep='\t', skiprows=2)
        out_dict[chr_name] = df

        if use_pbar:
            pbar.update(10)

    if use_pbar:
        pbar.close()

    return out_dict


def read_gff(path):
    """Load GFF file.
    TODO: Currently handles only our modified gff files; change to handle any gff

    Parameters
    ----------
    path: string
        Path to gff file
    Returns
    -------
    GFF as pandas DataFrame
    """

    col_names = ['seqname', 'source', 'feature', 'start', 'end',
                 'score', 'strand', 'frame', 'attribute']

    gff = pd.read_table(path, names=col_names, comment='#', engine='python')
    print('Reading gff file:\n{}'.format(path))

    return gff

def check_genome(input_str):
    """Determine reference genome (S288c or SK1), by comparing input to elements
    in the list of chromosome names/numbers in the reference genome annotations.

    Parameters
    ----------
    path: string
        Path to gff file (no default)
    Returns
    -------
    Reference genome as one of two strings: 'S288c' or 'SK1'"""
    S288C = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
             'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']
    SK1 = ['01', '02', '03', '04', '05', '06', '07', '08', '09',
           '10', '11', '12', '13', '14', '15', '16']

    if any('chr'+chr in input_str for chr in S288C):
        return 'S288C'
    elif any('chr'+chr in input_str for chr in SK1):
        return 'SK1'
    else:
        print("Error: cannot determine reference genome.")
