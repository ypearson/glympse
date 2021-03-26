import os
import sys
import pandas as pd
from collections import Counter
import pprint
import csv

pp = pprint.PrettyPrinter(indent=4)

def reads_diff(ref_read, test_read):
    ref_alt     = []
    ref_alt_sub = []
    new_variant = True
    for i in range(0,len(test_read)):
        if not ref_read[i] == test_read[i]:
            if new_variant:
                new_variant = False
                ref_alt_sub = [ref_read[i], test_read[i], i]
            else:
                ref_alt_sub[0] = ref_alt_sub[0]+ref_read[i]
                ref_alt_sub[1] = ref_alt_sub[1]+test_read[i]
        else:
            if not new_variant:
                ref_alt += [ref_alt_sub]
                new_variant = True
    if not new_variant: # check last one
        ref_alt += [ref_alt_sub]
    return ref_alt

if __name__ == '__main__':

    READS=0
    COUNTS=1
    NOISE_THRESHOLD = 5
    CSV_TITLE = ['position','ref','alt','alt_cov']
    variants_dict = {}

     # Add missing environment variable
    os.environ['PUZZLE_DATA'] = './puzzle'

    with open(os.path.join(os.environ['PUZZLE_DATA'], "EGFR_reference.txt"), 'r') as f:
        reference = f.read()

    # Remove newline character
    reference = reference[:-1]

    # Get reads from file
    df_reads = pd.read_csv(os.path.join(os.environ['PUZZLE_DATA'], "reads.csv"), index_col='read_id')

    # Get list of unique DNA positions
    positions_unique = df_reads['position'].unique().tolist()
    positions_unique.sort()

    for i in range(0,len(positions_unique)):

        # Get list of sequences per position,
        df_reads_per_position = df_reads.sort_values(['position','sequence']).loc[df_reads['position'] == positions_unique[i]]
        reads_per_position = df_reads_per_position['sequence'].values.tolist()

        # Build histogram of reads per position and extract peak frequencies given a threshold of NOISE_THRESHOLD
        reads_at_position_histo = Counter(reads_per_position)
        reads_at_position_filtered     = [x     for x, count in reads_at_position_histo.items() if count >= NOISE_THRESHOLD]
        reads_at_position_filtered_cnt = [count for x, count in reads_at_position_histo.items() if count >= NOISE_THRESHOLD]

        # Build dictionary of each position and peak reads
        variants_dict[positions_unique[i]] = [reads_at_position_filtered,reads_at_position_filtered_cnt]

    variants_discovered = {}
    variants_discovered_reduced = {}

    # Find all the variants_dict of filtered reads
    for p in positions_unique:

        for v in range(0,len(variants_dict[p][READS])):

            p0=p
            p1=p0+len(variants_dict[p][READS][v])

            reference_sub = reference[p0:p1]

            read_variant = reads_diff(reference_sub,variants_dict[p][READS][v])
            if read_variant:
                if p0 in variants_discovered:
                    variants_discovered[p0] += [variants_dict[p][COUNTS][v], read_variant]
                else:
                    variants_discovered[p0] = [variants_dict[p][COUNTS][v], read_variant]

    # Reduce the discovered variants_dict dictionary, write to CSV (**NOTE: running out of time, could do better on this section!)
    csv_strings = []
    for p in variants_discovered:

        count   = sum([x for x in variants_discovered[p] if isinstance(x,int) ])
        variant = variants_discovered[p][1][0]

        csv_strings += [[
        str(p+variants_discovered[p][1][0][2]),
        variants_discovered[p][1][0][0],
        variants_discovered[p][1][0][1], str(count)]]

    with open('result.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(CSV_TITLE)
        csvwriter.writerows(csv_strings)

