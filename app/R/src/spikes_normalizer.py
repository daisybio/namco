#!/usr/bin/env python3

"""Normalize an OTU Table with spike readcount and weight"""
import argparse
import logging
import os
import subprocess

import pandas as pd



def log_to_status_file(msg):
    with open('/var/log/shiny-server/status.txt', 'w') as stat_f_h:
        stat_f_h.write(msg)

log_to_status_file('Normalize OTU table for spike reads and weights')

parser = argparse.ArgumentParser()
parser.add_argument("otu_table_path", help="Path of original OTU Table. Output will be written to same folder", type=str)
parser.add_argument("spikes_stats_path", help="Path to spike stats file.", type=str)
args = parser.parse_args()

otu_table_path = args.otu_table_path
spikes_stats_path = args.spikes_stats_path



samples = {}
with open(spikes_stats_path, 'r') as stats_h:
    for line in stats_h:
        if line.startswith("#"):
            assert('#SampleID\tSpikeReads\tspikes_total_weight_in_g\tamount_spike' in line.strip())
        elif line == "\n":
            continue
        else:
            fields = line.strip().split("\t")
            sample_id = fields[0]
            spike_reads = int(fields[1])
            try:
                original_total_weight_in_g = float(fields[2])
            except ValueError:
                original_total_weight_in_g = 0
            amount = float(fields[3])
            samples[sample_id] = (spike_reads, original_total_weight_in_g, amount)


sum = 0
n = 0
for values in samples.values():
    count, _, amount = values
    if amount != 0:
        sum += count
        n += 1

mean = (sum / n)
logging.debug(f'mean={mean}')


otu_table = pd.read_csv(otu_table_path, delimiter="\t", index_col=0)
print(otu_table)
for file_id, values in samples.items():
    count, weight, amount = values
    thismean = mean
    factor = None
    try:
        if amount == 0 or count == 0:  # no spike sample
            # normalize to 10_000
            row_sum = otu_table[file_id].sum()
            otu_table[file_id] = otu_table[file_id] * (10000 / row_sum)
        else:
            # do spike normalization
            factor = 600 / (amount * 100)
            otu_table[file_id] = (otu_table[file_id] * thismean) / (count * weight * factor)
    except KeyError:
        print(f"{file_id} is in mapping file but not in OTU Table")


normalized_otu_path = os.path.join(os.path.dirname(otu_table_path), "SpikeNormalized-OTUs-Table.tab")
otu_table.to_csv(normalized_otu_path, sep="\t")
print(f'Wrote normalized OTU Table to {normalized_otu_path}')

# correct rights of output folders
subprocess.call(['chmod', '777', '-R', os.path.abspath(normalized_otu_path)])
