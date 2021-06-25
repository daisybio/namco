#!/usr/bin/env python3
import argparse
import collections
import gzip
import multiprocessing
import os
import shutil
import subprocess
import sys
from subprocess import Popen, PIPE
from typing import List, Dict

# Development of TUM CF Microbiome/NGS

#bowtie_index_dir = "/usr/local/bin/bowtie/spikesIndex"
bowtie_index_dir = "/tmp/spikesIndex"

def log_to_status_file(msg):
    with open('/var/log/shiny-server/status_rm_spikes.txt', 'w') as stat_f_h:
        stat_f_h.write(msg)

log_to_status_file('Starting spike removal')

ap = argparse.ArgumentParser()
ap.add_argument('ref', type=str, help='Fasta file with spikes reference sequences')
ap.add_argument('mapping', type=str, help='Mapping File')
ap.add_argument('inp', type=str, help='Input folder')
ap.add_argument('osamples', type=str, help='Output folder for samples')
ap.add_argument('ospikes', type=str, help='Output folder for spikes')
ap.add_argument('ostats', type=str, help='Output file for statistics')
ap.add_argument('omapping', type=str, help='Output file for reduced mapping file')
ap.add_argument('cores', type=int, help='Number of cores to use')
args = ap.parse_args()

spikes_ref_fa = args.ref  # fasta file for building an index and to align against; spike ins
mapping_file = args.mapping  # mapping file, listing relevant files
input_folder = args.inp  # fastq files
output_folder_samples = args.osamples  # output folder for samples
output_folder_spikes = args.ospikes  # output folder for spikes
output_file_spikes_counts = args.ostats  # output for stats (readcount of spikes)
output_file_reduced_mapping = args.omapping  # output for reduced mapping file
ncores = args.cores

os.makedirs(output_folder_samples, mode=0o777, exist_ok=True)
os.makedirs(output_folder_spikes, mode=0o777, exist_ok=True)

Cmd = List[str]


def call(cmd: Cmd) -> list:
    print("> Executing " + " ".join(cmd))
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout = []
    stderr = []
    with process.stdout as pipe:
        for line in iter(pipe.readline, b''):
            output = line.decode('utf-8')
            stdout.append(output)
    with process.stderr as pipe:
        for line in iter(pipe.readline, b''):
            output = line.decode('utf-8')
            stderr.append(output)
    return [stdout, stderr]


FastqPair = collections.namedtuple('FastqPair', ['R1', 'R2'])


def calc_spikes(r1_and_r2_file: FastqPair, mapping_lines, col_sample, col_weight, col_amount):
    fastq_R1 = r1_and_r2_file.R1
    fastq_R2 = r1_and_r2_file.R2

    # generate names for the output
    name_split = os.path.basename(fastq_R1).split("_")
    orig_sample_id = name_split[0]

    entry_line_in_mapping_file = mapping_lines[orig_sample_id]
    weight = entry_line_in_mapping_file[col_weight].strip()
    amount = entry_line_in_mapping_file[col_amount].strip()

    if amount != "0":
        name_split[0] += "_withoutSpikes"

    sample = '_'.join(name_split).replace("_R1_", "_R%_")
    fastq_unaligned = os.path.join(output_folder_samples, sample)
    fastq_aligned = os.path.join(output_folder_spikes,
                                 os.path.basename(fastq_R1).split("_R1_")[0].split("_")[0] + "_spikes_R%.fastq"
                                 )

    # check if we have a no spikes sample
    if amount == "0":
        # copy files to samples and create empty files in spikes
        shutil.copy(fastq_R1, fastq_unaligned.replace("%", "1"))
        shutil.copy(fastq_R2, fastq_unaligned.replace("%", "2"))
        with open(fastq_aligned.replace("%", "1"), "w"):
            pass
        with open(fastq_aligned.replace("%", "2"), "w"):
            pass
    else:
        # run bowtie2 to find spikes
        cmd = ["bowtie2", "-x", bowtie_index_dir, "-1", fastq_R1, "-2", fastq_R2, "--al-conc", fastq_aligned, "--un-conc",
               fastq_unaligned]
        stdout, stderr = call(cmd)
        # print(stdout)
        sys.stderr.write("\n".join(stderr))
        sys.stderr.flush()

    # count reads of spikes
    line_count = 0
    with open(fastq_aligned.replace("%", "1"), "r") as al:
        for _ in al:
            line_count += 1

    new_sample_id = sample.split("_")[0]
    with open(output_file_spikes_counts, 'a') as stats_h:
        stats_h.write(f'{new_sample_id}\t{str(line_count // 4)}\t{str(weight)}\t{str(amount)}\n')

    with open(output_file_reduced_mapping, 'a') as mapping_out_h:
        entry_line_in_mapping_file[col_sample] = new_sample_id
        mapping_out_h.write('\t'.join(entry_line_in_mapping_file).strip() + "\n")

    # os.remove(fastq_aligned.replace("%", "1"))
    # os.remove(fastq_aligned.replace("%", "2"))


def gunzip(file: str):
    new_filename = file.replace('.gz', '')
    with gzip.open(file, 'rb') as f_in:
        with open(new_filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return new_filename


def onexit():
    # correct rights of output folders
    subprocess.call(['chmod', '777', '-R', os.path.abspath(output_folder_samples)])
    subprocess.call(['chmod', '777', '-R', os.path.abspath(output_folder_spikes)])


# process mapping file to get valid IDs
valid_ids = []
ignored_ids = []  # e.g. entries with NA values
mapping_lines: Dict[str, List[str]] = {}  # from mapping file, full line for each relevant sample

index_of_sample_id_col = None
sample_col_name = "SampleID"
index_of_weight_col = None
weight_col_name = "total_weight_in_g"
index_of_amounts_col = None
amounts_col_name = "amount_spike"
with open(mapping_file, 'r') as mapping_file_h:
    header: str = next(mapping_file_h)

    if not header.startswith('#'):
        raise Exception('No header in mapping file')

    columns = header.split('\t')
    for i, col in enumerate(columns):
        colname = col.strip().lstrip('#')

        if colname == sample_col_name:
            index_of_sample_id_col = i
        elif colname == weight_col_name:
            index_of_weight_col = i
        elif colname == amounts_col_name:
            index_of_amounts_col = i

    if index_of_sample_id_col is None:
        onexit()
        raise Exception(f'Column {sample_col_name!r} missing in header')

    if index_of_weight_col is None:
        onexit()
        raise Exception(f'Column {weight_col_name!r} missing in header')

    if index_of_amounts_col is None:
        onexit()
        raise Exception(f'Column {amounts_col_name!r} missing in header')

    for line in mapping_file_h:
        if line.startswith('#'):
            continue

        fields = line.split('\t')
        sample_id = fields[index_of_sample_id_col].strip()
        total_weight_in_g = fields[index_of_weight_col].strip()
        amount = fields[index_of_amounts_col].strip()

        # warn if weight is not a valid number
        try:
            float(total_weight_in_g)
        except ValueError as e:
            print(
                f'Weight {total_weight_in_g!r} is not a valid floating point number for {sample_id!r}. Sample will be processed as no spikes',
                file=sys.stderr)
            fields[index_of_amounts_col] = "0"

        # if amount cannot be parsed to float change to 0
        try:
            float(amount)
        except ValueError as e:
            fields[index_of_amounts_col] = "0"

        valid_ids.append(sample_id)
        mapping_lines[sample_id] = fields

# building bowtie index from spikes file
spikes_ref_name = bowtie_index_dir
if not os.path.isdir(bowtie_index_dir):
    os.mkdir(bowtie_index_dir)
cmd = ["bowtie2-build", "-p", spikes_ref_fa, spikes_ref_name]
print(call(cmd))

log_to_status_file('Detecting spike sequences and writing filtered mapping file')
# write stats and reduced mapping header.
with open(output_file_spikes_counts, 'w') as stats_h:
    stats_h.write(f'#SampleID\tSpikeReads\tspikes_total_weight_in_g\tamount_spike\n')

with open(mapping_file, 'r') as mapping_file_h:
    header: str = next(mapping_file_h)
    assert (header.startswith('#'))
    with open(output_file_reduced_mapping, 'w') as mapping_out_h:
        mapping_out_h.write(header)

files_to_process = set()
for file in os.listdir(input_folder):
    if file.endswith(".fastq.gz"):
        sample_name, s1, l001, read, filenameend = file.split('_')
        if sample_name not in valid_ids:
            continue
        expected_paired_filename = f'{sample_name}_{s1}_{l001}_{"R1" if read == "R2" else "R2"}_{filenameend}'
        if not os.path.exists(os.path.join(input_folder, expected_paired_filename)):
            print(
                f'Paired file to {file} does not exist. Expected the file to be named {expected_paired_filename}. Ignoring pair...')
            continue
        if file.endswith(".gz") and not os.path.exists(os.path.join(input_folder, file.replace('.gz', ''))):
            file = gunzip(os.path.join(input_folder, file))
        if expected_paired_filename.endswith(".gz") and not os.path.exists(
                os.path.join(input_folder, expected_paired_filename.replace('.gz', ''))):
            expected_paired_filename = gunzip(os.path.join(input_folder, expected_paired_filename))

for file in os.listdir(input_folder):
    if file.endswith(".fastq"):
        sample_name, s1, l001, read, filenameend = file.split('_')
        if read == "R1":
            if sample_name not in valid_ids:
                print(f'Ignore file: {sample_name}', file=sys.stderr)
                continue
            expected_paired_filename = f'{sample_name}_{s1}_{l001}_{"R1" if read == "R2" else "R2"}_{filenameend}'
            if not os.path.exists(os.path.join(input_folder, expected_paired_filename)):
                print(f'Paired file to {file} does not exist. Expected the file to be named {expected_paired_filename}')
                continue
            files_to_process.add(FastqPair(os.path.join(input_folder, file),
                                           os.path.join(input_folder, expected_paired_filename)))

# calculate spikes on multiple cpu cores
#pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
pool = multiprocessing.Pool(processes=ncores)
child_processes = [
    pool.apply_async(calc_spikes,
                     (filepair, mapping_lines, index_of_sample_id_col, index_of_weight_col, index_of_amounts_col)) for
    filepair in
    files_to_process]
print([res.get() for res in child_processes])

os.rmdir(bowtie_index_dir)
onexit()
