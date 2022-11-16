"""
Used for testing the output of PhageFilter.

Usage:
    python benchmarking/test_output.py <simulated reads> <output file>

Example:
    python benchmarking/test_output.py examples/test_reads/simulated_reads.fa genomes_in_file.txt 
"""
import sys
import numpy as np




def get_true_maps(fasta_read_path):
    """
    returns a dictionary of true genomes and read counts
    """
    name2counts = {}
    with open(fasta_read_path) as read_file:
        line = read_file.readline()
        count = 0
        while line:
            if line[0] == ">":
                name = "_".join(line.strip(">").strip("\n").split("_")[0:-1])
                if name in name2counts:
                    name2counts[name] += 1
                else:
                    name2counts[name] = 1
            line = read_file.readline()
    return name2counts

def get_maps(output_path):
    """
    returns a dictionary of the output of PhageFilter.
    """
    name2counts = {}
    with open(output_path) as out_file:
        line = out_file.readline()
        count = 0
        while line:
            name, count = line.strip("\n").split(",")
            name2counts[name] = int(count)
            line = out_file.readline()
    return name2counts

def get_classification_metrics(true_map, out_map):
    """
    uses true map and PhageFilter map to obtain metrics
    of the classification accuracy.
    """
    TP, FP, FN = (0,0,0)
    for genome in out_map.keys():
        if genome in true_map.keys():
            TP += 1
        else:
            FP += 1
    for genome in true_map.keys():
        if genome not in out_map.keys():
            FN += 1
    # get recall and precision
    recall = TP / (TP + FN)
    precision = TP / (TP + FP)
    return recall, precision

def get_readcount_metrics(true_map, out_map):
    """
    uses true map and PhageFilter map to obtain metrics
    of the read count error.
    """
    absolute_difference = []
    for genome, count in out_map.items():
        if genome in true_map.keys():
            absolute_difference.append(abs(count - true_map[genome]))
    return absolute_difference

if __name__ == '__main__':
    fasta_read_path = sys.argv[1]
    result_read_path = sys.argv[2]

    # parse files
    true_map = get_true_maps(fasta_read_path)
    out_map = get_maps(result_read_path)

    # get classification metrics
    recall, precision = get_classification_metrics(true_map, out_map)
    print(f"recall: {recall}, precision: {precision}")

    # get read count/map metrics
    absolute_difference_array = get_readcount_metrics(true_map, out_map)
    max_error = max(absolute_difference_array)
    min_error = min(absolute_difference_array)
    average_error = np.average(absolute_difference_array)
    print(f"average read count error: {average_error}, max error: {max_error}, min error: {min_error}")
