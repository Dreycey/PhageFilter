"""
Benchmarking module for PhageFilter.

Use `python3 benchmarking/bench.py -h` for options. Further info and *examples* for each benchmark can be accessed using `python3 benchmarking/bench.py <benchmark_name> -h`.
"""
# standard libraries
from enum import Enum
import sys
from pathlib import Path
import argparse

# importing benchmarking tests
from bench.benchmarking_tests import (ParameterBenchmark, TaxonomicBenchmark, FilterBenchmark, PerformanceBenchmark,
                                      GenomeCountBenchmark, ReadLengthBenchmark, ThreadsBenchmarking, TunableDepthBenchmark,
                                      MemoryBenchmark)



class SubparserNames(Enum):
    performance_testing = "performance_testing"
    parameterization = "parameterization"
    genomecount = "genomecount"
    relative_performance = "relative_performance"
    filter_performance = "filter_performance"
    memory = "memory"
    readlength = "readlength"
    threads = "threads"
    depth = "depth"

def add_common_arguments(subparser):
    subparser.add_argument("-g", "--genome_dir", type=Path, help="Path to the genome directory", required=True)
    subparser.add_argument("-n", "--neg_genome_dir", type=Path, help="path to the contamination genome directory", required=True)
    subparser.add_argument("-c", "--config", type=Path, help="path to the configuration file", required=True)
    subparser.add_argument("-r", "--result_csv", type=Path, help="Path to the result CSV (output)", required=True)

def parseArgs(argv=None) -> argparse.Namespace:
    # Main parser description
    parser = argparse.ArgumentParser(description="Benchmarking module for PhageFilter. For detailed installation and dependency information, refer to the project documentation or README.")

    subparsers = parser.add_subparsers(help='Choose type of benchmarking', dest='sub_parser')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")

    # Example commands for each subparser
    examples = {
        SubparserNames.performance_testing.value: "python3 benchmarking/bench.py performance_testing -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_performance_benchmarking.csv",
        SubparserNames.genomecount.value: "python3 benchmarking/bench.py genomecount -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_genomes.csv",
        SubparserNames.threads.value: "python3 benchmarking/bench.py threads -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_threading.csv",
        SubparserNames.readlength.value: "python3 benchmarking/bench.py readlength -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_readlength.csv",
        SubparserNames.parameterization.value: "python3 benchmarking/bench.py parameterization -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_parameterization.csv",
        SubparserNames.relative_performance.value: "python3 benchmarking/bench.py relative_performance -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_relative_performance.csv",
        SubparserNames.filter_performance.value: "python3 benchmarking/bench.py filter_performance -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_filter_performance.csv",
        SubparserNames.memory.value: "python3 benchmarking/bench.py memory -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_filter_memory.csv",
        SubparserNames.depth.value: "python3 benchmarking/bench.py depth -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_depth.csv"
    }

    for name in SubparserNames:
        subparser = subparsers.add_parser(name.value, description=f"Example command:\n{examples[name.value]}")
        add_common_arguments(subparser)

    return parser.parse_args(argv)

def run_benchmark(args, benchmark_class, message):
    print(f"Performing {message} benchmarking...")
    benchmark_instance = benchmark_class()
    kwargs = {
        "pos_genome_path": args.genome_dir,
        "neg_genome_path": args.neg_genome_dir,
        "config": args.config,
        "result_csv": args.result_csv
    }
    benchmark_instance.run(**kwargs)


def main():
    # arguments
    args = parseArgs(sys.argv[1:])

    benchmarks = {
        SubparserNames.parameterization.value: (ParameterBenchmark, "parameterization"),
        SubparserNames.genomecount.value: (GenomeCountBenchmark, "genome count"),
        SubparserNames.relative_performance.value: (TaxonomicBenchmark, "relative performance"),
        SubparserNames.performance_testing.value: (PerformanceBenchmark, "performance testing"),
        SubparserNames.filter_performance.value: (FilterBenchmark, "filter performance"),
        SubparserNames.memory.value: (MemoryBenchmark, "filter memory & time"),
        SubparserNames.readlength.value: (ReadLengthBenchmark, "read length"),
        SubparserNames.threads.value: (ThreadsBenchmarking, "thread"),
        SubparserNames.depth.value: (TunableDepthBenchmark, "tree depth")
    }

    benchmark_class, message = benchmarks.get(args.sub_parser, (None, None))
    if benchmark_class:
        run_benchmark(args, benchmark_class, message)
    else:
        print(__doc__)


if __name__ == '__main__':
    main()
