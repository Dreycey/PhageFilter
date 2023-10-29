"""
Benchmarking module for PhageFilter.

Examples:
* performance benchmarking (when optimizing performance..)
```
python3 benchmarking/bench.py performance_testing -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_performance_benchmarking.csv
```

* running genome count benchmarking
```
python3 benchmarking/bench.py genomecount -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_genomes.csv
```

* perform threading tests
```
python3 benchmarking/bench.py threads -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_threading.csv
```

* read length benchmarking
```
python3 benchmarking/bench.py readlength -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_readlength.csv
```

* running parameterization benchmarking
```
python3 benchmarking/bench.py parameterization -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_parameterization.csv
```

* classification benchmarking
```
python3 benchmarking/bench.py relative_performance -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_relative_performance.csv
```

* filter performance benchmarking
```
python3 benchmarking/bench.py filter_performance -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_filter_performance.csv
```

* filter memory benchmarking
```
python3 benchmarking/bench.py memory -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_filter_memory.csv
```

* depth testing
```
python3 benchmarking/bench.py depth -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_depth.csv
```
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
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(help='Choose type of benchmarking', dest='sub_parser')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")

    # List of all subparser names
    subparser_names = [
        SubparserNames.performance_testing.value,
        SubparserNames.parameterization.value,
        SubparserNames.genomecount.value,
        SubparserNames.relative_performance.value,
        SubparserNames.filter_performance.value,
        SubparserNames.memory.value,
        SubparserNames.readlength.value,
        SubparserNames.threads.value,
        SubparserNames.depth.value
    ]

    for name in subparser_names:
        subparser = subparsers.add_parser(name)
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
