"""
Benchmarking module for PhageFilter.

Install / Dependencies:
    1. Yaml - https://pyyaml.org/wiki/PyYAMLDocumentation
    2. Kraken2 - https://github.com/DerrickWood/kraken2/wiki
    3. FastViromeExplorer - https://code.vt.edu/saima5/FastViromeExplorer
        3.A - install - javac -d bin src/*.java
        3.B - get SamTools - https://fastviromeexplorer.readthedocs.io/en/latest/
        3.C - get Kallisto - https://fastviromeexplorer.readthedocs.io/en/latest/
    4. For PhageFilter, make sure to build on release mode before running

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

# custom libraries
from bench.tools.phage_filter import PhageFilter
import bench.benchmarking_tests as bench_test

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
    subparser.add_argument("-db", "--database_name", default=Path("tree1/"), nargs='?', type=Path, help="path to the DB to use [default 'tree2/']", required=False)
    subparser.add_argument("-k", "--kmer_size", default=20, nargs='?', type=int, help="size of kmer to use for PhageFilter [default 20]", required=False)
    subparser.add_argument("--threads", default=4, nargs='?', type=int, help="number of threads to use [Default 4]", required=False)

def parseArgs(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(help='Choose type of benchmarking', dest='sub_parser')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")

    # if performance testing
    performance_testing_parser = subparsers.add_parser(SubparserNames.performance_testing.value)
    add_common_arguments(performance_testing_parser)

    # if parameterization
    parameterization_parser = subparsers.add_parser(SubparserNames.parameterization.value)
    add_common_arguments(parameterization_parser)

    # if genomecount
    genomecount_parser = subparsers.add_parser(SubparserNames.genomecount.value)
    add_common_arguments(genomecount_parser)

    # if relative_performance
    relative_performance_parser = subparsers.add_parser(SubparserNames.relative_performance.value)
    add_common_arguments(relative_performance_parser)

    # if filter_performance
    filter_performance_parser = subparsers.add_parser(SubparserNames.filter_performance.value)
    add_common_arguments(filter_performance_parser)

    # if memory & time bench
    memory_parser = subparsers.add_parser(SubparserNames.memory.value)
    add_common_arguments(memory_parser)

    # readlength test
    readlength_parser = subparsers.add_parser(SubparserNames.readlength.value)
    add_common_arguments(readlength_parser)

    # threads test
    threads_parser = subparsers.add_parser(SubparserNames.threads.value)
    add_common_arguments(threads_parser)

    # readlength test
    depth_parser = subparsers.add_parser(SubparserNames.depth.value)
    add_common_arguments(depth_parser)

    return parser.parse_args(argv)


def main():

    def performance_testing_action():
        print(f"Performing parameterization benchmarking...")
        PerformanceBenchmark().run(pos_genome_path=args.genome_dir,
                                    neg_genome_path=args.neg_genome_dir,
                                    config=args.config,
                                    result_csv=args.result_csv, 
                                    variation_count=10)
        
    def parameterization_action():
        print(f"Performing parameterization benchmarking...")
        ParameterBenchmark().run(pos_genome_path=args.genome_dir,
                                 neg_genome_path=args.neg_genome_dir,
                                 config=args.config, 
                                 result_csv=args.result_csv, 
                                 read_count=1000)

    def genomecount_action():
        print(f"Performing genome count benchmarking...")
        GenomeCountBenchmark().run(pos_genome_path=args.genome_dir,
                                    neg_genome_path=None,
                                    config=args.config,
                                    result_csv=args.result_csv)

    def readlength_action():
        print(f"Performing read length benchmarking...")
        ReadLengthBenchmark().run(pos_genome_path=args.genome_dir,
                                  neg_genome_path=args.neg_genome_dir,
                                  config=args.config,  
                                  result_csv=args.result_csv,
                                  variation_count = 2)

    def threads_action():
        print(f"Performing thread benchmarking...")
        ThreadsBenchmarking().run(pos_genome_path=args.genome_dir,
                                  neg_genome_path=args.neg_genome_dir,
                                  config=args.config, 
                                  result_csv=args.result_csv,
                                  variation_count = 3, 
                                  contamination_fraction = 0.5, 
                                  read_count=100000)
                                        
    
    def relative_performance_action():
        print(f"Performing relative performance benchmarking...")
        TaxonomicBenchmark().run(pos_genome_path=args.genome_dir,
                                 neg_genome_path=args.neg_genome_dir,
                                 config=args.config, 
                                 result_csv=args.result_csv,
                                 contamination_fraction = 0.5,
                                 read_count=100)

    def filter_performance_action():
        print(f"Performing filter performance benchmarking...")
        FilterBenchmark().run(pos_genome_path=args.genome_dir,
                              neg_genome_path=args.neg_genome_dir,
                              config=args.config,
                              result_csv=args.result_csv,
                              read_count=100)

    def memory_action():
        print(f"Performing filter memory & time benchmarking...")
        MemoryBenchmark().run(
            pos_genome_path=args.genome_dir, 
            neg_genome_path=args.neg_genome_dir,
            config=args.config,
            result_csv=args.result_csv, 
            contamination_fraction = 0.2,
            variation_count = 10
        )

    def depth_action():
        print(f"Performing tree depth benchmarking...")
        TunableDepthBenchmark().run(pos_genome_path=args.genome_dir, 
                                   neg_genome_path=args.neg_genome_dir, 
                                   config=args.config, 
                                   result_csv=args.result_csv, 
                                   read_count=100000,
                                   contamination_fraction=0.90
        )

    # arguments
    args = parseArgs(sys.argv[1:])

    # map subparsers to their corresponding functions and messages
    actions = {
        SubparserNames.parameterization.value: parameterization_action,
        SubparserNames.genomecount.value: genomecount_action,
        SubparserNames.relative_performance.value: relative_performance_action,
        SubparserNames.performance_testing.value: performance_testing_action,
        SubparserNames.filter_performance.value: filter_performance_action,
        SubparserNames.memory.value: memory_action,
        SubparserNames.readlength.value: readlength_action,
        SubparserNames.threads.value: threads_action,
        SubparserNames.depth.value: depth_action
    }

    # run benchmark type specified
    action = actions.get(args.sub_parser)
    if action:
        action()
    else:
        print(__doc__)

if __name__ == '__main__':
    main()
