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
* running genome count benchmarking
```
python benchmarking/bench.py genomecount -g examples/genomes/viral_genome_dir/ -r res_genomes.csv
```

* running parameterization benchmarking
```
python benchmarking/bench.py parameterization -g examples/genomes/viral_genome_dir/ -t examples/test_reads/ -r res_parameterization.csv
```

* relative performance benchmarking
```
python benchmarking/bench.py relative_performance -g examples/genomes/viral_genome_dir/ -r res_relative_performance.csv -t examples/test_reads/ -c benchmarking/config.yaml
```
"""
# standard libraries
from enum import Enum
import sys
from pathlib import Path
import argparse
# custom libraries
from bench.tools import PhageFilter
import bench.benchmarking_tests as bench_test




class SubparserNames(Enum):
    parameterization = "parameterization"
    genomecount = "genomecount"
    relative_performance = "relative_performance"

def add_common_arguments(subparser):
    subparser.add_argument("-g", "--genome_dir", type=Path, help="Path to the genome directory", required=True)
    subparser.add_argument("-r", "--result_csv", type=Path, help="Path to the result CSV (output)", required=True)
    subparser.add_argument("-db", "--database_name", default=Path("tree1/"), nargs='?', type=Path, help="path to the DB to use [default 'tree2/']", required=False)
    subparser.add_argument("-k", "--kmer_size", default=20, nargs='?', type=int, help="size of kmer to use [default 20]", required=False)
    subparser.add_argument("--threads", default=4, nargs='?', type=int, help="number of threads to use [Default 4]", required=False)

def parseArgs(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(help='Choose type of benchmarking', dest='sub_parser')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")

    # if parameterization
    parameterization_parser = subparsers.add_parser(SubparserNames.parameterization.value)
    parameterization_parser.add_argument("-t", "--test_directory", type=Path, help="path to the directory with simulated test reads.", required=True)
    add_common_arguments(parameterization_parser)

    # if genomecount
    genomecount_parser = subparsers.add_parser(SubparserNames.genomecount.value)
    add_common_arguments(genomecount_parser)

    # if relative_performance
    relative_performance_parser = subparsers.add_parser(SubparserNames.relative_performance.value)
    relative_performance_parser.add_argument("-t", "--test_directory", type=Path, help="path to the directory with simulated test reads.", required=True)
    relative_performance_parser.add_argument("-c", "--config", type=Path, help="path to the configuration file", required=True)
    add_common_arguments(relative_performance_parser)

    return parser.parse_args(argv)


def main():

    def parameterization_action():
        print(f"Performing parameterization benchmarking...")
        phagefilter = PhageFilter(kmer_size=args.kmer_size, filter_thresh=1.0)
        bench_test.benchtest_parameter_sweep(
            phagefilter, args.test_directory, args.database_name, args.genome_dir, args.result_csv)

    def genomecount_action():
        print(f"Performing genome count benchmarking...")
        phagefilter = PhageFilter(kmer_size=args.kmer_size, filter_thresh=1.0)
        bench_test.benchtest_genomecount(phagefilter, args.genome_dir,
                                                args.database_name, args.result_csv)

    def relative_performance_action():
        print(f"Performing relative performance benchmarking...")
        bench_test.benchtest_relative_performance(
            args.genome_dir, args.config, args.result_csv, args.test_directory)

    # arguments
    args = parseArgs(sys.argv[1:])

    # map subparsers to their corresponding functions and messages
    actions = {
        SubparserNames.parameterization.value: parameterization_action,
        SubparserNames.genomecount.value: genomecount_action,
        SubparserNames.relative_performance.value: relative_performance_action
    }

    # run benchmark type specified
    action = actions.get(args.sub_parser)
    if action:
        action()
    else:
        print(__doc__)

if __name__ == '__main__':
    main()
