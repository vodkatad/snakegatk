import json
from argparse import ArgumentParser


def main():
    args = get_args()
    with tsv as open(args.input, 'r'):
        header = info.readline().rstrip()
        fields = header.split("\t")
        for line in tsv.readlines():
            l = l.rstrip()
            values = l.split("\t")
            entry = { "id": values[2], "LAS_Validation": values[3], file: get_file(values[0], args.fastq) }
            if len(values) >= 4:
                #add fields


def get_args():
    parser = ArgumentParser(description="From a sample sheet in excel format to the json samples_map for las mdam")
    parser.add_argument("-o", "--output",
                        default="metadata.json",
                        help="path and name of output file. Directories must exist.")
    parser.add_argument("-f", "--fastq",
                        default="./fastq/",
                        help="the fastq dir")
    parser.add_argument("-i", "--input",
                        help="the tsv sample sheet")
    return parser.parse_args()


if __name__ == '__main__':
    main()
