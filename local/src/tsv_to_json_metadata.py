import json
import os
from argparse import ArgumentParser

# TODO: date format, get a metadata.json file to complete (honestly having the schema here and validating etc seems an overkill)

def main():
    args = get_args()
    with open(args.input, 'r') as tsv:
        header = tsv.readline().rstrip()
        fields = header.split("\t")
        with open(args.output, 'w') as json_meta:
            for line in tsv.readlines():
                line = line.rstrip()
                values = line.split("\t")
                entry = { "id": values[2], "LAS_Validation": values[3] == "TRUE", "files": get_files(values[0], args.fastq) }
                #if len(values) >= 4:
                for i in range(4, len(values)):
                    entry[fields[i]] = values[i]
                json_meta.write(json.dumps(entry) + "\n")        

def get_files(vendor_id, fq_dir):
    return [os.path.join(fq_dir, f) for f in [vendor_id+"S1_L004_R1_001.fastq.gz", vendor_id+"S1_L004_R2_001.fastq.gz"]]

def get_args():
    parser = ArgumentParser(description="From a sample sheet in excel format to the json samples_map for las mdam")
    parser.add_argument("-o", "--output",
                        default="metadata.json",
                        required=True,
                        help="path and name of output file. Directories must exist.")
    parser.add_argument("-f", "--fastq",
                        default="./fastq/",
                        help="the fastq dir")
    parser.add_argument("-i", "--input",
                        required=True,
                        help="the tsv sample sheet")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
