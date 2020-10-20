#!/usr/bin/env python3
import json
import os
import datetime
from argparse import ArgumentParser

# maybe-FIXME: it could be wiser to load an entire excel with all info on two sheets, 
# load everything then validate vs the schema rather than this "patching" of the
# manually generated first portion of the json, but honestly not urgent.

# TODO print stdout SAMPLES and SAMPLES_ORIG for conf.sk
def main():
    args = get_args()
    with open(args.input, 'r') as tsv:
        header = tsv.readline().rstrip()
        fields = header.split("\t")
        with open(args.json, 'a') as json_meta:
            entries = []
            samples = []
            orig_samples = []
            for line in tsv.readlines():
                line = line.rstrip()
                values = line.split("\t")
                id = values[2]#.replace('-','Y')
                entry = { "id": id, "LAS_Validation": values[3] == "TRUE", "files": get_files(values[0], args.fastq) }
                samples.append(values[2])
                orig_samples.append(values[0]+'_SA_L001')
                #if len(values) >= 4:
                for i in range(4, len(values)):
                    if fields[i] == "cloning_date":
                        if values[i] != "NA":
                            #entry[fields[i]] = str(datetime.datetime.strptime(values[i], '%Y/%m/%d').date())
                            try:
                                entry[fields[i]] = str(datetime.datetime.strptime(values[i], '%d-%m-%y').date())
                            except: # hurry and hate
                                entry[fields[i]] = None
                        else:
                            entry[fields[i]] = None
                    else:
                        entry[fields[i]] = values[i]
                entries.append(entry) # python triiiiicky indenting
                #json_meta.write(json.dumps(entry) + "\n")     
            print("SAMPLES_ORIG=" + str(orig_samples))
            print("SAMPLES=" + str(samples))
            jsons = json.dumps({"samples_map": entries}, indent=4)
            json_meta.write(jsons.strip('{}'))
            json_meta.write("}\n")

def get_files(vendor_id, fq_dir):
    return [os.path.join(fq_dir, f) for f in [vendor_id+"_SA_L001_R1_001.fastq.gz", vendor_id+"_SA_L001_R2_001.fastq.gz"]]

def get_args():
    parser = ArgumentParser(description="From a sample sheet in excel format to the json samples_map for las mdam")
    parser.add_argument("-f", "--fastq",
                        default="./fastq/",
                        help="the fastq dir")
    parser.add_argument("-i", "--input",
                        required=True,
                        help="the tsv sample sheet")
    parser.add_argument("-j", "--json",
                        required=True,
                        help="the json to be filled")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
