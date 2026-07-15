#!/usr/bin/env python3

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import yaml
from datetime import datetime
import argparse

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Validate and extract components from YAML configuration file""")

parser.add_argument("config",    help="YAML configuration file", type=str)
parser.add_argument("stats",     help="Statistics file", type=str)
parser.add_argument("--section", help="Section in YAML configuration file containing validation entries", type=str, default='validation.double_precision')

args = parser.parse_args()

f = open(args.config,'r')
yaml_document = f.read()
config = yaml.safe_load(yaml_document)
f.close()

class Stats:
    def __init__(self):
        self.entries = []
    
    @classmethod
    def load(self,doc):
        self = Stats()
        def parse_line(self,line):
            [time, index, name, avg_dec, avg_hex, min_dec, min_hex, max_dec, max_hex, non_missing] = line.split()
            entry = { 
                        'time': datetime.strptime(time, "%Y%m%d%H%M%S"),
                        'index': int(index),
                        'name': name,
                        'average' : float(avg_dec),
                        'average_hex' : avg_hex,
                        'minimum': float(min_dec),
                        'minimum_hex' : min_hex,
                        'maximum': float(max_dec),
                        'maximum_hex': max_hex,
                        'non_missing_points': int(non_missing)
                    }
            self.entries.append(entry)
        for line in doc.splitlines():
            line = line.strip()
            if line:
                if line[0] != '#':
                    parse_line(self,line)
        return self

    def filter(self,**kwargs):
        filtered = []
        for entry in self.entries:
            keep = True
            for (key,value) in kwargs.items():
                if key == 'time' and not isinstance(value,datetime):
                    value = datetime.strptime(str(value), "%Y%m%d%H%M%S")
                if key in entry:
                    if entry[key] != value:
                        keep = False
            if keep:
                filtered.append(entry)
        stats = Stats()
        stats.entries = filtered
        return stats

    def __getitem__(self,key):
        if len(self.entries) == 1:
            return self.entries[0][key]
        else:
            return [entry[key] for entry in self.entries]

    def __repr__(self):
        if not self.entries:
            return 'STATS_EMPTY'
        lines= '# ' + ' '.join(self.entries[0].keys())
        for entry in self.entries:
            lines +=  '\n' + ' '.join([str(value) for [key,value] in entry.items()])
        return lines
   

f = open(args.stats,'r')
stats = Stats.load(f.read())
f.close()

def extract_validation_list( config, keys ):
    try:
        if len(keys) == 1:
            return config[keys[0]]
        else:
            return extract_validation_list( config[keys[0]], keys[1:])
    except KeyError:
        print("WARNING: No validation section found with key \""+args.section+"\"")
        return []
    except TypeError:
        print("WARNING: No validation section found with key \""+args.section+"\"")
        return []


validation_passed = 0
validation_failed = 0
validation_bit_identical = 0
validation_happened = 0
validation_missing = 0

validation_list = extract_validation_list(config, args.section.split('.'))


for validate in validation_list:
    name = validate['name']
    time = validate['time']
    if 'average' in validate:
        norm_type = 'average'
    elif 'minimum' in validate:
        norm_type = 'minimum'
    elif 'maximum' in validate:
        norm_type = 'maximum'
    else:
        norm_type = None
    x_ref = validate[norm_type]
    hashes_ref = validate['hashes'] if 'hashes' in validate else []
    rtol = validate['relative_tolerance'] if 'relative_tolerance' in validate else 0

    validation_stats = stats.filter(name=name,time=time)
    if not validation_stats.entries:
        validation_missing += 1
        continue

    validation_happened += 1

    x = validation_stats[norm_type]
    hash  = validation_stats[norm_type+'_hex']

    diff = abs(x-x_ref)
    rdiff = diff/abs(x_ref)

    def format_float_scientific(val,precision):
        format_str = "<21."+str(precision-1)+"e"
        return format(val, format_str)

    print("Validating {norm: <8} {name: <6} at time {time} with reference {ref: <21} : {value: <21} (hash {hash})".format(
            norm=norm_type,
            name="'"+name+"'",
            time=time,
            ref=format_float_scientific(x_ref,precision=16),
            value=format_float_scientific(x,precision=16),
            hash=hash))

    if rdiff <= rtol:
        print("    SUCCESS: Relative difference (",rdiff,") below relative tolerance (",rtol,")")
        validation_passed += 1
        if hash in hashes_ref:
            print("             Exact match with hash ",hash)
            validation_bit_identical += 1
        print("\n")

    else:
        print("    FAILED:  Relative difference (",rdiff,") above relative tolerance (",rtol,")\n")
        validation_failed += 1

if validation_failed == 0:
    if validation_passed > 0:
        print("Validation PASSED:   ",str(validation_bit_identical)+"/"+str(validation_passed),"bit identical with recorded results")
        exit(0)
    if validation_passed == 0:
        if validation_missing > 0:
            print("WARNING: No validation entries match with recorded results")
            exit(0)
else:
    print("Validation FAILED:   ",validation_passed,"passes,", validation_failed, "failures")
    exit(1)
