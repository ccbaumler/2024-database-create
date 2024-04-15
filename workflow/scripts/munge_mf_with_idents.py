#! /usr/bin/env python

import sys
import argparse
import csv
import os.path
import time
import sourmash
from sourmash import manifest # also does 'import sourmash'


def main():
    p = argparse.ArgumentParser()

    # Add command-line arguments with default values from Snakemake
    p.add_argument('good_idents', nargs='?', help='list of valid identifiers')
    p.add_argument('bad_idents', nargs='?', help='list of suppressed and re-versioned identifiers')
    p.add_argument('old_mf', nargs='?', help='existing sourmash database manifest')
    p.add_argument('--report', nargs='?', help='details of removed etc., for humans')
    p.add_argument('-o', '--output', nargs='?', help='manifest cleansed of the impure')
    p.add_argument('-t', help='the Genbank assembly summary text file')
    p.add_argument('-s', help='the Genbank assembly summary history text file')

    args = p.parse_args()

    good_idents = set()
    good_idents_no_version = set()

    bad_idents = set()
    #bad_idents_no_version = set()

    bad_idents_dict = {}

    def get_suffix(name):
        ident = name.split(' ')[0]
        assert ident.startswith('GC')  
        suffix = ident[3:]
        return suffix

    def get_suffix_no_version(name):
        ident = name.split(' ')[0]
        suffix = get_suffix(ident)
        assert '.' in suffix
        return suffix.split('.')[0]

    with open(args.good_idents, 'rt') as fp:
        for line in fp:
            assert line.startswith('GC')
            line = line.strip()
            suffix = get_suffix(line)
            good_idents.add(suffix)

            suffix_no_version = get_suffix_no_version(line)
            good_idents_no_version.add(suffix_no_version)

    with open(args.bad_idents, 'rt') as fp:
        for line in fp:
            line = line.strip()
            assert line.startswith('GC')
            suffix = get_suffix(line)
            bad_idents.add(suffix)

            # Extract the version number from the suffix
            assembly, version = suffix.split('.')

            # Add the suffix and version as key-value pair to the dictionary
            if assembly not in bad_idents_dict:
                bad_idents_dict[assembly] = []
            
            bad_idents_dict[assembly].append(version)

    assert len(good_idents) == len(good_idents_no_version)

    assert len(bad_idents) == sum([len(k) for k in bad_idents_dict.values()])
    
    #print(len(bad_idents_dict.keys()))
    #print(len(bad_idents))
    
    #print(sum([len(x) for x in bad_idents_dict.values()]))
        
    old_version_dict = {k: v for k, v in bad_idents_dict.items() if len(v) > 1}
    #print(old_version_dict)

#
#    with open("bad_idents.txt", "w") as of:
#            of.write("\n".join(bad_idents))
#    with open("bad_idents_no_version.txt", "w") as of:
#                    of.write("\n".join(bad_idents_no_version))
#    with open("bad_idents_dict_keys.txt", "w") as of:
#        of.write("\n".join(bad_idents_dict.keys()))
#    with open("bad_idents_dict.csv", "w") as of:
#        for k in sorted( bad_idents_dict ):
#            for v in bad_idents_dict[k]:
#                of.write(k + '.' + v + '\n' )
#
    old_mf = manifest.BaseCollectionManifest.load_from_filename(args.old_mf)

    ## now go through and filter

    bad_list = []
    for k in bad_idents_dict:
        for v in bad_idents_dict[k]:
            bad_list.append( k + '.' + v )
    #print(bad_list)

    versioned_list = []
    for k in old_version_dict:
        for v in old_version_dict[k]:
            versioned_list.append( k + '.' + v )
    #print(versioned_list)

    removed_list = []
    updated_version_list = []
    keep_rows = []
    n_changed_version = 0
    
    bad = []
    suppressed_versioned = []
    

    for row in old_mf.rows:
        name = row['name']
        suffix = get_suffix(name)
        suffix_no_version = get_suffix_no_version(name)

        if suffix in good_idents:
            keep_rows.append(row)
        else:
            if suffix_no_version in good_idents_no_version:
                n_changed_version += 1
                updated_version_list.append(name)
            else:
                removed_list.append(name)
        if suffix in bad_list:
            bad.append(suffix)
    #print(bad)

    with open(args.s, "r", newline='') as fp:
        # skip header
        for x in range(2):
            next(fp)
        # create list from historical db    
        doc = [line.strip().split('\t') for line in fp]

    # Create a list of identifier and status
    for row in doc:
        name = row[0]
        suffix = get_suffix(name)

        if suffix in bad:
            suppressed_versioned.append(row[0::10])
        

    print(len(suppressed_versioned))
    
    n_removed = len(old_mf) - len(keep_rows)
    n_suspect_suspension = n_removed - n_changed_version
    new_mf = manifest.CollectionManifest(keep_rows)
    
    creat_time = time.ctime(os.path.getctime(args.t))
    mod_time = time.ctime(os.path.getmtime(args.t)) 


    print(f"\n\nFrom genome assemblies database:", file=sys.stderr)
    print(f"Loaded {len(good_idents)} identifiers from '{args.good_idents}'",
                          file=sys.stderr)
    print(f"(and loaded {len(good_idents_no_version)} identifiers without version number)", file=sys.stderr)
    print(f"File assembly database created on {creat_time}", file=sys.stderr)
    print(f"File assembly database last modified {mod_time}", file=sys.stderr)

    print(f"\nFrom '{args.old_mf}':", file=sys.stderr)
    print(f"Kept {len(keep_rows)} of {len(old_mf)} identifiers.",
            file=sys.stderr)
    
    print(f"\nFrom '{args.bad_idents}':", file=sys.stderr)
    print(f"Kept {len(bad_list)} of {len(bad_idents)} identifiers.", file=sys.stderr)
    
    print(f"\nNew manifest '{args.output}':", file=sys.stderr)
    print(f"Kept {len(new_mf)} identifiers.", file=sys.stderr)
    print(f"Removed {n_removed} total identifiers.",
            file=sys.stderr)
    print(f"Removed {n_changed_version} identifiers because of version change.",
            file=sys.stderr)
    print(f"Removed {n_suspect_suspension} identifiers because of suspected suspension of the genome.\n\n",
            file=sys.stderr)


    with open(args.output, 'w', newline='') as outfp:
        new_mf.write_to_csv(outfp, write_header=True)

    if args.report:
        with open(args.report, 'wt') as fp:
            print(f"From {len(old_mf)} in '{args.old_mf}':", file=fp)
            print(f"Kept {len(new_mf)} in '{args.output}.", file=fp)
            print(f"Removed {n_removed} total.", file=fp)
            print(f"Removed {n_suspect_suspension} identifiers because of suspected suspension of the genome.",
                  file=fp)
            print(f"Removed {n_changed_version} because of changed version.",
                  file=fp)
            print(f"---- {len(suppressed_versioned)} included into the bad list category.", file=fp)
            for item in suppressed_versioned:
                print("\n".join(str(i) for i in item), file=fp)
            print(f"---- {n_suspect_suspension} removed because presumed guilt ----", file=fp)
            print("\n".join(removed_list), file=fp)
            print(f"---- {n_removed} removed because version changed ----", file=fp)
            print("\n".join(updated_version_list), file=fp)


if __name__ == '__main__':
    sys.exit(main())
