#! /usr/bin/env python

# modified from https://github.com/michalbukowski/fetch-genomes/blob/main/fetch_genomes.py

import os
import argparse
import sys
import urllib.request
import hashlib

# Data that can be obtained from NCBI GenBank for a given genomic assembly,
# see parse_args function for more information.
assembly_formats = {
    'fna'  : 'genomic.fna.gz',
    'gbff' : 'genomic.gbff.gz',
    'gff'  : 'genomic.gff.gz',
    'rna'  : 'rna_from_genomic.fna.gz',
    'cds'  : 'cds_from_genomic.fna.gz',
    'prot' : 'translated_cds.faa.gz'
}

md5sums_fname = 'md5checksums.txt'


def fetch_genomes(urls, formats, output_dir):
    '''Fetches genomes by sending FTP requests via urllib. Arguments:
       urls       -- a list of URLs pointing to the genomes
       formats    -- formats of data to be retrieved
       output_dir -- a directory for the data to be saved to
    '''
    not_found = 0
    fetched   = 0
    existing  = 0 

    for url in urls:
        if url.startswith('https://'):
            url = 'ftp://' + url[8:]
        pos = url.rfind('/')
        asm_full_name = url[pos+1:]
        
        done = [False] * len(formats)
        for i, fmt in enumerate(formats):
            suffix  = assembly_formats[fmt]
            fnamein = f'{asm_full_name}_{suffix}'
            fpathout = f'{output_dir}/{fnamein}'
            if os.path.exists(fpathout):
                if os.path.isfile(fpathout):
                    done[i] = True
        if all(done):
            existing += len(formats)
            print(f'[INFO] All files requested for {asm_full_name} exist and are files, considered done')
            print(f'[INFO] Skipping {asm_full_name}, already fetched')
            continue
        print(f'\n[INFO] Fetching files for {asm_full_name}...')
        
        try:
            res = urllib.request.urlopen(url, timeout=60)
            lines = res.read().decode().rstrip().split('\n')
            flist = [ line.split()[-1] for line in lines ]
        except KeyboardInterrupt as e:
            raise e
        except:
            print(f'[ERROR] Cannot fetch file list from "{url}"')
            print(f'[WARNING] Skipping assembly {asm_full_name}...')
            continue
        print(f'[INFO] There is {len(flist)} files at "{url}"')

        # Fetch the file with MD5 sums for genome files, if unsuccessful, yield
        # a proper message and continue to next iteration/genome.
        full_path = f'{url}/{md5sums_fname}'
        try:
            res = urllib.request.urlopen(full_path, timeout=60)
            md5sums = res.read().decode().rstrip().split('\n')
        except KeyboardInterrupt as e:
            raise e
        except:
            print(f'[ERROR] Info on MD5 checksums cannot be fetched from "{full_path}"')
            print(f'[WARNING] Skipping assembly {asm_full_name}...')
            continue
        md5sums = [ line.split() for line in md5sums ]
        md5sums = { line[1].lstrip('./') : line[0] for line in md5sums }
        print(f'[INFO] MD5 checksums for {asm_full_name} successfully fetched')

        old_fetched = fetched
        for fmt in formats:
            suffix  = assembly_formats[fmt]
            fnamein = f'{asm_full_name}_{suffix}'
            fpathin = f'{url}/{fnamein}'
            
            if not fnamein in flist:
                print(f'[ERROR] No such file for {asm_full_name}: "{fpathin}"')
                print(f'[WARNING] Skipping {asm_full_name} assembly file: "{fpathin}"')
                not_found += 1
                continue
            
            fpathout = f'{output_dir}/{fnamein}'
            if os.path.exists(fpathout):
                if os.path.isfile(fpathout):
                    print(f'[INFO] The output path "{fpathout}" exists and is a file, considered done')
                    print(f'[INFO] Skipping {asm_full_name} assembly file: "{fpathin}", already fetched')
                    existing += 1
                else:
                    print(f'[ERROR] The output path "{fpathout}" exists and is not a file')
                    print(f'[WARNING] Skipping {asm_full_name} assembly file: "{fpathin}"')
                continue
            
            if not fnamein in md5sums:
                print(f'[ERROR] Cannot find MD5 checksum for {asm_full_name} assembly file: "{fpathin}"')
                print(f'[WARNING] Skipping {asm_full_name} assembly file: "{fpathin}"')
                continue
            
            try:
                res = urllib.request.urlopen(fpathin, timeout=60)
                content = res.read()
            except:
                print(f'[ERROR] {asm_full_name} assembly file cannot be fetched from: "{fpathin}"')
                print(f'[WARNING] Skipping {asm_full_name} assembly file: "{fpathin}"')
                continue
            print(f'[INFO] {asm_full_name} assembly file "{fpathin}" successfully fetched')
            
            md5sum = hashlib.md5(content).hexdigest()
            if md5sum == md5sums[fnamein]:
                print(f'[INFO] Correct MD5 checksum ({md5sum}) for {asm_full_name} assembly file: "{fpathin}"')
            else:
                print(f'[ERROR] Incorrect MD5 checksum ({md5sum}) ' + \
                      f'for {asm_full_name} assembly file ({md5sums[fnamein]}): "{fpathin}"')
                print(f'[WARNING] Skipping {asm_full_name} assembly file: "{fpathin}"')
                continue
            
            tmpfpathout = f'{output_dir}/.{fnamein}'
            try:
                with open(tmpfpathout, 'wb') as f:
                    f.write(content)
            except:
                print(f'[ERROR] Cannot save to "{fpathout}" the {asm_full_name} assembly file: "{fpathin}"')
                print(f'[WARNING] Skipping {asm_full_name} assembly file: "{fpathin}"')
                continue
            
            try:
                os.rename(tmpfpathout, fpathout)
            except:
                print(f'[ERROR] Cannot save to "{fpathout}" the {asm_full_name} assembly file: "{fpathin}"')
                print(f'[WARNING] Skipping {asm_full_name} assembly file: "{fpathin}"')
                os.remove(tmpfpathout)
            else:
                print(f'[INFO] {asm_full_name} assembly file "{fpathin}" successfully saved to "{fpathout}"')
                fetched += 1
    
    total = len(urls) * len(formats)
    left  = total-existing-not_found-fetched
    print(f'\n[INFO] Fetched {fetched} files out of {total} inferred ' + \
          f'(already existing: {existing}, not found on site: {not_found})')
    if left > 0:
        print(f'\n[WARNING] {left} files are still to be fetched')
    else:
        print('[INFO] All files have been successfully fetched')
    print('[INFO] Fetching genomes has been completed')

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='The input file containing only the https links to NCBI')
    p.add_argument('-f', '--format', type=str, nargs='+', default=['fna'], choices=assembly_formats.keys(), help='The file format you would like to acquire')
    p.add_argument('-o', '--output-dir', help='The location to send the downloaded files')

    args = p.parse_args()
    with open(args.input, "r") as fp:
        urls = [line.strip() for line in fp.readlines()]
    print(urls)
    fetch_genomes(urls = urls, formats = args.format, output_dir = args.output_dir)
    print('end')
if __name__ == '__main__':
    sys.exit(main())
# Example usage:
#urls = [
#    "https://example.com/genomes/genome1.zip",
#    "https://example.com/genomes/genome2.zip",
#    "https://example.com/genomes/genome3.zip"
#]
#formats = ["fasta", "gff"]
#output_dir = "output"
#fetch_genomes(urls, formats, output_dir)
#
