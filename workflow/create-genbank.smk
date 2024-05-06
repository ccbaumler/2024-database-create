###
# This workflow will create a genbank database for sourmash from the HPC cluster at UCDavis.
#
# To run:
# snakemake -s create-genbank.smk -j 6 --use-conda --rerun-incomplete --resources allowed_jobs=100
#
###
##

import time

#DATE = time.strftime("%Y%m%d")
DATE = "20240504"

DOMAINS=['plant',
         'vertebrate_mammalian',]

KSIZES=[21,31,51]

LOGS='logs'

wildcard_constraints:
    k = "\d+",
    D = "\w+",

rule all:
    input:
        expand("../dbs/genbank-{d}-{D}-k{k}.zip", d=DATE, D=DOMAINS, k=KSIZES),
        expand("data/lineages.{D}.csv", D=DOMAINS),

rule download_assembly_summary:
    output:
        good = 'data/assembly_summary.{D}.txt',
    shell: """
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.D}/assembly_summary.txt > {output.good}
    """

rule write_links_file:
    input:
        good = "data/assembly_summary.{D}.txt",
    output:
        links = "data/create.{d}-{D}.csv",
    run:
        url_pattern = r'https?://\S+'

        with open(str(output.links), 'wt') as fp:
            header = ["accession","name","url","organism_name","infraspecific_name","asm_name"]
            fp.write(','.join(header) + '\n')

            with open(str(input.good)) as tsvfile:
                lines = [line for line in tsvfile.readlines() if not line.startswith('#')]
                total = len(lines)
        
                for n, row in enumerate(lines):

                    row = row.strip().split('\t')
        
                    if n % 10 == 0:
                        print(f'...Writing {{output.links}}: Line {n} of {total}', end='\r', flush=True)
        
                    url = f'"{row[19]}"' if ',' in row[19] else row[19]
                    accession = f'"{row[0]}"' if ',' in row[0] else row[0]
                    organism_name = f'"{row[7]}"' if ',' in row[7] else row[7]
                    infraspecific_name = f'"{row[8]}"' if ',' in row[8] else row[8]
                    asm_name = f'"{row[15]}"' if ',' in row[15] else row[15]
        
                    elements = [accession, organism_name, infraspecific_name, asm_name]
                    elements = [e.strip('"') for e in elements if e != 'na']
                    name = ' '.join(elements)

                    if ',' in name:
                        name = f'"{name}"'
        
                    line = f"{accession},{name},{url},{organism_name},{infraspecific_name},{asm_name}\n"
                    fp.write(line)
        
                print(f'...Wrote {output.links}: Line {n+1} of {total}  ')


rule gather_sketch_db:
    input:
        links = "data/create.{d}-{D}.csv",
    output:
        failed = "data/create.{d}-{D}.failures.csv",
        db = temporary("../dbs/genbank-{d}-{D}.all.zip"),
    conda: "envs/directsketch.yaml",
    resources:
        allowed_jobs = 100
#        mem_mb = 8000,
#        time_min = 30,
    threads: 1
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        #k_list = lambda wildcards: f"k={','.join([f'{k}' for k in KSIZES])}",
    shell:'''
        sourmash scripts gbsketch {input.links} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled=1000,abund" -r 1 #2> {log}
    '''

rule extract_db:
    input:
        db = "../dbs/genbank-{d}-{D}.all.zip",
    output:
        db = "../dbs/genbank-{d}-{D}-k{k}.zip",
    conda: "envs/sourmash.yaml",
    shell:"""
        sourmash signature extract {input.db} -k {wildcards.k} --dna -o {output.db}
    """

### create a report with sourmash sig summarize for the databases... and sourmash compare(?)

# taxonomy rules, from https://github.com/ctb/2022-assembly-summary-to-lineages
rule download_ncbi_utils:
    output: "scripts/ncbi_taxdump_utils.py"
    shell:
        "curl -L https://raw.githubusercontent.com/ctb/2022-assembly-summary-to-lineages/main/ncbi_taxdump_utils.py > {output}"

rule download_taxscript:
    output: "scripts/make-lineage-csv.py"
    shell:
        "curl -L https://raw.githubusercontent.com/bluegenes/2022-assembly-summary-to-lineages/virus-tax/make-lineage-csv.py > {output}"

rule download_taxdump: # may need to restart this a couple times
    output:
        "taxdump/nodes.dmp",
        "taxdump/names.dmp"
    shell:
        "curl -L ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | (mkdir -p taxdump && cd taxdump && tar xzvf -)"

rule make_lineage_csv:
    input:
        "data/assembly_summary.{D}.txt",
        "taxdump/nodes.dmp",
        "taxdump/names.dmp",
        "scripts/make-lineage-csv.py",
        "scripts/ncbi_taxdump_utils.py",
    output:
        "data/lineages.{D}.csv"
    params:
        ictv_cmd = lambda w: " --ictv " if 'viral' in w.D else '',
    shell:
        "python scripts/make-lineage-csv.py taxdump/{{nodes.dmp,names.dmp}} {input[0]} -o {output} {params.ictv_cmd}"
