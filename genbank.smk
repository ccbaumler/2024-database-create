###
# This workflow will create a genbank database for sourmash from the HPC cluster at UCDavis.
# To run:
# snakemake -s genbank -j 1 --use-conda
#
# Note: the first run may take a few hours to build the wort sqlmf
# try linking an old sqlmf to skip some time and effort
# ln -s /group/ctbrowngrp/sourmash-db/wort-manifests/2023-05-31.wort.sqlmf $(pwd)/data/
###

wort-sigs is not an option anymore
need to rethink a more flexible approach without needing much resources
I think that I should take the existing genbank database and compare it to the 
new assembly reports (and historic ones) to create a fresh sparkly clean database
with ftpdirpaths piped directly into an updated database 

1. take clean database script
2. update to compare to existing database and
  2a. remove the bad genomes for historic
  2b. add the new latest genomes from summary
3. output a log of what has changed


import time

DATE = time.strftime("%Y%m%d")

DOMAINS=['archaea',
         'fungi',
         'protozoa',
         'bacteria',
         'viral']

TAG='latest'
DATABASES = ['/group/ctbrowngrp/irber/data/wort-data/wort-genomes/sigs']
KSIZES=[21,31,51]
LOGS='logs'

rule all:
    input:
        expand("data/{D}.missing.csv", D=DOMAINS),
        expand("data/{D}.manifest.csv", D=DOMAINS),
        expand("data/{D}.lineages.csv", D=DOMAINS),

rule build:
    input:
        expand("data/genbank-latest-{D}-k{k}.zip", D=DOMAINS, k=KSIZES)

rule check:
    input:
        expand("data/genbank-latest-{D}-k{k}.zip.check", D=DOMAINS, k=KSIZES)

rule tax:
    input:
        expand("data/{D}.lineages.csv", D=DOMAINS)

rule download_assembly_summary:
    output:
        'data/{D}.assembly_summary.txt'
    shell: """
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.D}/assembly_summary.txt > {output}
    """

rule genome_assembly_summary:
    input:
        'data/{D}.assembly_summary.txt'
    output:
        'data/{D}.assembly_summary.genome.txt'
    shell: """
        awk -F "\t" '$12=="Complete Genome" && $11=="latest" {print}' {input} > {output}
    """

rule make_idents:
    input:
        'data/{D}.assembly_summary.txt'
    output:
        "data/{D}.idents.csv"
    shell: """
        echo ident > {output}
        cut -f1 {input} | grep -v ^# >> {output}
    """

rule make_genome_idents:
    input:
        'data/{D}.assembly_summary.genome.txt'
    output:
        "data/{D}.idents.genome.csv"
    shell: """
        echo ident > {output}
        cut -f1 {input} | grep -v ^# >> {output}
    """

rule download_sigs_to_manifest:
    output: "scripts/sigs-to-manifest.py"
    shell:
        "curl -L https://raw.githubusercontent.com/sourmash-bio/database-examples/main/sigs-to-manifest.py > {output}"


rule make_database_check:
    input:
        script = "scripts/sigs-to-manifest.py",
    output:
        txt = f"data/{DATE}-wort-sigs.txt",
        manifest = f"data/{DATE}-wort-sigs.sqlmf",
    conda: "envs/sourmash.yaml"
    params:
        databases = DATABASES
    shell: """
        current_date=$(date +%Y%m%d)

        if [ -e data/*.sqlmf ]; then
            previous_manifests=$(find data/ -name "*.sqlmf" -print)
            closest_manifest=$(echo "$previous_manifests" | sort | head -n 1)

            echo "Closest previous manifest file to today ($current_date):"
            echo "$closest_manifest"

            closest_date=$(basename "$closest_manifest" | cut -d'-' -f1)

            if [ $current_date = $closest_date ]; then
                echo "Today's date matches previous manifest!"
                echo "Manifest will not be updated."
                exit 0
            fi

            echo "Creating text file of all wort signatures on Farm" 
            find {params.databases} -type f > {output.txt}

            echo "Updating the manifest file to current date with new content"
            python scripts/sigs-to-manifest.py \
                   --previous $closest_manifest --merge -F sql \
                   -o {output.manifest} {output.txt}

            chmod a-w $closest_manifest

        else
            find {params.databases} -type f > {output.txt}

            python scripts/sigs-to-manifest.py \
                 -F sql -o {output.manifest} {output.txt}
        fi
    """

rule picklist_check:
    input:
        database = f"data/{DATE}-wort-sigs.sqlmf",
        picklist = "data/{D}.idents.csv",
    output:
        missing = "data/{D}.missing.csv",
        manifest = "data/{D}.manifest.csv",
    conda: "envs/sourmash.yaml"
    log: f"{LOGS}/{{D}}.picklist-check.log"
    benchmark: f"{LOGS}/{{D}}.picklist-check.benchmark"
    params:
        databases = DATABASES
    shell: """
        sourmash sig check --picklist {input.picklist}:ident:ident \
            {input.database} --output-missing {output.missing} \
            --save-manifest {output.manifest} 2> {log}
        touch {output.missing}
    """

rule build_zip:
    input:
        databases = DATABASES,
        manifest = "data/{D}.manifest.csv",
    output:
        "data/genbank-latest-{D}-k{k}.zip"
    conda: "envs/sourmash.yaml"
    log: f"{LOGS}/{{D}}-{{k}}.build-zip.log"
    benchmark: f"{LOGS}/{{D}}-{{k}}.build-zip.benchmark"
    shell: """
        sourmash sig cat {input.manifest} -k {wildcards.k} -o {output} 2> {log}
    """

rule picklist_confirm:
    input:
        picklist = "data/{D}.idents.csv",
        zip = "data/genbank-latest-{D}-k{k}.zip",
    output:
        confirm = touch("genbank-latest-{D}-k{k}.zip.check")
    conda: "envs/sourmash.yaml"
    log: f"{LOGS}/{{D}}-k{{k}}.picklist-confirm.log"
    benchmark: f"{LOGS}/{{D}}-k{{k}}.picklist-confirm.benchmark"
    shell: """
        sourmash sig check --picklist {input.picklist}:ident:ident \
            {input.zip} --fail 2> {log}
    """

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
        "data/{D}.assembly_summary.txt",
        "taxdump/nodes.dmp",
        "taxdump/names.dmp",
        "scripts/make-lineage-csv.py",
        "scripts/ncbi_taxdump_utils.py",
    output:
        "data/{D}.lineages.csv"
    params:
        ictv_cmd = lambda w: " --ictv " if 'viral' in w.D else '',
    shell:
        "python scripts/make-lineage-csv.py taxdump/{{nodes.dmp,names.dmp}} {input[0]} -o {output} {params.ictv_cmd}"
