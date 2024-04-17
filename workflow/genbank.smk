###
# This workflow will create a genbank database for sourmash from the HPC cluster at UCDavis.
# To run:
# snakemake -s genbank -j 1 --use-conda
#
# Note: the first run may take a few hours to build the wort sqlmf
# try linking an old sqlmf to skip some time and effort
# ln -s /group/ctbrowngrp/sourmash-db/wort-manifests/2023-05-31.wort.sqlmf $(pwd)/data/
###

import time

DATE = time.strftime("%Y%m%d")

DOMAINS=['archaea',
         'fungi',
         'protozoa',
         'bacteria',
         'viral']

TAG='latest'

KSIZES=[21,31,51]
LOGS='logs'
DB_DATES = ['2022.03',]
DATABASES = [f'genbank-{date}-{domain}-k{ksize}.zip'
             for date in DB_DATES
             for domain in DOMAINS
             for ksize in KSIZES]

rule all:
    input:
        #expand("data/{D}.missing.csv", D=DOMAINS),
#        expand("genbank-{dates}-{D}-k{k}.zip", dates=DB_DATES, D=DOMAINS, k=KSIZES),
        expand("../dbs/genbank-{d}-{D}-k{k}.zip", d=DATE, D=DOMAINS, k=KSIZES),
        expand("data/{D}.{k}.mf.clean.csv", D=DOMAINS, k=KSIZES),
        expand("data/{D}.lineages.csv", D=DOMAINS),

rule build_genbank:
    input:
        expand("data/genbank-{d}-{D}-k{k}.zip", d=DATE, D=DOMAINS, k=KSIZES),

rule clean_genbank:
    input:
        expand("../dbs/genbank-{d}-{D}-k{k}.clean.zip", d=DATE, D=DOMAINS, k=KSIZES),

rule missing_genbank:
    input:
        expand("../dbs/genbank-{d}-{D}-k{k}.missing.zip", d=DATE, D=DOMAINS, k=KSIZES),

rule check_genbank:
    input:
        expand("data/genbank-{d}-{D}-k{k}.zip.check", d=DATE, D=DOMAINS, k=KSIZES),

rule tax_genbank:
    input:
        expand("data/{D}.lineages.csv", D=DOMAINS),

rule download_assembly_summary:
    output:
        good = 'data/{D}.assembly_summary.txt',
        bad = 'data/{D}.assembly_summary_historical.txt',
    shell: """
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.D}/assembly_summary.txt > {output.good}
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.D}/assembly_summary_historical.txt > {output.bad}
    """

rule get_ss_db:
    output:
        DATABASES
    conda: "envs/sourmash.yaml"
    params:
        db_date = DB_DATES
    shell: """
        for db_file in {output}; do
            echo "Checking if $db_file exists..."
            if [ -e /group/ctbrowngrp/sourmash-db/genbank-{params.db_date}/$db_file ]; then
    
                echo "$db_file exists!"
                echo "Linking existing file to $(pwd)"
    
                ln -s /group/ctbrowngrp/sourmash-db/genbank-{params.db_date}/$db_file $(pwd)/$db_file
            else
    
                echo "$db_file does not exist!"
                echo "Downloading file to $(pwd)"
    
                wget -O $db_file https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-{params.db_date}/$db_file
            fi
        done
    """

# Does this need to be updated because some scripts are not on main branches?
rule download_update_sourmash_dbs:
    output: "scripts/update_sourmash_dbs.py"
    shell: """
        curl -L "https://raw.githubusercontent.com/ccbaumler/2024-database-create/manifests/workflow/scripts/update_sourmash_dbs.py" > {output}
    """

rule manifest_manifest:
    input: DATABASES,
    output: "data/{D}.{k}.mf.csv",
    conda: "envs/sourmash.yaml",
    shell: """
        #processed_d=()
        for db_file in {input}; do
            if [[ $db_file == *{wildcards.D}*{wildcards.k}* ]]; then  #&& ! " ${{processed_d[@]}} " =~ [[:space:]]{wildcards.D}[[:space:]] ]]; then
                sourmash sig manifest -o {output} $db_file --no-rebuild
         #       processed_d+=({wildcards.D})
            fi
        done
    """

rule cleanse_manifest:
    input:
        script = "scripts/update_sourmash_dbs.py",
        good = "data/{D}.assembly_summary.txt",
        bad = "data/{D}.assembly_summary_historical.txt",
        manifest = "data/{D}.{k}.mf.csv",
    output:
        clean = "data/{D}.{k}.mf.clean.csv",
        reversion = "data/{D}.{k}.updated-versions.txt",
        report = "data/{D}.{k}.report.txt",
        missing = "data/{D}.{k}.missing-genomes.txt",
    conda: "envs/sourmash.yaml",
    shell: """
        {input.script} {input.manifest} -a {input.good} -b {input.bad} -o {output.clean} --updated-version {output.reversion} --report {output.report} --missing-genomes {output.missing}
    """

rule picklist_clean:
    input:
        clean = "data/{D}.{k}.mf.clean.csv",
        dbs = ["genbank-{0}-{{D}}-k{{k}}.zip".format(dates) for dates in DB_DATES],
    output:
        woohoo = "../dbs/genbank-{d}-{D}-k{k}.clean.zip",
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = 8000,
        time_min = 30
    #params:
    #    og_db = lambda dates, D, k: f"genbank-{date}-{D}-k{k}.zip" for date in DB_DATES,
    shell:'''
        echo "Cleaning {input.dbs}..."
        sourmash sig extract --picklist {input.clean}::manifest {input.dbs} -o {output.woohoo}
        echo "{input.dbs} cleaned and stored as {output.woohoo}"
    '''

rule gather_sketch_missing:
    input:
        reversion = "data/{D}.{k}.updated-versions.txt",
    output:
        woohoo = "../dbs/genbank-{d}-{D}-k{k}.clean.zip",
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = 8000,
        time_min = 30
    #params:
    #    og_db = lambda dates, D, k: f"genbank-{date}-{D}-k{k}.zip" for date in DB_DATES,
    shell:'''
        echo "Cleaning {input.dbs}..."
        sourmash sig extract --picklist {input.clean}::manifest {input.dbs} -o {output.woohoo}
        echo "{input.dbs} cleaned and stored as {output.woohoo}"
    '''



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
