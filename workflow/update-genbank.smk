###
# This workflow will create a genbank database for sourmash from the HPC cluster at UCDavis.
#
# To run:
# snakemake -s update-genbank.smk -j 151515151515151515151515151515 --use-conda --rerun-incomplete --resources allowed_jobs=100
#
# Note: the first run may take a few hours to build the wort sqlmf
# try linking an old sqlmf to skip some time and effort
# ln -s /group/ctbrowngrp/sourmash-db/wort-manifests/2023-05-31.wort.sqlmf $(pwd)/data/
###

import time

#DATE = time.strftime("%Y%m%d")
DATE="20240503",

DOMAINS=['archaea',
         'fungi',
         'protozoa',
         'bacteria',
         'viral']

#TAG='latest'

KSIZES=[21,31,51]
LOGS='logs'
DB_DATES = ['2022.03',]
DATABASES = [f'genbank-{date}-{domain}-k{ksize}.zip'
             for date in DB_DATES
             for domain in DOMAINS
             for ksize in KSIZES]

wildcard_constraints:
    k = "\d+",
    D = "\w+",

rule all:
    input:
        #expand("data/assembly_summary.{D}.txt", D=DOMAINS)
        #expand("data/{D}.missing.csv", D=DOMAINS),
        expand("../dbs/genbank-{d}-{D}-k{k}.clean.zip", d=DATE, D=DOMAINS, k=KSIZES),
        expand("../dbs/genbank-{d}-{D}-k{k}.update.zip", d=DATE, D=DOMAINS, k=KSIZES),
        expand("../dbs/genbank-{d}-{D}-k{k}.zip", d=DATE, D=DOMAINS, k=KSIZES),
        #expand("data/mf.{d}-{D}-k{k}.csv", d=DATE, D=DOMAINS, k=KSIZES),
        expand("data/lineages.{D}.csv", D=DOMAINS),
        #expand("versioned_sigs/{D}/", D=DOMAINS, k=KSIZES),

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
        expand("data/lineages.{D}.csv", D=DOMAINS),

rule download_assembly_summary:
    output:
        good = 'data/assembly_summary.{D}.txt',
        bad = 'data/assembly_summary_historical.{D}.txt',
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
    output: "data/mf.{d}-{D}-k{k}.csv",
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
        good = "data/assembly_summary.{D}.txt",
        bad = "data/assembly_summary_historical.{D}.txt",
        manifest = "data/mf.{d}-{D}-k{k}.csv",
    output:
        clean = "data/mf-clean.{d}-{D}-k{k}.csv",
        reversion = "data/updated-versions.{d}-{D}-k{k}.csv",
        report = "data/report.{d}-{D}-k{k}.txt",
        missing = "data/missing-genomes.{d}-{D}-k{k}.csv",
    conda: "envs/sourmash.yaml",
    shell: """
        {input.script} {input.manifest} -a {input.good} -b {input.bad} -o {output.clean} --updated-version {output.reversion} --report {output.report} --missing-genomes {output.missing}
    """

# checkpoint to let snakemake know not to continue until checkpoint is complete?
checkpoint picklist_clean_db:
    input:
        clean = "data/mf-clean.{d}-{D}-k{k}.csv",
        dbs = ["genbank-{0}-{{D}}-k{{k}}.zip".format(dates) for dates in DB_DATES],
    output:
        woohoo = protected("../dbs/genbank-{d}-{D}-k{k}.clean.zip"),
    conda: "envs/sourmash.yaml"
#    resources:
#        mem_mb = 8000,
#        time_min = 30
    shell:'''
        # Use `chmod +w ../dbs/*` on rerun
        if [ ! -e {output.woohoo} ]; then
            echo "Cleaning {input.dbs}..."
            sourmash sig extract --picklist {input.clean}::manifest {input.dbs} -o {output.woohoo}
            echo "{input.dbs} cleaned and stored as {output.woohoo}"
        fi
    '''

rule check_txt_reversioned:
    input:
        reversion = expand("data/updated-versions.{d}-{D}-k{k}.csv", d=DATE, D=DOMAINS, k=KSIZES),
        script = "scripts/check_txt_files.py",
    output:
        solo = "data/update.{d}-{D}.csv",
    shell: """
        files=""
        for file in "data/updated-versions.*{wildcards.d}-{wildcards.D}-k*.csv"; do
            files+=" $file"
        done
        {input.script} $files -o {output.solo}
    """

rule gather_sketch_reversioned:
    input:
        reversion = "data/update.{d}-{D}.csv",
    output:
        failed = "data/update.{d}-{D}.failures.csv",
        db = "../dbs/genbank-{d}-{D}.rever.zip"
    conda: "envs/directsketch.yaml"
    resources:
        allowed_jobs = 100
#        mem_mb = 8000,
#        time_min = 30,
    threads: 1
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        #k_list = lambda wildcards: f"k={','.join([f'{k}' for k in KSIZES])}",
    shell:'''
        sourmash scripts gbsketch {input.reversion} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled=1000,abund" -r 1 #2> {log}
    '''

# checkpoint to let snakemake know not to continue until checkpoint is complete?
checkpoint cat_to_clean_reversioned:
    input:
        dir = "../dbs/genbank-{d}-{D}.rever.zip",
        db = "../dbs/genbank-{d}-{D}-k{k}.clean.zip",
        missing = "data/update.{d}-{D}.failures.csv",
    output:
        woohoo = protected("../dbs/genbank-{d}-{D}-k{k}.update.zip"),
    conda: "envs/sourmash.yaml"
#    resources:
#        mem_mb = 8000,
#        time_min = 30
    shell: """
        sourmash sig cat {input.dir} {input.db} -k {wildcards.k} -o {output.woohoo}
    """

rule check_txt_missing:
    input:
        missing = expand("data/missing-genomes.{d}-{D}-k{k}.csv", d=DATE, D=DOMAINS, k=KSIZES),
        script = "scripts/check_txt_files.py",
    output:
        solo = "data/missing.{d}-{D}.csv",
    shell: """
        files=""
        for file in "data/missing-genomes.*{wildcards.d}-{wildcards.D}-k*.csv"; do
            files+=" $file"
        done
        {input.script} $files -o {output.solo}
    """

rule gather_sketch_missing:
    input:
         missing = "data/missing.{d}-{D}.csv",
    output:
        failed = "data/missing.{d}-{D}.failures.csv",
        db = "../dbs/genbank-{d}-{D}.miss.zip"
    conda: "envs/directsketch.yaml"
    resources:
        allowed_jobs = 100
#    resources:
#        mem_mb = 8000,
#        time_min = 30,
    threads: 1
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        #k_list = lambda wildcards: f"k={','.join([f'{k}' for k in KSIZES])}",
    shell:'''
        sourmash scripts gbsketch {input.missing} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled=1000,abund" -r 1 #2> {log}
    '''

checkpoint cat_to_clean_missing:
    input:
        dir = "../dbs/genbank-{d}-{D}.miss.zip",
        db = "../dbs/genbank-{d}-{D}-k{k}.clean.zip",
        missing = "data/missing.{d}-{D}.failures.csv",
    output:
        woohoo = protected("../dbs/genbank-{d}-{D}-k{k}.zip"),
    conda: "envs/sourmash.yaml"
#    resources:
#        mem_mb = 8000,
#        time_min = 30
    shell: """
        sourmash sig cat {input.dir} {input.db} -k {wildcards.k} -o {output.woohoo}
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
