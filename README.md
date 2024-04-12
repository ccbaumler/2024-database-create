I would like to create some databases for sourmash
- signature databases
- taxonomy databases

## 1. Get the data

First, make a `data` directory.

```
mkdir data && cd data
```

To download all the complete genomes of ncbi genbank, I will download genbank's `assembly_summary` file.
(View the file at https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt)

> WARNING: I can not find all the genome accessions within the base assembly summary. Don't know why... But, Go into each sub-directory and merge the files into one large and complete genome dataset.
> archaea/
> bacteria/
> fungi/
> invertebrate/
> plant/
> protozoa/
> unknown/
> vertebrate_mammalian/
> vertebrate_other/
> viral/                   

```
curl -O -L 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt'
```

I will filter the summary file for complete genomes and 
(Suggested here -> https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)

```
awk -F "\t" '$12=="Complete Genome" && $11=="latest" {print $20}' assembly_summary.txt > ftpdirpaths
```

> Column 11 == version_status and column 12 == assembly_level
> To get a breakdown of these columns use:
> `awk -F "\t" '{print $11}' assembly_summary_genbank.txt | sort | uniq` and `awk -F "\t" '{print $12}' assembly_summary_genbank.txt | sort | uniq`
> Column 11 contains only `latest` and -- its header -- `version_status`.
> Column 12 contains `Chromosome`, `Complete Genome`, `Scaffold`, and -- its header -- `assembly_level`

For a more robust analysis pipeline, we need to include the replaced and suspended genomes within genbank
(View the file at https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt)

```
curl -O -L 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt'
```

> Column 11 == version_status and column 12 == assembly_level
> To get a breakdown of these columns use:
> `awk -F "\t" '{print $11}' assembly_summary_genbank.txt | sort | uniq` and `awk -F "\t" '{print $12}' assembly_summary_genbank.txt | sort | uniq`
> Column 11 contains `replaced`, `suspended`, and -- its header -- `version_status`.
> Column 12 contains `Chromosome`, `Complete Genome`, `Scaffold`, and -- its header -- `assembly_level`


