###
# This workflow will create an allthebacteria database for sourmash from the HPC cluster at UCDavis.
#
# To run:
# snakemake -s allthebacteria.smk -j 1 --use-conda
#
###

RELEASES=[0.2,]

KSIZES=[21,31,51]
LOGS='logs'

# do not include the parent directory (--no-parent), host directory 'ftp.ebi.ac.uk' (--no-host-directory), or the long path of directories 'pub/databases/AllTheBacteria/Releases/0.1/' (--cut-dirs=5). Recursively download everything in a release (--recursive), continue the download if stopped (--continue), and do not re-download the same completed file (--no-clobber). Place all the files from the ftp server in a new directory and mkdir if directory does not exist (--directory-prefix)
#wget --no-parent --recursive --continue --no-clobber --no-host-directories --cut-dirs=5 --directory-prefix={VERSIONS}/ ftp://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/{VERSIONS}/

rule all:
    input:
        expand("data/r{r}-data-links.txt", r=RELEASES),
        #expand("allthebacteria-r{r}-metadata/sample_list.txt", r=RELEASES),
        expand("allthebacteria-r{r}-metadata/sylph.tsv", r=RELEASES),

rule make_links:
    input:
        script = "scripts/ftp_link_list.py",
    output:
        meta = "data/r{r}-metadata-links.txt",
        data = "data/r{r}-data-links.txt",
    shell: """
        {input.script} ftp.ebi.ac.uk pub/databases/AllTheBacteria/Releases/{wildcards.r}/metadata -s ftp -o {output.meta}

        {input.script} ftp.ebi.ac.uk pub/databases/AllTheBacteria/Releases/{wildcards.r}/assembly -s ftp -o {output.data}
    """

rule download_links:
    input:
        meta = "data/r{r}-metadata-links.txt",
        data = "data/r{r}-data-links.txt",
    output:
        tsv = "allthebacteria-r{r}-metadata/sample2species2file.tsv.gz",
        info = "allthebacteria-r{r}-metadata/ena_metadata.tsv.gz",
    shell: """
        # read the txt files, use only 1 link (-n), parallelize wget cmd upto 10 times (-P), path for downloads (-P), if interrupted continue the download (--continue), do not download anything that exists (--no-clobber), try up to 3 times per download (--tries), wait 1 sec per try (--wait), quiet output (-nv) 
        cat {input.data} | xargs -n 1 -P 10 wget -P allthebacteria-r{wildcards.r}-data/ --continue --no-clobber --tries=3 --wait=1 -nv && cat {input.meta} | xargs -n 1 -P 10 wget -P allthebacteria-r{wildcards.r}-metadata/ --continue --no-clobber --tries=3 --wait=1 -nv
    """

#checkpoint check:
#    input:
#        tsv = "allthebacteria-r{r}-metadata/sample_list.txt.gz",
#    output:
#        tsv = "allthebacteria-r{r}-metadata/sample_list.txt",
#    shell:"""
#        #https://stackoverflow.com/questions/20449543/shell-equality-operators-eq
#
#        tsv_count=$( gzip -cd {input.tsv} | awk -F"\\t" '{{print $3}}' | uniq | wc -l )
#
#        file_count=$( ls -l allthebacteria-r{wildcards.r}-data | wc -l )
#
#        # Both counts have an extra value (tsv_counts has 'header' and file_count has 'total')
#        result=$( (( tsv_count == file_count )); echo $? )
#
#        if [ $result -eq 0 ]; then
#            echo "The downloaded file count matches the expected file count from metadata"
#            echo "$file_count files downloads and $tsv_count files expected\n"
#            gzip -kvd {input.tsv}
#        else
#            echo "The downloaded file count DOES NOT matches the expected file count from metadata"
#            echo "$file_count files downloads and $tsv_count files expected\n"
#            echo "Stopping now!"
#        fi
#    """

rule sketch_sigs:
    input:
        #tsv = "allthebacteria-r{r}-metadata/sample_list.txt",
        info = "allthebacteria-r{r}-metadata/sylph.tsv.gz",
    output:
        info = "allthebacteria-r{r}-metadata/sylph.tsv"
    conda: "envs/branchwater.yaml"
    shell:"""
        gzip -kvd {input.info}

        tar_file_names=$(ls -1 allthebacteria-r{wildcards.r}-data)
        seq_dir='allthebacteria-seqs'
        sig_dir='allthebacteria-sigs'

        if [ ! -d "$seq_dir" ]; then
          echo "Making $seq_dir..." ; echo ''

          mkdir $seq_dir

        else
          echo "Using existing '$seq_dir'" ; echo ''

        fi

        if [ ! -d "$sig_dir" ]; then
          echo "Making $sig_dir..." ; echo ''

          mkdir $sig_dir

        else
          echo "Using existing '$sig_dir'" ; echo ''

        fi

        for file in $tar_file_names; do
          echo "Extracting all files from $file ..." ; echo ''

          #apparently, this is optimized for xz files and will use all available CPUs to extract
          pv allthebacteria-r{wildcards.r}-data/$file | tar --use-compress-program="xz -T0 -q" --skip-old-files -xf - -C $seq_dir

          echo "Extraction complete!" ; echo ''
          dir_name="${{file%%.*}}"

          if [ -d "$seq_dir/$dir_name" ]; then
            echo "Building 'manysketch.csv' for '$seq_dir/$dir_name...'" ; echo ''
            echo name,genome_filename,protein_filename > $seq_dir/$dir_name/manysketch.csv

            num_files=$(ls -1 $seq_dir/$dir_name/* | wc -l)

            # find value for 5% increment (bc is for floating point math)
            increment=$( bc <<< "$num_files * 20 / 100" )  # 5% of total files
            target_line=0

            line_count=0
            echo "working"
            ls -1 $seq_dir/$dir_name/* | while read filepath; do
              filename_ext=$(basename "$filepath")
              filename="${{filename_ext%%.*}}"

              #grab the name hidden in the metadata and remove any commas'
              name=$(awk -F'\\t' -v fname="$filename" '$3 == fname {{print $17, $3; exit}}' {output.info} | sed 's/,//g')
              echo "$name,$filepath," >> $seq_dir/$dir_name/manysketch.csv

              #count iteration of lines
              line_count=$(($line_count + 1))
          
              # line count and target percent check
              if [ $line_count -ge $target_line ]; then
                echo "Progress: $line_count/$num_files lines written"
          
                #Update the target for the next 5% increment
                target_line=$(($target_line + $increment))
              fi
            done

            echo "'$seq_dir/$dir_name/manysketch.csv' built. Sketching all seq files in '$seq_dir/$dir_name'!!!" ; echo ''
            sourmash scripts manysketch $seq_dir/$dir_name/manysketch.csv -p k=21,k=31,k=51,scaled=1000,abund -o $sig_dir/$dir_name.zip
            echo "Sketching completed at '$sig_dir/$dir_name.zip'! Removing sequence files at '$seq_dir/$dir_name'"  ; echo ''

            rm -r $seq_dir/$dir_name 

          else
            echo "'$dir_name' is not a directory." ; echo ''

          fi

        done
   """
