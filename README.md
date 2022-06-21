# Hybrid_Assembly_Waffles
Hybrid Assembly using NML HPC Waffles

Required Packages (must be made into conda environemtns as seen in the script): Nanostat Porechop Filtlong Flye Quast Fastp FastQC Unicycler bwa samtools pilon mummer snippy vcftools tabix

It is designed to have a preparation file that will look something like this:

#!/bin/bash cat text.txt |while read line; do echo "Submitting job for file $line" sbatch hybridassemblywaffles.sh "$line" done

That allows you to create a list of filenames in "text.txt" that corresponds to a variable screated with each line of the .txt file. In the hybridassembly.sh file each line will be "$line". This prep file will submit the sbatch for each name in each line of text.txt and feed the line in as "$line" in the hybridassembly.sh file. I usually have everything in one directory and have subdirectories named for each isolate with the Long and Short reads in that subdirectory.

(for example:
/"$line"/"$line"_(SR_1/LR_ALL).fastq)

Replace you HPC/slurm parameters at the top of the page.
First the script does the QC, trimming and filtering of the Nanopore long reads then assembles these reads using Flye. Then it does the QC, trimming and filtering of the Illumina short reads. The Nanopore and Illumina read are then assembled into a hybrid genome using Unicycler and the quality of the assembly is determined using Quast. Next, a folder is made to hold the polishing files and we make that file that current directory. Then the polishing steps begin with 15 pilon polishing steps. Then the Illumina reads are mapped to the polished genome and Quast is run on the genome. The assembly is then renamed to $line"_complete_hybrid_assembly.fasta.
There may be additional snps still in the genome, however, removing these snps by snippy has bugs on Waffles with the tabix program so that part is commented out. Most of the time, when I need this step I swich HPCs or complete the snippy steps on my computer than return the files to Waffles for the pilon polishing.

I hope to eventually make this a Nextflow pipeline but it works for now.
