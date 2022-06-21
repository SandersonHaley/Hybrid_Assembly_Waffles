#!/bin/bash
#SBATCH --account=hpc_p_white
#SBATCH --constraint=skylake
#SBATCH --time=21-00:00:00
#SBATCH --mem=64G

line=$1
module load python
source $HOME/nanostat/bin/activate
#Check quality of raw long reads with Nanostat
#NanoStat -o "$line"/nanostat_results_unprocessed -p "$line" -n "$line"_nanostat_results --fastq "$line"/*LR_ALL.fastq

module load StdEnv/2020 gcc/9.3.0
module load porechop/0.2.4 
#Trim adaptors from long reads
#porechop -i "$line"/"$line"_LR_ALL.fastq -o "$line"/"$line"_LR_ALL_trimmed.fastq --verbosity 1


module load nixpkgs/16.09
module load nixpkgs/16.09
module load gcc/5.4.0
module load intel/2016.4
module load filtlong
#Filter out bad and too short long reads 
#filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 "$line"/"$line"_LR_ALL_trimmed.fastq | gzip > "$line"/"$line"_LR_ALL_trimmed_filtered.fastq.gz


module load python
source $HOME/nanostat/bin/activate
#Get new quality of the remaining reads used for long read assembly
#NanoStat -o "$line"/nanostat_results -p "$line" -n "$line"_nanostat_results --fastq "$line"/*.fq.gz


module load StdEnv/2020
module load gcc/9.3.0
module load flye/2.9
#Assemble long reads using Flye
#flye --nano-raw "$line"/"$line"_LR_ALL_trimmed_filtered.fastq.gz -g 6m -o "$line"/flye_output -t 15 --plasmid -i 5


module load StdEnv/2020
module load gcc/9.3.0
module load quast/5.0.2
#Get quality of long read assembly with Quast
#quast.py "$line"/flye_output/assembly.fasta -o "$line"/flye_output/quast

#Move on to preparing short reads for hybrid assembly


module load gentoo/2020
module load fastp
#Trim and remove adapters using fastp from short reads
#fastp -i "$line"/"$line"_1.fq.gz -o "$line"/"$line"_SR_1.trimmed.fastq -I "$line"/"$line"_2.fq.gz -O "$line"/"$line"_SR_2.trimmed.fastq --detect_adapter_for_pe


module load StdEnv/2020
module load gcc/9.3.0
module load nixpkgs/16.09
module load fastqc/0.11.9
#check quality of short reads being used for hybrid assembly
#fastqc "$line"/"$line"_SR_*.trimmed.fastq.gz


module load python
source $HOME/unicycler/bin/activate
module load nixpkgs/16.09
module load gcc/9.3.0
module load StdEnv/2020
module load spades/3.14.1
module load racon/1.4.13
module load bowtie2/2.4.1
module load blast+/2.12.0
module load java/1.8
module load pilon
module load samtools/1.13

#Hybrid assembly using unicycler with the long read assembly and the trimmed short reads
#unicycler -t 15 -o "$line"/"$line"_unicycler_out_trimmed -1 "$line"/"$line"_SR_1.trimmed.fastq -2 "$line"/"$line"_SR_2.trimmed.fastq --mode bold -l "$line"/flye_output/assembly.fasta 

#If assembly with Unicycler is not successful you can use my other work flow to polish Flye assemblies with short reads since Flye closes the most genomes and polishing with short reads will deal with many of the errors in the sequence


module load StdEnv/2020
module load gcc/9.3.0
module load quast/5.0.2
#Check quality of the hybrid assembly if successful
#quast.py "$line"/"$line"_unicycler_out_trimmed/assembly.fasta -o "$line"/"$line"_unicycler_out_trimmed/quast_results

#Make a directory to do polishing in and move to that directory
#mkdir "$line"/"$line"_polishing
cd "$line"/"$line"_polishing




module load StdEnv/2020
module load bwa/0.7.17
module load gcc/9.3.0
module load samtools/1.13
module load pilon/1.24
module load mummer/4.0.0beta2


#Polish hybrid assembly with raw short reads using Pilon 15 times

#bwa index ~/Hybrid_Assembly_Test/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta 
#bwa mem  -v 2 -M -t 15 ~/Hybrid_Assembly_Test/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T ~/Hybrid_Assembly_Test/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta -F 3844 -q 60 | samtools sort --reference ~/Hybrid_Assembly_Test/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome ~/Hybrid_Assembly_Test/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta --frags "$line".pilon.bam --output "$line".pilon.0 


#nucmer -p "$line".pilon.0.nucmer ~/Hybrid_Assembly_Test/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta "$line".pilon.0.fasta
#show-snps -C -T -r "$line".pilon.0.nucmer.delta > "$line".pilon.0.nucmer.delta.log


#bwa index "$line".pilon.0.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.0.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.0.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.0.fasta > "$line".pilon.bam


#samtools index "$line".pilon.bam

#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.0.fasta --frags "$line".pilon.bam --output "$line".pilon.1 


#nucmer -p "$line".pilon.1.nucmer "$line".pilon.0.fasta "$line".pilon.1.fasta
#show-snps -C -T -r "$line".pilon.1.nucmer.delta > "$line".pilon.1.nucmer.delta.log



#bwa index "$line".pilon.1.fasta
#bwa mem  -v 2 -M -t 15 "$line".pilon.1.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.1.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.1.fasta > "$line".pilon.bam 

#samtools index "$line".pilon.bam

#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.1.fasta --frags "$line".pilon.bam --output "$line".pilon.2


#nucmer -p "$line".pilon.2.nucmer "$line".pilon.1.fasta "$line".pilon.2.fasta
#show-snps -C -T -r "$line".pilon.2.nucmer.delta > "$line".pilon.2.nucmer.delta.log

#bwa index "$line".pilon.2.fasta
#bwa mem  -v 2 -M -t 15 "$line".pilon.2.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.2.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.2.fasta > "$line".pilon.bam
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.2.fasta --frags "$line".pilon.bam --output "$line".pilon.3

#nucmer -p "$line".pilon.3.nucmer "$line".pilon.2.fasta "$line".pilon.3.fasta
#show-snps -C -T -r "$line".pilon.3.nucmer.delta > "$line".pilon.3.nucmer.delta.log


#bwa index "$line".pilon.3.fasta
#bwa mem  -v 2 -M -t 15 "$line".pilon.3.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.3.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.3.fasta > "$line".pilon.bam
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.3.fasta --frags "$line".pilon.bam --output "$line".pilon.4


#nucmer -p "$line".pilon.4.nucmer "$line".pilon.3.fasta "$line".pilon.4.fasta
#show-snps -C -T -r "$line".pilon.4.nucmer.delta > "$line".pilon.4.nucmer.delta.log


#bwa index "$line".pilon.4.fasta
#bwa mem  -v 2 -M -t 15 "$line".pilon.4.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.4.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.4.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam



#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.4.fasta --frags "$line".pilon.bam --output "$line".pilon.5 


#nucmer -p "$line".pilon.5.nucmer "$line".pilon.4.fasta "$line".pilon.5.fasta
#show-snps -C -T -r "$line".pilon.5.nucmer.delta > "$line".pilon.5.nucmer.delta.log


#bwa index "$line".pilon.5.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.5.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.5.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.5.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam

#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.5.fasta --frags "$line".pilon.bam --output "$line".pilon.6 


#nucmer -p "$line".pilon.6.nucmer "$line".pilon.5.fasta "$line".pilon.6.fasta
#show-snps -C -T -r "$line".pilon.6.nucmer.delta > "$line".pilon.6.nucmer.delta.log


#bwa index "$line".pilon.6.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.6.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.6.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.6.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam

#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.6.fasta --frags "$line".pilon.bam --output "$line".pilon.7 

#nucmer -p "$line".pilon.7.nucmer "$line".pilon.6.fasta "$line".pilon.7.fasta
#show-snps -C -T -r "$line".pilon.7.nucmer.delta > "$line".pilon.7.nucmer.delta.log


#bwa index "$line".pilon.7.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.7.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.7.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.7.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.7.fasta --frags "$line".pilon.bam --output "$line".pilon.8 


#nucmer -p "$line".pilon.8.nucmer "$line".pilon.7.fasta "$line".pilon.8.fasta
#show-snps -C -T -r "$line".pilon.8.nucmer.delta > "$line".pilon.8.nucmer.delta.log



#bwa index "$line".pilon.8.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.8.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.8.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.8.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.8.fasta --frags "$line".pilon.bam --output "$line".pilon.9 


#nucmer -p "$line".pilon.9.nucmer "$line".pilon.8.fasta "$line".pilon.9.fasta
#show-snps -C -T -r "$line".pilon.9.nucmer.delta > "$line".pilon.9.nucmer.delta.log


#bwa index "$line".pilon.9.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.9.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.9.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.9.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.9.fasta --frags "$line".pilon.bam --output "$line".pilon.10 


#nucmer -p "$line".pilon.10.nucmer "$line".pilon.9.fasta "$line".pilon.10.fasta
#show-snps -C -T -r "$line".pilon.10.nucmer.delta > "$line".pilon.10.nucmer.delta.log


#bwa index "$line".pilon.10.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.10.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.10.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.10.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.10.fasta --frags "$line".pilon.bam --output "$line".pilon.11 

#nucmer -p "$line".pilon.11.nucmer "$line".pilon.10.fasta "$line".pilon.11.fasta
#show-snps -C -T -r "$line".pilon.11.nucmer.delta > "$line".pilon.11.nucmer.delta.log


#bwa index "$line".pilon.11.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.11.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.11.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.11.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.11.fasta --frags "$line".pilon.bam --output "$line".pilon.12 


#nucmer -p "$line".pilon.12.nucmer "$line".pilon.11.fasta "$line".pilon.12.fasta
#show-snps -C -T -r "$line".pilon.12.nucmer.delta > "$line".pilon.12.nucmer.delta.log

#bwa index "$line".pilon.12.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.12.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.12.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.12.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.12.fasta --frags "$line".pilon.bam --output "$line".pilon.13

#nucmer -p "$line".pilon.13.nucmer "$line".pilon.12.fasta "$line".pilon.13.fasta
#show-snps -C -T -r "$line".pilon.13.nucmer.delta > "$line".pilon.13.nucmer.delta.log



#bwa index "$line".pilon.13.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon.13.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon.13.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.13.fasta > "$line".pilon.bam 
#samtools index "$line".pilon.bam

#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon.13.fasta --frags "$line".pilon.bam --output "$line".pilon.14 

#nucmer -p "$line".pilon.14.nucmer "$line".pilon.13.fasta "$line".pilon.14.fasta
#show-snps -C -T -r "$line".pilon.14.nucmer.delta > "$line".pilon.14.nucmer.delta.log

#Most assemblies will contain no snippy after the first 12 rounds of pilon. If there still is SNPs then additional polishing with SNIPPY and then java -Xmx30G -jar $EBROOTPILON/pilon.jar is required

#Run snippy and do corrections on the assembly (need tabix and vcf-consensus). You use the last fasta generated in the Pilon polishing as the reference. 

module load gentoo/2020
module load snippy/4.6.0
module load StdEnv/2020
module load vcftools/0.1.16
#snippy --cpus 8 --outdir "$line"_snippy --ref "$line".pilon.14.fasta --R1 ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz --R2 ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz --force
#cp ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/snps.vcf ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf 
#bgzip -c ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf > ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf.gz

module load nixpkgs/16.09
module load intel/2016.4
module load tabix/0.2.6

#tabix -f -p vcf ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf.gz

module load StdEnv/2020
module load vcftools/0.1.16

#cat "$line".pilon.14.fasta| vcf-consensus ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf.gz > ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa
 

#Now polished the corrected.fa from snippy throguh 15 more Pilon iterations
module load StdEnv/2020
module load bwa/0.7.17
module load gcc/9.3.0
module load samtools/1.13
module load pilon/1.24
module load mummer/4.0.0beta2


#bwa index ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa
#bwa mem  -v 2 -M -t 15 ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa -F 3844 -q 60 | samtools sort --reference ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam



#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa --frags "$line".pilon2.bam --output "$line".pilon2.0 


#nucmer -p "$line".pilon2.0.nucmer ~/Hybrid_Assembly_Test/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa "$line".pilon2.0.fasta
#show-snps -C -T -r "$line".pilon2.0.nucmer.delta > "$line".pilon2.0.nucmer.delta.log


#bwa index "$line".pilon2.0.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.0.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.0.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.0.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.0.fasta --frags "$line".pilon2.bam --output "$line".pilon2.1 


#nucmer -p "$line".pilon2.1.nucmer "$line".pilon2.0.fasta "$line".pilon2.1.fasta
#show-snps -C -T -r "$line".pilon2.1.nucmer.delta > "$line".pilon2.1.nucmer.delta.log


#bwa index "$line".pilon2.1.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.1.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.1.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.1.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
 
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.1.fasta --frags "$line".pilon2.bam --output "$line".pilon2.2 



#nucmer -p "$line".pilon2.2.nucmer "$line".pilon2.1.fasta "$line".pilon2.2.fasta
#show-snps -C -T -r "$line".pilon2.2.nucmer.delta > "$line".pilon2.2.nucmer.delta.log


#bwa index "$line".pilon2.2.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.2.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.2.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.2.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.2.fasta --frags "$line".pilon2.bam --output "$line".pilon2.3 


#nucmer -p "$line".pilon2.3.nucmer "$line".pilon2.2.fasta "$line".pilon2.3.fasta
#show-snps -C -T -r "$line".pilon2.3.nucmer.delta > "$line".pilon2.3.nucmer.delta.log


#bwa index "$line".pilon2.3.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.3.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.3.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.3.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam



#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.3.fasta --frags "$line".pilon2.bam --output "$line".pilon2.4 


#nucmer -p "$line".pilon2.4.nucmer "$line".pilon2.3.fasta "$line".pilon2.4.fasta
#show-snps -C -T -r "$line".pilon2.4.nucmer.delta > "$line".pilon2.4.nucmer.delta.log

#bwa index "$line".pilon2.4.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.4.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.4.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.4.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam


#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.4.fasta --frags "$line".pilon2.bam --output "$line".pilon2.5 


#nucmer -p "$line".pilon2.5.nucmer "$line".pilon2.4.fasta "$line".pilon2.5.fasta
#show-snps -C -T -r "$line".pilon2.5.nucmer.delta > "$line".pilon2.5.nucmer.delta.log
#bwa index "$line".pilon2.5.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.5.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.5.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.5.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.5.fasta --frags "$line".pilon2.bam --output "$line".pilon2.6 

#nucmer -p "$line".pilon2.6.nucmer "$line".pilon2.5.fasta "$line".pilon2.6.fasta
#show-snps -C -T -r "$line".pilon2.6.nucmer.delta > "$line".pilon2.6.nucmer.delta.log
#bwa index "$line".pilon2.6.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.6.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.6.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.6.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.6.fasta --frags "$line".pilon2.bam --output "$line".pilon2.7 

#nucmer -p "$line".pilon2.7.nucmer "$line".pilon2.6.fasta "$line".pilon2.7.fasta
#show-snps -C -T -r "$line".pilon2.7.nucmer.delta > "$line".pilon2.7.nucmer.delta.log
#bwa index "$line".pilon2.7.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.7.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.7.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.7.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.7.fasta --frags "$line".pilon2.bam --output "$line".pilon2.8 

#nucmer -p "$line".pilon2.8.nucmer "$line".pilon2.7.fasta "$line".pilon2.8.fasta
#show-snps -C -T -r "$line".pilon2.8.nucmer.delta > "$line".pilon2.8.nucmer.delta.log
#bwa index "$line".pilon2.8.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.8.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.8.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.8.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.8.fasta --frags "$line".pilon2.bam --output "$line".pilon2.9 

#nucmer -p "$line".pilon2.9.nucmer "$line".pilon2.8.fasta "$line".pilon2.9.fasta
#show-snps -C -T -r "$line".pilon2.9.nucmer.delta > "$line".pilon2.9.nucmer.delta.log
#bwa index "$line".pilon2.9.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.9.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.9.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.9.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.9.fasta --frags "$line".pilon2.bam --output "$line".pilon2.10 

#nucmer -p "$line".pilon2.10.nucmer "$line".pilon2.9.fasta "$line".pilon2.10.fasta
#show-snps -C -T -r "$line".pilon2.10.nucmer.delta > "$line".pilon2.10.nucmer.delta.log
#bwa index "$line".pilon2.10.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.10.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.10.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.10.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.10.fasta --frags "$line".pilon2.bam --output "$line".pilon2.11 

#nucmer -p "$line".pilon2.11.nucmer "$line".pilon2.10.fasta "$line".pilon2.11.fasta
#show-snps -C -T -r "$line".pilon2.11.nucmer.delta > "$line".pilon2.11.nucmer.delta.log
#bwa index "$line".pilon2.11.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.11.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.11.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.11.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.11.fasta --frags "$line".pilon2.bam --output "$line".pilon2.12 

#nucmer -p "$line".pilon2.12.nucmer "$line".pilon2.11.fasta "$line".pilon2.12.fasta
#show-snps -C -T -r "$line".pilon2.12.nucmer.delta > "$line".pilon2.12.nucmer.delta.log
#bwa index "$line".pilon2.12.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.12.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.12.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.12.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.12.fasta --frags "$line".pilon2.bam --output "$line".pilon2.13 

#nucmer -p "$line".pilon2.13.nucmer "$line".pilon2.12.fasta "$line".pilon2.13.fasta
#show-snps -C -T -r "$line".pilon2.13.nucmer.delta > "$line".pilon2.13.nucmer.delta.log
#bwa index "$line".pilon2.13.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.13.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz | samtools view -u -T "$line".pilon2.13.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.13.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#java -Xmx30G -jar $EBROOTPILON/pilon.jar --genome "$line".pilon2.13.fasta --frags "$line".pilon2.bam --output "$line".pilon2.14 

#nucmer -p "$line".pilon2.14.nucmer "$line".pilon2.13.fasta "$line".pilon2.14.fasta
#show-snps -C -T -r "$line".pilon2.14.nucmer.delta > "$line".pilon2.14.nucmer.delta.log 
#Most of the genomes should now have no snps. Check before and after snippy to see if that remaining SNPs in the genomes with them are just repeated showing up. If this is the case, you cannot resolve them. If they are not the same repeat the snippy polishing and rounds of pilon. 



module load StdEnv/2020
module load bwa/0.7.17
module load gcc/9.3.0
module load samtools/1.13

#Map the short reads to the assembly
bwa index "$line".pilon2.14.fasta
bwa mem -v 2 -M -t 15 "$line".pilon2.14.fasta ~/Hybrid_Assembly_Test/"$line"/"$line"_1.fq.gz ~/Hybrid_Assembly_Test/"$line"/"$line"_2.fq.gz|
samtools view -u -T "$line".pilon2.14.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.14.fasta >> "$line"_illumina_mapped.bam


module load StdEnv/2020
module load gcc/9.3.0
module load quast/5.0.2

#get quality of the complete polished genome
quast.py "$line".pilon2.14.fasta -o "$line"_quast_output
mv "$line".pilon2.14.fasta "$line"_complete_hybrid_assembly.fasta
