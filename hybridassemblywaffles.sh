#!/bin/bash
#SBATCH -p ExternalResearch
#SBATCH -c 1
#SBATCH --time=0
#SBATCH --mem=150G

line=$1
source activate longreadprep
#Check quality of raw long reads with Nanostat
NanoStat -o "$line"/nanostat_results_unprocessed -p "$line" -n "$line"_nanostat_results --fastq "$line"/*LR_ALL.fastq


#Trim adaptors from long reads
porechop -i "$line"/"$line"_LR_ALL.fastq -o "$line"/"$line"_LR_ALL_trimmed.fastq --verbosity 1



#Filter out bad and too short long reads 
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 "$line"/"$line"_LR_ALL_trimmed.fastq | gzip > "$line"/"$line"_LR_ALL_trimmed_filtered.fastq.gz

#Get new quality of the remaining reads used for long read assembly
NanoStat -o "$line"/nanostat_results -p "$line" -n "$line"_nanostat_results --fastq "$line"/*.fq.gz

#Assemble long reads using Flye
flye --nano-raw "$line"/"$line"_LR_ALL_trimmed_filtered.fastq.gz -g 6m -o "$line"/flye_output -t 15 --plasmid -i 5


source activate quast
#Get quality of long read assembly with Quast
quast.py "$line"/flye_output/assembly.fasta -o "$line"/flye_output/quast

#Move on to preparing short reads for hybrid assembly

source activate hybridassembly
#Trim and remove adapters using fastp from short reads
fastp -i "$line"/"$line"_SR_1.fq.gz -o "$line"/"$line"_SR_1.trimmed.fastq -I "$line"/"$line"_SR_2.fq.gz -O "$line"/"$line"_SR_2.trimmed.fastq --detect_adapter_for_pe


#check quality of short reads being used for hybrid assembly
fastqc "$line"/"$line"_SR_*.trimmed.fastq.gz


#Hybrid assembly using unicycler with the long read assembly and the trimmed short reads
unicycler -t 15 -o "$line"/"$line"_unicycler_out_trimmed -1 "$line"/"$line"_SR_1.trimmed.fastq -2 "$line"/"$line"_SR_2.trimmed.fastq --mode bold -l "$line"/flye_output/assembly.fasta 

#If assembly with Unicycler is not successful you can use my other work flow to polish Flye assemblies with short reads since Flye closes the most genomes and polishing with short reads will deal with many of the errors in the sequence


source activate quast
#Check quality of the hybrid assembly if successful
quast.py "$line"/"$line"_unicycler_out_trimmed/assembly.fasta -o "$line"/"$line"_unicycler_out_trimmed/quast_results

#Make a directory to do polishing in and move to that directory
mkdir "$line"/"$line"_polishing
cd "$line"/"$line"_polishing

source activate assemblypolishing

#Polish hybrid assembly with raw short reads using Pilon 15 times

bwa index ~/6isolatehybridredo/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta 
bwa mem  -v 2 -M -t 15 ~/6isolatehybridredo/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T ~/6isolatehybridredo/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta -F 3844 -q 60 | samtools sort --reference ~/6isolatehybridredo/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam


pilon --genome ~/6isolatehybridredo/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta --frags "$line".pilon.bam --output "$line".pilon.0 


nucmer -p "$line".pilon.0.nucmer ~/6isolatehybridredo/"$line"/"$line"_unicycler_out_trimmed/assembly.fasta "$line".pilon.0.fasta
show-snps -C -T -r "$line".pilon.0.nucmer.delta > "$line".pilon.0.nucmer.delta.log


bwa index "$line".pilon.0.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.0.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.0.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.0.fasta > "$line".pilon.bam


samtools index "$line".pilon.bam

pilon --genome "$line".pilon.0.fasta --frags "$line".pilon.bam --output "$line".pilon.1 


nucmer -p "$line".pilon.1.nucmer "$line".pilon.0.fasta "$line".pilon.1.fasta
show-snps -C -T -r "$line".pilon.1.nucmer.delta > "$line".pilon.1.nucmer.delta.log



bwa index "$line".pilon.1.fasta
bwa mem  -v 2 -M -t 15 "$line".pilon.1.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.1.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.1.fasta > "$line".pilon.bam 

samtools index "$line".pilon.bam

pilon --genome "$line".pilon.1.fasta --frags "$line".pilon.bam --output "$line".pilon.2


nucmer -p "$line".pilon.2.nucmer "$line".pilon.1.fasta "$line".pilon.2.fasta
show-snps -C -T -r "$line".pilon.2.nucmer.delta > "$line".pilon.2.nucmer.delta.log

bwa index "$line".pilon.2.fasta
bwa mem  -v 2 -M -t 15 "$line".pilon.2.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.2.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.2.fasta > "$line".pilon.bam
samtools index "$line".pilon.bam


pilon --genome "$line".pilon.2.fasta --frags "$line".pilon.bam --output "$line".pilon.3

nucmer -p "$line".pilon.3.nucmer "$line".pilon.2.fasta "$line".pilon.3.fasta
show-snps -C -T -r "$line".pilon.3.nucmer.delta > "$line".pilon.3.nucmer.delta.log


bwa index "$line".pilon.3.fasta
bwa mem  -v 2 -M -t 15 "$line".pilon.3.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.3.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.3.fasta > "$line".pilon.bam
samtools index "$line".pilon.bam


pilon --genome "$line".pilon.3.fasta --frags "$line".pilon.bam --output "$line".pilon.4


nucmer -p "$line".pilon.4.nucmer "$line".pilon.3.fasta "$line".pilon.4.fasta
show-snps -C -T -r "$line".pilon.4.nucmer.delta > "$line".pilon.4.nucmer.delta.log


bwa index "$line".pilon.4.fasta
bwa mem  -v 2 -M -t 15 "$line".pilon.4.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.4.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.4.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam



pilon --genome "$line".pilon.4.fasta --frags "$line".pilon.bam --output "$line".pilon.5 


nucmer -p "$line".pilon.5.nucmer "$line".pilon.4.fasta "$line".pilon.5.fasta
show-snps -C -T -r "$line".pilon.5.nucmer.delta > "$line".pilon.5.nucmer.delta.log


bwa index "$line".pilon.5.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.5.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.5.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.5.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam

pilon --genome "$line".pilon.5.fasta --frags "$line".pilon.bam --output "$line".pilon.6 


nucmer -p "$line".pilon.6.nucmer "$line".pilon.5.fasta "$line".pilon.6.fasta
show-snps -C -T -r "$line".pilon.6.nucmer.delta > "$line".pilon.6.nucmer.delta.log


bwa index "$line".pilon.6.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.6.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.6.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.6.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam

pilon --genome "$line".pilon.6.fasta --frags "$line".pilon.bam --output "$line".pilon.7 

nucmer -p "$line".pilon.7.nucmer "$line".pilon.6.fasta "$line".pilon.7.fasta
show-snps -C -T -r "$line".pilon.7.nucmer.delta > "$line".pilon.7.nucmer.delta.log


bwa index "$line".pilon.7.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.7.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.7.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.7.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam


pilon --genome "$line".pilon.7.fasta --frags "$line".pilon.bam --output "$line".pilon.8 


nucmer -p "$line".pilon.8.nucmer "$line".pilon.7.fasta "$line".pilon.8.fasta
show-snps -C -T -r "$line".pilon.8.nucmer.delta > "$line".pilon.8.nucmer.delta.log



bwa index "$line".pilon.8.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.8.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.8.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.8.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam


pilon --genome "$line".pilon.8.fasta --frags "$line".pilon.bam --output "$line".pilon.9 


nucmer -p "$line".pilon.9.nucmer "$line".pilon.8.fasta "$line".pilon.9.fasta
show-snps -C -T -r "$line".pilon.9.nucmer.delta > "$line".pilon.9.nucmer.delta.log


bwa index "$line".pilon.9.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.9.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.9.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.9.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam


pilon --genome "$line".pilon.9.fasta --frags "$line".pilon.bam --output "$line".pilon.10 


nucmer -p "$line".pilon.10.nucmer "$line".pilon.9.fasta "$line".pilon.10.fasta
show-snps -C -T -r "$line".pilon.10.nucmer.delta > "$line".pilon.10.nucmer.delta.log


bwa index "$line".pilon.10.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.10.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.10.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.10.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam


pilon --genome "$line".pilon.10.fasta --frags "$line".pilon.bam --output "$line".pilon.11 

nucmer -p "$line".pilon.11.nucmer "$line".pilon.10.fasta "$line".pilon.11.fasta
show-snps -C -T -r "$line".pilon.11.nucmer.delta > "$line".pilon.11.nucmer.delta.log


bwa index "$line".pilon.11.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.11.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.11.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.11.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam


pilon --genome "$line".pilon.11.fasta --frags "$line".pilon.bam --output "$line".pilon.12 


nucmer -p "$line".pilon.12.nucmer "$line".pilon.11.fasta "$line".pilon.12.fasta
show-snps -C -T -r "$line".pilon.12.nucmer.delta > "$line".pilon.12.nucmer.delta.log

bwa index "$line".pilon.12.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.12.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.12.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.12.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam


pilon --genome "$line".pilon.12.fasta --frags "$line".pilon.bam --output "$line".pilon.13

nucmer -p "$line".pilon.13.nucmer "$line".pilon.12.fasta "$line".pilon.13.fasta
show-snps -C -T -r "$line".pilon.13.nucmer.delta > "$line".pilon.13.nucmer.delta.log



bwa index "$line".pilon.13.fasta 
bwa mem  -v 2 -M -t 15 "$line".pilon.13.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon.13.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon.13.fasta > "$line".pilon.bam 
samtools index "$line".pilon.bam

pilon --genome "$line".pilon.13.fasta --frags "$line".pilon.bam --output "$line".pilon.14 

nucmer -p "$line".pilon.14.nucmer "$line".pilon.13.fasta "$line".pilon.14.fasta
show-snps -C -T -r "$line".pilon.14.nucmer.delta > "$line".pilon.14.nucmer.delta.log

#Most assemblies will contain no snippy after the first 12 rounds of pilon. If there still is SNPs then additional polishing with SNIPPY and then java -Xmx30G -jar $EBROOTPILON/pilon.jar is required

#Run snippy and do corrections on the assembly (need tabix and vcf-consensus). You use the last fasta generated in the Pilon polishing as the reference. 

source activate snippy
snippy --cpus 8 --outdir "$line"_snippy --ref "$line".pilon.14.fasta --R1 ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz --R2 ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz --force
source activate vcftools
cp ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/snps.vcf ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf 
source activate assemblypolishing
bgzip -c ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf > ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf.gz

#source activate tabix

#tabix -f -p vcf ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf.gz

#source activate vcftools

#cat "$line".pilon.14.fasta| vcf-consensus ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/corrections.vcf.gz > ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa
 

#Now polished the corrected.fa from snippy throguh 15 more Pilon iterations
#source activate assemblypolishing


#bwa index ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa
#bwa mem  -v 2 -M -t 15 ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa -F 3844 -q 60 | samtools sort --reference ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam



#pilon --genome ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa --frags "$line".pilon2.bam --output "$line".pilon2.0 


#nucmer -p "$line".pilon2.0.nucmer ~/6isolatehybridredo/"$line"/"$line"_polishing/"$line"_snippy/"$line"_corrected.fa "$line".pilon2.0.fasta
#show-snps -C -T -r "$line".pilon2.0.nucmer.delta > "$line".pilon2.0.nucmer.delta.log


#bwa index "$line".pilon2.0.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.0.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.0.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.0.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam


#pilon --genome "$line".pilon2.0.fasta --frags "$line".pilon2.bam --output "$line".pilon2.1 


#nucmer -p "$line".pilon2.1.nucmer "$line".pilon2.0.fasta "$line".pilon2.1.fasta
#show-snps -C -T -r "$line".pilon2.1.nucmer.delta > "$line".pilon2.1.nucmer.delta.log


#bwa index "$line".pilon2.1.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.1.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.1.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.1.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
 
#pilon --genome "$line".pilon2.1.fasta --frags "$line".pilon2.bam --output "$line".pilon2.2 



#nucmer -p "$line".pilon2.2.nucmer "$line".pilon2.1.fasta "$line".pilon2.2.fasta
#show-snps -C -T -r "$line".pilon2.2.nucmer.delta > "$line".pilon2.2.nucmer.delta.log


#bwa index "$line".pilon2.2.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.2.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.2.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.2.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam


#pilon --genome "$line".pilon2.2.fasta --frags "$line".pilon2.bam --output "$line".pilon2.3 


#nucmer -p "$line".pilon2.3.nucmer "$line".pilon2.2.fasta "$line".pilon2.3.fasta
#show-snps -C -T -r "$line".pilon2.3.nucmer.delta > "$line".pilon2.3.nucmer.delta.log


#bwa index "$line".pilon2.3.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.3.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.3.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.3.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam



#pilon --genome "$line".pilon2.3.fasta --frags "$line".pilon2.bam --output "$line".pilon2.4 


#nucmer -p "$line".pilon2.4.nucmer "$line".pilon2.3.fasta "$line".pilon2.4.fasta
#show-snps -C -T -r "$line".pilon2.4.nucmer.delta > "$line".pilon2.4.nucmer.delta.log

#bwa index "$line".pilon2.4.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.4.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.4.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.4.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam


#pilon --genome "$line".pilon2.4.fasta --frags "$line".pilon2.bam --output "$line".pilon2.5 


#nucmer -p "$line".pilon2.5.nucmer "$line".pilon2.4.fasta "$line".pilon2.5.fasta
#show-snps -C -T -r "$line".pilon2.5.nucmer.delta > "$line".pilon2.5.nucmer.delta.log
#bwa index "$line".pilon2.5.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.5.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.5.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.5.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.5.fasta --frags "$line".pilon2.bam --output "$line".pilon2.6 

#nucmer -p "$line".pilon2.6.nucmer "$line".pilon2.5.fasta "$line".pilon2.6.fasta
#show-snps -C -T -r "$line".pilon2.6.nucmer.delta > "$line".pilon2.6.nucmer.delta.log
#bwa index "$line".pilon2.6.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.6.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.6.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.6.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.6.fasta --frags "$line".pilon2.bam --output "$line".pilon2.7 

#nucmer -p "$line".pilon2.7.nucmer "$line".pilon2.6.fasta "$line".pilon2.7.fasta
#show-snps -C -T -r "$line".pilon2.7.nucmer.delta > "$line".pilon2.7.nucmer.delta.log
#bwa index "$line".pilon2.7.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.7.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.7.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.7.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.7.fasta --frags "$line".pilon2.bam --output "$line".pilon2.8 

#nucmer -p "$line".pilon2.8.nucmer "$line".pilon2.7.fasta "$line".pilon2.8.fasta
#show-snps -C -T -r "$line".pilon2.8.nucmer.delta > "$line".pilon2.8.nucmer.delta.log
#bwa index "$line".pilon2.8.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.8.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.8.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.8.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.8.fasta --frags "$line".pilon2.bam --output "$line".pilon2.9 

#nucmer -p "$line".pilon2.9.nucmer "$line".pilon2.8.fasta "$line".pilon2.9.fasta
#show-snps -C -T -r "$line".pilon2.9.nucmer.delta > "$line".pilon2.9.nucmer.delta.log
#bwa index "$line".pilon2.9.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.9.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.9.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.9.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.9.fasta --frags "$line".pilon2.bam --output "$line".pilon2.10 

#nucmer -p "$line".pilon2.10.nucmer "$line".pilon2.9.fasta "$line".pilon2.10.fasta
#show-snps -C -T -r "$line".pilon2.10.nucmer.delta > "$line".pilon2.10.nucmer.delta.log
#bwa index "$line".pilon2.10.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.10.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.10.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.10.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.10.fasta --frags "$line".pilon2.bam --output "$line".pilon2.11 

#nucmer -p "$line".pilon2.11.nucmer "$line".pilon2.10.fasta "$line".pilon2.11.fasta
#show-snps -C -T -r "$line".pilon2.11.nucmer.delta > "$line".pilon2.11.nucmer.delta.log
#bwa index "$line".pilon2.11.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.11.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.11.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.11.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.11.fasta --frags "$line".pilon2.bam --output "$line".pilon2.12 

#nucmer -p "$line".pilon2.12.nucmer "$line".pilon2.11.fasta "$line".pilon2.12.fasta
#show-snps -C -T -r "$line".pilon2.12.nucmer.delta > "$line".pilon2.12.nucmer.delta.log
#bwa index "$line".pilon2.12.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.12.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.12.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.12.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.12.fasta --frags "$line".pilon2.bam --output "$line".pilon2.13 

#nucmer -p "$line".pilon2.13.nucmer "$line".pilon2.12.fasta "$line".pilon2.13.fasta
#show-snps -C -T -r "$line".pilon2.13.nucmer.delta > "$line".pilon2.13.nucmer.delta.log
#bwa index "$line".pilon2.13.fasta 
#bwa mem  -v 2 -M -t 15 "$line".pilon2.13.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz | samtools view -u -T "$line".pilon2.13.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.13.fasta > "$line".pilon2.bam 
#samtools index "$line".pilon2.bam
#pilon --genome "$line".pilon2.13.fasta --frags "$line".pilon2.bam --output "$line".pilon2.14 

#nucmer -p "$line".pilon2.14.nucmer "$line".pilon2.13.fasta "$line".pilon2.14.fasta
#show-snps -C -T -r "$line".pilon2.14.nucmer.delta > "$line".pilon2.14.nucmer.delta.log 
#Most of the genomes should now have no snps. Check before and after snippy to see if that remaining SNPs in the genomes with them are just repeated showing up. If this is the case, you cannot resolve them. If they are not the same repeat the snippy polishing and rounds of pilon. 




#Map the short reads to the assembly
bwa index "$line".pilon.14.fasta
bwa mem -v 2 -M -t 15 "$line".pilon.14.fasta ~/6isolatehybridredo/"$line"/"$line"_SR_1.fq.gz ~/6isolatehybridredo/"$line"/"$line"_SR_2.fq.gz|
samtools view -u -T "$line".pilon.14.fasta -F 3844 -q 60 | samtools sort --reference "$line".pilon2.14.fasta >> "$line"_illumina_mapped.bam



source activate quast
#get quality of the complete polished genome
quast.py "$line".pilon.14.fasta -o "$line"_quast_output
mv "$line".pilon.14.fasta "$line"_complete_hybrid_assembly.fasta
