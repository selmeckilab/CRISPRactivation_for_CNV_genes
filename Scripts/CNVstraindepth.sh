#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --job-name="CNVstrains"
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pvzande@umn.edu
#SBATCH -p agsmall,aglarge
#SBATCH -A selmecki
#SBATCH --array=

#Change to correct directory
cd /projects/standard/selmecki/pvzande/CNVstrains
samples=(`cat "CNVstrains.txt"`) #A list of CNV strain names to plot depth. Adjust array numbers for number of strains
strain=${samples[${SLURM_ARRAY_TASK_ID}]} #designate a sample for each arrayed job

#File Name Entry - adjust file paths and names as necessary. 
read1paths=(`cat "CNVstrainsR1paths.txt"`)
read1=${read1paths[${SLURM_ARRAY_TASK_ID}]}
read2paths=(`cat "CNVstrainsR2paths.txt"`)
read2=${read2paths[${SLURM_ARRAY_TASK_ID}]}

#Reference fasta - path to reference genome for mapping
reference_fasta="/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta"


#Load modules
module load fastqc
module load trimmomatic/0.39
module load bwa/0.7.17
module use /home/selmecki/shared/software/modulefiles.local
module load samtools
module load gatk


#Check sequencing QC for all files enclosed in the directory with FastQC
fastqc ./*.fastq.gz

# Trim using trimmomatic
java -jar $TRIMMOMATIC/trimmomatic.jar PE -threads 8 -phred33 -trimlog ./${strain}.trimlog ${read1} ${read2} ${strain}_trimpairedr1.fastq.gz ${strain}_trimunpairedr1.fastq.gz ${strain}_trimpairedr2.fastq.gz ${strain}_trimunpairedr2.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

# Align to reference fasta using BWA MEM
bwa mem -t 128 -R "@RG\tID:${strain}\tPL:ILLUMINA\tPM:HiSeq\tSM:${strain}" ${reference_fasta} ./${strain}_trimpairedr1.fastq.gz ./${strain}_trimpairedr2.fastq.gz > ${strain}_trimmed_bwa.sam

# Samtools sort, index, remove duplicates, and reindex, and generate depth file
samtools view -bS -o ${strain}_trimmed_bwa.bam ${strain}_trimmed_bwa.sam
samtools flagstat ${strain}_trimmed_bwa.bam > ${strain}_trimmed_bwa.stdout
samtools sort ${strain}_trimmed_bwa.bam -o ${strain}_trimmed_bwa_sorted.bam
samtools index ${strain}_trimmed_bwa_sorted.bam
samtools rmdup ${strain}_trimmed_bwa_sorted.bam ${strain}_trimmed_bwa_sorted_rmdup.bam
samtools index ${strain}_trimmed_bwa_sorted_rmdup.bam
samtools depth -aa ${strain}_trimmed_bwa_sorted_rmdup.bam > ${strain}_depth.txt

