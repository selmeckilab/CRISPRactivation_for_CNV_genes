#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --job-name="assembleguidesMiseq"
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=%x_%u_%A_%a.out
#SBATCH --error=%x_%u_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pvzande@umn.edu
#SBATCH -p msismall
#SBATCH -A selmecki
#SBATCH --array=0-2

cd /projects/standard/selmecki/pvzande/CRISPRa/Selmecki_Project_002

module load pear/0.9.11
module load vsearch/2.3.4

files=(`cat "R1filenames.txt"`)
file=${files[${SLURM_ARRAY_TASK_ID}]}
name=$(echo ${file} | awk '{sub("_R1_001.fastq.gz", ""); print}')

fastq_path=/projects/standard/selmecki/pvzande/CRISPRa/Selmecki_Project_002/${name}
fastq_R1=${fastq_path}_R1_001.fastq.gz
fastq_R2=${fastq_path}_R2_001.fastq.gz
pear_out=./pear_output/${name}
        
pear -f ${fastq_R1} -r ${fastq_R2} -o ${pear_out} -j 16 -u 0.01 -q 20 -m 350 -n 100

assembled_reads=${pear_out}.assembled.fastq
filter_out=./vsearch_trim/${name}_trimmed.fastq
filter_rc=./vsearch_trim/${name}_trimmed_rc.fastq

#vsearch --fastx_filter ${assembled_reads} --fastq_stripleft 134 --fastaout ${filter_out} #This version does not have the stripright command yet. 
vsearch --fastx_revcomp ${assembled_reads} --fastaout ${filter_rc}
vsearch --fastx_filter ${filter_rc} --fastq_stripleft 80 --fastaout ${filter_out} 

aggregate_fasta=./vsearch_aggregate/${name}.aggregate.fasta
vsearch --derep_fulllength ${filter_out} --relabel seq --output ${aggregate_fasta} --sizeout

#Counting - note the Fguides.txt file has been modified to include 8 nucleotide overhangs on each side.
#It ended up being easier to just count from the filtered fastq file rather than the aggregated fasta, so I am just using that. 
module load parallel

export file_tocount=${filter_out}

FASTAOUTPUT=counts2/${name}counts.csv

# Define the function to search for barcode occurrences
count_barcode_occurrences() {
  barcode="$1"
  fastq_files=$(echo ${file_tocount})
  count=$(grep -c "$barcode" "$fastq_files")
  echo "$barcode,$count"
}

export -f count_barcode_occurrences

# Use parallel processing to search for barcode occurrences
echo "Barcode,Count" > "$FASTAOUTPUT"
cat Fguides.txt | parallel count_barcode_occurrences {} >> "$FASTAOUTPUT"

