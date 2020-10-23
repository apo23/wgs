#!/usr/bin/bash 

#script meant to run in personal folder

# $1 all_genome.fasta
# $2 read 1
# $3 read 2

module load bowtie2
module load samtools
module load r
module load megahit
module load prodigal
module load blast

bowtie2-build ../databases/$1 $1

bowtie2 -p 5 -x $1 -1 ../fastq/$2 -2 ../fastq/$3 -S output.sam

samtools view -@ 5 -b output.sam > output.bam

samtools sort -@ 12 -l 1 output.bam > output_sorted.bam

samtools index -@ 12 output_sorted.bam 

samtools idxstats -@ 12 output_sorted.bam > output_sorted_counts

grep ">" ../databases/$1|cut -f 2 -d ">" >association.tsv


echo """data = read.table('output_sorted_counts', sep = '\t', header = F)
        data_names = read.csv('association.tsv', sep = ' ', header = F)
        data_names = data_names[-which(data_names[,3] == ''),]
        true_names = paste(data_names[,2],data_names[,3])
        pie(data[,3], true_names)""" > pie_plot.rscript

Rscript pie_plot.rscript


#scp aroubert@core.cluster.france-bioinformatique.fr:/shared/projects/uparis_m2_bi_2020/metagenomique_2/aroubert/Rplots.pdf .
# Commande pour copier Rplots.pdf au local pwd pour avoir le chemin vers actuel





megahit -1 ../fastq/$2 -2 ../fastq/$3 --k-list 21 -t 12 -o megahit_result

prodigal -d nuc_file.fasta -i megahit_result/final.contigs.fa > prodigal_output.gbk

sed "s:>:*\n>:g" nuc_file.fasta | sed -n "/partial=00/,/*/p"|grep -v "*" > genes_full.fna

blastn -db ../databases/resfinder.fna -query genes_full.fna -out blastn_results.out -outfmt "6 qseqid sseqid evalue pident qcovhsp" -perc_identity 0.8 -qcov_hsp_perc 0.8 -evalue 0.003 -num_threads 12
