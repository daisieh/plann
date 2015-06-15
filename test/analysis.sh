#!/bin/bash

SAMPLES=*.gb

# extract the fasta sequences:
for GB in $SAMPLES
do
	echo "making $GB.fasta"
	perl $REPOS/phylogenomics/converting/gb_to_fasta.pl -input $GB
done;

for GB in $SAMPLES
do
	for FASTA in $SAMPLES
	do
		echo "comparing $FASTA.fasta to $GB"
		perl $REPOS/plann/plann.pl -ref $GB -fasta $FASTA.fasta -out $GB.$FASTA -organism "Populus balsamifera" -sample "bals"
	done;
done;

head -n 2 *.results.txt > all_results.txt

echo "Running tbl2asn"
tbl2asn -Vbv -p . -t plann.sbt

# library(ggplot2)
# data <- read.table("~/Documents/DRAFTS/plann/results.tsv")
# ggplot(data) + geom_point(shape=3,color="black",aes(x=1:5,y=P_balsamifera)) + geom_point(shape=3,color="red",aes(x=1:5,y=P_trichocarpa)) + geom_point(shape=3,color="green",aes(x=1:5,y=P_fremontii)) + geom_point(shape=3,color="orange",aes(x=1:5,y=Salix_interior)) + geom_point(shape=3,color="blue",aes(x=1:5,y=Jatropha_curcas)) +  scale_x_discrete(labels=row.names(data)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
