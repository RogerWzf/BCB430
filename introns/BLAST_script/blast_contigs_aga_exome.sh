#!bin/bash

cd /media/Storage4/Zhifan/Introns

/media/Storage4/mattpreston/Misc\ Programs/CreateBlastDB/CreateBlastDB –r Anolis_RSO6ex3_5_DeNovo_All_Trinity.fasta -o outputFile 

blastn –db outputFile -query Anolis_RSO6ex3_5_Trinity_NR.fasta -out contigs_against_exome.xml -outfmt 5 –num_threads 20 –max_target_seqs 5 –evalue 0.00001
