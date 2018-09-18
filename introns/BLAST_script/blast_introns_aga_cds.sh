#!bin/bash

cd /media/Storage4/Zhifan/Introns

/media/Storage4/mattpreston/Misc\ Programs/CreateBlastDB/CreateBlastDB –r Anolis_carolinensis.AnoCar2.0.cds.all.fa -o AnoliscdsDB 

blastn –db AnoliscdsDB -query result.fasta -out introns_aga_cds.xml -outfmt 5 –num_threads 20 –max_target_seqs 5 –evalue 0.00001
