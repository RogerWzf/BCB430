#!bin/bash

cd /media/Storage4/Zhifan/Introns

/media/Storage4/mattpreston/Misc\ Programs/CreateBlastDB/CreateBlastDB –r RRHintrons.fasta -o RRHintronsDB 

blastn –db RRHintronsDB -query RRH_putative_introns.fasta -out RRHintrons_aga_putative.xml -outfmt 5 –num_threads 20 –max_target_seqs 5 –evalue 0.00001

blastn –db RRHintronsDB -query RRH_verify_introns.fasta -out RRHintrons_aga_verify.xml -outfmt 5 –num_threads 20 –max_target_seqs 5 –evalue 0.00001