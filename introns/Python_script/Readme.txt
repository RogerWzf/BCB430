sequencing.py is a python script used for getting putative intron sequences. 
    It reads in a blast result file (.xml) and a contigs file (.fasta) and return a fasta file.
verify.py is a python script uded for verify the putative intron sequences. 
    It reads in a BLAST result file (.xml) and a putative introns fasta file 
    and returns a fasta file contain verified introns.
find_introns_for_protein.py is a python script used to find introns sequences for 
    certain protein in the verified_intron files.
    It reads in a fasta file and returns a fasta file.