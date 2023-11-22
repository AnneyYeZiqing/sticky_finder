Before using Sticky finder or Sticky finder counter, please make sure to download fasta.py and expSS_score.py


Sticky finder:

A python program that searches for the complementary sequences with a minimum length of 6 nts (default) between identical or non-identical transcripts. 

Usage: modify the relevant parameters (i.e. RNA sequence of interest, RNA secondary structure in dot-bracket format, and name of the output files) inside sticky_finder_v3.py, then run sticky_finder_v3.py in a python shell (3.7 and above).

Output format: it produces two files. File1 in the script contains inverted repeats and palindromes. File2 has sense-antisense sequences. 



Sticky finder counter:

A python program that counts the number of complementary sequences with a minimum length of 6 nts (default) and a user-defined %GC content, with the fasta file as an input.

Usage: sticky_finder_v3_counter.py fasta_file_randonmized_seqs GC%threshold of the complementary sequences

Output format: transcript name,gene name,sequence length,total # of palindromes,total # of inverted repeats,total number of sense-antisense sequences