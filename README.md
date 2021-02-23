# HPVChimera_Gb
HPVChimera_Gb

This script aims to check putative human papillomavirus (HPV) nucleotide sequences for chimera presence.
Requirements for the script: a) Fasta file with contigs´ sequences b) HPVdatabase c) Nucleotide-Nucleotide BLAST 2.10.1+
You will find all files needed for the database (updated November 17th, 2020), however, you can create it yourself with the following command:
makeblastdb -in HPV_COMPLETE_GENOMES.fasta -parse_seqids -blastdb_version 5  -title "HPVcomplete" -dbtype nucl

This script requires three parameters as input.
    1.	Fasta file with your contigs´s sequences (Full Path)
    2.	Workspace folder (Full Path)
    3.	HPVdatabase (HPV_COMPLETE_GENOMES.fasta) (Full Path)
Example

./HPVChimera_Gb.sh /PATH/contigssequences.fasta /PATH/$workspacefolder /PATH/HPVdatabase.fasta

After being executed, you will see the message "Check result in: /$workspacefolder/chimera/result.txt"

Summary of steps:

The putative HPV nucleotide sequences are compared to a database of known HPV sequences with BLAST. The database (HPVdatabase) comprises HPV sequences from all oficially established HPV types present in the International HPV Reference Center database (hpvcenter.se) together with all sequences from non-oficially established HPV types found in the PaVe database (https://pave.niaid.nih.gov/). Steps that are performed:
1.	Check if the contig´s sequence shows at least 85% sequence identity to any of the HPV types present in the database. We decided to set the cut-off at 85% identity, to give a slight margin for the 90% homology within the HPV L1 gene required for 2 sequences to be considered the same HPV type.
2.	Check if the contig´s sequence shows at least 60% aligning coverage to the top hit HPV type.
3.	Contigs´s sequences are divided in 3 equal segments. Chimeras are reported if a) at least one of the segments doesn´t show at least > 85% similarity to the top HPV hit aligned for the total sequence or b) any of the 3 segments shows an alignment coverage <70%.


Example of output: (HPVchimera_results.txt)

    Name of 1st contig  NO_chimera  HPV85 
    Name of 2nd contig  Chimera Contig sequence does not show any similarity to any HPV type given.
    Name of 3rd contig  Chimera The sequence aligns to 2 different HPV types with same bitscore.
    Name of 4th contig  Chimera Contig sequence does not show >85% sequence identity to any HPV type given.
    Name of 5th contig  Chimera At least one contig segment sequence does not show any similarity to any HPV type given.
    Name of 6th contig  Chimera At least one sequence segment aligns to 2 different HPV types with same bitscore.
    Name of 7th contig  Chimera At least one contig segment does not show >85% sequence identity to any HPV type given.
    Name of 8th contig  Chimera The top hit alignment does not cover >70% of contigs segment sequence.
