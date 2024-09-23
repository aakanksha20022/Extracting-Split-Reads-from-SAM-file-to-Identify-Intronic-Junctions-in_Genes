Extracting-Split-Reads-from-SAM-file-to-Identify-Intronic-Junctions-in-Genes
---
This script extracts split reads from an alignment file to identify the locations of introns in genes. To run the script, two input files are required:
1. A .sam file, a special format to hold alignments
2. A tab-separated file. The file has a header and contains three columns. The first is a gene id. The second is a transcript id. The third is the location of the gene. The location is in the
   format TGME49_chrVIII:6,793,066..6,795,596(-). This string includes the name of the chromosome where the gene is encoded (TGME49_chrVIII), the start position (6793066), the end position
   (6795596) and the strand (-).

USAGE:
python3 myScript.py mySamFile.sam myInputTable.txt

The output is a tab-delimited file called 2875662.txt containing:
1. Gene ID
2. Junction start
3. Junction end
4. Number of reads supporting the junction
