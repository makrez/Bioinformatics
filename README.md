# Bioinformatics Scripts

Collection of scripts used for bioinformatics analyses.

## extract_fasta_from_genbank

Python program that let's the user easily extract fasta sequences based on features.

Requires python3 and the following packages: [biopython](https://www.biopython.org) and [fire](https://github.com/google/python-fire).

Example usage:

`extract_fasta.py --file A_lyrata.gbk --type gene --gene rbcL`
