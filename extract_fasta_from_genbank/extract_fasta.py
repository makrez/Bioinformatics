#!/usr/bin/python3
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import fire

def extract_fasta(file,type,gene):
    """
    Program for extracting fasta sequences from genbank files.
    Arguments:
     - file: file in genbank format
     - type: Feature type (e.g. rRNA, CDS, gene, tRNA, ...)
     - gene: Name of the gene / rRNA, CDS etc. (e.g. rbcL). If set to 'all', all
       entries of the particular type will be returned.
    Usage: extract_fasta.py --file test.gbk --type gene --gene rbcL
    """

    # Check if feature type exists
    parser = SeqIO.parse(f'{file}', 'genbank')
    unique_feature_types = []
    for record in parser:
        for feature in record.features:
            unique_feature_types.append(feature.type)
    features_unique = set(unique_feature_types)

    try:
        assert f'{type}' in features_unique
    except:
        raise ValueError("Feature type not in genbank file. "+
        "Change the --type flag. " +
        "Available feature types in this file are: " +
        ', '.join(features_unique))

    # Check if gene exists
    genes=[]
    for record in SeqIO.parse(f'{file}', 'genbank'):
        for feature in record.features:
            if feature.type == 'gene':
                genes.append(feature.qualifiers.get('gene'))
    genes_flat = set([item for sublist in genes for item in sublist])

    try:
        if f'{gene}' != 'all':
            assert f'{gene}' in genes_flat
    except:
        raise ValueError("Gene not found in genbank file."+
        " Change the --gene flag. "+
        "Availabe genes in this genbank file are: " +
        ', '.join(genes_flat))


    for record in SeqIO.parse(f'{file}', 'genbank'):

        for feature in record.features:

            if feature.type == f'{type}':
                if f'{gene}' == 'all':
                    print(">" + record.id + "_"  +
                            ''.join(dict(feature.qualifiers)[f'gene']) +
                            "_" + "_Start_" + str(feature.location.start) )
                    print(feature.location.extract(record).seq)

                elif ''.join(dict(feature.qualifiers)['gene']) == f'{gene}':
                    print(">" + record.id + "_"  +
                            ''.join(dict(feature.qualifiers)[f'gene']) +
                            "_" + str(feature.location.start) )
                    print(feature.location.extract(record).seq)

if __name__ == '__main__':
    fire.Fire(extract_fasta)
