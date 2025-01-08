#!/usr/bin/env python3

import os
import sys
import argparse
from Bio import SeqIO
import numpy as np
import Bio.Data.CodonTable 

def main():
    # set arguments
    description = "use to identify putative pyl-containing proteins in code 15 annotated genomes"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-s", "--source", required=True, help="directory code 15 prokka annotated genomes")
    parser.add_argument("-o", "--out", default='./output/', help="output directory, otherwise default 'output'")
    args = parser.parse_args()

    source = args.source
    out = args.out

    ext_dir = os.path.join(out, "genome_extension")
    ext_csv = os.path.join(out, "extension_lengths.csv")

    os.makedirs(ext_dir, exist_ok=True)

    # recode tag to pyl
    new_code = Bio.Data.CodonTable.generic_by_id[15]
    new_code.forward_table['TAG'] = 'O'

    extension_lengths = {}

    # get extension lengths
    def get_extension(genome, path):
        extension_lengths[genome] = []

        with open(os.path.join(ext_dir, genome + '_out.tsv'), 'w') as g:
            g.write('genome\tgene\tannotation\tlen_readthrough\tlen_tag\n')

            for record in SeqIO.parse(path, 'genbank'):
                for feature in record.features:
                	# only use coding sequences, otherwise will return tRNAs
                    if feature.type == "CDS":  
                        seq = feature.location.extract(record.seq)
                        tag_pos = None

                        for i in range(0, len(seq) - 2, 3):
                            codon = seq[i: i + 3]
                            if codon == 'TAG':
                                tag_pos = i / 3  
                                break

                        if tag_pos is not None:
                            len_readthrough = len(seq) // 3
                            extension_length = len_readthrough - tag_pos
                            extension_lengths[genome].append(extension_length)

                            gene = str(feature.qualifiers.get('locus_tag', [''])[0])[16:]  
                            annotation = feature.qualifiers.get('product', [''])[0]  
                            g.write(f"{genome}\t{gene}\t{annotation}\t{len_readthrough}\t{tag_pos}\n")

    # apply to source
    for item in os.listdir(source):
        item_path = os.path.join(source, item)
        
        if os.path.isdir(item_path):
            genome_path = os.path.join(item_path, item)
            if os.path.isfile(genome_path + '.gbk'):  
                get_extension(item, genome_path + '.gbk')

    # write to output
    with open(ext_csv, 'w') as f:
        f.write('genome,total,average_ext,median_ext,stdev,min,max\n')
        for genome, lengths in extension_lengths.items():
            f.write(f"{genome},{len(lengths)},{np.average(lengths)},{np.median(lengths)},{np.std(lengths)},{min(lengths)},{max(lengths)}\n")


if __name__ == "__main__":
    main()
