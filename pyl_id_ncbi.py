#!/usr/bin/env python3

import os
import sys
import argparse
from Bio import SeqIO
import numpy as np
import Bio.Data.CodonTable


    # set arguments
    description = "Identify putative pyl-containing proteins in code 15 annotated genomes"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-s", "--source", required=True,
        help="Directory containing GenBank files (.gbk)"
    )
    parser.add_argument(
        "-o", "--out", default='./output/',
        help="Output directory, default is './output/'"
    )
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
                    if feature.type == "CDS":
                        seq = feature.location.extract(record.seq)

                        tag_pos = None
                        for i in range(0, len(seq) - 2, 3):
                            codon = str(seq[i: i + 3]).upper()  # convert to string
                            if codon == 'TAG':
                                tag_pos = i // 3  # use integer division
                                break

                        if tag_pos is not None:
                            len_readthrough = len(seq) // 3
                            extension_length = len_readthrough - tag_pos
                            extension_lengths[genome].append(extension_length)

                            # Handle different possible qualifiers
                            gene = feature.qualifiers.get('locus_tag', [''])[0] \
                                or feature.qualifiers.get('gene', [''])[0] \
                                or 'unknown'

                            annotation = feature.qualifiers.get('product', [''])[0]

                            g.write(
                                f"{genome}\t{gene}\t{annotation}\t{len_readthrough}\t{tag_pos}\n"
                            )

    # apply to source
    for filename in os.listdir(source):
        if filename.endswith('.gbk'):
            genome = os.path.splitext(filename)[0]
            genome_path = os.path.join(source, filename)
            get_extension(genome, genome_path)

    # write to output
    with open(ext_csv, 'w') as f:
        f.write('genome,total,average_ext,median_ext,stdev,min,max\n')
        for genome, lengths in extension_lengths.items():
            if lengths:
                f.write(
                    f"{genome},{len(lengths)},{np.average(lengths):.2f},{np.median(lengths):.2f},"
                    f"{np.std(lengths):.2f},{min(lengths)},{max(lengths)}\n"
                )
            else:
                f.write(f"{genome},0,NA,NA,NA,NA,NA\n")


if __name__ == "__main__":
    main()
