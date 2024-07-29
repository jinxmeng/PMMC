#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-05-24, 12:20:10
# modified date: 2024-05-24, 12:20:10

import sys
from Bio import SeqIO

if len(sys.argv) == 1:
    sys.exit("Usage: extract_fasta.py [in_f *.gbk] [out_f *fa]")

def extract_biosynthetic_cds(gbk, out_f):
    name = gbk.split("/")[-1].replace(".gbk", "")
    out_f = open(out_f, "w")
    for record in SeqIO.parse(gbk, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "gene_kind" in feature.qualifiers and feature.qualifiers["gene_kind"][0] == "biosynthetic":
                    cds = feature.extract(record.seq)
                    gene_kind = "biosynthetic"
                    location = feature.location
                    out_f.write(f">{name}|{record.id}|{gene_kind}|{location}\n{cds}\n")

def main(gbk, out_f):
    extract_biosynthetic_cds(gbk, out_f)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

