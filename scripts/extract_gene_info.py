#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-05-26, 15:30:58
# modified date: 2024-05-26, 15:30:58

import sys, re
from Bio import SeqIO

if len(sys.argv) != 3:
    sys.exit("Usage: extract_gene_info.py [*gbk] [*fa] [*gff] [out_f]")

def parse_gbk(gbk):
    x = {}
    x["gene"] = {}
    x["genome"] = re.findall(r"(\S+?)_", gbk)[0]
    for record in SeqIO.parse(gbk, "genbank"):
        x["contigs"] = record.id
        for feature in record.features:
            if feature.type == "CDS" and "gene_kind" in feature.qualifiers:
                tag = feature.qualifiers["locus_tag"][0]
                x["gene"][tag] = [re.findall(r"_(\d+)", tag)[0], feature.qualifiers["gene_kind"][0]]
    return x
            
def main(gbk):
    gene_info = parse_gbk(gbk)


if __name__ == "__main__":
    main()


