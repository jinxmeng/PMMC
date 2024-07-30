#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-03-12, 23:14:49
# modified date: 2024-03-12, 23:14:49

'''
Jinxin Meng, 20240312
Generate the TPM for each KO by summaring the gene abundance.
Required files:
1. KO annotation with field gene name and feature name with "\t" or "\s" delimited
2. gene rc and len provided by samtools coverage [rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq]
'''

import sys

if len(sys.argv) != 4:
    sys.exit("Usage: calcu_tpm_for_KO.py [in_f <gene ID|feature ID>] [cvg file generated from samtools coverage] [out_f]")

def main(in_f, cvg, out_f):
    KO = {}
    with open(in_f, "r") as f:
        for i in f:
            l = i.strip().split()
            KO[l[0]] = l[1]

    total = 0
    with open(cvg, "r") as f:
        f.readline()
        for i in f:
            l = i.strip().split()
            if int(l[3]) == 0:
                continue
            # if float(l[5]) < 10:
            #    continue
            total += int(l[3])/int(l[2])
    
    res = {}
    with open(cvg, "r") as f:
        f.readline()
        for i in f:
            l = i.strip().split()
            if int(l[3]) == 0:
                continue
            # if float(l[5]) < 10:
            #     continue
            tpm = (int(l[3])/int(l[2])/total)*1000000
            if l[0] not in KO:
                continue
            if KO[l[0]] not in res:
                res[KO[l[0]]] = 0
            res[KO[l[0]]] += tpm
    out_f = open(out_f, "w")
    for k, v in res.items():
        out_f.write(k + "\t" + str(round(v, 4)) + "\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])

