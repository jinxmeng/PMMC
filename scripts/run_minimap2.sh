#!/usr/bin/bash

( [ $# -ne 4 ] ) &&  { echo -e "Usage: $0 [fq1] [fq2] [*.mmi] [out_prefix]" && exit 2; }

( [ -e $4.sort.bam ] ) && { echo -e "Skip sample: ${4##*/} .." && exit 0; }

minimap2 --MD -t 32 -ax sr $3 $1 $2 2>>$4.log |\
    perl -e 'while(<>){if(/^@/){print "$_"; next};chomp;@l=split/\t/;;next if ($l[1] & 0x4) != 0;$ref_length=0;$match_counts=0;$cov_length=0;while($l[5]=~/(\d+)[M=XID]/g){$cov_length+=$1};while($l[5]=~/(\d+)[MDN=X]/g){$ref_length+=$1};foreach $k(@l[11..$#l]){if($k=~/MD:Z:/){while($k=~/(\d+)/g){$match_counts+=$1}}};$identity=$match_counts/$cov_length*100;print "$_\n" if $identity>=95}' |\
    samtools view -@ 32 -bS | samtools sort -@ 32 -o $4.sort.bam -
