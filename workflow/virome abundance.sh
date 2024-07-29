minimap2 -d final_vOTUs.mmi final_vOTUs.fa -I20g
cat clean_fq.list | parallel -j 5 --colsep="\t" run_short_mapping.sh -i {2},{3} -x final_vOTUs.mmi -o sort.bam/{1} -b -t minimap2 -p 24

cat ../sample_name | parallel -j 5 samtools coverage sort.bam/{}.sort.bam -o sort.bam/{}.cvg
le clean_fq.list | cut -f1 | grep CNR | parallel -j 5 samtools coverage sort.bam/{}.sort.bam -o sort.bam/{}.cvg

combine_file_zy_folder_allsample.py -D sort.bam/ -suffix .cvg -o final_vOTUs.rc -t 1 -n 1 -v 4
combine_file_zy_folder_allsample.py -D sort.bam/ -suffix .cvg -o final_vOTUs.cvg -t 1 -n 1 -v 6
profile_filter_by_cvg.py final_vOTUs.rc final_vOTUs.cvg 50 final_vOTUs.rc2 &
profile_rc2tpm.pl final_vOTUs.rc2 final_vOTUs.len final_vOTUs.tpm

# map rate
ls *cvg | parallel -j 5 -q perl -e '$n=0;open I, "$ARGV[0]";while(<I>){chomp;next if (/rname/); @s=split/\t/;$n+=$s[3];} print "$ARGV[1]\t$n\n"' {} {/.} > ../map_rate/stat.mapped.reads.counts.tsv
csvtk join -t -f 1 stat.mapped.reads.counts.tsv stat.total.reads.counts.tsv | perl -lne '@s=split/\t/;$n=$s[1]/$s[2];print "$_\t$n"' > stat.mapped.rate.tsv
