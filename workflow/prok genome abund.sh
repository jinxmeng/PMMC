# abundance
le ../06.taxa/final_genomospecies.list | perl -e 'while(<>){chomp;open I, "$_"; $_=~/cies\/(\S+).fa/;$x=$1;while(<I>){chomp;if(/>(\S+)/){print ">$x\.$1\n"}else{print "$_\n"}}}' > genomospecies.fa &
grep '>' genomospecies.fa | sed 's/>//' | perl -ne 'chomp;if(/(\S+)\.(\S+)/){print "$1.$2\t$1.fa\n"}'  > genomospecies.stb &

minimap2 -d genomospecies.mmi genomospecies.fa -I20g
cat ../clean_fq_list | parallel -j 5 --colsep="\t" run_short_mapping.sh -i {2},{3} -x genomospecies.mmi -o sort_bam/{1} -b -t minimap2 -p 24
cat ../17.public_wp/CNP0000824.fq.list | parallel -j 5 --colsep="\t" run_short_mapping.sh -i {2},{3} -x genomospecies.mmi -o sort_bam/{1} -b -t minimap2 -p 24

cat ../sample_name | parallel -j 5 run_coverm.sh sort_bam/{}.sort.bam genomospecies.list cvm.out/{}
cat ../17.public_wp/CNP0000824.sample_name | parallel -j 5 run_coverm.sh sort_bam/{}.sort.bam genomospecies.list cvm.out/{}
cut -f1-2 cvm.out/PN10RA.cvm | sed '1,2d' > genomospecies.len
combine_file_zy_folder_allsample.py -D cvm.out -suffix .cvm -n 1 -v 3 --skip 2 -o genomospecies.rc
combine_file_zy_folder_allsample.py -D cvm.out -suffix .cvm -n 1 -v 4 --skip 2 -o genomospecies.cvg
combine_file_zy_folder_allsample.py -D cvm.out -suffix .cvm -n 1 -v 6 --skip 2 -o genomospecies.tpm
combine_file_zy_folder_allsample.py -D cvm.out -suffix .cvm -n 1 -v 8 --skip 2 -o genomospecies.ab

# stat
ls *cvm | parallel -j 5 -q perl -e '$n=0;open I, "$ARGV[0]";while(<I>){chomp;next if (/Genome/ or /unmapped/); @s=split/\t/;$n+=$s[2];} print "$ARGV[1]\t$n\n"' {} {/.} > ../map_rate/stat.mapped.reads.counts.tsv
le ../../02.assembly/stat/stat_contigs.tsv | cut -f1,3 > stat.total.reads.counts.tsv
le stat.mapped.reads.counts.tsv| grep CNR | cut -f1 | grep -w -f - ../../17.public_wp/CNP0000824_fq_info | cut -f1,3 >> stat.total.reads.counts.tsv
csvtk join -t -f 1 stat* | sed '1d' | perl -lne '@s=split/\t/;$n=$s[1]/$s[2];print "$_\t$n"' > stat.mapped.rate.tsv


# inStrain profile sort_bam/PN1RA.sort.bam genomospecies.fa -o inStrain.out/PN1RA.IS -p 32 -s genomospecies.stb
# cat ../sample_name | parallel -j 2 samtools index sort_bam/{}.sort.bam -@ 32 -b -o sort_bam/{}.sort.bam.bai
# cat ../sample_name | parallel -j 4 inStrain profile sort_bam/{}.sort.bam genomospecies.fa -o inStrain.out/{}.IS -p 32 -s genomospecies.stb
