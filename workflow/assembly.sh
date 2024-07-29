cat ../sample_name | parallel -j 5 run_megahit.sh ../01.cleandata/{}_dehost_1.fq.gz,../01.cleandata/{}_dehost_2.fq.gz 21,41,61,81,101,121,141 {}
# ls *fa | parallel -j 20 so  {} \> assembly_so/{/.}.so
# parse_so.pl assembly_so/ assembly_so.merge.tsv
# perl -e 'open I, "find assembly_log/*log |"; open O, ">pig_contigs_info"; print O "sample\tcontigs_count\ttotal_len\tmin_len\tmax_len\tavg_len\tN50_len\n";while(<I>){chomp;$_=~/\/(.*).log/;$a=$1; open IN, "<$_"; $b=-1; while(<IN>){chomp;if(/final contigs/){$b=$.} if($. eq $b+1){$_=~/- (\d+) contigs, total (\d+) bp, min (\d+) bp, max (\d+) bp, avg (\d+) bp, N50 (\d+) bp/; print O "$a\t$1\t$2\t$3\t$4\t$5\t$6\n"} }}'
# 20231119
cat ../sample_name2 | while read i;do run_megahit.sh ../01.cleandata/${i}_rmhost_1.fq.gz,../01.cleandata/${i}_rmhost_2.fq.gz 21,41,61,81,101,121,141 ${i};done
cat ../sample_name2 | parallel -j 4 run_megahit.sh ../01.cleandata/rmhost_fq/{}_rmhost_1.fq.gz,../01.cleandata/rmhost_fq/{}_rmhost_2.fq.gz 21,41,61,81,101,121,141 {} &
# 20231129
cat ../sample_name3 | parallel -j 4 run_megahit.sh ../01.cleandata/rmhost_fq/{}_rmhost_1.fq.gz,../01.cleandata/rmhost_fq/{}_rmhost_2.fq.gz 21,41,61,81,101,121,141 {} &
# 20231130
run_megahit.sh ../01.cleandata/rmhost_fq/PN6RC_rmhost_1.fq.gz,../01.cleandata/rmhost_fq/PN6RC_rmhost_1.fq.gz 21,41,61,81,101,121,141 PN6RC
# 20231201
cat ../sample_name3 | parallel -j 4 run_megahit.sh ../01.cleandata/{}_rmhost_1.fq.gz,../01.cleandata/{}_rmhost_2.fq.gz 21,41,61,81,101,121,141 {} &
# 20231204
cat ../sample_name4 | parallel -j 4 run_megahit.sh ../01.cleandata/{}_rmhost_1.fq.gz,../01.cleandata/{}_rmhost_2.fq.gz 21,41,61,81,101,121,141 {} &
# 20231222
cut -f1 ../sample_list | parallel -j 4 run_megahit.sh ../01.cleandata/{}_rmhost_1.fq.gz,../01.cleandata/{}_rmhost_2.fq.gz 21,41,61,81,101,121,141 {} &

# statistics
grep 'pe,\|contigs, total' *log | perl -ne 'chomp;if(/contigs/){$_=~/- (\d+) contigs, total (\d+) bp, min (\d+) bp, max (\d+) bp, avg (\d+) bp, N50 (\d+) bp/;print "$1\t$2\t$3\t$4\t$5\t$6\n"}else{$_=~/(.*).log:.*: pe, (\d+) reads, (\d+) max/;print "$1\tpe\t$2\t$3\n"}' | awk '{if(NR%2==1){ORS="\t";print}else{ORS="\n";print}}' | csvtk add-header -t -n sample,seq_strategy,input_reads,reads_max_len,contigs,total_len,min_len,max_len,avg_len,N50_len > ../stat_contigs.tsv
# stat > 2000
ls *.fa | parallel -j 10 -q perl -e '%h;$x=$ARGV[0];while(<>){chomp;if(/>/){$_=~/len=(\d+)$/; $h{a}++ unless $1<2000;}} print "$x\t$h{a}\n"' {} >> ../stat_contigs_m2000.tsv

parallel -j 3 --colsep='\t' seqkit seq -g -m 1500 {2} -o contigs_m1500/{1}.m1500.fa :::: ../contigs_fa_list
parallel -j 3 --colsep='\t' seqkit seq -g -m 5000 {2} -o contigs_m5000/{1}.m5000.fa :::: ../contigs_fa_list
grep TP ../../contigs_fa_list | parallel -j 5 --colsep='\t' -q perl -e 'open I, "seqkit seq --quiet -g -m 500 $ARGV[1] |";open O, ">$ARGV[0].m500.fa";while(<I>){chomp;if(/>(\S+)/){print O ">$ARGV[0]_$1\n"}else{print O "$_\n"}} close O' {1} {2}
ls *fa | parallel -j 10 samtools faidx {}

# contigs_m1500
ls *fa | while read i;do echo "seqkit replace -p \" .*\" -r \"\" $i -o ${i}2";done | parallel -j 2 &
sourmash sketch dna *fa 2> sourmash.out.log
sourmash compare --csv contigs_m1500.sourmash.csv -o contigs_m1500.sourmash.mat contigs_m1500/*sig -p 16 -q
sourmash plot contigs_m1500.sourmash.mat --pdf

ls *fa | parallel -j 2 seqkit stat {} >> ../stat.contigs.m500.tsv2
le stat.contigs.m500.tsv | sed 's/,//g;1d' | perl -e '$x;while(<>){chomp;@s=split/\s+/;$x+=$s[3]}; print "$x\n"'
