# flow_run_dbcan4.sh gene.faa run_dbcan4
diamond blastp -d /share/data1/database/run_dbcan4/db_20240122/CAZy.dmnd --outfmt 6 --min-score 60 --query-cover 50 --max-target-seqs 5 -p 32 -q geneset.faa -o CAZyme.btp > /dev/null 2>/dev/null
parse_CAZy_blast_for_geneset.pl CAZyme.btp CAZyme

# strain genomes cgc pul substrate prediction
cat ../06.taxa/final_strains.list | parallel -j 10 run_dbcan {} meta --tools diamond -c cluster --cgc_substrate --out_dir dbcan4.out/{/.} --dia_cpu 10 --hmm_cpu 10 --tf_cpu 10 --stp_cpu 10 --db_dir /share/data1/database/run_dbcan4/db_20240122/
cut -f2 ../06.taxa/add_genome/fa.filepath | parallel -j 8 run_dbcan {} meta --tools diamond -c cluster --cgc_substrate --out_dir dbcan4.out/{/.} --dia_cpu 10 --hmm_cpu 10 --tf_cpu 10 --stp_cpu 10 --db_dir /share/data1/database/run_dbcan4/db_20240122/

# CAZyme
le total_genomes.name | parallel -j 30 cut -f1,5 dbcan4.out/{}/overview.txt \| sed \'1d\' \>  dbcan4.CAZyme.count/{}.CAZyme.out &
ls * | parallel -j 20 -q perl -e 'open I, "$ARGV[0]";open O, ">$ARGV[1]";%h;while(<I>){chomp;@s=split/\t/;@s2=split/\+/, $s[1];$len=$#s2+1; for(@s2){$_=~s/_.*//g;$h{$_}+=1/$len}} for (keys %h){print O "$_\t$h{$_}\n"}' {} {}2 
combine_file_zy_folder_allsample.py -D dbcan4.CAZyme.count -o dbcan4.CAZyme.count.profile -suffix .CAZyme.out2


wc -l dbcan4.stat/*cgc.list | awk '{x=$2;sub(".cgc.list", "", $x);sub("dbcan4.stat/", "", $x);print $2"\t"$1}' > dbcan4.cgc.tsv
# 获取 cgc 列表和序列相应的位置
cat total_genomes.name | parallel -j 80 -q perl -e '%start;%end; open I, "dbcan4.out/$ARGV[0]/cgc.out"; while(<I>){chomp;@s=split/\t/;next if /\+\+\+\+\+/; $cgc="$s[5]\|$s[4]"; $start{$cgc}=$s[6] if ! exists $start{$cgc}; $end{$cgc}=$s[7] if $end{$cgc} < $s[7]} open O, ">dbcan4.stat/$ARGV[0].cgc.list"; for (keys %start){$x=$_;$x=~s/\|.*//;print O "$ARGV[0]\t$_\t$x\t$start{$_}\t$end{$_}\n"}' {} &
# 生成提取序列的id
ls ../total_genomes.name | parallel --plus -j 60 -q perl -e 'open I, "$ARGV[0]"; open O, ">$ARGV[1]"; while(<I>){chomp;@s=split/\t/;print O "$s[2]:$s[3]-$s[4]\n"}' {} {/.}.extract.info
# 提取序列 
cat ../total_genomes.filepath | parallel -j 10 --colsep="\t" samtools faidx {2} -r {1}.cgc.extract.info -o {1}.cgc.extract.fa
cat dbcan4.stat/*.cgc.extract.fa > dbcan4.cgc.extract.fa
# mmseqs
mmseqs easy-cluster dbcan4.cgc.extract.fa dbcan4.cgc.extract.mmseqs /share/data1/mjx/tmp/ --cluster-mode 2 --cov-mode 1 --min-seq-id 0.95 -c 0.9 --kmer-per-seq-scale 0.8 --threads 112
# abundance
minimap2 -d dbcan4.cgc.extract.mmseqs_rep_seq.mmi dbcan4.cgc.extract.mmseqs_rep_seq.fasta
le clean_fq.list | parallel -j 5 --colsep="\t" run_minimap2.sh {2} {3} dbcan4.cgc.extract.mmseqs_rep_seq.mmi dbcan4.cgc.profile/{1}
cut -f1 clean_fq.list | parallel -j 6 samtools coverage -d 0 dbcan4.cgc.profile/{}.sort.bam -o dbcan4.cgc.profile/{}.cvg
combine_file_zy_folder_allsample.py -D dbcan4.cgc.profile/ -suffix .cvg -o dbcan4.cgc.profile.rc -t 1 -n 1 -v 4
combine_file_zy_folder_allsample.py -D dbcan4.cgc.profile/ -suffix .cvg -o dbcan4.cgc.profile.cvg -t 1 -n 1 -v 6
profile_filter_by_cvg.py dbcan4.cgc.profile.rc dbcan4.cgc.profile.cvg 50 dbcan4.cgc.profile.rc2 &


# puls
wc -l dbcan4.stat/*PULs.extract.info | awk '{x=$2;sub(".PULs.extract.info", "", $x);sub("dbcan4.stat/", "", $x);print $2"\t"$1}' > dbcan4.puls.count.stat &
cat ../total_genomes.filepath | parallel -j 80 -q perl -e '%h; open I, "dbcan4.out/$ARGV[0]/substrate.out";while(<I>){chomp;@s=split/\t/;$h{$s[1]}++ unless /^#/}; open O, ">dbcan4.stat/$ARGV[0].PULs.stat";for (sort keys %h){print O "$_\t$h{$_}\n"}' {}
cat dbcan4.stat/*substrate.out | grep -v PULID | cut -f2| sort | uniq -c | sort -rnk1 > dbcan4.puls.count
# PULs 序列列表和PULs的对应关系
cat ../total_genomes.name | parallel -j 20 cp ../dbcan4.out/{}/substrate.out {}.substrate.out
cat ../total_genomes.name | parallel -j 30 -q perl -e '%cgc;open I, "$ARGV[0]"; while(<I>){chomp;@s=split/\t/;$cgc{$s[1]}="$s[2]:$s[3]-$s[4]"} open I, "$ARGV[1]"; open O, ">$ARGV[2]";while(<I>){chomp;next if /^#/;@s=split/\t/; print O "$s[1]\t$cgc{$s[0]}\n"}' {}.cgc.list {}.substrate.out {}.PULs.extract.info
# 提取序列
le ../total_genomes.name | parallel -j 20 -q perl -e '%pul;open I, "$ARGV[0]"; while(<I>){chomp;@s=split/\t/;$pul{$s[1]}=$s[0];} open I, "$ARGV[1]"; open O, ">$ARGV[2]";while(<I>){chomp;if(/>(\S+)/){if(exists $pul{$1}){$x=1; print O ">$1|$pul{$1}\n"}else{$x=0}} if($x==1 & !/>/){print O "$_\n"}}' {}.PULs.extract.info {}.cgc.extract.fa {}.PULs.extract.fa
cat ../total_genomes.name | parallel -j 20 seqkit replace -p \"\^\" -r \"{}_\" -o {}.PULs.extract2.fa {}.PULs.extract.fa
cat dbcan4.stat/*PULs.extract2.fa > dbcan4.PULs.extract.fa
# mmseqs
mmseqs easy-cluster dbcan4.PULs.extract.fa dbcan4.PULs.extract.mmseqs /share/data1/mjx/tmp/ --cluster-mode 2 --cov-mode 1 --min-seq-id 0.95 -c 0.9 --kmer-per-seq-scale 0.8 --threads 112
# abundance
minimap2 -d dbcan4.PULs.extract.mmseqs_rep_seq.mmi dbcan4.PULs.extract.mmseqs_rep_seq.fasta
le clean_fq.list | parallel -j 5 --colsep="\t" run_minimap2.sh {2} {3} dbcan4.PULs.extract.mmseqs_rep_seq.mmi dbcan4.PULs.profile/{1}
cut -f1 clean_fq.list | parallel -j 6 samtools coverage -d 0 dbcan4.PULs.profile/{}.sort.bam -o dbcan4.PULs.profile/{}.cvg
combine_file_zy_folder_allsample.py -D dbcan4.PULs.profile/ -suffix .cvg -o dbcan4.PULs.profile.rc -t 1 -n 1 -v 4
combine_file_zy_folder_allsample.py -D dbcan4.PULs.profile/ -suffix .cvg -o dbcan4.PULs.profile.cvg -t 1 -n 1 -v 6
profile_filter_by_cvg.py dbcan4.PULs.profile.rc dbcan4.PULs.profile.cvg 50 dbcan4.PULs.profile.rc2 &
profile_filter_by_cvg.py dbcan4.PULs.profile.rc dbcan4.PULs.profile.cvg 80 dbcan4.PULs.profile.rc3 &
cat total_genomes.name | parallel -j 50 cat dbcan4.stat/{}.PULs.extract.info \| cut -f1 \| sort -u \| xargs echo {} | sed 's/ /\t/' > dbcan4.puls.uniq.tsv
le dbcan4.puls.uniq.tsv | perl -e '%h;open I, "final_strains.taxa2";while(<I>){chomp;@s=split/\t/;$h{$s[0]}=$s[6]}; %t; %n; while(<>){chomp;@s=split/\t/;@s2=split/ /, $s[1]; for(@s2){$t{$h{$s[0]}}{$_}++; $n{$_}++};} @head=sort keys %t; @n=sort keys %n; print "name\t".join("\t", @head)."\n"; for $i (@n){@line=();push @line, $i; for $j (@head) { if (exists $t{$j}{$i}){push @line, $t{$j}{$i}}else{push @line, 0;}} print join("\t", @line)."\n"}' > dbcan4.puls.uniq.tsv2


# substrate 统计
cat ../06.taxa/final_strains.name | parallel -j 80 -q perl -e '%h; open I, "dbcan4.out/$ARGV[0]/substrate.out";while(<I>){chomp;@s=split/\t/;$h{$s[2]}++ unless /^#/}; open O, ">dbcan4.stat/$ARGV[0].substrates.stat";for (sort keys %h){print O "$_\t$h{$_}\n"}' {}
cat ../06.taxa/add_genome/fa.name | parallel -j 20 -q perl -e '%h; open I, "dbcan4.out/$ARGV[0]/substrate.out";while(<I>){chomp;@s=split/\t/;$h{$s[2]}++ unless /^#/}; open O, ">dbcan4.stat/$ARGV[0].substrates.stat";for (sort keys %h){print O "$_\t$h{$_}\n"}' {}
# cat total_genomes.name | while read i;do cat dbcan4.stat/$i.substrate.out | sed '1d' | cut -f3 | sort -u | wc -l | xargs echo -e "$i\t";done > dbcan4.substrate.class.count &
cat total_genomes.name | while read i;do cat dbcan4.stat/$i.substrate.out | sed '1d' | cut -f3 | sort -u | perl -e '@a;while(<>){chomp;push @a, $_;}print "'$i'\t".join(",", @a)."\n"';done > dbcan4.substrate.class.count2
le dbcan4.substrate.class.count2 | perl -e '%h;open I, "final_strains.taxa2";while(<I>){chomp;@s=split/\t/;$h{$s[0]}=$s[6]};%t; while(<>){chomp;@s=split/\t/;@s2=split/,/, $s[1]; for(@s2){$t{$h{$s[0]}}{$_}++}} for (keys %t){@n=keys $t{$_};$c=scalar @n;print "$_\t$c\t".join(",",@n)."\n"}' > dbcan4.substrate.class.count3
cat dbcan4.stat/*substrate.out  | grep -v PULID | cut -f3 | sort -u | wc -l
le dbcan4.substrate.class.count2 | perl -e '%h;open I, "final_strains.taxa2";while(<I>){chomp;@s=split/\t/;$h{$s[0]}=$s[6]}; %t; %n; while(<>){chomp;@s=split/\t/;@s2=split/,/, $s[1]; for(@s2){$t{$h{$s[0]}}{$_}++; $n{$_}++};} @head=sort keys %t; @n=sort keys %n; print "name\t".join("\t", @head)."\n"; for $i (@n){@line=();push @line, $i; for $j (@head) { if (exists $t{$j}{$i}){push @line, $t{$j}{$i}}else{push @line, 0;}} print join("\t", @line)."\n"}' > dbcan4.substrate.class.count4
cat dbcan4.out/*/substrate.out | cut -f2,3 | awk '$1!="PULID"' | sort --parallel=20 -u > dbcan4.PULs.substrate.tsv


# checkm coverage . -x fasta PN1RA.ckm.cvg PULs_sort_bam/PN1RA.sort.bam -t 16
# cat ../sample_name | while read i;do samtools coverage PULs_sort_bam/$i.sort.bam > PULs_sam_cvg/$i.cvg;done
# combine_file_zy_folder_allsample.py -D PULs_sam_cvg/ -suffix .cvg -o PULs.rc -t 1 -n 1 -v 4
# combine_file_zy_folder_allsample.py -D PULs_sam_cvg/ -suffix .cvg -o PULs.cvg -t 1 -n 1 -v 6 
# profile_filter_by_cvg.py PULs.rc PULs.cvg 50 PULs.rc.filter
# perl -ne 'chomp;print "$_" if /name/;$_=~/\|(.*)$/;print "$1\n"' PULs.filter.rc > PULs.filter.rc2
# profile_summary_by_group.py PULs.rc.filter PULs.summary.rc
# seqkit fx2tab -n -l dbcan4_out.PULs.extract.fasta | awk '{print $1"\t1\t"$2}' > dbcan4_out.PULs.extract.bed
# cat ../sample_name | while read i;do samtools index PULs_sort_bam/$i.sort.bam && bedtools coverage -a dbcan4_out.PULs.extract.bed -b PULs_sort_bam/$i.sort.bam > PULs_bed_cvg/$i.cvg;done
# cat ../sample_name | while read i;do echo "samtools index PULs_sort_bam/$i.sort.bam && bedtools coverage -a dbcan4_out.PULs.extract.bed -b PULs_sort_bam/$i.sort.bam > PULs_bed_cvg/$i.cvg";done > r2.sh
# 合并数据
# combine_file_zy_folder_allsample.py -D dbcan4_out/ -suffix .substrates.stat -s 0 -o dbcan4_out.substrates.stat
# combine_file_zy_folder_allsample.py -D dbcan4_out/ -suffix .PULs.stat -o dbcan4_out.PULs.stat
# wc -l *cgc.list | awk '{if($2!="total"){gsub(/.cgc.list/,"");print $2"\t"$1}}' | sort -rnk2 > ../dbcan4_out.cgc.stat

