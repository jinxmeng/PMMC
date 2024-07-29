# 对基因组进行BGCs预测
cat ../06.taxa/final_strains.list | parallel -j 30 antismash --taxon bacteria -c 4 --genefinding-tool prodigal --output-dir BGCs.out/{/.} {}
cat ../06.taxa/add_genome/fa.filepath | parallel -j 30 --colsep="\t" antismash --taxon bacteria -c 4 --genefinding-tool prodigal --output-dir BGCs.out/{1} {2}
\ls -d * | parallel -j 80 parse_antismash.py {}/index.html {}.out
cut -f1 ../../06.taxa/add_genome/fa.filepath | parallel -j 80 parse_antismash.py {}/index.html {}.out
ls *out | perl -ne 'chomp;$x=$_;$x=~s/.out//;open I, "<$_";while(<I>){chomp;if(!/contigs/){print "$x\t$_\n"}}' > ../BGCs.out.stat

# 如果BGCs所在的contigs小于5K，则过滤；和MIBIG比较
ls */*region*gbk | perl -ne 'chomp;$x=$_;$x=~s/\//_/;print "cp $_ ../BGCs.gbk/$x\n"' | sh 
cat BGCs.gbk.list | xargs grep '^LOCUS' |  perl -lne '$_=~/BGCs.gbk\/(\S+?)\.(region\d+.gbk).*\s+(\d+) bp/;print "$1\t$1\.$2\t$3"' > BGCs.gbk.len
bigscape -i BGCs.gbk.len_filter -o bigscape.out --pfam_dir /share/data1/database/pfam/release36 -c 100 --mibig

# To prevent sampling biases in quantitative analysis (taxonomic and functional compositions of GCCs/GCFs, GCF and GCC distances to reference databases as well as GCF 
# metagenomic abundances), the 73,864 BGCs were further dereplicated by retaining only the longest BGC per GCF per species, resulting in a total of 30,182 BGCs.
# 确定新颖性
bigslice --query BGCs.gbk.m5k.rep --query_name jinxmeng_run0006 --n_ranks 3 /share/data1/database/bigslice/full_run_result/ --run_id 0006
extract_data_from_bigscape.py -i /share/data1/database/bigslice/full_run_result/reports/1/data.db -tn bgc -o bigslice.out.bgc
extract_data_from_bigscape.py -i /share/data1/database/bigslice/full_run_result/reports/1/data.db -tn gcf_membership -o bigslice.out.gcf_membership
csvtk join --left-join -t -f "2;1" bigslice.out.gcf_membership bigslice.out.bgc > bigslice.out.metadata
le bigslice.out.metadata | perl -e '%h;while(<>){chomp;next if /gcf_id/;@s=split/\t/;$h{$s[1]}{$s[2]}=$_;} for (keys %h){@x=sort keys %{$h{$_}};print "$h{$_}{$x[-1]}\n"}' | awk '$3>900' > bigslice.out.metadata.d_m900

# 聚类分析
# cat BGCs.gbk.m5k.rep.list | perl -lne '$_=~/(\S+?)\.re/;print "mkdir -p bigslice.in/dataset_1/$1"' | sort -u | parallel -j 30
# cat BGCs.gbk.m5k.rep.list | perl -lne '$_=~/(\S+?)\.re/;print "cp BGCs.gbk.m5k/$_\.gbk bigslice.in/dataset_1/$1/"' | parallel -j 30
# for i in *.gbk; do sed -i 's/Version      :: False/Version      :: 6.1.1/' $i; done # 修改MIBiG数据库的格式，才能被bigslice读入，进行聚类分析
# bigslice -i bigslice.in bigslice.out # 聚类分析
# bash bigslice.out/start_server.sh 33001

# 定量
cat BGCs.gbk.m5k.rep.filepath | parallel -j 60 extract_biosynthetic_cds.py {} {.}.fa &
cat BGCs.gbk.m5k.rep/*fa | seqkit replace -p ".*\/" -r "" > BGCs.gbk.m5k.rep.biosyn.fa
minimap2 -d BGCs.gbk.m5k.rep.biosyn.mmi BGCs.gbk.m5k.rep.biosyn.fa
cat ../clean_fq_list | parallel -j 8 --colsep="\t" run_minimap2.sh {2} {3} BGCs.gbk.m5k.rep.biosyn.mmi profile/{1} &
cat ../17.public_wp/CNP0000824.fq.list | parallel -j 8 --colsep="\t" run_minimap2.sh {2} {3} BGCs.gbk.m5k.rep.biosyn.mmi profile/{1} &
le clean_fq.list | parallel -j 5 --colsep="\t" samtools coverage -d 0 -o profile/{1}.cvg profile/{1}.sort.bam
combine_file_zy_folder_allsample.py -D profile/ -suffix .cvg -o profile.rc -t 1 -n 1 -v 4
combine_file_zy_folder_allsample.py -D profile/ -suffix .cvg -o profile.cvg -t 1 -n 1 -v 6
profile_filter_by_cvg.py profile.rc profile.cvg 50 profile.rc2
profile_filter.py profile.rc2 profile.rc3

# seqkit fx2tab -n -l BGCs.gbk.m5k.rep.biosyn.fa -o BGCs.gbk.m5k.rep.biosyn.len
# profile_rc2tpm.pl profile.rc3 BGCs.gbk.m5k.rep.biosyn.len profile.rc3.tpm

# perl -ne 'chomp;print "$_" if /name/;$_=~/\|(.*)$/;print "$1\n"' PULs.filter.rc > PULs.filter.rc2
profile_summary_by_group.py PULs.rc.filter PULs.summary.rc

