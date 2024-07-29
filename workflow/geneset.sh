# gene calling
cat ../sample_name | parallel -j 20 flow_prodigal.sh /share/data1/mjx/proj/04.black_pig_metagenome_20230529/02.assembly/contigs_m500/{}.m500.fa prodigal/{}.m500 5 >/dev/null 2>/dev/null
ls *fa | parallel -j 10 prodigal -f gff -p single -q -i {/.}.fa -a {/.}.faa -d {/.}.ffn -o {/.}.gff

# merge primary geneset
ls /share/data1/mjx/proj/04.black_pig_metagenome_20230529/10.geneset/prodigal/*ffn > gene.ffn.list
perl -e 'while(<>){chomp; open I, "seqkit seq -g -m 100 $_ |"; while(<I>){chomp;if(/^(\S+)/){print "$1\n"}else{print "$_\n"}}}' gene.ffn.list > gene.ffn
mmseqs easy-cluster gene.ffn cluster . --cluster-mode 2 --cov-mode 1 --min-seq-id 0.95 -c 0.9 --kmer-per-seq-scale 0.8 --threads 112
ls /share/data1/mjx/proj/04.black_pig_metagenome_20230529/10.geneset/prodigal/*faa > gene.faa.list
perl -e 'while(<>){chomp; open I, "$_"; while(<I>){chomp;if(/^(\S+)/){print "$1\n"}else{print "$_\n"}}}' gene.faa.list > gene.faa
cp mmseqs/cluster_rep_seq.fasta geneset.ffn
# seqkit replace -p "\s" -r "" geneset.ffn -o geneset.ffn2
seqkit fx2tab -n -l geneset.ffn -o geneset.len
cut -f1 geneset.len > geneset.id
seqkit grep -f geneset.id mmseqs/gene.faa -o geneset.faa

# taxa
seqkit split --by-part 4 --by-part-prefix geneset. geneset.faa
for i in {1..4};do diamond blastp -d /share/data1/database/ncbi/nr_20240207/nr.dmnd --outfmt 6 --evalue 0.00001 --max-target-seqs 10 -p 100 -q ../geneset.faa.split/geneset.00$i.faa -o nr.00$i.btp >/dev/null 2>&1 ;done
perl -ne 'chomp;@s=split/\t/;if($s[0] ne $a){print "$s[0]\t$s[1]\n";$a=$s[0]}' nr.*.btp > nr.combind.btp
sed 's/\..*$//' nr.combind.btp > nr.tsv
perl -lne '$_=~/(\d+)$/;print "$_\t$1"' /share/data1/database/ncbi/nr_20240207/prot.accession2taxid.FULL.chunk.list | parallel -j 50 --colsep="\t" -q perl -e '%h;open I, "$ARGV[0]";while(<I>){chomp;@s=split/\t/;$h{$s[0]}=$s[1]} open I, "$ARGV[1]"; open O, ">$ARGV[2]"; while(<I>){chomp;@s=split/\t/;print O "$_\t$h{$s[1]}\n" if exists $h{$s[1]}}' {1} nr.tsv nr.taxid.chunk_{2}
cat nr.taxid.chunk_* > nr.taxid
perl -e '%h;open I, "/share/data1/database/ncbi/taxonomy_20231110/nr.taxid.table";while(<I>){chomp;@s=split/\t/;$h{$s[0]}=$s[1]} open I, "nr.taxid"; while(<I>){chomp;@s=split/\t/;if(exists $h{$s[2]}){print "$_\t$h{$s[2]}\n"}else{print "$_\tUnknown\n"}}' > nr.tax
perl -lne '$_=~/supk__(.*?);/;print "$1"' nr.tax | awk '$1!=""' | sort --parallel=20 | uniq -c | awk '{print $2"\t"$1}' > nr.stat.supk.tsv

# abundance
minimap2 -d geneset.mmi geneset.ffn -I200g &
cat ../clean_fq_list | parallel -j 4 --colsep="\t" run_short_mapping.sh -i {2},{3} -x geneset.mmi -o sort_bam/{1} -b -t minimap2 -p 32 
# cat ../17.public_wp/CNP0000824.fq.list | parallel -j 4 --colsep="\t" run_short_mapping.sh -i {2},{3} -x geneset.mmi -o sort_bam/{1} -b -t minimap2 -p 32 &
cat ../sample_name | parallel -j 5 samtools coverage abundance/{}.sort.bam -o abundance/{}.cvg
combine_file_zy_folder_allsample.py -D abundance -suffix .cvg -o geneset.rc -t 1 -n 1 -v 4
combine_file_zy_folder_allsample.py -D abundance -suffix .cvg -o geneset.cvg -t 1 -n 1 -v 6
profile_filter_by_cvg.py dbcan4.cgc.profile.rc dbcan4.cgc.profile.cvg 50 dbcan4.cgc.profile.rc2 &

# rarefaction
convert_cluster2matrix.pl ../mmseqs/cluster_cluster.tsv gene.profile
perl -ne 'chomp;$_=~/(\S+)\s(\S+?)_/;print "$2\t$1\n"' ../mmseqs/cluster_cluster.tsv > gene.cluster.tsv
perl -ne 'chomp;@s=split/\s/;open O, ">>gene.cluster.tsv.split/$s[0].gene.list";print O "$s[1]\n"' gene.cluster.tsv
le gene.cluster.tsv.split.filepath| perl -lne '$_=~/lit\/(\S+).gene/;print "$1\t$_"' | csvtk join -t -H --left-join -f 1 - ../../sample_group  | cut -f2,9 > gene.cluster.tsv.split.by_breed.filepath
run_subsample.py gene.cluster.tsv.split.filepath gene.cluster.subsample.tsv
run_subsample_by_group.py gene.cluster.tsv.split.by_breed.filepath gene.cluster.tsv.split.by_breed.subsample.tsv &
run_subsample_by_group.py gene.cluster.tsv.split.by_region.filepath gene.cluster.tsv.split.by_region.subsample.tsv &

ls *cvg | sed 's/.cvg//g'|parallel -j 10 awk \'\{if\(\$6\>50 \&\& \$6\!="coverage"\)\{print \$1\}\}\' {}.cvg \> ../rare2/gene.cluster.tsv.split2/{}.gene.list
ls *cvg | parallel -j 10 -q perl -e 'open I, "$ARGV[0].cvg";open O, ">../rare2/gene.cluster.tsv.split/$ARGV[0].gene.list"; while(<I>){chomp;@s=split/\t/;print O "$s[0]\n" if ($s[5]>50)}' {/.} &

# contrast
perl -lane 'if(/>(\S+)/){print ">ChenC_2021\|$1"}else{print "$_"}' /share/data1/database/catalog_geneset/pig_PIGC_ChenCongying_2021_NC/PIGC90_17237052_cds.fa > ChenC_2021.fa
perl -lane 'if(/>(\S+)/){print ">HuJ_2024\|$1"}else{print "$_"}' /share/data1/database/catalog_geneset/pig_HuJun_2024_ISME/GeneCatalog.fa > HuJ_2024.fa
perl -lane 'if(/>(\S+)/){print ">XiaoL_2016\|$1"}else{print "$_"}' /share/data1/database/catalog_geneset/pig_LiangXiao_2016_NM/287sample_7.7M.GeneSet.fa > XiaoL_2016.fa
perl -lane 'if(/>(\S+)/){print ">LiM_2024\|$1"}else{print "$_"}' ../geneset.ffn > LiM_2024.fa
mmseqs easy-cluster LiM_2024.fa XiaoL_2016.fa HuJ_2024.fa ChenC_2021.fa cluster . --cluster-mode 2 --cov-mode 1 -v 0 --min-seq-id 0.95 -c 0.9 --kmer-per-seq-scale 0.8 --threads 112
run_cluster2tab.py cluster_cluster.tsv > cluster.tsv
/share/data1/software/miniconda3/envs/rbase_4.3/bin/Rscript mjx.R
cut --complement -f 1 cluster.tsv | tail -n+2 | sort --parallel=80 | uniq -c > cluster.pattern.tsv
cat cluster.pattern.tsv | perl -lpe 's/^\s//g;s/\s+/\t/g' | csvtk add-header -t -n count,LiM_2024,XiaoL_2016,ChenC_2021,HuJ_2024 | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' > cluster.pattern.stat.tsv
cut -f1 cluster.tsv | tail -n+2 | cut -d "|" -f1 | sort --parallel=20 | uniq -c > cluster.gene.rep.tsv

# map rate
cat ../clean_fq_list | parallel -j 4 --colsep="\t" run_short_mapping.sh -i {2},{3} -x geneset.mmi -o sort_bam/{1} -b -t minimap2 -p 3
cat ../public_wp_CNP0000824.sample_name | parallel -j 2 samtools flagstat {}.sort.bam \> {}.flagstat 

# function
cat *gene.list | sort --parallel=20 -u > known.gene.list
perl -e 'open I, "known.gene.list";%h;while(<I>){chomp;$h{$_}=1} open I, "../geneset.id";while(<I>){chomp;print "$_\n" unless exists $h{$_}}' > unknown.gene.list
le nr/nr.taxid | perl -e '%h;while(<>){chomp;@s=split/\t/;$h{$s[0]}=1} open I, "function/unknown.gene.list"; while(<I>){chomp;print "$_\n" unless exists $h{$_}}' > geneset.without.any.annotations.list &

# gene count
ls *gff | parallel -j 5 -q perl -e 'open I, "$ARGV[0]";$x=0;while(<I>){chomp; $x+=1 if /Prodigal_v2.6.1\sCDS/;} print "$ARGV[0]\t$x\n"' {} > ../stat.prodigal.gene.count.tsv 
ls *gff | parallel -j 5 -q perl -e 'open I, "$ARGV[0]";$x=0;while(<I>){chomp; $x+=1 if /partial=00/;} print "$ARGV[0]\t$x\n"' {} > ../stat.prodigal.gene.count.complete.tsv
ls *gff | parallel -j 3 -q perl -e '$x=$ARGV[0];open I, "$x";$len=0; $n=0; while(<I>){chomp;next if /^#/; @s=split/\s/;$len+=$s[4]-$s[3]+1; $n++};$avg=$len/$n;print "$x\t$len\t$n\t$avg\n"' {} >> ../stat.prodigal.gene.len.tsv &
