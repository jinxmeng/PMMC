aenv gtdbtk
gtdbtk classify_wf --genome_dir ../05.dRep/total_genomes_cls_95_fa --out_dir classify_4045 -x fa --cpus 80 --pplacer_cpus 24 --skip_ani_screen

# 95聚类的结果进行物种注释
# 如果是unclassified，就找到这个代表基因组的簇，换别的基因组看看有没有物种注释结果
# 有的话合并进来
# 代表基因组并且有物种注释的基因组作为最终的基因组集合
# 株水平也按照代表基因组过滤

cut -f1 classify_4019_summary.tsv | perl -e '%h;while(<>){chomp;$h{$_}=1};open I, "../05.dRep/total_genomes.list";while(<I>){chomp;@s=split/\//;$s[-1]=~s/.fa//;$x=$s[-1];$x=~s/vae_/vae\./;print "$s[-1]\t$x\t$_\n" if exists $h{$s[-1]}}' > final_genomospecies_rename.list
grep PN final_genomospecies_rename.list | parallel -j 1 --colsep="\t" -q perl -e 'open I, "$ARGV[0]";open O, ">$ARGV[1]";while(<I>){chomp;if(/>.*(k\S+)/){print O ">$1\n"}else{print O "$_\n"}}' {3} final_genomospecies/{2}.fa
grep -v PN final_genomospecies_rename.list | parallel -j 3 --colsep="\t" -q perl -e 'open I, "$ARGV[0]";open O, ">$ARGV[1]";while(<I>){chomp;if(/>NODE_(\d+)/){print O ">NODE_$1\n"}else{print O "$_\n"}}' {3} final_genomospecies/{2}.fa
cat final_genomospecies.list | parallel -j 10 -q perl -e '$x=$ARGV[1];open I, "final_genomospecies/$ARGV[0]";open O, ">final_genomospecies2/$ARGV[0]";while(<I>){chomp;if(/>(\S+)/){print O ">$x\.$1\n"}else{print O "$_\n"}}' {/} {/.}
seqkit stat *fa > ../final_genomospecies.stat

cut -f1 classify_4019_summary.tsv | perl -e '%h;while(<>){chomp;$h{$_}=1}; open I, "../05.dRep/total_genomes_cls_95_res";while(<I>){chomp;@s=split/\t/;print "$_\n" if exists $h{$s[2]}}' > final_strains.cls
grep PN9RH.bin.188 classify_unclassified.list.cls >> final_strains.cls
cat final_strains_cls | perl -ne 'chomp;@s=split/\t/;@s2=split/,/,$s[3];for (@s2){print "$_\n"}' | perl -e '%h;while(<>){chomp;$h{$_}=1} open I, "../05.dRep/total_genomes_cls_99.list";while(<I>){chomp;@s=split/\//;$s[-1]=~s/.fa//;$x=$s[-1];$x=~s/vae_/vae./;print "$s[-1]\t$x\t$_\n" if exists $h{$s[-1]}}' > final_strains_rename.list
grep PN final_strains_rename.list | parallel -j 10 --colsep="\t" -q perl -e 'open I, "$ARGV[0]";open O, ">$ARGV[1]";while(<I>){chomp;if(/>.*(k\S+)/){print O ">$1\n"}else{print O "$_\n"}}' {3} final_strains/{2}.fa
grep -v PN final_strains_rename.list | parallel -j 3 --colsep="\t" -q perl -e 'open I, "$ARGV[0]";open O, ">$ARGV[1]";while(<I>){chomp;if(/>NODE_(\d+)/){print O ">NODE_$1\n"}else{print O "$_\n"}}' {3} final_strains/{2}.fa

# contrast to other catalogs
dRep compare cls_95 -g fa.list -d -pa 0.9 -sa 0.95 -nc 0.30 -cm larger -p 50 --S_algorithm fastANI
