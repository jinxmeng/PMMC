# 20231214 metabat2 single coverage binning
cat file_path_list | parallel -j 5 --colsep="\t" flow_metabat2.sh {2},{3} {4} {1} metabat2_single
cat ../../sample_name | parallel -j 10 sh contigs2bin.sh bin_fa/{}* \> contigs2bin/{}.tsv &

# metabinner
# cut -f1 ../sample_name | parallel -j 2 /share/data1/mjx/bin/flow_metabin.sh /share/data1/mjx/proj/04.black_pig_metagenome_20230529/02.assembly/contigs_m1500_fna/{}.contigs_m1500.fna metabinner/depth_res/{}.depth metabinner_res/{}

# vamb mcvg
# 根据距离拆分成14个类分别进行分箱
sourmash-->dist-->hclust-->cutree
for i in `seq 1 14`;do awk '$2=="'"$i"'"' contigs_m1500.sourmash.hclust.cutree.tsv | cut -d. -f1 | grep -w -f - ../../file_path_list > cluster$i.list;done
ls -d cluster* | parallel -j 10 mkdir {/.}.contigs
for i in `seq 1 14`;do cut -f1-2 cluster$i.list | awk '{print $0"\tcluster"'"$i"'".contigs"}';done > contigs_rename.list
# 改名字
cat contigs_rename.list | parallel -j 2 --colsep="\t" -q perl -e '$s="|";open I, "$ARGV[0]"; open O, ">$ARGV[1]/$ARGV[2].fa"; while(<I>){chomp;if(/>(\S+)/){print O ">$ARGV[2]$s$1\n"}else{print O "$_\n"}}' {2} {3} {1}
# concatenate.py
for i in `seq 1 14`;do concatenate.py cluster$i.m1500.fa cluster$i.contigs/*fa -m 1500 --keepnames --nozip;done
# mkdir
ls -d cluster*contigs | parallel -j 10 mkdir {/.}.sort_bam
# minimap2 index
for i in `seq 1 14`; do minimap2 -d cluster$i.m1500.mmi cluster$i.m1500.fa -I100g;done
# run_minimap2
for i in `seq 1 14`;do cat cluster$i.list | parallel -j 10 --colsep="\t" echo run_minimap2.sh {3} {4} cluster$i.m1500.mmi cluster$i.sort_bam/{1};done > r1.cmd.sh
cat r1.cmd.sh | parallel -j 3 {}
# vamb
# for i in `seq 1 14`;do echo "vamb -o \"|\" --outdir cluster$i.vamb.out --fasta cluster$i.m1500.fa --bamfiles cluster$i.sort_bam/*.sort.bam --seed 2024 -p 56 --minfasta 200000" ;done
vamb -o "|" --outdir cluster1.vamb.out --fasta cluster1.m1500.fa --bamfiles cluster1.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster2.vamb.out --fasta cluster2.m1500.fa --bamfiles cluster2.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster3.vamb.out --fasta cluster3.m1500.fa --bamfiles cluster3.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster4.vamb.out --fasta cluster4.m1500.fa --bamfiles cluster4.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster5.vamb.out --fasta cluster5.m1500.fa --bamfiles cluster5.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster6.vamb.out --fasta cluster6.m1500.fa --bamfiles cluster6.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster7.vamb.out --fasta cluster7.m1500.fa --bamfiles cluster7.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster8.vamb.out --fasta cluster8.m1500.fa --bamfiles cluster8.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster9.vamb.out --fasta cluster9.m1500.fa --bamfiles cluster9.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster10.vamb.out --fasta cluster10.m1500.fa --bamfiles cluster10.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster11.vamb.out --fasta cluster11.m1500.fa --bamfiles cluster11.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster12.vamb.out --fasta cluster12.m1500.fa --bamfiles cluster12.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster13.vamb.out --fasta cluster13.m1500.fa --bamfiles cluster13.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000
vamb -o "|" --outdir cluster14.vamb.out --fasta cluster14.m1500.fa --bamfiles cluster14.sort_bam/*.sort.bam --seed 2024 -p 112 --minfasta 200000

# vamb single
# concatenate.py PN10RA.m1500.fa /share/data1/mjx/proj/04.black_pig_metagenome_20230529/02.assembly/contigs/PN10RA.fa -m 1500
# ls -d * | parallel -j 3 sh ../../contigs2bin.sh {}/* \> {}.tsv

# 对每个样本两种分析结果进行简单的去冗余
cat sample_name | parallel -j 3 dRep dereplicate dRep.out/{} -p 40 -g merge_fa.list/{}.tsv --ignoreGenomeQuality --S_algorithm fastANI -pa 0.9 -sa 0.99 -nc 0.3 -cm larger --skip_plots  &
ls -d * | while read i;do ls /share/data1/mjx/proj/04.black_pig_metagenome_20230529/03.binning/dRep/dRep.out/$i/dereplicated_genomes/*fa;done | parallel -j 30 ln -s {} ../bins/
ls -d * | while read i;do ls /share/data1/mjx/proj/04.black_pig_metagenome_20230529/03.binning/dRep/dRep.out/$i/dereplicated_genomes/*fna;done | perl -ne 'chomp;$_=~/out\/(\S+?)\/de\S+\/(\S+).fna/;print "ln -s $_ ../bins/$1.$2.fa\n"' | parallel -j 30 {} &
ls | grep -v -f ../sample_name | while read i;do cat $i;done | perl -ne 'chomp;$_=~/bins\/(\S+)\/(\S+).fna/;print "ln -s $_ ../bins/$1.$2.fa\n"' | sh
