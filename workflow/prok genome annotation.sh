# protein
cat ../06.taxa/final_strains.list | parallel -j 10 prodigal -q -f gff -p meta -a faa/{/.}.faa -d ffn/{/.}.ffn -o gff/{/.}.gff -i {}

# rRNA
cat ../06.taxa/final_strains.list | parallel -j 5 cmsearch -Z 30 --hmmonly --cut_ga --noali --cpu 2 --tblout rRNA/{/.}.tblout -o rRNA/{/.}.res /share/data1/database/Rfam/Rfam-v14.9/Rfam.rRNA.cm {} &
perl -e '@files=@ARGV; @class=("5S_rRNA", "SSU_rRNA_archaea", "SSU_rRNA_bacteria", "LSU_rRNA_archaea", "LSU_rRNA_bacteria"); print "name\t".join("\t", @class)."\n"; for (@files){open I, "$_"; $n=$_; %h=(); @l=(); while(<I>){chomp; next if /^#/; @s=split/\s+/; $h{$s[2]}++}; for (@class){if(exists $h{$_}){push @l, $h{$_};}else{push @l, 0;}}; print "$n\t".join("\t", @l)."\n"}' *tblout > ../rRNA.stat
perl -e '@files=@ARGV; @class=("5S_rRNA", "SSU_rRNA_archaea", "SSU_rRNA_bacteria", "LSU_rRNA_archaea", "LSU_rRNA_bacteria"); print "name\t5S\t16S\t23S\n"; for (@files){open I, "$_"; $n=$_; %h=(); @l=(); while(<I>){chomp; next if /^#/; @s=split/\s+/; $h{$s[2]}++}; for (@class){if(exists $h{$_}){push @l, $h{$_};}else{push @l, 0;}}; print "$n\t$l[0]\t".($l[1]+$l[2])."\t".($l[3]+$l[4])."\n"}' *tblout > ../rRNA.stat2

# tRNA
cat ../06.taxa/final_strains.list | parallel -j 5 tRNAscan-SE -B -L -q -o tRNA/{/.}.out {} &
perl -e '@files=@ARGV;print "name\ttRNA\n";for (@files){open I, "$_"; $n=$_; $count=0;while(<I>){chomp;next if /Sequence/ or /Name/ or /----/;@s=split/\s+/;$count+=1}; print "$n\t$count\n"}' *out > ../tRNA.stat 

csvtk join --left-join -f "1" -t gff.cds.stat rRNA.stat2 tRNA.stat > gene.anno.res

