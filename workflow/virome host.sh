# cat ../06.taxa/final_genomospecies.list | parallel -j 30 minced -minNR 2 {} genomospecies.csp/{/.}.csp genomospecies.csp/{/.}.gff
cut -f2 ../06.taxa/total_genomes.filepath | parallel -j 30 minced -minNR 2 {} genomostrain.csp/{/.}.csp genomostrain.csp/{/.}.gff 

# ls genomospecies.csp/*csp | parallel -j 20 -q perl -e '$x=$ARGV[0];open I, "$ARGV[1]"; open O, ">$ARGV[2]"; $n=1; %h; while(<>){chomp;next unless /^\S+\s+\S+\s+(\S+)\s+\[/;$s=$1;unless (exists $h{$s}){print O ">$x|crisp_$n\n$s\n";$n++; $h{$s}++}}' {/.} {} {}.fa
ls genomostrain.csp/*csp | parallel -j 20 -q perl -e '$x=$ARGV[0];open I, "$ARGV[1]"; open O, ">$ARGV[2]"; $n=1; %h; while(<>){chomp;next unless /^\S+\s+\S+\s+(\S+)\s+\[/;$s=$1;unless (exists $h{$s}){print O ">$x|crisp_$n\n$s\n";$n++; $h{$s}++}}' {/.} {} {}.fa
# cat genomospecies.csp/*fa > genomospecies.csp.fa
cat genomostrain.csp/*fa > genomostrain.csp.fa

# makeblastdb -in genomospecies.csp.fa -dbtype nucl -out genomospecies.csp
makeblastdb -in genomostrain.csp.fa -dbtype nucl -out genomostrain.csp

# crisp
blastn -query ../20.vir_cluster/final_vOTUs.fa -db genomostrain.csp -evalue 1e-2 -out final_vOTUs.csp.btn -outfmt 6 -num_alignments 99999 -num_threads 80 -word_size 8
filter_blast -i final_vOTUs.csp.btn -o final_vOTUs.csp.btn.f --evalue 1e-5 --score 45 --tops 20
le final_vOTUs.csp.btn | perl -e 'while(<>){chomp;@s=split /\s+/;($a,$b)=($s[0],$s[1]);$b=~s/\|crisp.*$//; push @{$g{$a}},$b unless exists $h{"$a $b"}; $h{"$a $b"}++;} for(sort keys %g){@a=@{$g{$_}};print "$_\t".(join ",",@a)."\n";}' > final_vOTUs.csp.list

# genome map
cut -f1 final_vOTUs.csp.list | perl -e '%h;while(<>){chomp;$h{$_}=1} open I, "../20.vir_cluster/final_vOTUs.name";while(<I>){chomp;print "$_\n" unless exists $h{$_}}'
le final_vOTUs.name | perl -lne 'if(/(\S+\|\S+_\d+?)_fragment.*/){print "$_\t$1"}elsif(/(\S+\|\S+_\d+?)_\d$/){print "$_\t$1"}else{print "$_\t$_"}'

# blast strains
blastn -query ../20.vir_cluster/final_vOTUs.fa -db /share/data1/mjx/proj/04.black_pig_metagenome_20230529/06.taxa/final_strains -evalue 1e-2 -out final_vOTUs.genome.btn -outfmt 6 -num_alignments 999999 -num_threads 112
connect_blast final_vOTUs.genome.btn final_vOTUs.genome.btn.conn 1
le ../20.vir_cluster/final_vOTUs.len | sort -rnk2 > final_vOTUs.len.sort
filter_blast -i final_vOTUs.genome.btn.conn -o final_vOTUs.genome.btn.conn.f --evalue 1e-10 --qfile final_vOTUs.len.sort --qper 30 --identity 90
le final_vOTUs.genome.btn.conn.f | perl -e 'while(<>){chomp;@s=split /\s+/;($a,$b)=($s[0],$s[1]);$b=~s/[\|:].*$//; push @{$g{$a}},$b unless exists $h{"$a $b"}; $h{"$a $b"}++;} for(sort keys %g){@a=@{$g{$_}};print "$_\t".(join ",",@a)."\n";}' > final_vOTUs.genome.list
cat final_vOTUs.csp.list final_vOTUs.genome.list | sort | perl -ne 'chomp;@a = split /[,\t]/; foreach $x (@a){print "$a[0]\t$x\n"}' | awk '$1!=$2' | sort -u > final_vOTUs.host
csvtk join --left-join -H -t -f "2;1" final_vOTUs.host ../06.taxa/final_strains.taxa2 > final_vOTUs.host.taxa

le final_vOTUs.host.taxa | cut -f1,9 | perl -e '%h;while(<>){chomp;@s=split/\t/;$h{$s[0]}{$s[1]}++};for $i(keys %h){$x=scalar keys $h{$i};for $j(keys $h{$i}){if($x>1){print "$i\t$j\t$h{$i}{$j}\tmulti-genera\n"}else{print "$i\t$j\t$h{$i}{$j}\tsingleton\n" } } }' > final_vOTUs.host.taxa2
le final_vOTUs.host.taxa | cut -f1,8 | perl -e '%h;while(<>){chomp;@s=split/\t/;$h{$s[0]}{$s[1]}++};for $i(keys %h){$x=scalar keys $h{$i};for $j(keys $h{$i}){if($x>1){print "$i\t$j\t$h{$i}{$j}\tmulti-families\n"}else{print "$i\t$j\t$h{$i}{$j}\tsingleton\n" } } }' > final_vOTUs.host.taxa3
