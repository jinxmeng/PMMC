makeblastdb -in virus.fa -dbtype nucl -out virus -hash_index
blastn -query virus.fa -db virus -outfmt '6 std qlen slen' -max_target_seqs 999999 -out virus.btn -num_threads 80
anicalc.py -i virus.btn -o virus.ani
aniclust.py --fna virus.fa --ani virus.ani --out virus.cls --min_ani 95 --min_tcov 85 --min_qcov 0
cut -f1 virus.cls | seqkit grep -f - virus.fa -o virus.cls.fa &

VIBRANT_run.py -virome -i virus.cls.fa -f nucl -folder VIBRANT.out -t 80 -l 5000 -no_plot
perl -e '%h;open I, "$ARGV[0]";while(<I>){chomp;$h{$_}=1} open I, "$ARGV[1]";$x=1;while(<I>){chomp;@s=split/\t/;if(/contig_id/ & $x==1){print "$_\n";$x=0}; if(exists $h{$s[0]}){print "$_\n"}}' vOTUs.name ../19.phages/virus.ckv.qs.tsv > vOTUs.ckv.qs.part01
perl -e '%h;open I, "vOTUs.ckv.qs.part01";while(<I>){chomp;@s=split/\t/;$h{$s[0]}=1} open I, "vOTUs.name";while(<I>){chomp;if(! exists $h{$_}){print "$_\n"}}' > vOTUs.name.for.ckv
seqkit grep -f vOTUs.name.for.ckv vOTUs.fa > vOTUs.name.for.ckv.fa
checkv end_to_end -t 80 vOTUs.name.for.ckv.fa vOTUs.name.for.ckv.out --quiet
awk '$1!="contig_id"' vOTUs.ckv.qs.part02 | cat vOTUs.ckv.qs.part01 - > vOTUs.ckv.qs

awk '$1=="contig_id" || $8=="Complete" || $8=="High-quality" || $8=="Medium-quality"' vOTUs.ckv.qs > final_vOTUs.qs
tail -n+2 final_vOTUs.qs | cut -f1 > final_vOTUs.name
seqkit grep -f final_vOTUs.name vOTUs.fa > final_vOTUs.fa
le final_vOTUs.name | perl -e '%h;while(<>){chomp;$h{$_}=1} open I, "vOTUs.lysogenic.lytic";while(<I>){chomp;@s=split/\t/;print "$s[0]\t$s[1]\n" if exists $h{$s[0]}}' > final_vOTUs.lysogenic.lytic
flow_prodigal.sh final_vOTUs.fa final_vOTUs 112

# compare_to_PVD
cat /share/data1/database/catalog_bacteriophages/PVD/final_PVD.vOTUs.rename.fa ../final_vOTUs.fa > total_vOTUs.fa
makeblastdb -in total_vOTUs.fa -dbtype nucl -out total_vOTUs -hash_index
blastn -query total_vOTUs.fa -db total_vOTUs -outfmt '6 std qlen slen' -max_target_seqs 999999 -out total_vOTUs.btn -num_threads 80
anicalc.py -i total_vOTUs.btn -o total_vOTUs.ani
aniclust.py --fna total_vOTUs.fa --ani total_vOTUs.ani --out total_vOTUs.cls --min_ani 95 --min_tcov 85 --min_qcov 0

# PVD-specific vOTU: 45919 
# Bpig-specific vOTU: 8625
# shared vOTU: 2252

# phylogenetic
aenv ViPTreeGen
ViPTreeGen --ncpus 80 ../final_vOTUs.fa phylo.out &
# source activate /share/data2/guorc/Software/conda/r4.1
# /share/data2/guorc/Software/ViPTreeGen-master/ViPTreeGen --ncpus 80 final_vOTUs.fa phylo.out

le final_vOTUs.rename.tsv | perl -e '%h;while(<>){chomp;@s=split/\t/;$h{$s[0]}=$s[1]}; open I, "final_vOTUs.fa";while(<I>){chomp;if(/>(\S+)/){$x=$h{$1};print ">$x\n"}else{print "$_\n"}}' > final_vOTUs.fna

