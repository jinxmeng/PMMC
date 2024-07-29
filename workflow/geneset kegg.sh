diamond blastp -d /share/data1/database/KEGG/KEGG20230401.dmnd --outfmt 6 --min-score 60 --query-cover 50 --max-target-seqs 10 -p 112 -q geneset.faa -o kegg.btp >/dev/null 2>&1
perl -ne 'chomp;@s=split/\t/;if($s[0] ne $a){$s[1]=~/\|(.*)/;print "$s[0]\t$1\n";$a=$s[0]}' kegg.btp > kegg.tsv

# 根据聚类结果，分配每个样本KO
get_profile_from_cluster.py kegg.tsv ../10.geneset/mmseqs/cluster_cluster.tsv kegg.count.profile

# KO丰度
cat ../sample_name | parallel -j 5 calcu_tpm_for_KO.py kegg.tsv ../10.geneset/abundance/{}.cvg KOs.ab/{}.tpm
combine_file_zy_folder_allsample.py -D KOs.ab -suffix .tpm -o KOs.tpm
kegg_convert_levels.py -i KOs.tpm --lineage F -o KOs.tpm

# gsea
le /share/data1/database/KEGG/KO_level_A_B_C_D_Description | perl -e 'open I, "pathway.id";%h;while(<I>){chomp;$h{$_}=1};while(<>){chomp;@s=split/\t/;$s[5]=~/(.*)\s\[\S+:(\S+?)\]/;print "$2:$1\t$s[6]\n" if exists $h{$s[4]}}' > pathway_library.tsv
../calcu_gsea_v2.py ../KOs.tpm sample_group group_pair pathway_library.tsv gsea &

cat /share/data1/database/KEGG/KEGG20230401/KEGG_API/database/KO_Pathway.txt | grep map | sed 's/ko://;s/path://' | csvtk join -t -H --left-join -f "1" KOs.list - | awk '{if($2!=""){print $2}}' | sort -u | wc -l

