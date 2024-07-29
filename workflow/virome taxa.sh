# genomad
genomad end-to-end --quiet --threads 80 --splits 32 ../20.vir_cluster/final_vOTUs.fa genomad.out /share/data1/database/genomad/genomad_db/

# vcontact2
vcontact2_gene2genome --proteins ../20.vir_cluster/final_vOTUs.faa --output final_vOTUs.g2g.csv --source-type Prodigal-FAA
vcontact2 --db 'ProkaryoticViralRefSeq211-Merged' --threads 112 --raw-proteins ../20.vir_cluster/final_vOTUs.faa --proteins-fp final_vOTUs.g2g.csv --output-dir vcontact2.out

# grorc flow
diamond blastp --db /share/data1/database/virus_tax/virus_guorc/tot.faa.dmnd --query final_vOTUs.faa --threads 112 --id 30 --subject-cover 50 --query-cover 50 --min-score 50 --outfmt 6 --out final_vOTUs.btp --quiet &
grep '^>' final_vOTUs.faa | sed 's/^>//g' | perl -pne 's/_(\d+) #.*//' | sort | uniq -c | awk '{print $2"\t"$1}' > final_vOTUs.gn
vctg_stat.v2.pl final_vOTUs.btp /share/data1/database/virus_tax/virus_guorc/tot.tax final_vOTUs.gn final_vOTUs.btp.tax
awk '$2!="NA"' final_vOTUs.btp.tax > final_vOTUs.btp.tax.f
msort -k 1,rn3,rn6 final_vOTUs.btp.tax.f | perl -ne 'chomp;@s=split /\s+/;next if exists $h{$s[0]}; $pct=$s[2]/$s[4]; next unless ($pct>=0.2 or $s[2]>=10) and $s[5]>=30; printf "$s[0]\t$s[2]/$s[4]\t%.2f\t$s[1]\n",$s[5];$h{$s[0]}=1;' > final_vOTUs.tax.family
