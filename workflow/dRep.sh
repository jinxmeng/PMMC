dRep compare bins_cls_99 -g bins.list -d -pa 0.9 -sa 0.99 -nc 0.30 -cm larger -p 50 --S_algorithm fastANI
dRep compare bins_cls_95 -g bins_cls_99.list -d -pa 0.9 -sa 0.99 -nc 0.30 -cm larger -p 50 --S_algorithm fastANI

# plus isolates
dRep compare total_genomes_cls_99 -g total_genomes.list -d -pa 0.9 -sa 0.99 -nc 0.30 -cm larger -p 50 --S_algorithm fastANI
parse_dRep_by_ckm2.pl total_genomes_cls_99/data_tables/Cdb.csv total_genomes_quality_report.tsv total_genomes_cls_99_res
cut -f3 total_genomes_cls_99_res | perl -e '%h;while(<>){chomp;$h{$_}=1};open I, "total_genomes.list";while(<I>){chomp;@s=split/\//;$s[-1]=~s/.fa//;print "$_\n" if exists $h{$s[-1]}}' > total_genomes_cls_99.list
dRep compare total_genomes_cls_95 -g total_genomes_cls_99.list -d -pa 0.9 -sa 0.95 -nc 0.30 -cm larger -p 50 --S_algorithm fastANI
