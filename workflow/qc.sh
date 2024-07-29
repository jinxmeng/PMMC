# 20230530
cat ../sample_name | parallel -j 10 run_fastp_rmhost.sh ../00.Data/{}_1.fq.gz,../00.Data/{}_2.fq.gz pig {}
cat ../sample_name2 | parallel -j 10 run_fastp_rmhost.sh ../00.data/stage2_20230928/rawdata/{}_1.fq.gz,../00.data/stage2_20230928/rawdata/{}_2.fq.gz pig {} 
# 20231119
cat ../00.data/stage3_20231118/X102SC23102967-Z02-J002_01/sample_list | parallel -j 1 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat ../00.data/stage3_20231118/X102SC23102967-Z02-J002_02/sample_list | parallel -j 6 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat ../00.data/stage3_20231118/X102SC23102967-Z02-J002_03/sample_list | parallel -j 6 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat ../00.data/stage3_20231118/X102SC23102967-Z02-J002_04/sample_list | parallel -j 2 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat ../00.data/stage3_20231118/X102SC23102967-Z02-J002_05/sample_list | parallel -j 2 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat ../00.data/stage3_20231118/X102SC23102967-Z02-J001/sample_list | parallel -j 2 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat ../00.data/stage3_20231118/X101SC23083595-Z01-J008/sample_list | parallel -j 2 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat ../00.data/stage3_20231118/X101SC23083595-Z01-J007/sample_list | parallel -j 2 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
run_fastp_rmhost.sh ../01.cleandata/rmhost_fq/PN6RC_rmhost_1.fq.gz,../01.cleandata/rmhost_fq/PN6RC_rmhost_2.fq.gz pig PN6RC
# 20231129
cat ../00.data/sample_list | parallel -j 6 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
# 20231201
cat /share/data1/mjx/proj/04.black_pig_metagenome_20230529/00.data/stage7_20231201/sample_list | parallel -j 2 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat /share/data1/mjx/proj/04.black_pig_metagenome_20230529/00.data/stage8_20231201/X102SC23102967-Z02-J003/sample_list | parallel -j 4 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
cat /share/data1/mjx/proj/04.black_pig_metagenome_20230529/00.data/stage8_20231201/X101SC23083595-Z01-J004/sample_list | parallel -j 4 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
# 20231202
cat /share/data1/mjx/proj/04.black_pig_metagenome_20230529/00.data/stage8_20231201/X101SC23083595-Z01-J004/sample_list2 | parallel -j 4 --colsep="\t" run_fastp_rmhost.sh {2},{3} pig {1}
# 20231222
cat ../sample_list | parallel -j 4 --colsep='\t' run_fastp_rmhost.sh {2},{3} pig {1}

# reads before fastp
grep -A 2 "before filtering" *.log | perl -ne 'if(/(\S+?)_.*Read/){print "$1\t"}elsif(/reads:\s(\d+)$/){print "$1\t"}elsif(/bases:\s(\d+)$/){print "$1\n"}' | awk 'BEGIN{print"sample\tr1_reads\tr1_bases\tr2_reads\tr2_bases"}{if(NR%2==1){ORS="\t";print}else{ORS="\n";print $2"\t"$3}}' > ../stat_raw_reads.tsv
# reads after fastp
grep -A 2 "after filtering" *.log | perl -ne 'if(/(\S+?)_.*Read/){print "$1\t"}elsif(/reads:\s(\d+)$/){print "$1\t"}elsif(/bases:\s(\d+)$/){print "$1\n"}' | awk 'BEGIN{print"sample\tr1_reads\tr1_bases\tr2_reads\tr2_bases"}{if(NR%2==1){ORS="\t";print}else{ORS="\n";print $2"\t"$3}}' > ../stat_clean_reads.tsv
# map2host rate 
grep overall *.log | perl -ne '$_=~/(\S+?)_.*:(\S+)%/;print "$1\t$2\n"' | sort -rnk2 | csvtk add-header -t -n "sample,map2host(%)" > ../stat_map2host_rate.tsv
# unmap2host reads
grep 'mates make up the pairs; of these' *log | perl -ne '$_=~/(\S+?)_.*log:\s+(\d+?)\s/;print "$1\t$2\n"' |  csvtk add-header -t -n "sample,unmap_reads" > ../stat_unmap2host_reads.tsv

# 20m
le clean_fq_list | parallel -j 2 --colsep="\t" seqkit head -n 20000000 {2} -o clean_fq_20m/{1}_rmhost_20m_1.fq.gz
le clean_fq_list | parallel -j 2 --colsep="\t" seqkit head -n 20000000 {3} -o clean_fq_20m/{1}_rmhost_20m_2.fq.gz

