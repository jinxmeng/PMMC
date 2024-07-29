ls *list | parallel -j 2 run_contigs_rename.pl {} {/.}
$ ls *fa | cut -d. -f1 | while read i;do echo run_virFinder_SOP.sh $i.m5000.fa 5000 $i 240;done
run_virFinder_SOP.sh PN10.m5000.fa 5000 PN10 240
run_virFinder_SOP.sh PN11.m5000.fa 5000 PN11 240
run_virFinder_SOP.sh PN12.m5000.fa 5000 PN12 240
run_virFinder_SOP.sh PN13.m5000.fa 5000 PN13 240
run_virFinder_SOP.sh PN14.m5000.fa 5000 PN14 240
run_virFinder_SOP.sh PN15.m5000.fa 5000 PN15 240
run_virFinder_SOP.sh PN16.m5000.fa 5000 PN16 240
run_virFinder_SOP.sh PN1.m5000.fa 5000 PN1 240
run_virFinder_SOP.sh PN2.m5000.fa 5000 PN2 240
run_virFinder_SOP.sh PN3.m5000.fa 5000 PN3 240
run_virFinder_SOP.sh PN4.m5000.fa 5000 PN4 240
run_virFinder_SOP.sh PN5.m5000.fa 5000 PN5 240
run_virFinder_SOP.sh PN6.m5000.fa 5000 PN6 240
run_virFinder_SOP.sh PN7.m5000.fa 5000 PN7 240
run_virFinder_SOP.sh PN8.m5000.fa 5000 PN8 240
run_virFinder_SOP.sh PN9.m5000.fa 5000 PN9 240
# run_virFinder_SOP.sh TP.m5000.fa 5000 TP 240
run_virFinder_SOP.sh isolate.m5000.fa 5000 isolate 240

cat PN*.phages/virus.fa isolate.phages/virus.fa > virus.fa &
cat PN*phages/ckv.out/quality_summary.tsv isolate.phages/ckv.out/quality_summary.tsv > ckv_quality_summary.tsv
cat PN*phages/vs2.out/final-viral-score.tsv isolate.phages/vs2.out/final-viral-score.tsv > virus.vs2.score.tsv
