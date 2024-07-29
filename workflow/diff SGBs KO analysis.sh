diamond blastp -d /share/data1/database/KEGG/KEGG20230401.dmnd --outfmt 6 --min-score 60 --query-cover 50 --max-target-seqs 10 -p 112 -q faa.fa -o KO.btp >/dev/null 2>&1
perl -ne 'chomp;@s=split/\t/;if($s[0] ne $a){$s[1]=~/\|(.*)/;print "$s[0]\t$1\n";$a=$s[0]}' KO.btp > KO.tsv
 grep -f xx2 KO.tsv | perl -lne '$_=~/(\S+)\.\d+\s(K\d+)/;print "$1\t$2"' | csvtk join --left-join -t -f 1 - species.taxa | perl -lne '$_=~/\S+\s(K\d+)\s.*p__(\S+?);/;print "$1\t$2"' | sed 's/_\w$//g' | sort | uniq -c | perl -lne '$_=~/(\d+)\s+(K\d+)\s+(\S+)/;print "$2\t$3\t$1"'  > map00900.taxa.count
