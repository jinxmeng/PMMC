diamond blastp -d /share/data1/database/CARD/CARD_v3.2.7.dmnd --outfmt 6 --query-cover 80 --id 80 --max-target-seqs 5 -p 32 -q ../10.geneset/geneset.faa -o ARGs.btp > /dev/null 2>/dev/null
perl -ne 'chomp;@s=split/\t/;if($s[0] ne $a){$s[1]=~/(ARO:\d+)/;print "$s[0]\t$1\n";$a=$s[0]}' ARGs.btp > ARGs.tsv

