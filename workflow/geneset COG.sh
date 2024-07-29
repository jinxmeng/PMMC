emapper.py --cpu 112 -i ../10.geneset/geneset.faa --itype proteins --output_dir . --output COG
cut -f1,7 COG.emapper.annotations | grep -v '^#\|-' > COG.annotation.tsv
