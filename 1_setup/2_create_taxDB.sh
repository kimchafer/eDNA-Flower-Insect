
mkdir -p /storage2/jihoonkim/eDNA/ncbi
cd /storage2/jihoonkim/eDNA/ncbi
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvzf taxdump.tar.gz

# run gb_to_fasta.py

mmseqs createdb /storage2/jihoonkim/eDNA/ncbi/cox1/cox1_sequences_with_taxid.fasta \
/storage2/jihoonkim/eDNA/mmseqsDB/cox1_taxDB

mmseqs createtaxdb /storage2/jihoonkim/eDNA/mmseqsDB/cox1_taxDB \
/storage2/jihoonkim/eDNA/mmseqsDB/mmseqs_taxDB \
--ncbi-tax-dump /storage2/jihoonkim/eDNA/ncbi/taxonomy











mmseqs createtaxdb /storage2/jihoonkim/eDNA/ncbi/cox1/cox1_sequences_with_taxid.fasta \
/storage2/jihoonkim/eDNA/mmseqsDB/mmseqs_taxDB \
--ncbi-tax-dump /storage2/jihoonkim/eDNA/ncbi/taxonomy





