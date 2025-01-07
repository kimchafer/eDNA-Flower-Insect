
mkdir -p /storage2/jihoonkim/eDNA/ncbi
cd /storage2/jihoonkim/eDNA/ncbi
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvzf taxdump.tar.gz

mmseqs createtaxdb /storage2/jihoonkim/eDNA/mmseqsDB/cox1/cox1_refDB \
/storage2/jihoonkim/eDNA/mmseqsDB/taxonomy/mmseqs_taxDB \
--ncbi-tax-dump /storage2/jihoonkim/eDNA/ncbi/taxonomy







