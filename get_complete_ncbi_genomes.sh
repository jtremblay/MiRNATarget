## Bacteria

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > assembly_summary_complete_genomes.txt
#awk -F '\t' '{print $20}' assembly_summary_selected.txt > assembly_summary_complete_genomes_selected.txt
mkdir -p GbBac

for next in $(cat assembly_summary_complete_genomes.txt); do
   filename="${next##*/}"
   if [ ! -f GbBac/${filename}_genomic.fna.gz ]; then
      echo GbBac/"$filename"_genomic.fna.gz not found!
      wget -P GbBac "$next"/"$filename"_genomic.fna.gz;
   fi
done
gunzip GbBac/*.gz
#cat GbBac/*.fna > all_complete_Gb_bac.fasta

## Archaea
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt -O assembly_summary_arch.txt
#awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary_arch.txt > assembly_summary_arch_complete_genomes.txt
#mkdir GbArch
#for next in $(cat assembly_summary_arch_complete_genomes.txt); do 
#   filename="${next##*/}"
#   wget -P GbArch "$next"/"$filename"_genomic.fna.gz;
#done
#gunzip GbArch/*.gz
##cat GbArch/*.fna > all_complete_Gb_arch.fasta
