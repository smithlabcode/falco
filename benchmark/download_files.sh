files="SRR10124060 SRR10143153 SRR3897196 SRR9624732 SRR1853178 SRR6387347
SRR891268 SRR1772703 SRR9878537 SRR6059706"
mkdir -p tests/fastq
for i in $files
do
  echo "Downloading ${i}..."
  fastq-dump --outdir tests/fastq --skip-technical --readids \
               --read-filter pass --dumpbase --split-3 --clip \
               ${i}
done

