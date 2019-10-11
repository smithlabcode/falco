files="SRR10124060 SRR10143153 SRR3897196 SRR9624732 SRR1853178 SRR6387347
SRR891268 SRR1772703 SRR9878537"
for i in $files
do
  echo "Downloading ${i}..."
  fastq-dump --outdir  --gzip --skip-technical --readids \
               --read-filter pass --dumpbase --split-3 --clip \
               ${i}
done

