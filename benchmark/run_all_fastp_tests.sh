# change this directory to any directory with a set of fastq/fastq.gz/sam/bam
# files to run the benchmarking in all of them
mkdir -p outs/fastp
for i in `ls tests/fastq`
do
 # Remove the "-p" flag to run without overrepresentation
 a=`echo $(basename $i) | sed 's/.fastq//g'`
 mkdir -p outs/fastp/${a}
 echo "[$(date) - fastp] $a"
 time fastp -V -A -G -Q -L -w 1 -i tests/fastq/${a}.fastq \
       -h outs/fastp/${a}/fastp_report.html \
       -j outs/fastp/${a}/fastp_data.json \
       1>outs/fastp/${a}output_${a}.txt \
       2>outs/fastp/${a}/error_${a}.txt
done
