# change this directory to any directory with a set of fastq/fastq.gz/sam/bam
# files to run the benchmarking in all of them
dirs="tests/fastq"
for i in `ls ${dirs}`
do
 # Remove the "-p" flag to run without overrepresentation
 echo $i;
 time fastp -V -A -G -Q -L -p -i ${dirs}/${i} \
       -h outs/fastp/${i}.html \
       -j outs/fastp/${i}.json \
       1>outs/fastp/output_${i}.txt \
       2>outs/fastp/error_${i}.txt
done
