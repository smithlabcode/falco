mkdir -p outs/fastqc
for i in `ls tests/fastq`
do
  a=`echo $(basename $i) | sed 's/.fastq//g'`
  echo "[$(date) - fastqc] $a"
  mkdir -p outs/fastqc/${a}
  time fastqc -o outs/fastqc/${a} tests/fastq/${a}.fastq \
       1>outs/fastqc/${a}/${a}.output \
       2>outs/fastqc/${a}/${a}.error
done
