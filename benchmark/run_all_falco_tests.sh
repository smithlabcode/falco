mkdir -p outs/falco
for i in `ls tests/fastq`
do
  a=`echo $(basename $i) | sed 's/.fastq//g'`
  echo "[$(date) - falco] $a"
  mkdir -p outs/falco/${a}
  time falco -o outs/falco/${a} tests/fastq/${a}.fastq \
       1>outs/falco/${a}/${a}.output  \
       2>outs/falco/${a}/${a}.error
done
