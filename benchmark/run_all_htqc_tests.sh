mkdir -p outs/htqc
for i in tests/fastq/test*.fastq
do
  a=$(basename $i)
  echo $a
  time ht-stat -S -i $i -o outs/htqc \
       1>outs/falco/output_${a} \
       2>outs/falco/error_${a}
done
