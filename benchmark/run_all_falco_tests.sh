mkdir -p outs/falco
for i in tests/fastq/test*.fastq
do
  a=$(basename $i)
  echo $a
  time falco -o outs/falco -i $i \
       1>outs/falco/output_${a} \
       2>outs/falco/error_${a}
done
