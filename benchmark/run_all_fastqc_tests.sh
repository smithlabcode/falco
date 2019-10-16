mkdir -p outs/fastqc
for i in tests/fastq/*
do
  a=$(basename $i)
  echo $a
  time fastqc -o outs/fastqc $i \
       1>outs/fastqc/output_${a} \
       2>outs/fastqc/error_${a}
done
