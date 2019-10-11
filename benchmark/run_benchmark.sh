# change this directory to any directory with a set of fastq/fastq.gz/sam/bam
# files to run the benchmarking in all of them
dirs="tests/fastq"
for i in `ls ${dirs}`
do
  # falco
  # Change the Configuration/limits.txt file and set overrepresented ignore 1
  # and duplication ignore 1 to disable overrepresentation
  echo "running falco on ${i}"
  time falco -o benchmark/outs/falco ${dirs}/${i} \
       1>benchmark/outs/falco/output_${i}.txt \
       2>benchmark/outs/falco/error_${i}.txt
  echo "running fastp on ${i}"

  # Remove the "-p" flag to run without overrepresentation
  time fastp -V -A -G -i ${dirs}/${i} \
       -h benchmark/outs/fastp/${i}.html \
       -j benchmark/outs/fastp/${i}.json \
       1>benchmark/outs/falco/output_${i}.txt \
       2>benchmark/outs/falco/error_${i}.txt

  # Change the Configuration the same way as above to disable the
  # overrepresented module on FastQC
  echo "running fastqc on ${i}"
  time fastqc -o benchmark/outs/fastqc benchmark/tests/${i} \
       1>benchmark/outs/fastqc/output_${i}.txt \
       2>benchmark/outs/fastqc/error_${i}.txt
done
