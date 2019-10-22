#Trim adapters and sort reads by inline barcodes
#Adaptor + Index 12: GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG

mkdir ~/demultiplexed
mkdir ~/raw
mv ~/filtered/QualFilt_PPA_i3.fastq ~/raw/QualFilt_PPA_i3.fastq
mv ~/raw/QualFilt_PPA_i2.fastq ~/filtered/QualFilt_PPA_i2.fastq

/Users/user/documents/stacks/bin/process_radtags -p ~/raw -b ~/Index03_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

