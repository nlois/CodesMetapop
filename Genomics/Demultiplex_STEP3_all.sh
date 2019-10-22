#Trim adapters and sort reads by inline barcodes

#Adapter + Index 1: GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
mkdir ~/demultiplexed
mkdir ~/raw_i1
mv ~/filtered/QualFilt_PPA_i1.fastq ~/raw_i1/QualFilt_PPA_i1.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i1 -b ~/Index01_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

#Adapter + Index 2: GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
mkdir ~/raw_i2
mv ~/filtered/QualFilt_PPA_i2.fastq ~/raw_i2/QualFilt_PPA_i2.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i2 -b ~/Index02_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

#Adapter + Index 3: GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
mkdir ~/raw_i3
mv ~/filtered/QualFilt_PPA_i3.fastq ~/raw_i3/QualFilt_PPA_i3.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i3 -b ~/Index03_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

#Adapter + Index 4: GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
mkdir ~/raw_i4
mv ~/filtered/QualFilt_PPA_i4.fastq ~/raw_i4/QualFilt_PPA_i4.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i4 -b ~/Index04_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

#Adapter + Index 5: GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
mkdir ~/raw_i5
mv ~/filtered/QualFilt_PPA_i5.fastq ~/raw_i5/QualFilt_PPA_i5.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i5 -b ~/Index05_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

#Adapter + Index 6: GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
mkdir ~/raw_i6
mv ~/filtered/QualFilt_PPA_i6.fastq ~/raw_i6/QualFilt_PPA_i6.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i6 -b ~/Index06_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

#Adapter + Index 7: GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
mkdir ~/raw_i7
mv ~/filtered/QualFilt_PPA_i7.fastq ~/raw_i7/QualFilt_PPA_i7.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i7 -b ~/Index07_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

#Adapter + Index 12: GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
mkdir ~/raw_i12
mv ~/filtered/QualFilt_PPA_i12.fastq ~/raw_i12/QualFilt_PPA_i12.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i12 -b ~/Index12_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
