#Trim adapters and sort reads by inline barcodes
#Adaptor + Index 1: CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTG
mkdir ~/raw_i2
mv ~/filtered/QualFilt_PPA_i2.fastq ~/raw_i2/QualFilt_PPA_i2.fastq
/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i2 -b ~/Index02_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTG --adapter_mm 1 --filter_illumina

