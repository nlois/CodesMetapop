#Trim adapters and sort reads by inline barcodes
#Adaptor + Index 1: CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTG
mkdir ~/demultiplexed
mkdir ~/raw_i1
mv ~/filtered/QualFilt_PPA_i1.fastq ~/raw_i1/QualFilt_PPA_i1.fastq

/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i1 -b ~/Index01_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTG --adapter_mm 1 --filter_illumina

