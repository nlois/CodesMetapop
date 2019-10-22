#Trim adapters and sort reads by inline barcodes
#Adaptor + Index 3: CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCAGACGTGTG
mkdir ~/raw_i3
mv ~/filtered/QualFilt_PPA_i3.fastq ~/raw_i3/QualFilt_PPA_i3.fastq

/Users/user/documents/stacks/bin/process_radtags -p ~/raw_i3 -b ~/Index03_Barcodes.txt -o ~/demultiplexed -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCAGACGTGTG --adapter_mm 1 --filter_illumina

