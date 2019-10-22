#ALL_andDEAD
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./popGen_dataSET_andDEAD.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --min_maf 0.01 &
# 4705 SNPs; loci 3882 (includes invariant)
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./popGen_dataSET_andDEAD.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --write_single_snp --min_maf 0.01 &
#2645 SNPs; loci 3882 (includes invariant)

#ALL
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./popGen_dataSET.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --min_maf 0.01 &
#4975 SNPs; 4061 loci (includes invariant)
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./popGen_dataSET.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --write_single_snp --min_maf 0.01 &
#2779; 4061 loci (includes invariant)

#NORTH
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./NORTH.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --min_maf 0.01 &
#4572 SNPs; 4228 loci (includes invariant)
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./NORTH.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --write_single_snp --min_maf 0.01 &
#2770 SNPs; 4228 loci (includes invariant)loci (includes invariant)

#SOUTH
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./SOUTH.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --min_maf 0.01 &
#4918 SNPs; 3953 loci (includes invariant)
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./SOUTH.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --write_single_snp --min_maf 0.01 &
#2718 SNPs; 3953 loci (includes invariant)

#NORTH_andDEAD
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./NORTH_andDEAD.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --min_maf 0.01 &
#4140 SNPs; 3658 loci (includes invariant)
/programs/stacks-1.44/bin/populations -b 1 -P ./ -M ./NORTH_andDEAD.txt -r 0.8 -p 1 -m 20 -t 15 --structure --vcf --write_single_snp --min_maf 0.01 &
#2440 SNPs; 3658 loci (includes invariant)



bcftools query -l  ...