mkdir ~/filtered

# Cortar segun resultados de FastQC
/users/user/documents/FastXToolkit/bin/fastx_trimmer -Q33 -l 90 -i ~/data/PPA_i1.fastq -o ~/filtered/Filt_PPA_i1.fastq

# Filtrar reads de mala calidad (Todas arriba de 20 de score filtre)
/users/user/documents/FastXToolkit/bin/fastq_quality_filter -Q33 -v -q 20 -p 100 -i ~/data/PPA_i1.fastq -o ~/filtered/QualFilt_PPA_i1.fastq
/users/user/documents/FastXToolkit/bin/fastq_quality_filter -Q33 -v -q 20 -p 100 -i ~/data/PPA_i2.fastq -o ~/filtered/QualFilt_PPA_i2.fastq
/users/user/documents/FastXToolkit/bin/fastq_quality_filter -Q33 -v -q 20 -p 100 -i ~/data/PPA_i3.fastq -o ~/filtered/QualFilt_PPA_i3.fastq
/users/user/documents/FastXToolkit/bin/fastq_quality_filter -Q33 -v -q 20 -p 100 -i ~/data/PPA_i4.fastq -o ~/filtered/QualFilt_PPA_i4.fastq
/users/user/documents/FastXToolkit/bin/fastq_quality_filter -Q33 -v -q 20 -p 100 -i ~/data/PPA_i5.fastq -o ~/filtered/QualFilt_PPA_i5.fastq
/users/user/documents/FastXToolkit/bin/fastq_quality_filter -Q33 -v -q 20 -p 100 -i ~/data/PPA_i6.fastq -o ~/filtered/QualFilt_PPA_i6.fastq
/users/user/documents/FastXToolkit/bin/fastq_quality_filter -Q33 -v -q 20 -p 100 -i ~/data/PPA_i7.fastq -o ~/filtered/QualFilt_PPA_i7.fastq
/users/user/documents/FastXToolkit/bin/fastq_quality_filter -Q33 -v -q 20 -p 100 -i ~/data/PPA_i12.fastq -o ~/filtered/QualFilt_PPA_i12.fastq

