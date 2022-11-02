#We here start with the trimmed fastq files for the 347 selected files. 
#Before going in, we need to upload the fasta file that will be used as the reference genome. 
cd ~/Labbet/Bioinformatics/Saved_transcriptome_files/HomoS_GRCh38.98_cdna_ncrna_ERCC
sftp -q jakob-sens2021562@bianca-sftp.uppmax.uu.se:jakob-sens2021562
put Homo_sapiens.GRCh38.fa.gz
put Homo_sapiens.GRCh38.98_original.gtf.gz

ssh -A jakob-sens2021562@bianca.uppmax.uu.se
cd /proj/sens2021562
mkdir 220929_LGI1_B_velocyto

cd /proj/sens2021562/220929_LGI1_B_velocyto

#Step one now is to align the data. We started by using Kallisto, but that was not compatible with the downstream analyses, so we
#had to go over to Star. 

#See https://nbisweden.github.io/workshop-RNAseq/1911/lab_kallisto.html for the current workflow. 
#mkdir /proj/sens2021562/fasta_gtf_etc
#cd /proj/sens2021562/fasta_gtf_etc
#mv /proj/sens2021562/nobackup/wharf/jakob/jakob-sens2021562/Homo_sapiens.GRCh38.98_original.gtf.gz Homo_sapiens.GRCh38.98_original.gtf.gz
#gunzip Homo_sapiens.GRCh38.98_original.gtf.gz
#So now we have placed the fasta and gtf files in a useful space.

cd /proj/sens2021562/220929_LGI1_B_velocyto/selected_fastq

cd /proj/sens2021562/fasta_gtf_etc
nano star_genome.sh
#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 05:00:00 
#SBATCH -J star_genome # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END
module load bioinfo-tools
module load star/2.7.9a

STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles /sw/data/reference/Homo_sapiens/GRCh38/concat/Homo_sapiens.GRCh38.dna.concat.fa --runThreadN 16

sbatch star_genome.sh
jobinfo -u jakob

#And here comes the real analyis

cd /proj/sens2021562/220929_LGI1_B_velocyto
mkdir star
cd /proj/sens2021562/220929_LGI1_B_velocyto/run_folder
nano star.sh
#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 05:00:00 
#SBATCH -J star # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END
cd /proj/sens2021562/220929_LGI1_B_velocyto/selected_fastq
module load bioinfo-tools
module load star/2.7.9a

arr1=($(echo *_R1_001_val_1.fq.gz))
arr2=($(echo *_R2_001_val_2.fq.gz))

for i in {0..114}; 
do
STAR  --runThreadN 16 --outSAMstrandField intronMotif --readFilesCommand zcat --genomeDir /proj/sens2021562/fasta_gtf_etc --outSAMtype BAM SortedByCoordinate \
--readFilesIn ${arr1[$i]} ${arr2[$i]} \
--outFileNamePrefix /proj/sens2021562/220929_LGI1_B_velocyto/star2/${arr1[$i]:0:14}/
done

nano star2.sh
#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 05:00:00 
#SBATCH -J star # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END
cd /proj/sens2021562/220929_LGI1_B_velocyto/selected_fastq
module load bioinfo-tools
module load star/2.7.9a

arr1=($(echo *_R1_001_val_1.fq.gz))
arr2=($(echo *_R2_001_val_2.fq.gz))

for i in {115..230}; 
do
STAR  --runThreadN  --outSAMstrandField intronMotif --readFilesCommand zcat --genomeDir /proj/sens2021562/fasta_gtf_etc --outSAMtype BAM SortedByCoordinate \
--readFilesIn ${arr1[$i]} ${arr2[$i]} \
--outFileNamePrefix /proj/sens2021562/220929_LGI1_B_velocyto/star3/${arr1[$i]:0:14}/
done

nano star3.sh
#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 05:00:00 
#SBATCH -J star # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END
cd /proj/sens2021562/220929_LGI1_B_velocyto/selected_fastq
module load bioinfo-tools
module load star/2.7.9a

arr1=($(echo *_R1_001_val_1.fq.gz))
arr2=($(echo *_R2_001_val_2.fq.gz))

for i in {231..347}; 
do
STAR  --runThreadN  --outSAMstrandField intronMotif --readFilesCommand zcat --genomeDir /proj/sens2021562/fasta_gtf_etc --outSAMtype BAM SortedByCoordinate \
--readFilesIn ${arr1[$i]} ${arr2[$i]} \
--outFileNamePrefix /proj/sens2021562/220929_LGI1_B_velocyto/star3/${arr1[$i]:0:14}/
done

sbatch star.sh
sbatch star2.sh
sbatch star3.sh
jobinfo -u jakob

cd /proj/sens2021562/220929_LGI1_B_velocyto
mkdir velocyto
cd /proj/sens2021562/220929_LGI1_B_velocyto/run_folder
nano velocyto.sh

#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 10:00:00 
#SBATCH -J velocyto # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END

module load bioinfo-tools
module load velocyto/0.17.17

cd /proj/sens2021562/220929_LGI1_B_velocyto/star

velocyto run-smartseq2 -o /proj/sens2021562/220929_LGI1_B_velocyto/velocyto \
JR*/Aligned.sortedByCoord.out.bam \
/proj/sens2021562/fasta_gtf_etc/Homo_sapiens.GRCh38.98_original.gtf

sbatch velocyto.sh
jobinfo -u jakob


cp -r /proj/sens2021562/220929_LGI1_B_velocyto/velocyto \
/proj/sens2021562/nobackup/wharf/jakob/jakob-sens2021562/velocyto_out

cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/Velocity
sftp -q jakob-sens2021562@bianca-sftp.uppmax.uu.se:jakob-sens2021562
get -r velocyto_out







cd /proj/sens2021562/220929_LGI1_B_velocyto/velocyto
mkdir JR1166
cd /proj/sens2021562/220929_LGI1_B_velocyto/run_folder
nano velocyto.sh

#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 20:00:00 
#SBATCH -J velocyto_1166 # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END

module load bioinfo-tools
module load velocyto/0.17.17

cd /proj/sens2021562/220929_LGI1_B_velocyto/star

velocyto run-smartseq2 -o /proj/sens2021562/220929_LGI1_B_velocyto/velocyto/JR1166 \
JR1166_1*/*.out.bam \
/proj/sens2021562/fasta_gtf_etc/Homo_sapiens.GRCh38.98_original.gtf

sbatch velocyto.sh

cd /proj/sens2021562/220929_LGI1_B_velocyto/velocyto
mkdir JR1227
cd /proj/sens2021562/220929_LGI1_B_velocyto/run_folder
nano velocyto1227.sh

#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 20:00:00 
#SBATCH -J velocyto_1227 # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END

module load bioinfo-tools
module load velocyto/0.17.17

cd /proj/sens2021562/220929_LGI1_B_velocyto/star

velocyto run-smartseq2 -o /proj/sens2021562/220929_LGI1_B_velocyto/velocyto/JR1227 \
JR1227_1*/*.out.bam \
/proj/sens2021562/fasta_gtf_etc/Homo_sapiens.GRCh38.98_original.gtf

sbatch velocyto1227.sh

cd /proj/sens2021562/220929_LGI1_B_velocyto/velocyto
mkdir JR1284
cd /proj/sens2021562/220929_LGI1_B_velocyto/run_folder
nano velocyto1284.sh

#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 20:00:00 
#SBATCH -J velocyto_1282 # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END

module load bioinfo-tools
module load velocyto/0.17.17

cd /proj/sens2021562/220929_LGI1_B_velocyto/star

velocyto run-smartseq2 -o /proj/sens2021562/220929_LGI1_B_velocyto/velocyto/JR1284 \
JR1284_1*/*.out.bam \
/proj/sens2021562/fasta_gtf_etc/Homo_sapiens.GRCh38.98_original.gtf

sbatch velocyto1284.sh

jobinfo -u jakob


#We needed to change the names of the files to get sensical names. 
cd /proj/sens2021562/220929_LGI1_B_velocyto/star
R
bamList <- list.files(".", recursive = TRUE, pattern = ".bam")
sapply(bamList[4:length(bamList)], function(x) file.rename(x, paste0(substr(x, 1, 15), gsub("|/Aligned.sortedByCoord|", "", x))))