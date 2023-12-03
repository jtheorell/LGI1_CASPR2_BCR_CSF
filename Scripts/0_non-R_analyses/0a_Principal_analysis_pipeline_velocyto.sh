ssh -A jakob-sens2021562@bianca.uppmax.uu.se

#Step one now is to align the data. We started by using Kallisto, but that was not compatible with the downstream analyses, so we
#had to go over to Star. We will have to run this in much parallel. 
#FIrst, how many files are there? 
cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/selected_fastq
ls | wc -l
#1048, so 524 total cells. we divide this in ten parts, and run these in parallel. We will only keep one of the scripts here, 
#as the only difference is that we change: 
#for i in {0..52};
#for i in {53..104};
#for i in {105..156};
#for i in {157..208};
#for i in {209..260};
#for i in {261..312};
#for i in {313..364};
#for i in {365..417};
#for i in {418..470};
#for i in {471..523};

cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto
mkdir star
mkdir run_folder
cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/run_folder
nano star0.sh
#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 04:00:00 
#SBATCH -J star0 # Name of the job
#SBATCH --mail-user=jakob.theorell@ki.se
#SBATCH --mail-type=END
cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/selected_fastq
module load bioinfo-tools
module load star/2.7.9a

arr1=($(echo *_R1_001_val_1.fq.gz))
arr2=($(echo *_R2_001_val_2.fq.gz))

for i in {0..52};
do
STAR  --runThreadN 16 --outSAMstrandField intronMotif --readFilesCommand zcat --genomeDir /proj/sens2021562/fasta_gtf_etc --outSAMtype BAM SortedByCoordinate \
--readFilesIn ${arr1[$i]} ${arr2[$i]} \
--outFileNamePrefix /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/star/${arr1[$i]:0:14}/
done

sbatch star0.sh #
sbatch star53.sh #
sbatch star105.sh #
sbatch star157.sh #
sbatch star209.sh #
sbatch star261.sh #
sbatch star313.sh #
sbatch star365.sh #
sbatch star418.sh #
sbatch star471.sh #

jobinfo -u jakob

#We needed to change the names of the files to get sensical names. 
cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/star
R
bamList <- list.files(".", recursive = TRUE, pattern = "Coord.out.bam")
sapply(bamList[1:length(bamList)], function(x) file.rename(x, paste0(substr(x, 1, 15), gsub("|/Aligned.sortedByCoord|", "", x))))


#And here comes a large set of velocyto runs. 
cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto
mkdir velocyto

cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/velocyto
mkdir JR1166
cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/run_folder
nano velocyto1166.sh

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

cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/star

velocyto run-smartseq2 -o /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/velocyto/JR1166 \
JR1166_1*/*.out.bam \
/proj/sens2021562/fasta_gtf_etc/Homo_sapiens.GRCh38.98_original.gtf

sbatch velocyto1166.sh

mkdir /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/velocyto/JR1227

cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/run_folder
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

cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/star

velocyto run-smartseq2 -o /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/velocyto/JR1227 \
JR1227_1*/*.out.bam \
/proj/sens2021562/fasta_gtf_etc/Homo_sapiens.GRCh38.98_original.gtf

sbatch velocyto1227.sh

mkdir /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/velocyto/JR1284

cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/run_folder
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

cd /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/star

velocyto run-smartseq2 -o /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/velocyto/JR1284 \
JR1284_1*/*.out.bam \
/proj/sens2021562/fasta_gtf_etc/Homo_sapiens.GRCh38.98_original.gtf

sbatch velocyto1284.sh

jobinfo -u jakob

#And here the cells are exported

cp -r /proj/sens2021562/nobackup/221215_LGI1_B_velocyto/velocyto \
/proj/sens2021562/nobackup/wharf/jakob/jakob-sens2021562/velocyto_out

cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/Velocity
sftp -q jakob-sens2021562@bianca-sftp.uppmax.uu.se:jakob-sens2021562
get -r velocyto_out


cp -r /proj/sens2021562/220929_LGI1_B_velocyto/raw_fastq \
/proj/sens2021562/nobackup/wharf/jakob/jakob-sens2021562/220929_LGI1_B_velocyto_raw_fastq_real