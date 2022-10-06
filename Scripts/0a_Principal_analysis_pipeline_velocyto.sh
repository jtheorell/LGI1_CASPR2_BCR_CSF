



#We here start with the trimmed fastq files
ssh -A jakob-sens2021562@bianca.uppmax.uu.se
cd /proj/sens2021562
mkdir 220929_LGI1_B_velocyto

cd /proj/sens2021562/220816_LGI1_tracer_reanalysis

#Now, we get over to doing the actual jobs
mkdir raw_fastq/200129
mkdir run_folder
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis/run_folder
nano 200129_bcl2fastq.sh
#And add
#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 02:00:00 # Two hours  
#SBATCH -J 200129_bcl2fastq # Name of the job
# go to some directory
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis
# load software modules
module load bioinfo-tools
module load bcl2fastq/2.20.0
bcl2fastq --runfolder-dir \
/proj/sens2021562/220816_LGI1_tracer_reanalysis/bcl/200128_NB501183_0771_AH5CHHBGXF \
-p 16 --output-dir \
/proj/sens2021562/220816_LGI1_tracer_reanalysis/raw_fastq/200129 \
--no-lane-splitting

sbatch 200129_bcl2fastq.sh
jobinfo -u jakob

#As JR0022 and JR0051 were unsuccessful, we will delete all the directories containing their data. 
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis/raw_fastq/200129/LGI1_CSF
rm -rf JR0022*
rm -rf JR0051*

###########################
#TRIM GALORE
###########################

cd /proj/sens2021562/220816_LGI1_tracer_reanalysis
mkdir trimmed_fastq/200129
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis/run_folder
nano 200129_trimGalore.sh

#!/bin/bash -l
#SBATCH -A sens2021562 # Project name
#SBATCH -p core # Asking for cores (as opposed to multiple nodes)
#SBATCH -n 16 # Number of cores
#SBATCH -t 01:00:00 # One hour  
#SBATCH -J 200129_trimgalore # Name of the job
# go to some directory
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis/raw_fastq/200129/LGI1_CSF
# load software modules
module load gnuparallel/20180822
module load bioinfo-tools
module load TrimGalore/0.6.1
parallel --xapply -j 16 trim_galore \
--nextera --paired --fastqc \
-o /proj/sens2021562/220816_LGI1_tracer_reanalysis/trimmed_fastq/200129 \
::: *R1_001.fastq.gz ::: *R2_001.fastq.gz

sbatch 200129_trimGalore.sh
jobinfo -u jakob

###########################
#Multiqc
###########################

#Now for some multiqc
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis
mkdir -p Posttrim_fastqc_result/200129

module load bioinfo-tools
module load MultiQC/1.12

multiqc --outdir Posttrim_fastqc_result/200129 trimmed_fastq

#Running this live does not seem like a great idea: takes about 20 minutes. 
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis
cp Posttrim_fastqc_result/200129/multiqc_report.html /proj/sens2021562/nobackup/wharf/jakob/jakob-sens2021562/200129_multiqc_report.html

#Then on the local system: 
cd ~/Labbet/2022/220815_JR1227_TraCeR_reanalysis/Data/MultiQC_reports

sftp -q jakob-sens2021562@bianca-sftp.uppmax.uu.se:jakob-sens2021562
get multiqc_report.html

###########################
#Combining files from run 1 and 2.
###########################

#Now, what we are going to do is to combine all files from run 1 and run 2 with identical names. 
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis/trimmed_fastq
R
filelist1 <- list.files("200129")
filelist2 <- list.files("200311")
identical(filelist1, filelist2)
#FALSE
#This is because the second experiment contains data also from JR0147. But for the length of fileList1, they are identical: 
identical(filelist1, filelist2[1:length(filelist1)])
quit()

#Now, we create an array of all the names we are going to concatenate. 
cd /proj/sens2021562/220816_LGI1_tracer_reanalysis/trimmed_fastq/200129
arr1=($(echo *.fq.gz))

cd /proj/sens2021562/220816_LGI1_tracer_reanalysis/trimmed_fastq
#And then we do it
for i in "${!arr1[@]}"; 
do
echo ${arr1[$i]}
cat 200129/${arr1[$i]} 200311/${arr1[$i]} > 200129-0311/${arr1[$i]} 
done

#And with that, we leave this dataset and go on to the combined one. 
arr1=($(echo JR1227*))
for i in {70..90}; 
do
stat ${arr1[$i]}
done