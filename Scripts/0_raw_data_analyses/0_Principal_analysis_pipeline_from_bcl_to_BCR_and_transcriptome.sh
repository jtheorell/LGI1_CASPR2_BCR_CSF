#All computing here described was conoducted using the 
#Computational Biology Research Group cloud computing service.
#I have left out all information about creating and changing 
#directories, and similar information, 
#as the setup on another system would differ anyway. 

############################################
##bcl2fastq
############################################
#Much obliged to the following authors: 
#http://bioinformatics.cvr.ac.uk/blog/how-to-demultiplex-illumina-data-and-generate-fastq-files-using-bcl2fastq/

#The following specific information was added before submitting as a batch job:

module add bcl2fastq/2.20.0.422
bcl2fastq --runfolder-dir \
<PATH TO DIR> \
-p 16 --output-dir \
<PATH TO OUTPUT DIR> \
--no-lane-splitting \
--sample-sheet <PATH TO SAMPLE SHEET>

############################################
#First fastqc
############################################

#The following specific information was added before submitting as a batch job:

module add fastqc/0.11.8
arr1=($(echo *001.fastq.gz))
fastqc -t 16 \
-o <PATH TO OUTPUT DIR> \
${arr1[*]}

############################################
##trim-galore
############################################

#The following specific information was added before submitting as a batch job:

module add trim_galore/0.5.0
arr1=($(echo *R1_001.fastq.gz))
arr2=($(echo *R2_001.fastq.gz))
for i in "${!arr1[@]}";
do
trim_galore --nextera --paired --fastqc \
-o <PATH TO OUTPUT DIR> \
${arr1[$i]} ${arr2[$i]}
done

############################################
##Multiqc
############################################
#Here, multiqc was used for both the pre- and post trimming data, using version 0.9
#and standard settings. 

############################################
##Salmon
############################################

#Before this, an transcriptome index file containing the HomoS_GRCh38.98_cdna_ncrna
#information as well as information about ERCC was created and creatively called
#HomoS_GRCh38.98_cdna_ncrna_ERCC_index

#The following specific information was added before submitting as a batch job:

#And enter
module add salmon
arr1=($(echo *R1_001_val_1.fq.gz))
arr2=($(echo *R2_001_val_2.fq.gz))
for i in "${!arr1[@]}"; 
do
salmon quant --seqBias --validateMappings --rangeFactorizationBins 4 \
-i HomoS_GRCh38.98_cdna_ncrna_ERCC_index \
  -l A -1 ${arr1[$i]} -2 ${arr2[$i]} -p 30 \
  -o <PATH TO OUTPUT DIRECTORY>/${arr1[$i]:0:12}/
done

#This whole directory was then exported and imported to R using tximeta, see analysis downstream. 

############################################
##Multiqc for salmon data
############################################

############################################
##Bracer
############################################

#The custom trnascriptome was added for the analysis, but as the TPM numbers for the
#BCRs only were used for relative measurements internal to each cell downstream,
#it was uncessesarily complicated to do this, and the step is therefore omitted here
#as it would not change the results. 

#The following specific information was added before submitting as a batch job:

module add bracer
arr1=($(echo *R1_001_val_1.fq.gz))
arr2=($(echo *R2_001_val_2.fq.gz))
for i in "${!arr1[@]}"; 
do
bracer assemble --resume_with_existing_files -p 16 \
${arr1[$i]:0:12} <PATH TO OUTPUT DIRECTORY> ${arr1[$i]} ${arr2[$i]} 
done

#Then, summarisation was conducted interactively:
module add bracer
bracer summarise --include_multiplets <PATH TO OUTPUT DIRECTORY>

