# Reduced responsiveness of general stress response genes contributes to the higher toxicity of a pesticide under warming 
# Vienna Delnat, Janne Swaegers, Jana Asselman and Robby Stoks
# Command line used on the VSC (Vlaams Supercomputer Centrum) using Putty


###### Bioconda ######

# Download Miniconda from link and install on DATA server
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Activate channels 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels r
conda config --add channels conda-forge

# Update conda and install bioconda (bwa)
conda update -n base -c defaults conda
conda install bwa

# Packages added in environments
conda create -y -n denovoV3_env python=3.6.7 trinity=2.8.4 trinotate=3.1.1 transdecoder=5.5.0 samtools=1.9 salmon=0.12.0 cd-hit=4.6.8 igv=2.4.17 r=3.5.1 bioconductor-deseq2=1.22.1 bioconductor-edger=3.24.1 bioconductor-tximport=1.10.0 bioconductor-qvalue=2.14.0 r-fastcluster=1.1.25 bioconductor-goseq=1.34.1

conda create -y -n pythonV2_env python=2.7 hisat2=2.1.0 busco=3.0.2

conda create -y -n referenceV3_env python=3.6.7 stringtie=1.3.4 samtools=1.9 gffcompare=0.10.6 cd-hit=4.6.8 igv=2.4.17 r=3.5.1 bioconductor-deseq2=1.22.1 bioconductor-edger=3.24.1 bioconductor-ballgown=2.14.0 bioconductor-tximport=1.10.0 bioconductor-qvalue=2.14.0 r-fastcluster=1.1.25


###### Combine all left (forward) and right (reverse) reads in two separate files using cat ######

cat *_1.fq.gz > /data/reads.ALL.left.fq.gz

cat *_2.fq.gz > /data/reads.ALL.right.fq.gz


###### Denovo assembly using Trinity ######

source activate denovoV3_env

Trinity --seqType fq --left /data/reads.ALL.left.fq.gz --right /data/reads.ALL.right.fq.gz --CPU 10 --max_memory 60G


###### Statistics of Trinity (N50) ######

conda activate denovoV3_env
TrinityStats.pl /Denovo/trinity/trinity_out_dir/Trinity20190213.fasta


###### Quality check using ExN50 instead of N50 ######

source activate denovoV3_env

contig_ExN50_statistic.pl /Denovo/countSalmonCDHIT/Trinity.isoform.TMM.EXPR.matrix /Denovo/trinityCDHIT/trinity_out_dir/Trinity20190213.cdhit.fasta > ExN50.stats

# Show ExN50 statistics
cat ExN50.stats | column -t


###### Assembly improvement using CDHIT (removal of duplicates/redundancies) ######

source activate denovoV3_env

cd-hit-est -i /Denovo/trinity/trinity_out_dir/Trinity20190213.fasta -o /Denovo/trinity/trinity_out_dir/Trinity20190213.cdhit.fasta -c 0.95 -n 10 -T 15  


###### Quality check of fasta file and comparison with other transcriptomes BUSCO ######

# Diptera dataset
tar -xvzf diptera_odb9.tar.gz

source activate pythonV2_env

run_busco -m transcriptome -i /Denovo/trinityCDHIT/trinity_out_dir/Trinity20190213.cdhit.fasta -o Trinity.culexCDHIT.filter -l /data/diptera_odb9 -sp aedes


###### Map samples against de novo assembly ######

source activate pythonV2_env

hisat2-build /Denovo/trinityCDHIT/trinity_out_dir/Trinity20190213.cdhit.fasta /Denovo/alignmentHisatCDHIT/TrinityCDHITIndexed

hisat2 -p 15 -x /Denovo/alignmentHisatCDHIT/TrinityCDHITIndexed -1 /data/reads.ALL.left.fq.gz -2 /data/reads.ALL.right.fq.gz -S /Denovo/alignmentHisatCDHIT/alignmentCDHIT20190318.sam


###### Convert sam to bam and sort bam file using Samtools ######

source activate denovoV3_env

/Denovo/alignmentHisatCDHIT/alignmentCDHIT20190318.bam /Denovo/alignmentHisatCDHIT/alignmentCDHIT20190318.sam
samtools view -b -S -o 

samtools sort /Denovo/alignmentHisatCDHIT/alignmentCDHIT20190318.bam -o /Denovo/alignmentHisatCDHIT/alignmentCDHIT20190318.sorted.bam

samtools index /Denovo/alignmentHisatCDHIT/alignmentCDHIT20190318.sorted.bam


###### Count data as input for R using salmon ######

source activate denovoV3_env

salmon index -t /Denovo/trinityCDHIT/trinity_out_dir/Trinity20190213.cdhit.fasta -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientSolventControl1_1.fq.gz -2 /FASTQ/AmbientSolventControl1_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientSolventControl1

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientSolventControl2_1.fq.gz -2 /FASTQ/AmbientSolventControl2_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientSolventControl2

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientSolventControl3_1.fq.gz -2 /FASTQ/AmbientSolventControl3_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientSolventControl3

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientSolventControl4_1.fq.gz -2 /FASTQ/AmbientSolventControl4_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientSolventControl4

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientSolventControl5_1.fq.gz -2 /FASTQ/AmbientSolventControl5_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientSolventControl5

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientSolventControl6_1.fq.gz -2 /FASTQ/AmbientSolventControl6_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientSolventControl6

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientLowChlorpyrifos1_1.fq.gz -2 /FASTQ/AmbientLowChlorpyrifos1_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientLowChlorpyrifos1

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientLowChlorpyrifos2_1.fq.gz -2 /FASTQ/AmbientLowChlorpyrifos2_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientLowChlorpyrifos2

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientLowChlorpyrifos3_1.fq.gz -2 /FASTQ/AmbientLowChlorpyrifos3_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientLowChlorpyrifos3

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientLowChlorpyrifos4_1.fq.gz -2 /FASTQ/AmbientLowChlorpyrifos4_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientLowChlorpyrifos4

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientLowChlorpyrifos5_1.fq.gz -2 /FASTQ/AmbientLowChlorpyrifos5_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientLowChlorpyrifos5

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientLowChlorpyrifos6_1.fq.gz -2 /FASTQ/AmbientLowChlorpyrifos6_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientLowChlorpyrifos6

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientHighChlorpyrifos1_1.fq.gz -2 /FASTQ/AmbientHighChlorpyrifos1_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientHighChlorpyrifos1

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientHighChlorpyrifos2_1.fq.gz -2 /FASTQ/AmbientHighChlorpyrifos2_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientHighChlorpyrifos2

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientHighChlorpyrifos3_1.fq.gz -2 /FASTQ/AmbientHighChlorpyrifos3_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientHighChlorpyrifos3

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientHighChlorpyrifos4_1.fq.gz -2 /FASTQ/AmbientHighChlorpyrifos4_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientHighChlorpyrifos4

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientHighChlorpyrifos5_1.fq.gz -2 /FASTQ/AmbientHighChlorpyrifos5_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientHighChlorpyrifos5

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/AmbientHighChlorpyrifos6_1.fq.gz -2 /FASTQ/AmbientHighChlorpyrifos6_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/AmbientHighChlorpyrifos6

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingSolventControl1_1.fq.gz -2 /FASTQ/WarmingSolventControl1_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingSolventControl1

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingSolventControl2_1.fq.gz -2 /FASTQ/WarmingSolventControl2_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingSolventControl2

salmon quant -p 12 --validateMappings --gcBias -l IU -i ./Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 ./FASTQ/WarmingSolventControl3_1.fq.gz -2 /FASTQ/WarmingSolventControl3_2.fq.gz -o ./Denovo/countSalmonCDHIT/validatemappings/WarmingSolventControl3

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingSolventControl4_1.fq.gz -2 /FASTQ/WarmingSolventControl4_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingSolventControl4

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingSolventControl5_1.fq.gz -2 /FASTQ/WarmingSolventControl5_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingSolventControl5

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingSolventControl6_1.fq.gz -2 /FASTQ/WarmingSolventControl6_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingSolventControl6

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingLowChlorpyrifos1_1.fq.gz -2 /FASTQ/WarmingLowChlorpyrifos1_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingLowChlorpyrifos1

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingLowChlorpyrifos2_1.fq.gz -2 /FASTQ/WarmingLowChlorpyrifos2_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingLowChlorpyrifos2

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingLowChlorpyrifos3_1.fq.gz -2 /FASTQ/WarmingLowChlorpyrifos3_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingLowChlorpyrifos3

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingLowChlorpyrifos4_1.fq.gz -2 /FASTQ/WarmingLowChlorpyrifos4_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingLowChlorpyrifos4

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingLowChlorpyrifos5_1.fq.gz -2 /FASTQ/WarmingLowChlorpyrifos5_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingLowChlorpyrifos5

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingLowChlorpyrifos6_1.fq.gz -2 /FASTQ/WarmingLowChlorpyrifos6_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingLowChlorpyrifos6

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingHighChlorpyrifos1_1.fq.gz -2 /FASTQ/WarmingHighChlorpyrifos1_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingHighChlorpyrifos1

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingHighChlorpyrifos2_1.fq.gz -2 /FASTQ/WarmingHighChlorpyrifos2_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingHighChlorpyrifos2

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingHighChlorpyrifos3_1.fq.gz -2 /FASTQ/WarmingHighChlorpyrifos3_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingHighChlorpyrifos3

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingHighChlorpyrifos4_1.fq.gz -2 /FASTQ/WarmingHighChlorpyrifos4_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingHighChlorpyrifos4

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingHighChlorpyrifos5_1.fq.gz -2 /FASTQ/WarmingHighChlorpyrifos5_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingHighChlorpyrifos5

salmon quant -p 12 --validateMappings --gcBias -l IU -i /Denovo/countSalmonCDHIT/validatemappings/salmon_index -1 /FASTQ/WarmingHighChlorpyrifos6_1.fq.gz -2 /FASTQ/WarmingHighChlorpyrifos6_2.fq.gz -o /Denovo/countSalmonCDHIT/validatemappings/WarmingHighChlorpyrifos6

