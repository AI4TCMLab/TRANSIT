# Pig single-cell RNA-seq data analysis
########## QC 1 ###########
genomedir="/home/tsh/software/cellranger/refdata-gex-mm10-2020-A/" ## Path to the reference genome index (Note: mm10 used)
datadir="/home/tsh/E-MTAB-7895/MF17010_GT17-12037_S7/0day_batch1/" ## Path to the raw data directory
sample='MF17010_GT17-12037' ## Sample identifier
#date
#for s in $sample
#do
#date
cellranger count --id=0day_batch_cellrange_out \	## --id: Specifies the output directory name
         --transcriptome=/home/tsh/software/cellranger/refdata-gex-mm10-2020-A \	## --transcriptome: Path to the CellRanger-compatible reference genome
         --fastqs=/home/tsh/E-MTAB-7895/MF17010_GT17-12037_S7/0day_batch1 \	## --fastqs: Directory containing FASTQ files (mkfastq or custom)
         --sample=MF17010_GT17-12037 \			## --sample: Must match the FASTQ filename prefix for software identification
         --nosecondary				## --nosecondary: Generates expression matrix only; skips secondary analysis (PCA/Clustering)

cellranger count --id=0day_batch2_cellrange_out  --transcriptome=/home/tsh/software/cellranger/refdata-gex-mm10-2020-A  --fastqs=/home/tsh/E-MTAB-7895/MF17010_GT17-12037_S7/0day_batch2   --sample=MF17010_GT17-12037  --nosecondary

date
wait
done
exit