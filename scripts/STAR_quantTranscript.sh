#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --account=project_XXXXX
#SBATCH --mem=140G
#SBATCH --time=4:00:00

### calculate genome indexes based on the Arabidopsis genome. 
###You don't need to recalculte the indexes if they already exist.

#STAR --runMode genomeGenerate --genomeDir star_index  \
#--genomeFastaFiles TAIR10_chr_all.fas \
#--sjdbGTFfile /projappl/project_2009761/star/atRTD3_TS_21Feb22_transfix.gtf \
#-runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 48 --limitGenomeGenerateRAM 
140000000000

##for each fastq file, quantiy the expression.
cd /scratch/project_XXXX/star_genome/sra/Botrytis
files=*fastq
for f in $files
do
output_f=`echo $f |  sed 's/.fastq/_/g'`;
echo $f
echo $output_f
#quant
STAR --genomeDir ../../star_index --readFilesIn $f --outReadsUnmapped unmapped.txt \
--outFileNamePrefix ../../star_align/Botrytis/$output_f --runThreadN 15  --quantMode 
TranscriptomeSAM GeneCounts
done
