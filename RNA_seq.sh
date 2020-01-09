#RNA分析的一整套流程
#启动RNA分析环境
source activate rna
##shell中使用到for循环时，是通过IFS同样也是来定义分隔符。
IFS=$'\n' 
##原始数据文件夹
RAWDATA=/home/sda1/workspace/MetalRNA/00.rawdata/
ls 00.rawdata/*1.fq.gz > 1
ls 00.rawdata/*2.fq.gz > 2
paste 1 2 > sample.txt
mkdir 01.fastqc
step1:check quality of sequence reads
for i in `ls 00.rawdata/*`
do
    {
    nohup fastqc $i -o 01.fastqc/ -t 30
    } &
done && wait

multiqc 01.fastqc/ -o 01.fastqc/

sleep 5s

##step2:filter bad quality reads and remove adaptors
mkdir 02.cleandata
DIR=02.cleandata/
TruSeq3=/home/zeng/miniconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa

for id in `cat sample.txt`
do 
         fq1=$(echo ${id} | cut -f1)
         fq2=$(echo ${id} | cut -f2)
         out1=$(basename ${fq1})
         out2=$(basename ${fq2})
         {
         nohup trimmomatic\
                 PE \
                 $fq1 $fq2 \
                 ${DIR}${out1}_paired.fq.gz ${DIR}${out1}_unpaired.fq.gz ${DIR}${out2}_paired.fq.gz ${DIR}${out2}_unpaired.fq.gz \
                 ILLUMINACLIP:$TruSeq3:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36
         } & 
done && wait
echo "finish filter"

sleep 5s

######mapgenome
mkdir /home/sda1/workspace/MetalRNA/04.mapgenome
genome=/home/sda1/workspace/MetalRNA/03.reference/Hordeum_vulgare.IBSC_v2.dna.toplevel.fa
genome_prefix=/home/sda1/workspace/MetalRNA/03.reference/Hordeum_vulgare.IBSC_v2
cdir=/home/sda1/workspace/MetalRNA/04.mapgenome/
hisat2-build -p 70 $genome $genome_prefix 1>/home/sda1/workspace/MetalRNA/03.reference/hisat2-build.log 2>/home/sda1/workspace/MetalRNA/03.reference/hisat2-build.err
cd /home/sda1/workspace/MetalRNA/04.mapgenome/
for id in `cat clean.txt`
do
         cid=$(echo $id | cut -f1)
         cfq1=$(echo $id | cut -f2)
         cfq2=$(echo $id | cut -f3)
         {
         nohup hisat2 --new-summary --mm -p 4 -x $genome_prefix -1 ${cfq1} -2 ${cfq2} -S ${cdir}${cid}.sam 2> ${cdir}${cid}.err
         } &
done && wait
sleep 5s
##排序

for id in `ls ${cdir}*.sam`
do      
	samtools sort -m 1G -o ${id}.bam -@ 2 ${id} &&
        { 
	samtools index -c ${id}.bam
        } &
done && wait
#
sleep 5s
###转录组reads
#
gtf1=/home/sda1/workspace/MetalRNA/03.reference/Hordeum_vulgare.IBSC_v2.44.gtf
threads=60
cd /home/sda1/workspace/MetalRNA/04.mapgenome/ 
featureCounts -T $threads -a $gtf1 -o counts.txt -t exon *.bam
#
###转录组矩阵定量ls
text=/home/sda1/workspace/MetalRNA/clean.txt
matrix1=/home/sda1/workspace/Cai_drought/04.mapgenome/step3.featureCounts2matrix.R
input1=$(wc -l $text | awk '{print $1}')
input2="counts.txt"
Rscript $matrix1 $input1 $input2 
#
###标准化
#
TRINITY_HOME=/home/zeng/miniconda3/envs/rna/opt/trinity-2.6.6
genes_TPM_EXPR_matrix="genes.TPM.not_cross_norm"
$TRINITY_HOME/util/support_scripts/run_TMM_scale_matrix.pl --matrix $genes_TPM_EXPR_matrix > genes.TMM.EXPR.matrix
#
#
#
#DEG
genes_counts_matrix=/home/sda1/workspace/MetalRNA/04.mapgenome/genes.counts.matrix
samples_file=/home/sda1/workspace/MetalRNA/sample10.txt
TRINITY_HOME=/home/zeng/miniconda3/envs/rna/opt/trinity-2.6.6
contras_file=/home/sda1/workspace/MetalRNA/sample11.txt
perl $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
         --matrix $genes_counts_matrix \
         --method DESeq2 \
         --samples_file $samples_file \
         --contrasts $contras_file
sleep 2s

相关性
mkdir /home/sda1/workspace/MetalRNA/05.Ptr
cd /home/sda1/workspace/MetalRNA/05.Ptr/
TRINITY_HOME=/home/zeng/miniconda3/envs/rna/opt/trinity-2.6.6/
PtR=$TRINITY_HOME/Analysis/DifferentialExpression/PtR
samples_file=/home/sda1/workspace/MetalRNA/sample10.txt
genes_TMM_EXPR_matrix=/home/sda1/workspace/MetalRNA/04.mapgenome/genes.TMM.EXPR.matrix
gene=AT3G19900

$PtR --matrix $genes_TMM_EXPR_matrix --samples $samples_file --compare_replicates
$PtR --matrix $genes_TMM_EXPR_matrix --samples $samples_file --sample_cor_matrix
$PtR --matrix $genes_TMM_EXPR_matrix --samples $samples_file --indiv_gene_cor $gene --top_cor_gene_count 10 --min_rowSums 3
$PtR --matrix $genes_TMM_EXPR_matrix --samples $samples_file --prin_comp 3 --add_prin_comp_heatmaps 10



#需要注意的几项
#（1）sample10.txt的格式一定要注意，注意名称和genes.counts.matrix中的名称匹配
#（2）genes.counts.matrix是点不是__横线
#（3）最后要注意 EOF 还是少用吧，会出各种bug。
