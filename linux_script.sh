Linux命令行
##服务器管理
conda creat -n #创建conda环境
mount -o remount /dev/mapper/ME4084-01  ##挂载硬盘
edquota liujianfeng #修改磁盘配额  edquota -p quser1 quser2 
useradd -d "/storage2/zhangfan" -m zhangfan   #创建账户 
passwd zhangfan #设置密码
usermod -d /home/Data/data7/data2_user/youzhen -m youzhen  #给用户更换文件夹
chown root:root testfile  ##一次性更改用户和组
usermod -d /data0/home/mypic -m mypic ##移动用户位置
/home/Software/anaconda3/bin/conda init bash #小服务器激活conda环境   conda activate base
cat /proc/cpuinfo| grep "processor"| wc -l  #查看线程数
vi /opt/maui/maui.cfg    service  maui restart  ##修改使用核数
#服务器联网
vncserver :2c
vncserver -kill :2
#载入环境
module load anaconda
####检查服务器问题
fuser -cu data1 ##查看当前文件夹内的进程
umount -v data1  ##umount命令用于卸载已经加载的文件系统
mount /dev/sdl data2 ##挂载文件系统
ifconfig  ###查找服务器ip
vncserver :2 ##服务器联网 vncserver -kill :2
losetup -d /dev/loop1 ##卸载设备
fuser -ck /dev ##杀死指定文件内所有进程
ps auxf|grep 'gatk.sh'|grep -v grep|awk '{print $2}'|xargs kill -9   ##杀死脚本 之后再kill ID
ps -aux | sort -k4,4nr | head -n 8    ###展示最消耗内存的前8个程序
split -l 50 list -d -a 2 split_file   ###分割文件，每个文件50行，分割文件名为split_filex
find . ! -name 'file.txt' -type f -exec rm -f {} +    ###只保留文件夹中的file.txt

#Linux系统下比较两个文件并删除相同部分
awk '{print $0}' file1 file2 |sort|uniq -u
for i in `ls`; do mv -f $i `echo $i | sed 's/reheadered/reheader/'`; done ##文件批量改名
du -h --max-depth=1   ##查看文件夹大小
find . -name "*" -type f -size 72c  | xargs -n 1 rm -f   ##删除特定大小的文件
###寻找多个文件的交集（例子中为6）
find /path/to/files -name "*bed" -type f -exec bash -c ' sort $1 | uniq ' _ {} \; | sort | uniq -c | awk '{if($1==6){print $2}}'
#对blast以后的结果过滤
cat file.txt |sort -k 12nr |awk '!a[$1]++{print $0}' >file1.txt
##txt文件转为csv
cat test.txt | tr "old value" "new value" >test.csv
###awk删除列
cat file |awk ' { $5=null;print $0 }'

#从组装结果中提取未注释到的转录子
cat merge.gtf |grep -vi ensue |awk '$3~/transcript/{print $0}' > unannotation.gtf 
gffread -w unannotation.fa -g Ref.fa unannotation.gtf
OR
fasta_formatter -i out.fa -o sss.fa -w 1000000000
cat sss.fa |grep 'M' -A1 > unannotation.fa

##fastq转fasta
sed -n '1~4s/^@/>/p;2~4p' INFILE.fastq > OUTFILE.fasta

##显示fasta序列长度
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"%":$0 }'  data.fa >data2.fa
awk -F"%" '{print $1"\t"length($2)}'  data2.fa >data3.txt

# 1. 如果你要用前两个，就
# module load Anaconda/3.7.2
# source activate bioconda3
# 就可以使用了
# 2. 如果你要使用CPAT，就先输入 workon python_3.7.2，就可以使用了
# !!!   如果你提交PBS脚本，在脚本文件里不要写conda，module，workon这样的命令，通过前面的命令获得软件的绝对路径，脚本里直接写软件的绝对路径

# 获取转录子长度
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' sss.fa |awk '{print $1"\t"length($3)}' > long.txt

#gtf处理—提取并计数有多少类feature 
cat Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.96.gtf |grep -v ^# |awk '{print $3}' |sort |uniq -c
#Gtf处理—筛选特定行—第一列为染色体+XY的行
cat ***.gtf | awk '$1 ~ /[0~9]/ || $1 ~ /[X|Y]/ '  > chr.txt 
#Gtf处理—查看每条染色体的基因数
cat ***.gtf | awk '$1 ~ /[0~9]/ || $1 ~ /[X|Y]/' |awk '$3 == "gene"' |awk '{print  $1 }'  |sort |uniq -c 
#gtf处理—计算基因的长度
cat Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.96.gtf |awk '$3=="gene" {split($10,x,";");name = x[1];gsub("\"","", name);print name,$5-$4+1}' > len_gene_CAU.txt
#gtf处理-计算CDS长度
cat Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.96.gtf |awk '$3=="CDS" {len=$5-$4 +1;size += len;print $10, size}' > len_cds_CAU.txt
# gtf处理-根据基因名列表提取gtf文件
gene=($(cat gene.txt));for gene in ${gene[@]};do grep "\"${gene}\"" 123.gtf  >> result.gtf ;done


###将文件 BLM.txt 分成若干个小文件，每个文件2482行(-l 2482)，文件前缀为BLM_ ，系数不是字母而是数字（-d），后缀系数为四位数（-a 4）
split -l 50 ../BLM/BLM.txt -d -a 4 BLM_


#######fst 脚本#########
bedtools merge -i layer_low.bed -c 5 -o min > layer_low_merge.bed
bedtools merge -i layer_high.bed -c 5 -o max > layer_high_merge.bed
bedtools intersect -wa -wb -a pasa.bed -b pekin_low_merge.bed pekin_high_merge.bed layer_low_merge.bed layer_high_merge.bed -sorted -filenames > result.txt
#获取基因
awk '$8=="gene"' result.txt > gene.txt

cat gene.txt |awk '{print $4,$11,$15}' > genename.txt
cat genename.txt |awk '{print $1}' > name.txt
sed -i "s/.TU./.model./g" name.txt
sed -i 's/$/&\\>/g' name.txt
cat Swissprot.result |grep -f name.txt > anno.txt



##vcftools 筛选特定范围内的snp
vcftools --vcf genotype_id.vcf --chr A1 --from-bp 33000000 --to-bp 33781370 --recode --out A1_analysis_pos_bp
##vcftools filter
vcftools --gzvcf greenegg_snp.vcf.gz --max-missing 0.5 --maf 0.05 --minDP 10  --minQ 30 --recode --recode-INFO-all --out green_maf
##vcftools提取个体
vcftools --gzvcf in.vcf.gz --recode  --keep id  --out out.vcf
##vcf去除indel
vcftools --remove-indels --recode --recode-INFO-all --vcf raw.vcf --stdout >raw.snp.vcf
##vcf只保留indel
vcftools --keep-only-indels  --recode --recode-INFO-all --vcf raw.vcf --stdout >raw.indel.vcf
##提取特定位点的vcf文件
vcftools --gzvcf 627_rm_contig.vcf.recode.vcf.gz --positions ../anno/HIGH_position_formatted --recode --out specific_position.vcf
##vcf 2 bed
sed -e 's/chr//' file.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2,$2+1,$4"/"$5,"+"}}'  ##snp
sed -e 's/Hic_asm_//' input.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2,$2+length($5),$4"/"$5,"+"}}' > test.bed  ##sv
vcf2bed < input.vcf > output.vcf

##使用过vep注释
module load ensembl-vep/98.2  ##学院
source activate vep   ##小服务器
bedtools sort -i unsort.gff > sort.gffcc
bgzip sort.gff
tabix sort.gff.gz
nohup vep --custom CAU_duck_sort.gff.gz,,gff --fasta /home/Data/data2/zhangfan/Genome/Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.dna.toplevel.fa --gff ../CAU_duck_sort.gff.gz --input_file elovl6l.vcf -o test --fork 2 &
condaa annotation  ###大服务器
vep --custom ~/zhangfan/Genome/Anas104.sort.gtf.gz,,gtf --fasta /storage-01/poultrylab1/zhangfan/Genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa --gtf ~/zhangfan/Genome/Anas104.sort.gtf.gz --input_file sig.vcf.recode.vcf -o vep --fork 50 &

####提取基因组指定位置序列
bedtools getfasta  -fi <fasta> -bed

##bcftools 过滤无突变位点
bcftools filter -e 'N_PASS(GT="RR")  > 0' mallard20.recode.vcf > test
bcftools filter -e 'N_PASS(GT="RR")  > 1' mallard20.recode.vcf > test ##多于1个有变异则剔除
##提取特定bam片段
samtools view -hb  green.bam Hic_asm_1:104398442-104488442  > green_target.region.bam

#找两个snp的交集
bcftools isec dir/0000.vcf.gz yeya2yeya_samtools_filter.vcf.gz -p dir1

##找到当前文件夹下所有以vcf结尾的文件名
find . -type f -name "*.vcf" > doc.txt
find . -name *.txt -exec cat {} \;> MD5_all_files

##计算覆盖度
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' ../../../1_duck_eggshell_color/yeya_header.fa |awk '{print $1"\t"length($2)}' > genome.fasta.tag
bedtools makewindows -g genome.fasta.tag -w 100000  -i winnum>genome.windows_100k.bed
bedtools coverage -mean -sorted -bed -g genome.fasta.tag  -a genome.windows_100k.bed -b 1_GNAQ.bam > 1.100kbin.txt

##获得第四列小于0.00001的行数
cat feed_days_map | awk -F" " '$4< 0.00001{print $0}'|wc -l

###vlookup 脚本
#!/bin/bash
awk 'NR==FNR{a[$1]=$0;next}NR>FNR{ print $0"\t" (($1 in a)?a[$1]:"NA")}' $1 $2
./vlookup aa bb > result ##bb是list，aa是待比对的文件

awk 'NR==FNR{a[$1]=$2}NR>FNR{print $0,a[$1]}' 2.txt  1.txt

###查找多个文件的交集
awk ' FNR == 1 { b++ } { a[$0]++ } END { for (i in a) { if (a[i] == b) { print i } } } ' 1.txt 2.txt 3.txt

##计算LDdecay
#单群体：产生“Out.Prefix.png” and “Out.Prefix.pdf”
PopLDdecay -InVCF ALLchr.vcf.gz -OutStat LDDecay.stat.gz 
Plot_OnePop.pl -inFile LDDecay.stat.gz -output Out.Prefix
#多个群体
PopLDdecay -InVCF In.vcf.gz -OutStat wild.stat.gz -SubPop wildName.list 
PopLDdecay -InVCF In.vcf.gz -OutStat cul.stat.gz -SubPop culName.list
Plot_MultiPop.pl -inList multi.list -output OutputPrefix

##R语言提取子集
JFK <- hflights[which(hflights$Dest == 'JFK'), c('TaxiIn','TaxiOut')]

###picard 构建参考基因组interval
java -jar /storage-01/poultrylab1/liguangsheng/SVE/src/picard-tools-2.5.0/picard.jar  ScatterIntervalsByNs   R=Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa  OT=BOTH   O=output.interval_list

###gcta 计算kinship
gcta64 --bfile test --make-grm-gz --make-grm-alg 1 --out kinship --autosome-num 40 --thread-num 10



##Rstudio-severse
systemctl status rstudio-server  ##查看rstudio运行状态
rstudio-server stop
rstudio-server start
##rstudio配置文件
 /etc/rstudio/rserver.conf
 /etc/rstudio/rsession.conf
# 开放端口
iptables -I INPUT -p tcp --dport 8787 -j ACCEPT
iptables -I OUTPUT -p udp --dport 8787 -j ACCEPT
# 保存防火墙规则, /etc/iptables/rules.v4
sudo netfilter-persistent save

##提取包含特定碱基的reads
seqkit grep --by-seq --max-mismatch 1 --pattern "ATCGAAG" test.fq
