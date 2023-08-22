### 甲基化

表观遗传修饰是指对基因组功能的相关修饰，**通过一系列生物学修饰改变基因的活性而不是DNA的核苷酸序列影响基因的表达。**对基因组功能的相关修饰主要包括对**DNA、RNA、以及组蛋白**等的修饰，这些修饰改变了染色质的局部电化学特性和构象，从而调节基因的转录活性。

组蛋白是染色质的重要组成部分，主要分为H2A、H2B、H3、H4，与DNA缠绕可形成核小体。**组蛋白修饰**是在组蛋白N末端的氨基酸残基上发生的共价修饰，主要包括甲基化、乙酰化、泛素化、磷酸化、羰基化、糖基化等。

**DNA甲基化**是表观遗传学领域一个重要的研究方向，真核生物中最常见的DNA修饰非**5-甲基胞嘧啶（5mC）**莫属了，然而在原核生物中最常见的DNA修饰方式则为**N6-methyladenine (6mA)**，即腺嘌呤第6位氮原子甲基化修饰。

**人类是真核生物**。人的参考基因组约30亿碱基，上面不到1%是 CpG位点，可以被甲基化，也就是说不到3千万个 CpG 位点。这些 CpG 位点中，大约 60~80% 被甲基化。主要是而启动子等特殊区域存在 未被甲基化的CpG 岛，那些区域的CpG 位点比较富集。目前研究表明，肿瘤细胞的甲基化水平平均是低于正常细胞的。

亚硫酸盐是甲基化探测的“金标准”，不管是芯片或者甲基化测序，都要先对DNA样品进行亚硫酸盐处理，使非甲基化的C变成U，而甲基化的C保持不变，从而在后续的测序或者杂交后区分出来。

甲基化检测方法多达上百种，哪怕是基于NGS的测序技术，也有BS-Seq、MeDIP-Seq、RRBS-Seq、WGBS、MBD-Seq、SMRT 等，

    WGBS（Whole Genome Bisulfite Sequence）全基因组甲基化测序，利用重亚硫酸氢盐使DNA中未发生甲基化的胞嘧啶(C)脱氨基转变成尿嘧啶(U)，而甲基化的胞嘧啶保持不变，然后通过PCR将U变为A，仅有甲基化的C可以成功保留，最后通过测序就可判断CpG位点是否发生甲基化。
    
    简化甲基化测序 (Reduced representation bisulfite sequencing, RRBS)是一种准确、高效、经济的DNA甲基化研究方法，通过酶切 (Msp I) 富集启动子及CpG岛区域，并进行Bisulfite测序，同时实现DNA甲基化状态检测的高分辨率和测序数据的高利用率。作为一种高性价比的甲基化研究方法，简化甲基化测序在大规模临床样本的研究中具有广泛的应用前景。
    
    oxBS-Seq（oxidativ ebisulfite sequencing）化学氧化结合重亚硫酸盐测序，在哺乳动物中，5mC 可以在TET酶的作用下转换成 5hmC，而传统的 WGBS 方法不能区分 5mC 和 5hmC。oxBS-Seq 技术先将 5hmC 氧化为甲酰基修饰（5fC），进而被重亚硫酸盐处理转换为U，从而排除了 5hmC 的干扰，实现 DNA 甲基化的精准检测。
    
    RRHP（Reduced Representation 5-Hydroxymethylcytosine Profiling）DNA 羟甲基化，即 DNA 甲基化中的5-甲基胞嘧啶易发生氧化形成5-羟甲基胞嘧啶。Msp I 酶切完之后，末端修复，加上接头，然后 β-GT 反应将 5hmC 进行糖基化保护，无法被 Msp I 酶识别；然后进行再次酶切，与常规的 C 以及 5mC 相连的 P5 接头均被切下来，最后仅有含 5hmC 的片段含有两端接头，可以被扩增后建库测序。
    
    MeDIP（Methylated DNA Immunoprecipitation Sequencing） 甲基化 DNA 免疫共沉淀测序，先通过与5mC特异性结合的抗体加入到变性的基因组DNA片段中，富集胞嘧啶甲基化的基因组片断，然后对富集的片段进行高通量测序。



可以使用minfi包的read.metharray.exp函数读取，你前面下载的该数据集的RAW.tar 里面的各个样本的idat文件，就被批量加载到R里面啦。（必须注意的的是， 下面代码里面**GSE68777/idat文件夹**里面有idat文件

也可以使用The Chip Analysis Methylation Pipeline，但是不得不说，每次安装 ChAMP  都得脱一层皮，它的依赖包实在是太多了。其中一个ChAMPdata_2.18.0.tar.gz就是680M文件。首先可以读取R包自带的芯片的idat原始文件，如下

```
# BiocManager::install("ChAMP",ask = F,update = F)
library("ChAMP")
testDir=system.file("extdata",package="ChAMPdata")
myLoad <- champ.load(testDir,arraytype="450K")
```

这个数据集就在ChAMPdata_2.18.0.tar.gz就是680M文件，很可怕！champ.load函数作用于 testDir 这个文件夹，我们就去看看这个 testDir 文件夹里面的内容，发现除了idat芯片原始数据之外，还有一个csv文件。

其实目前甲基化信号值矩阵差异分析的R包非常多，比如  IMA, minfi, methylumi, RnBeads and  wateRmelon，以及我们一直强推的ChAMP！不过我没有那么多时间去一一解读啦，我相信minfi或者champ就够用了。你的目的其实是做完甲基化信号值矩阵有**3个层次的差异分析**：DMP，DMR，DMB：

- DMP代表找出Differential Methylation Probe（差异化CpG位点）
- DMR代表找出Differential Methylation Region（差异化CpG区域）
- Block代表Differential Methylation Block（更大范围的差异化region区域）

methylKit是一个R软件包，用于分析和注释通过高通量亚硫酸氢盐测序获得的DNA甲基化信息。该软件包旨在处理RRBS及其变体的测序数据。但是，如果提供正确的输入格式，它可能会处理全基因组的亚硫酸氢盐测序数据。

脊椎动物的DNA甲基化通常发生在CpG二核苷酸处，但是某些组织中（如胚胎干细胞）非CPG的位点也会被甲基化。DNA甲基化可以作为基因调控的表观遗传控制机制。甲基化可能会阻碍转录因子的结合和/或甲基化的碱基可能会被甲基结合域蛋白结合，而后者会募集染色质重塑因子。在这两种情况下，受调节基因的转录都会发生。另外，异常的DNA甲基化模式已经与许多人类恶性肿瘤相关联并且可以以可预测的方式使用。在恶性组织中，与正常组织相比，DNA被低甲基化或高甲基化。高和低甲基化位点的位置为许多疾病提供了独特的特征。

该软件主要使用`methylKit`对象存储甲基化数据，这些对象在分析过程中可以与`GRanges`对象无缝衔接，作者提供了函数`as(object,"GRanges")`可以将其转化为`GRanges`对象进行其他分析。

In this unit we will demonstrate how to read idat files from the illumina 450K DNA methylation array. We make use the the Bioconductor minfi package [cite 24478339]

RRBS甲基化分析流程

【基本流程】Bismark—Methylkit—Genomation—GO/KEGG

 fastqc 质控

trim-galore 去接头

二、比对基因组

    安装 bismark

从github下载

git clone https://github.com/FelixKrueger/Bismark.git
cd Bismark
make

添加环境变量

echo export PATH = .../software/Bismark:$PATH >> ~/.bashrc

bismark --version


          Bismark - Bisulfite Mapper and Methylation Caller.
    
                       Bismark Version: v0.22.1
        Copyright 2010-19 Felix Krueger, Babraham Bioinformatics
              www.bioinformatics.babraham.ac.uk/projects/
                https://github.com/FelixKrueger/Bismark

安装成功

- 构建参考基因组

cd Genome
mkdir hg19
cd hg19
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz &

构建甲基化比对基因组文件

bismark_genome_preparation hg19

- 排序

首先按照比对的基因组坐标进行排序

```shell
sambamba sort \
-m 8GB \
--tmpdir tmp \
-t 5 \
-o output_sorted.bam input_file.bam
```

- 去重

去除多重比对、重复、未比对上的reads

```shell
sambamba markdup \
--overflow-list-size 1000000 \
--tmpdir tmp \
-t 5 \
input_sorted.bam \
output.bam \
2> MarkDup_report.txt
```

四、提取甲基化信息

MethylDackel extract \
path_to_genome_file/Bismark/hg19/hg19.fa \
demo.bam \
--opref path_to_save_folder/sample

五、将CpG位点保存为`GR`文件

由于测序是区分正负链的，而在分析的时候不区分，所以需要合并正负链的信息。
 还需要将与基因组CpG位点不匹配的位点去除，因此需要load一个全基因组CpG位点文件。

寻找差异甲基化区域（DMR）-- DSS 包

`DMRfinder` 是一款用于WGBS的C位点提取、鉴定以及DMR分析的软件，具体过程是：Bismark 软件比对后，用 DMRfinder 提取CpG位点，然后识别差异甲基化区域（call DMR）并与DSS包的call DMR结果做比较。

RRBSdata

An RRBS data set with 12 samples and 10,000 simulated DMRs

methylSig is a method for testing for differential methylated cytosines (DMCs) or regions (DMRs) in whole-genome bisulfite sequencing (WGBS) or reduced representation bisulfite sequencing (RRBS) experiments. methylSig uses a beta-binomial model to test for significant differences between groups of samples. Several options exist for either site-specific or sliding window tests, combining strands, and for variance estimation.

Methylation calls output by either MethylDackel or Bismark can be read by the bsseq::read.bismark() function from the bsseq R/Bioconductor package.

This function accepts bedGraphs from MethylDackel and either the coverage or genome-wide cytosine reports from Bismark. Options to consider when reading data are:

    colData, a data.frame or DataFrame whose rows are samples and columns are phenotype data. The row ordering should match the ordering of files in files. This matrix will be needed for downstream differential methylation testing.
    strandCollapse, a logical (TRUE/FALSE) indicating whether or not to collapse +/- CpG data onto the + strand. Note, this can only be TRUE when the input type is the genome-wide cytosine report from Bismark. MethylDackel has an option to destrand data when methylation calls are made so that the output is already destranded. In this case, strandCollapse should be FALSE.

In this manual, we will show how to use the methylKit package. methylKit is an R package for analysis and annotation of DNA methylation information obtained by high-throughput bisulfite sequencing. The package is designed to deal with sequencing data from RRBS and its variants. But, it can potentially handle whole-genome bisulfite sequencing data if proper input format is provided.

We start by reading in the methylation call data from bisulfite sequencing with methRead function. Reading in the data this way will return a methylRawList object which stores methylation information per sample for each covered base. By default methRead requires a minimum coverage of 10 reads per base to ensure good quality of the data and a high confidence methylation percentage.

The methylation call files are basically text files that contain percent methylation score per base. Such input files may be obtained from AMP pipeline developed for aligning RRBS reads or from processBismarkAln function. However, “cytosineReport” and “coverage” files from Bismark aligner can be read in to methylKit as well.

https://github.com/ben-laufer/DMRichR



