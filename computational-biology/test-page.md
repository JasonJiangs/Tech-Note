# Test Page

<div>

<img src="https://s3-us-west-2.amazonaws.com/secure.notion-static.com/373d2201-6ff3-4440-8485-de94aeb4cfa6/Untitled.png" alt="Untitled">

 

<figure><img src="../.gitbook/assets/Untitled.png" alt=""><figcaption></figcaption></figure>

</div>

In general, based on the multiqc result, which integrates will the other fastqc results together to generate a summary. We can see that the quality check mostly fails at `Per Base Sequence Content`, `Per Sequence GC Content`, `Sequence Duplication Levels`, and `Overrepresented Sequence`. Although there are failed cases in the results, we do not decide to perform pre-process in this case, since the decision should be guided by the context of specific experiment and the nature of the data.

In our cases, the context is using RNA-seq to do differential expression analysis and ChIP-seq to explore an unknown transcription factor. So

Quality of sequence reads. You can use a general QC package called “fastqc” to help you check different reads quality measurement.

```bash
# before
fastqc -f fastq -o /home/pfq7pm/test_proj/result/fastqc/ /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/*.fastq.gz
# after
fastqc -f bam -o /scratch/pfq7pm/test_proj/fastqc_after/ /scratch/pfq7pm/test_proj/bowtie2/*.bam
fastqc -f bam -o /scratch/pfq7pm/test_proj/fastqc_after/ /scratch/pfq7pm/test_proj/hisat2/*.bam
```

[fastqc Result Analysis](https://www.notion.so/fastqc-Result-Analysis-49638c14da1a48f8878d3c456fbdcceb)

[fastqc analysis after alignment](https://www.notion.so/fastqc-analysis-after-alignment-54c715e15d7a4946bafaf3ecebbebac2)

#### 2. Mapping

Mappability (percentage of sequence reads that can be mapped to the corresponding genome) of both RNA-seq and ChIP-seq data. Such information is available in the output of the mapping software (bowtie2 for ChIP-seq and hisat2 for RNA-seq)

To extract the mappability information from the output, use the **`samtools flagstat`** command on the resulting SAM file. This command will generate a summary report of the alignment statistics.

```bash
#!/bin/bash
#SBATCH --job-name=samtools
#SBATCH --cpus-per-task=16
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH -A zanglab
#SBATCH --output="/home/pfq7pm/test_proj/log/slurm_log/**flagstat**.out"
#SBATCH --mem=64G

module load gcc/9.2.0
module load samtools

for i in LNCaP_DHT_RNA_rep1 LNCaP_DHT_RNA_rep2 LNCaP_Veh_RNA_rep1 LNCaP_Veh_RNA_rep2
do
samtools flagstat /scratch/pfq7pm/test_proj/hisat2/"$i"_batch.sam
done

for i in LNCaP_DHT_AR_1 LNCaP_DHT_AR_2 LNCaP_Veh_AR_1 LNCaP_Veh_AR_2
do
samtools flagstat /scratch/pfq7pm/test_proj/bowtie2/"$i".sam
done
```

[Output for Samtools flagstat](https://www.notion.so/Output-for-Samtools-flagstat-15bf55733bbd4299ad3537a942da5625)

By order `mapped rate` / `properly pared rate` / `singletons rate`:

1. LNCaP\_DHT\_RNA\_rep1: 96.14% / 90.90% / 2.82%
2. LNCaP\_DHT\_RNA\_rep2: 96.14% / 90.72% / 2.96%
3. LNCaP\_Veh\_RNA\_rep1: 95.85% / 90.00% / 3.17%
4. LNCaP\_Veh\_RNA\_rep2: 96.03% / 90.53% / 2.92%
5. LNCaP\_DHT\_AR\_1: 86.67% / NA / NA
6. LNCaP\_DHT\_AR\_2: 78.08% / NA / NA
7. LNCaP\_Veh\_AR\_1: 88.54% / NA / NA
8. LNCaP\_Veh\_AR\_2: 88.95% / NA / NA

**Check Log of Hisat2 for RNA-seq data:**

**In order of** LNCaP\_DHT\_RNA\_rep1 LNCaP\_DHT\_RNA\_rep2 LNCaP\_Veh\_RNA\_rep1 LNCaP\_Veh\_RNA\_rep2

```
Time loading forward index: 00:00:17
Time loading reference: 00:00:02
Multiseed full-index search: 00:37:09
190669387 reads; of these:
  190669387 (100.00%) were paired; of these:
    17348542 (9.10%) aligned concordantly 0 times
    145862611 (76.50%) aligned concordantly exactly 1 time
    27458234 (14.40%) aligned concordantly >1 times
    ----
    17348542 pairs aligned concordantly 0 times; of these:
      940160 (5.42%) aligned discordantly 1 time
    ----
    16408382 pairs aligned 0 times concordantly or discordantly; of these:
      32816764 mates make up the pairs; of these:
        19244958 (58.64%) aligned 0 times
        10900759 (33.22%) aligned exactly 1 time
        2671047 (8.14%) aligned >1 times
94.95% overall alignment rate
Time searching: 00:37:12
Overall time: 00:37:31
Time loading forward index: 00:00:03
Time loading reference: 00:00:01
Multiseed full-index search: 00:33:58
186962133 reads; of these:
  186962133 (100.00%) were paired; of these:
    17341035 (9.28%) aligned concordantly 0 times
    144137299 (77.09%) aligned concordantly exactly 1 time
    25483799 (13.63%) aligned concordantly >1 times
    ----
    17341035 pairs aligned concordantly 0 times; of these:
      921762 (5.32%) aligned discordantly 1 time
    ----
    16419273 pairs aligned 0 times concordantly or discordantly; of these:
      32838546 mates make up the pairs; of these:
        18940173 (57.68%) aligned 0 times
        11161265 (33.99%) aligned exactly 1 time
        2737108 (8.34%) aligned >1 times
94.93% overall alignment rate
Time searching: 00:34:00
Overall time: 00:34:03
Time loading forward index: 00:00:04
Time loading reference: 00:00:00
Multiseed full-index search: 00:33:35
175175817 reads; of these:
  175175817 (100.00%) were paired; of these:
    17509310 (10.00%) aligned concordantly 0 times
    130857051 (74.70%) aligned concordantly exactly 1 time
    26809456 (15.30%) aligned concordantly >1 times
    ----
    17509310 pairs aligned concordantly 0 times; of these:
      884661 (5.05%) aligned discordantly 1 time
    ----
    16624649 pairs aligned 0 times concordantly or discordantly; of these:
      33249298 mates make up the pairs; of these:
        19265626 (57.94%) aligned 0 times
        10977090 (33.01%) aligned exactly 1 time
        3006582 (9.04%) aligned >1 times
94.50% overall alignment rate
Time searching: 00:33:36
Overall time: 00:33:40
Time loading forward index: 00:00:03
Time loading reference: 00:00:01
Multiseed full-index search: 00:29:56
179684279 reads; of these:
  179684279 (100.00%) were paired; of these:
    17009142 (9.47%) aligned concordantly 0 times
    136308009 (75.86%) aligned concordantly exactly 1 time
    26367128 (14.67%) aligned concordantly >1 times
    ----
    17009142 pairs aligned concordantly 0 times; of these:
      907204 (5.33%) aligned discordantly 1 time
    ----
    16101938 pairs aligned 0 times concordantly or discordantly; of these:
      32203876 mates make up the pairs; of these:
        18809388 (58.41%) aligned 0 times
        10528399 (32.69%) aligned exactly 1 time
        2866089 (8.90%) aligned >1 times
94.77% overall alignment rate
Time searching: 00:29:58
Overall time: 00:30:01
```

The alignment rate for each:

1. LNCaP\_DHT\_RNA\_rep1: 94.95%
2. LNCaP\_DHT\_RNA\_rep2: 94.93%
3. LNCaP\_Veh\_RNA\_rep1: 94.50%
4. LNCaP\_Veh\_RNA\_rep2: 94.77%

**Check Log of Bowtie2 for ChIP-seq data:**

```
52644155 reads; of these:
  52644155 (100.00%) were unpaired; of these:
    7015750 (13.33%) aligned 0 times
    31021896 (58.93%) aligned exactly 1 time
    14606509 (27.75%) aligned >1 times
86.67% overall alignment rate
50161537 reads; of these:
  50161537 (100.00%) were unpaired; of these:
    10997425 (21.92%) aligned 0 times
    26652956 (53.13%) aligned exactly 1 time
    12511156 (24.94%) aligned >1 times
78.08% overall alignment rate
46205142 reads; of these:
  46205142 (100.00%) were unpaired; of these:
    5296219 (11.46%) aligned 0 times
    27746572 (60.05%) aligned exactly 1 time
    13162351 (28.49%) aligned >1 times
88.54% overall alignment rate
47745017 reads; of these:
  47745017 (100.00%) were unpaired; of these:
    5274770 (11.05%) aligned 0 times
    28739080 (60.19%) aligned exactly 1 time
    13731167 (28.76%) aligned >1 times
88.95% overall alignment rate
```

The alignment rate of each:

1. LNCaP\_DHT\_AR\_1: 86.67%
2. LNCaP\_DHT\_AR\_2: 78.08%
3. LNCaP\_Veh\_AR\_1: 88.54%
4. LNCaP\_Veh\_AR\_2: 88.95%

Comparison between`samtools flagstat` and `Hisat2 / Bowtie2` , the mapping rates are different in RNA-seq data, possible reasons:

1. **Differences in how they classify and count reads:** Hisat2's reported alignment rate includes both reads that mapped uniquely and those that mapped to multiple locations. **`samtools flagstat`**, on the other hand, also provides statistics on secondary alignments and supplementary alignments. These are alignments for reads that have multiple primary alignment records.
2. **Filtering of alignments:** Hisat2's reported alignment rate may consider all alignments found, including those that were of low quality or that didn't meet certain thresholds set during the alignment process. **`samtools flagstat`** works on the output BAM file, which may have had low-quality alignments or certain types of alignments filtered out.
3. **Handling of paired-end data:** Hisat2 and **`samtools flagstat`** may handle paired-end reads differently, especially in the case of orphan reads (where one read of a pair aligns and the other doesn't) or improperly paired reads.

#### 3. Replicate correlation

Replicate correlation (reproducibility). Replicate correlation is one of the most important QC measurements. You are required to design different measurements to check this for RNA-seq and ChIP-seq data (for example, expression correlation for RNA-seq and peak overlap for ChIP-seq).

**Check for RNA-seq (Replicate correlation), do gene and transcript sample respectively (gene\_count\_matrix / transcript\_count\_matrix):**

1.  Use R to generate correlation coefficient between each part of replicates. Do normalization and take Spearman correlation (rank-based measures) as an example.

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/8eb58759-01eb-45a5-b9d8-91e0c3754679/Untitled.png)

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/5f76a46c-80ec-4287-aa07-efc935e6c6ee/Untitled.png)
2.  Visualization the Correlation with `heatmap()` and `pheatmap()`.

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/a081d4ca-402f-42de-b685-06956298a8d9/Untitled.png)

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/165e6c07-f099-42e8-b258-1da391fafba3/Untitled.png)

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/4e226697-df78-4cbd-9259-7bb88299004c/Untitled.png)

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/d56d1f10-6641-4957-9ee1-d613c5e731f0/Untitled.png)

1.  Create a scatter plot for each pair of replicates to visualize the correlation with `ggplot2`.

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/4b597ff0-66d7-4c47-a43e-2fb6105759d2/Untitled.png)

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/e6fe6ed4-90b8-4be7-99ac-45d97f1478e3/Untitled.png)

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/c842be65-5854-4253-9cd5-de3da322e58e/Untitled.png)

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/d98f395e-74ce-4298-bb5d-2e2ef813e16d/Untitled.png)

> Between gene\_count\_matrix and transcript\_count\_matrix:
>
> 1. **Gene Count Matrix**: In a gene count matrix, the rows represent genes, and the values in the matrix represent the number of sequencing reads that have been mapped to each gene. This gives an overall picture of the expression level of each gene, but it doesn't provide information about the expression of individual transcripts (also called isoforms) that may be produced from each gene.
> 2. **Transcript Count Matrix**: In a transcript count matrix, the rows represent individual transcripts rather than whole genes. This allows for the examination of alternative splicing events, since it can show the expression level of each individual transcript that is produced from a gene. Therefore, a transcript count matrix provides a more detailed view of gene expression, at the cost of being more complex to analyze.
>
> When performing replicate correlation analysis as described above, we are measuring the similarity in expression patterns between different samples (typically biological replicates) in RNA-seq data. This can give several valuable pieces of information:
>
> 1. **Quality Control**: If the replicates of the same condition are not highly correlated, it may indicate a problem with the experiment or the sequencing. For example, there could have been a technical error, or the biological samples might not actually be as similar as assumed.
> 2. **Reproducibility**: High correlation between replicates indicates that your results are reproducible. This is important for ensuring the validity of your findings.
> 3. **Expression Patterns**: By examining the correlation between different conditions (rather than replicates), you can gain insight into the similarities and differences in gene expression patterns between those conditions.

**Check for ChIP-seq Data (Peak overlap and correlation of Peak Intensity)**

1.  **Peak overlap**: The overlap of detected peaks between replicates can be used to measure reproducibility. Overlaps can be quantified using the Jaccard index or other set similarity metrics.

    1. **Call Peak:** Using MACS2 to identify regions of the genome that have been enriched in ChIP-seq experiment. The BED files in the output contains the genomic coordinates of each peak.
    2. **Identify overlapping peaks**: `bedtool intersect` allows to find overlapping regions between two sets of genomic coordinates. For example, find peaks that are represent in both replicate 1 and replicate 2.

    ```bash
    bedtools intersect -a replicate1_peaks.bed -b replicate2_peaks.bed -wa > overlapping_peaks.bed
    ```

    The "-wa" flag tells Bedtools to write the original entry in file A (specified by "-a") for each overlap with B. Without this flag, Bedtools intersect will only output the portions of the entries in A that overlap with B.

    So, if you have a region in file A (replicate1\_peaks.bed) that partially overlaps with a region in file B (replicate2\_peaks.bed), using the "-wa" flag will result in the output of the entire original region from file A, not just the overlapping portion.

    In other words, "-wa" ensures that the whole peak from the "-a" file is written to the output when an overlap with a peak in the "-b" file is detected. This can be helpful in understanding the complete genomic region of interest from file A that is also present in file B.

    c. **Quantify Overlap**: Once the set of overlapping peaks are found, we can calculate the Jaccard Coefficient with bedtools.

    \$$ J(A,B) = \frac{|A \cap B|}{|A \cup B|} \$$

    ```bash
    -bash-4.2$bedtools jaccard -a LNCaP_DHT_AR_1_PK_summits.bed -b LNCaP_DHT_AR_2_PK_summits.bed 
    intersection   union     jaccard    n_intersections
        89         18678    0.00476496      89

    intersection   union     jaccard    n_intersections
      1166839     2371753    0.491973       5400
    ```

    The jaccard index is approximately 0.0048, which is quite low and suggest that there is very little overlap between peak sets from the two replicates.

    > **Issue**: the calculated jaccard coefficient should be approaching to 1 to suggest that the replicates are similar, which means they have a large overlap of peaks. However, the results are extremely low. The result indicates several possibilities: 1. High biological or technical variability between replicates; 2. low quality of one or both datasets; 3. inappropriate peak calling parameters. **Solution**: The potential problem is that lots of false positive cases are recorded in the multiple hypothesis testing that MACS2 performs to call peaks. So, a lower q-value means a stricter threshold, which will result in fewer peaks being called, but with higher confidence. Make an adjustment to a lower q = 0.01 instead of the default value 0.05.
    >
    > intersection union jaccard n\_intersections 74 14567 0.00507998 74
    >
    > Still unsolve the problem, we then switch to use narrowPeak file for the calculation.
    >
    > The **`NAME_summits.bed`** and **`NAME_peaks.narrowPeak`** files contain slightly different information about the peaks called by MACS2.
    >
    > 1. **`NAME_summits.bed`**: This file contains the precise location of the peak summits, which are the single base pair positions within each peak region that have the highest pileup of reads. This file will only contain one summit per peak region, even if the region contains multiple local maxima.
    > 2. **`NAME_peaks.narrowPeak`**: This file contains broader peak regions, which are the areas around the peak summits that have a significant pileup of reads. Each entry in this file represents a contiguous region of the genome where the read pileup is above the threshold determined by MACS2.
    >
    > The discrepancy in the Jaccard index between these two files could be due to the differences in their contents.
    >
    > The **`NAME_summits.bed`** file represents the very highest-confidence points within each peak, so there will be fewer entries in this file compared to the **`NAME_peaks.narrowPeak`** file. If the peak summits in the two replicates are not exactly aligned, which can often be the case due to the inherent variability and noise in ChIP-seq data, the Jaccard index calculated from the **`NAME_summits.bed`** files might be low.
    >
    > On the other hand, the **`NAME_peaks.narrowPeak`** file contains broader regions that have a significant pileup of reads, so there is more opportunity for overlap between the replicates, which can result in a higher Jaccard index.
    >
    > Therefore, it's not uncommon to see a higher Jaccard index when using **`NAME_peaks.narrowPeak`** files compared to **`NAME_summits.bed`** files.
    >
    > intersection union jaccard n\_intersections 1166839 2371753 0.491973 5400 In this case, it is a pretty high Jaccard coefficient value which refers a good identical characteristics between replicates. This means that about 49.2% of the total regions identified in the two replicates overlap, suggesting a substantial agreement between the two datasets. Although there isn't a universally accepted threshold for a "good" Jaccard Index in this context, a higher value generally indicates better reproducibility. Given that the index ranges from 0 (no overlap) to 1 (complete overlap), a value of 0.492 can be considered relatively good, indicating that the peak sets from your two replicates are largely consistent with each other.
2.  **Correlation of Peak Intensities**: Similar to expression correlation in RNA-seq, you can calculate the correlation of peak intensities between replicates.

    Pearson Correlation: 0.846`LNCaP_DHT_AR_2_PK_summits.bed`

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/4ade7541-e186-4789-976d-ea8b9e2de5dc/Untitled.png)

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/b4e8f227-3274-49d2-88ac-e81b66e0ba68/Untitled.png)

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/eae0c329-f413-4369-b06c-a33a8b2f7d55/Untitled.png)

#### 4. Duplicate rate in ChIP-seq

Find this information in your macs2 output and try to understand what this measurement can tell us.

In this case, using Python to operate on the BED files, and both duplicate rate is 0%.

> Understanding the duplicate rate can provide insights into the quality of your ChIP-seq experiment. Low duplicate rates suggest good library complexity, meaning there is a wide variety of different DNA fragments in your sample, which is generally a good thing. High duplicate rates could indicate over-amplification during PCR or a low amount of starting material, which might limit the number of different DNA fragments that can be sequenced. This could potentially lead to biased results in downstream analyses. However, it's worth noting that some level of duplicates in ChIP-seq data is expected due to the nature of protein-DNA interactions and the way ChIP-seq works.

#### 5. Other QC measurements

Other potential QC measurement you think is important for this study.

**For RNA-seq:**

1. **Then we use PCA to reduce the dimensionality of the data and to visualize the variation between samples. (use `gene_fpkm_matrix`)**
2.  reduce the dimension to 4x2

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/e2e89f61-e4ba-4ebc-8615-088474f0320c/Untitled.png)
3.  use log2 transformation to standardize the data.

    [Why gene expression data should be log2 transformed?](https://www.biostars.org/p/242573/)

    1. **Variance Stabilization**: Gene expression data often have a property where the variance increases with the mean. That is, highly expressed genes have a greater absolute variability than lowly expressed genes. Log transformation can help stabilize this variance across the range of values.
    2. **Normalizing Distribution**: Many statistical techniques assume that the data follow a Gaussian (or "normal") distribution. While gene expression data are typically not normally distributed, the log transformation can make their distribution more symmetric and closer to a Gaussian distribution.
    3. **Scaling Down Large Values**: Log transformation reduces the impact of very large expression values that can disproportionately affect downstream analyses, especially those involving Euclidean distances (like hierarchical clustering or PCA).
    4. **Interpretability**: Log2-transformed values are easier to interpret. For example, a doubling of expression corresponds to an increase of 1 in log2 scale, a quadrupling corresponds to an increase of 2, and so on.
    5. **Comparing Expression Levels**: When comparing expression levels across samples, raw counts can be misleading due to differences in sequencing depth, RNA composition, etc. Log transformation (usually after some form of normalization) can make these comparisons more meaningful.
4.  then get the visualization

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/88cf2ff9-13e2-49e7-8a41-52742d886052/Untitled.png)
5. find the bottom 10 genes loadings and top 10 genes loadings: from the figure shown above, we can see:

* **`LNCaP_DHT_RNA_rep1_batch`** and **`LNCaP_DHT_RNA_rep2_batch`** samples have positive values for PC1. This indicates that these two samples, which are both treated with DHT, are similar to each other in terms of gene expression profiles and are different from the Vehicle-treated samples.
* **`LNCaP_Veh_RNA_rep1_batch`** and **`LNCaP_Veh_RNA_rep2_batch`** samples have negative values for PC1. This means that these two samples, which are both Vehicle-treated, are similar to each other and are different from the DHT-treated samples.
*   However, when take a look at PC2, we can see that take the negative of LNCaP\_Veh\_RNA\_rep1\_batch is almost equal to LNCaP\_Veh\_RNA\_rep2\_batch. To gain more insight into what's driving the variation along PC2, we could look at the loadings of PC2 (the weights used to calculate PC2 from the original variables), see which genes have the highest loadings, and check whether these genes have any biological functions or pathways in common. This might give we some clues about what's causing the symmetric pattern along PC2.

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/71667de2-f1bc-453e-841d-0c5af5580191/Untitled.png)

    ![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/503af5f5-3fdd-43bd-bc6f-74dc50d6a466/Untitled.png)

1. **Clustering: by using clustering algorithms (heirarchical clustering and K-means) to group genes or samples based on their expression patterns, which can be visualized as a heat map.**

> **Issues in Clustering** stack overflow: reduce the size of data by filtering out lowly expressed genes by:

```python
df = pd.read_csv('../data/rnaseq/gene_fpkm_matrix.csv', index_col=0)
# Log transform data
df_log2 = np.log2(df + 1)
# filter out lowly expressed genes
# Calculate average expression level for each gene
average_expression = df_log2.mean(axis=1)
# Set threshold
threshold = 6
# Remove lowly expressed genes
df_filtered = df_log2[average_expression > threshold]
```

**For ChIP-seq**

1. **Irreproducible Discovery Rate (IDR)**: IDR is a specific measure developed for assessing reproducibility in high-throughput experiments like ChIP-seq. It assesses the consistency of rank orders of peaks between replicates.
2.  **Consensus Peak Sets**: Create a consensus set of peaks from multiple replicates and then check the percentage of consensus peaks that are found in each individual replicate.

    For example, if you have 3 replicates, **`len(replicate_data) // 2 + 1`** would be **`3 // 2 + 1 = 1 + 1 = 2`**. This means that a peak must be present in at least 2 of the 3 replicates to be included in the consensus set.

    If you have 4 replicates, **`len(replicate_data) // 2 + 1`** would be **`4 // 2 + 1 = 2 + 1 = 3`**. So a peak must be present in at least 3 of the 4 replicates to be included in the consensus set.

    **Result**: of course 100% for both of the replicates since there are only two.
