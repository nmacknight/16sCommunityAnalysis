# Qiime2
### Adapted to the White Plague Microbiome Dataset :microbe:	

**Metadata**

The metadata file includes the sample ID and the Golay Barcode at its most fundamental. Some examples of Metadata files have had other details that is relevant to the experimental design but thats not necessary to perform the Qiime2 Pipeline.

**Sequencer Files**

The Qiime2 tutorials will have you download a file called emp-paired-end-sequences.zip. What we want is what is packaged in that .zip which I *think* we are getting from the sequencer:

1.) barcodes (barcodes.fastq.gz)

2.) forward reads (forward.fastq.gz)

3.) reverese reads. (reverse.fastq.gz)

You can see how to import these files here: https://docs.qiime2.org/2022.11/tutorials/importing/

> Files from Sequencers come in different formats. You can either ask your sequencer or go through Qiime2's importing tutorial to see which style aligns with your data.

In the tutorials, they have you mkdir emp-paired-end-sequences and then import the three files above into that directory. 

From here, you should be able to run the qiime2 protocols.

**Practicing the Qiime2 Pipeline with White Plague 2017 Microbiome Data.**


**Update your PATH before Demultiplexing, DADA2**

```
/yourpathtoimportedfiles/
```

## 1. DEMULTIPLEXING

**Tutorial Example**

  *example: demux*
  
  qiime demux emp-paired \
    --i-seqs sequences.qza \
    --m-barcodes-file sample-metadata.tsv \
    --m-barcodes-column barcode-sequence \
    --o-per-sample-sequences demux.qza \
    --o-error-correction-details demux-details.qza
    
> Additionally, you may need to add "--p-no-golay-error-correction \" to this. This is because the golay barcodes used by MRDNA were 8 nucleotides long, while Qiime2 is expecting them to be 12. This added line ignores the error.

**Customized for this project:**
```
mkdir DemultiplexedSeqs
qiime demux emp-paired \
  --i-seqs ./emp-paired-end-sequences/emp-paired-end-sequences.qza \
  --p-no-golay-error-correction \
  --m-barcodes-file 101018NM2-sample-metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences ./DemultiplexedSeqs/demux.qza \
  --o-error-correction-details ./DemultiplexedSeqs/demux-details.qza
```
> RUN TIME: 35 minutes (as emp-paired)


**Summary of the demultiplexing results**

qiime demux summarize \
  --i-data ./DemultiplexedSeqs/demux.qza \
  --o-visualization ./DemultiplexedSeqs/demux.qzv
> RUN TIME: (2 Min.)

qiime tools view ./DemultiplexedSeqs/demux.qzv


## 2. DADA2

**Tutorial Example:**

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 120 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza

**Customized for this project:**

```
mkdir DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./DemultiplexedSeqs/demux.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 220 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 226 \
  --o-representative-sequences ./DADA2/rep-seqs-dada2.qza \
  --o-table ./DADA2/table-dada2.qza \
  --o-denoising-stats ./DADA2/stats-dada2.qza
```
> RUN TIME: (25 minutes (315mb demux.qza)), (40+minutes (463mb demux.qza))

```
qiime metadata tabulate \
  --m-input-file ./DADA2/stats-dada2.qza \
  --o-visualization ./DADA2/stats-dada2.qzv
```
> RUN TIME: <1 min.

```
qiime feature-table summarize \
  --i-table ./DADA2/table-dada2.qza \
  --o-visualization ./DADA2/table-dada2.qzv \
  --m-sample-metadata-file 101018NM2-sample-metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data ./DADA2/rep-seqs-dada2.qza \
  --o-visualization ./DADA2/rep-seqs-dada2.qzv
 ```
> Run Time < 1 min.


**Repeat the above commands for each Sequencing Run. The below command will then merge the runs.**


## 3. Merging Sequencing Runs

**Tutorial Example**

qiime feature-table merge \
  --i-tables table-1.qza \
  --i-tables table-2.qza \
  --o-merged-table table.qza
qiime feature-table merge-seqs \
  --i-data rep-seqs-1.qza \
  --i-data rep-seqs-2.qza \
  --o-merged-data rep-seqs.qza
  
**Customized for this project:**

```
cd ../
qiime feature-table merge \
  --i-tables ./091018NMillcus515F-qiime2/table.qza \
  --i-tables ./101018NM1-qiime2/DADA2/table-dada2.qza \
  --i-tables ./101018NM2-qiime2/DADA2/table-dada2.qza \
  --o-merged-table merged-table.qza
qiime feature-table merge-seqs \
  --i-data ./091018NMillcus515F-qiime2/rep-seqs.qza \
  --i-data ./101018NM1-qiime2/DADA2/rep-seqs-dada2.qza \
  --i-data ./101018NM2-qiime2/DADA2/rep-seqs-dada2.qza \
  --o-merged-data rep-seqs.qza
```

> Run Time <1 min each.

Summary of merged artifact. Note how the sample-metadata has been concatenated from the other sequencing run metadata files.  

**Tutorial Example:**

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv
  
**Customized for this project:**

```
qiime feature-table summarize \
  --i-table merged-table.qza \
  --o-visualization merged-table.qzv \
  --m-sample-metadata-file ../../sample-metadata.tsv
```

> Run Time < 1min.

**View**

```
qiime tools view merged-table.qzv
```

## 4. GENERATE A TREE FOR PHYLOGENETIC DIVERSITY ANALYSES

**Tutorial Example**

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  
**Customized for this project:**

```
mkdir PhylogeneticTree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment ./PhylogeneticTree/aligned-rep-seqs.qza \
  --o-masked-alignment ./PhylogeneticTree/masked-aligned-rep-seqs.qza \
  --o-tree ./PhylogeneticTree/unrooted-tree.qza \
  --o-rooted-tree ./PhylogeneticTree/rooted-tree.qza
```
> RunTime: 15 min. (617 kb reps.seqs.qza)
  
> I skipped the Rarefying step in the Qiime2 Tutorial. My intention is to do Total Sum Scaling (TSM) normalization. 

## 5. TAXONOMIC CLASSIFICATION

**Classifiers**

> A "classifier" is basically a reference database catered to your dataset. Qiime2 provides some pretrained classifiers which you should use only if you know they are appropriate for your dataset. Fortunately, the V4 region (515F to 806R) region of the 16s rRNA gene is incredibly commonly studied so there are pretained datasets appropriate for our use. 
> For more details: https://docs.qiime2.org/2022.11/tutorials/feature-classifier/
> In this tutorial, I downloaded the most recent Greengenes 13_8 reference datasets "gg_13_8_otus.tar.gz" in: /Users/nickmacknight/Desktop/UMiami_NOAA/qiime2-tutorials/WhitePlague-tutorial/
> While my White plague example here used Greeengenes, consider Silva which is more recent ---> silva-138-99-515-806-nb-classifier.qza
> Opening this zipped file you will see the 99_otus.fasta and 99_otu_taxonomy.txt. I already have the reps_seq.qza from earlier in the analysis so we are set to make our own classifier.
> The example has 85 because the sequences need to have only 85% similarity for it to cluster, the sole purpose is that this is less computationally intensive so its better for a tutorial as it takes less time. But for REAL research questions we want 99 for 99% sequence similarity before clustering into a representative sequence.

**Tutorial Example:**

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 85_otus.fasta \
  --output-path 85_otus.qza
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 85_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza

**Customized for this project:**

```
mkdir ../Taxonomy
cd ../Taxonomy
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ./gg_13_8_otus/rep_set/99_otus.fasta \
  --output-path 99_otus.qza
```

>Run Time: 2 min.

```
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ./gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza
```

**Extract Reference Reads**

**Tutorial Example:**

qiime feature-classifier extract-reads \
  --i-sequences 99_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 120 \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs.qza
  
  
**Customized for this project:**

The F and R primers were the same as example. I adjusted the p-trunc-len to match what we chopped off earlier in the analysis.

```
qiime feature-classifier extract-reads \
  --i-sequences 99_otus.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-trunc-len 220 \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs.qza
```

> RUN TIME: 20 minutes.

**Train the Classifier**
**Tutorial Example:**

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza
> RUN TIME: 8 minutes

**Test the Classifier**
**Tutorial Example:**

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

**Customized for this project:**

```
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads ./../Seqs/rep-seqs.qza \
  --o-classification taxonomy.qza
```

> Run Time: 9 minutes

Filter those missing features out of your feature table. This helps with an error below during taxa bar plots.
https://forum.qiime2.org/t/valueerror-feature-ids-found-in-the-table-are-missing-from-the-taxonomy/5232/7

```
qiime feature-table filter-features \
  --i-table ./../Seqs/merged-table.qza \
  --m-metadata-file taxonomy.qza \
  --o-filtered-table id-filtered-table.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
  
qiime tools view taxonomy.qzv
```

**Tutorial Example**

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

**Customized for this project:**
```
qiime taxa barplot \
  --i-table ./id-filtered-table.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file ./../sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```
> Run Time: < 1 min.

```
qiime tools view taxa-bar-plots.qzv
```
