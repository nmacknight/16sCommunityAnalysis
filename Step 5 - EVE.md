# Expression Variance and Evolution (EVE) Model

Designed for quantitative data, the EVE model has been utilized in RNAseq analysis. This pipeline applies the EVE Model to microbial abundances rather than gene counts. When I first used the EVE model, it had just been published and the code wasnt even organized yet for public use. So I had some conversations with Lars Gronvold, the author of the github code for EVE. He managed to organize EVE for public use and through that conversation, discussed the use of the EVE model with microbial data. Because the model was designed for quantitative data, microbial abundance counts were considered appropriate. As many researchers who work with 16s or rnaseq know, some tools are shared between 16s and rnaseq.

By applying the EVE model, researchers will be now able to identify lineage-specific and phylosymbiotic bacteria in a formal, scalable, and replicable manner.

This analysis is a new approach to microbial ecology research that I hope to bring attention to. Personally, I see it as a new "tenet" in microbial ecology when compariing host species. Other tenets which are common in coral microbial ecology compare i.) differential abundance and ii.) diversity, few look at iii.) networks, and it requires a lot of repeat sampling/independent studies and time (about a decade of research) to anecdotally classify bacteria as "species-specific".

This model allows for formal and scalable identification of these species-specific bacteria and the implications of these co-evolved bacteria allow researchers to address emerging questions and historical questions we havent previously had the accurate tools to address including:

(I use the terms "species-specific" and "lineage-specific" synonymously. I also use the term "co-evolved" bacteria as a parental term that includes both lineage-specific and phylosymbiotic bacteria.)

- Identification of the Core Microbiome (lineage-specific bacteria are empirically supported bacteria of the core microbiome)
- Host-Microbe Co-evolution (Some coral lineages will notably possess more/less of these lineage-specific bacteria)
- Integral Members of Microbial Networks (Central/Hub bacteria with a lineage-speific signature)
- Preventative Mediators of Dysbiosis (Co-evolved bacteria when they are retained/enriched during disturbances)
- The Microbiome's Evolutionary Constraints (The Pathobiome - Lineage-specific bacteria associated with susceptibility)
- Evolved Holobiont Dependence Differences among Hosts (Is the microbiome equally relevant at mitiginating disturbances? Not always or equally across coral lineages and this helps us evaluate that holobiont dependence.)

With the EVE model, we can investigate these topics, whereas before the methods to do so have been missing or not satisfying the consitency we aim for.

Important! > To use the EVE Model, you need a dataste of quantitative counts (e.g. bacteria abundances or gene counts) and the phylogenetic distance between host species. I will explain how to obtain the phylogenetic distances in the pipeline.

> Here is the cool thing, you are not giving EVE your metadata, only which sample belongs to which host species (which is a bit of metadata really), the point is that EVE does not know what your treatment is or the response or any dependent variables besides the count matrix and phylogenetic tree.
> To better visualize lineage-specific bacteria, this will likely be a bacteria that within a species, there will be little variation in abundance, but the abundance between species is quite distinct or divergent. This is what makes the bacteria lineage-specific, its not responding to the treatment in its abundance and is divergent in its abundance relative to other species.
> By contrast, a highly variable bacteria will have a lot of variation in its abundance within a coral species, but between coral species, the abundance will have a similar average abundance. This means there isnt a lineage-specific relationship between highly variable bacteria and their coral host and perhaps this variation is the result of the treatment conditions.
> Finally, we are able to identify phylosymbiotic bacteria. This is a stricter classifcation of bacteria from the pool of lineage-specific bacteria. A phylosymbiotic bacteria will have an abundance pattern that recapitulates the host phylogenetic distance. I perform an additional Pegal's Lambda to determine signficance of this but the EVE output, specifically the low alpha value, points to candidates for phylosymbiosis.
> Perhaps you can now imagine the exciting questions that can be formally addressed now that those classifications are empirically supported and can be done at scale.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
#newest version of eve. 
devtools::install_gitlab("sandve-lab/evemodel")

library(evemodel)
```

### Input data

> File "Bacteria_EVE_input.txt" is a count matrix of bacteria abundances with rows as bacteria IDs and columns as sample IDs.

```{r}
# White Plague Microbiome. read tab separated table
setwd("~/Desktop/White Plague /Microbiome/Data/Analysis/Annotated")
bacteriaTbl <- read.delim("Bacteria_EVE_input.txt")
head(bacteriaTbl)
colnames(bacteriaTbl)
#logtransform bacterial data
dim(bacteriaTbl)
```

The table needs to be converted to a matrix:
```{r}
# The first column in the table is the Bacteria ID
# which will be the rownames of the matrix.

#White Plague 7 Species Microbiome:
bacteriaMat <- as.matrix(bacteriaTbl[,-1])
rownames(bacteriaMat) <- bacteriaTbl$species
summary(log10(bacteriaMat[,1:68]+1)) #seeing what kind of transformation is ideal. EVE recommends log transformed
bacteriaMat <- log10(bacteriaMat[,1:68]+1) #applying the log transform 
```

### Species phylogeny

The species phylogeny is in Newick format (from Orthofinder) which can be read using the `read.tree` function in the `ape` package:
> You will need transcriptomic or genomic data to infer a phylogenetic tree. If you do not have this specific to your project on hand, you can source host lineage transcriptomic or genomic references from NCBI and import that into Orthofinder or MEGA as another popular tool to create a phylogenetic tree. The Newick format is the numbers and parenthesis you see below, it is the visual phylogenetic tree written out.

```{r}
library(ape)

#White Plague 7 Species Bacteria Tree:
speciesTree <-read.tree(text="((Ppor:0.0952845,Past:0.101606):0.790892,((Cnat:0.122831,(Mcav:0.0936041,(Ofav:0.174805,Oann:0.174805):.3):0.307649):0.658835,Ssid:0.154932):0.790892);")

# plot the species tree
plot(speciesTree)
add.scale.bar(x=0,y=7,length = 0.1)
```

### Mapping species in tree to the columns in the expression matrix

The `evemodel` methods needs to know which species each column in the expression matrix belongs to. This is done by creating a vector of species names corresponding to the tip labels in the species tree.

```{r}
# the species names in the tree is given by the tip.labels
speciesTree$tip.label

#bacteria dataset
colnames(bacteriaMat)

# remove the trailing number so that we get a vector with the species for each column
colSpecies <- sub("_.*$","",colnames(bacteriaMat))

colSpecies
```


## Running the beta shared test

The beta shared test, (a.k.a. phylogenetic ANOVA), can detect bacteria with increased or decreased ratios of expression divergence to diversity (this ratio is the beta parameter). The model can be used for purposes such as identifying bacteria with high expression divergence between species as candidates for expression level adaptation (lineage-specific), and bacteria with high expression diversity within species as candidates for expression level conservation and or plasticity (highly variable).

This works by finding a shared beta that gives the maximum likelihood across all bacterias and comparing that to the model where the beta is fitted to each individual bacteria.

> Basically, the shared beta can be imagined as the average variation across your dataset, and then individual bacteria are examined to determine if their pattern of variation is signficantly (lineage-specific or highly variable).


```{r runTest, cache=TRUE}
#Bacteria results
res <- betaSharedTest(tree = speciesTree, bacteria.data = bacteriaMat, colSpecies = colSpecies)
res$sharedBeta
log(res$sharedBeta)
# 1.54493
```
### Results: LRT

The log likelihood ratio between the individual and shared beta fit indicates whether the individual beta was a better fit, i.e. the bacteria has an increased or decreased ratios of expression divergence to diversity. The log likelihood ratio test statistic is given `LRT` in the returned result and should follow a chi-squared distribution with one degree of freedom.

```{r}
# plot likelihood ratio test statistic histogram
hist(res$LRT,freq = F)

# Plot the chi-squared distribution with one degree of freedom
x = seq(0.5,10,length.out = 100)
y = dchisq(x,df = 1)
lines(x,y,col="red")
```

P-value can then be calculated using:

```{r}
pval = pchisq(res$LRT,df = 1,lower.tail = F)
```

### Result: fitted parameters

The shared beta:

```{r}
res$sharedBeta
```

The parameters for each gene given shared beta:

```{r}
head(res$sharedBetaRes$par)
```

The parameters for each gene given individual beta:

```{r}
head(res$indivBetaRes$par)
```

Note that, with the exception of theta, the fitted parameters are rather unstable

Combine LRT, beta, theta, sigma2, alpha:
```{r}
head(cbind(res$LRT,res$indivBetaRes$par,pval))
a <- cbind(res$LRT,res$indivBetaRes$par,pval,bacteriaMat[,1:68])
#genes
rownames(a) <- exprTbl$Entry
#bacteria
rownames(a) <- rownames(bacteriaMat)
```

Give Column Name to LRT column

```{r}
colnames(a)
colnames(a)[colnames(a)==""] <- "LRT"
head(a)
```

Export results

```{r}
write.csv(a, file="R_EVE_results_WP_bacteria_7species.csv")
```

Visualize LRT v Beta by volcano plot for gene data
```{r}
Bac_EVE <- read.csv("R_EVE_results_WP_bacteria_7species.csv")
colnames(Bac_EVE)
colnames(Shared_EVE)[colnames(Shared_EVE)=="X"] <- "Entry"
head(Bac_EVE)

library(tidyverse)
Bac_EVE_sig <- Bac_EVE %>% filter(pval <= 0.1)
write.csv(Bac_EVE_sig, file="Signficant EVE Bacteria.csv")

# Make a basic volcano plot
with(Bac_EVE, plot(log(beta),LRT, pch=20, main="EVE", xlim=c(-5,5),ylim=c(-1,110)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(Bac_EVE, LRT>3.4), points(log(beta), LRT, pch=20, col="mediumpurple1"))
with(subset(Bac_EVE, LRT<3.4 & abs(log(beta))>1.5), points(log(beta), LRT, pch=20, col="pink1"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(Bac_EVE_sig, LRT>20 & abs(log(beta))>0), textxy(log(beta), LRT, labs=Bacteria, cex=0.7, offset=0.6))
with(subset(Bac_EVE_sig, LRT<5 & abs(log(beta))>2), textxy(log(beta), LRT, labs=Bacteria, cex=0.7, offset=0.6))
```

### Isolating Lineage Specific "LS" and Highly Variable "HV" Bacteria using the sharedbeta metric.
```{r,results='hide',tidy=TRUE}
# P-value can then be calculated using:
pval = pchisq(res$LRT,df = 1,lower.tail = F)

# The shared beta:
res$sharedBeta
4.687645

log(4.687645)

#Combine LRT, beta, theta, sigma2, alpha:
head(cbind(res$LRT,res$indivBetaEVE$par,pval))
colnames(bacteriaMat)
dim(bacteriaMat)
# [1]  371   76
a <- cbind(res$LRT,res$indivBetaRes$par,pval,bacteriaMat)
#rownames(a) <- NodesignTbl.t$Orthogroup

# Give Column Name to LRT column
colnames(a)
colnames(a)[colnames(a)==""] <- "LRT"
head(a)
EVE_results <- as.data.frame(a)
EVE_results <- tibble::rownames_to_column(EVE_results, "Orthogroup")

EVE_results$type <- ifelse(EVE_results$beta<log(res$sharedBeta),"Lineage Specific","Highly Variable")
EVE_results$significant <- ifelse(EVE_results$pval<=0.1,"Significant","Not Significant")
EVE_results$category <- ifelse(EVE_results$significant == "Significant",EVE_results$type,"NS")
colnames(EVE_results)
EVE_results <- EVE_results %>% relocate(type, .after = pval)
EVE_results <- EVE_results %>% relocate(significant, .after = type)
EVE_results <- EVE_results %>% relocate(category, .after = significant)

Shared_EVE <- EVE_results
Shared_EVE_sig <- Shared_EVE %>% filter(pval <= 0.1)
Shared_EVE_LS <- Shared_EVE_sig %>% filter(log(beta)< (log(res$sharedBeta)))
dim(Shared_EVE_LS)  #108
Shared_EVE_HV <- Shared_EVE_sig %>% filter(log(beta)> (log(res$sharedBeta)))
dim(Shared_EVE_HV)  #160

```

# Figure Creation

### Volcano Plot - Figure 1A, Supp Data 1
```{r EVE_volcano,fig.height=3,fig.width=3,dpi=300}
p <- ggplot(data = Shared_EVE,
       aes(x=log(beta),y=LRT,color=category))+
  geom_point()+
  geom_vline(xintercept = log(res$sharedBeta), col="black", linetype="dashed")+
  annotate("text", x = 3, y = -4, label = "Shared beta", size = 3)+
  geom_hline(yintercept = -log(0.1), col="black",linetype = "dashed")+
  annotate("text", x = -2, y = 4.5, label = "p = 0.1", size = 3)+
  scale_color_manual(values=c("pink1","mediumpurple1","black"))+
  #xlim(-5,5)+
  theme_light()
p + labs(color = "EVE Gene Category",title = "EVE Genes Volcano Plot",subtitle = "108 Lineage-Specific and 160 Highly Variable Genes") +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9))




p <- ggplot() +
  # First, plot black dots (NS) in the background
  geom_point(data = Shared_EVE %>% filter(category == "NS"), 
             aes(x = log(beta), y = LRT), 
             color = "black", alpha = 0.6, size = 1.5) +
  
  # Then, plot purple dots (Lineage Specific) in the middle layer
  geom_point(data = Shared_EVE %>% filter(category == "Lineage Specific"), 
             aes(x = log(beta), y = LRT), 
             color = "mediumpurple1", alpha = 0.8, size = 1.8) +
  
  # Finally, plot pink dots (Highly Variable) on top so they are most visible
  geom_point(data = Shared_EVE %>% filter(category == "Highly Variable"), 
             aes(x = log(beta), y = LRT), 
             color = "pink1", size = 2) +

  # Add dashed lines for significance and shared beta
  geom_vline(xintercept = log(res$sharedBeta), col="black", linetype="dashed") +
  annotate("text", x = 3, y = -4, label = "Shared beta", size = 3) +
  geom_hline(yintercept = -log(0.1), col="black", linetype="dashed") +
  annotate("text", x = -2, y = 4.5, label = "p = 0.1", size = 3) +

  # Customize color legend manually
  scale_color_manual(values=c("pink1", "mediumpurple1", "black")) +

  # Apply theme and labels
  theme_light() +
  labs(color = "EVE Gene Category", 
       title = "EVE Genes Volcano Plot",
       subtitle = "108 Lineage-Specific and 160 Highly Variable Genes") +
  theme(plot.title = element_text(face = "bold"), 
        plot.subtitle = element_text(size = 9))

p  # Display the plot

```

### PCA of LS and HV Bacteria
```{r}

# LS
sig_LS_Bac <- Shared_EVE_LS

sig_LS_Bac <- t((sig_LS_Bac[,c(1,12:length(sig_LS_Bac))]))

colnames(sig_LS_Bac) <- sig_LS_Bac[1,]
sig_LS_Bac. <- sig_LS_Bac[-1,]
sig_LS_Bac. <- cbind(SampleID = rownames(sig_LS_Bac.), sig_LS_Bac.)

sig_LS_Bac_meta <- merge(metadata,sig_LS_Bac., by="SampleID")

write.csv(sig_LS_Bac_meta, file="path/sig_LS_Bac_meta.csv") # converted to numeric. 
sig_LS_Bac_meta <- read.csv("path/sig_LS_Bac_meta.csv")

# HV
sig_HEV_Bac <- Shared_EVE_HV

sig_HEV_Bac <- t((sig_HEV_Bac[,c(1,12:length(sig_HEV_Bac))]))

colnames(sig_HEV_Bac) <- sig_HEV_Bac[1,]
sig_HEV_Bac. <- sig_HEV_Bac[-1,]
sig_HEV_Bac. <- cbind(SampleID = rownames(sig_HEV_Bac.), sig_HEV_Bac.)

sig_HEV_Bac_meta <- merge(metadata,sig_HEV_Bac., by="SampleID")

write.csv(sig_HEV_Bac_meta, file="path/sig_HEV_Bac_meta.csv") # converted to numeric. 
sig_HEV_Bac_meta <- read.csv("path/sig_HEV_Bac_meta.csv")


#Lineage-Specific
LS_PCA_input <- sig_LS_Bac_meta[,11:36]
pca_res_LS <- prcomp(LS_PCA_input, scale. = TRUE)
LSPCA <- autoplot(pca_res_LS, data = sig_LS_Bac_meta, colour = 'Coral.Species',shape = 'Outcome', frame = TRUE, frame.type = 'norm')
LSPCA <- LSPCA + theme(legend.position="bottom")+
          scale_fill_manual(breaks = c("Ofav","Mcav","Acer","Past"), 
                            values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
          scale_color_manual(breaks = c("Ofav","Mcav","Acer","Past"), 
                            values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
LSPCA

#Phylosymbiotic
Phylosymbiotic_Bacteria_Transposed <- read.csv("Phylosymbiotic Bacteria_Transposed.csv")
PS_PCA_input <- Phylosymbiotic_Bacteria_Transposed[,11:43]
pca_res_PS <- prcomp(PS_PCA_input, scale. = TRUE)
PSPCA <- autoplot(pca_res_PS, data = Phylosymbiotic_Bacteria_Transposed, colour = 'Species',shape = 'Status', frame = TRUE, frame.type = 'norm')
PSPCA <- PSPCA + theme(legend.position="bottom")+
          scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                            values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
          scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                            values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
PSPCA

#Highly Variable
HEV_PCA_input <- sig_HEV_Bac_meta[,12:666]
pca_res_HEV <- prcomp(HEV_PCA_input, scale. = TRUE)
HEVPCA <- autoplot(pca_res_HEV, data = sig_HEV_Bac_meta, colour = 'Coral.Species',shape = 'Outcome',frame = TRUE, frame.type = 'norm')
HEVPCA <- HEVPCA + theme(legend.position="bottom")+
          scale_fill_manual(breaks = c("Ofav","Mcav","Acer","Past"), 
                            values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
          scale_color_manual(breaks = c("Ofav","Mcav","Acer","Past"), 
                            values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
HEVPCA 

library(patchwork)
HEVPCA
HEVPCA | LSPCA | PSPCA 

violin+(LRTPCA/betaPCA)
LRTPCA|violin|betaPCA

Results: LSPCA separates by species, organized by phylogeny as expected. Porites and Siderastraea overlap a lot. The orbicella are divergent. Montastraea is more similar to O.annularis. Cnat seems highly variable per usual. PC1 loading is quite high (46.73%)
Discussion: The takeaway from HEVPCA is that the species cluster on top of each other as expected showing low host influence.
```


### Histogram - Supp Fig 1
```{r}

Bac_EVE_sig$Frequency <- rowSums(Bac_EVE_sig[, 12:113] > 0, na.rm = TRUE)

#Histogram
#HEV_Bac <- read.csv("Putative HEV EVE Bacteria.csv")
#Bacteria_EVE_input <- read.csv("Bacteria_EVE_input.csv") #all bacteria
HEV_Bac_Freq30 <- Bac_EVE_sig %>% filter(Frequency > 30) #keep only bacteria with a frequency greater than 10. 
HEV_Bac_Freq20 <- Bac_EVE_sig %>% filter(Frequency > 20) #keep only bacteria with a frequency greater than 10. 
HEV_Bac_Freq10 <- Bac_EVE_sig %>% filter(Frequency > 10) #keep only bacteria with a frequency greater than 10. 
HEV_Bac_Freq5 <- Bac_EVE_sig %>% filter(Frequency > 5) #keep only bacteria with a frequency greater than 5.
HEV_Bac_Freq3 <- Bac_EVE_sig %>% filter(Frequency > 3) #keep only bacteria with a frequency greater than 3. 

hist(Bac_EVE_sig$Frequency,freq = F)

# Plot the chi-squared distribution with one degree of freedom
x = seq(0.2,10,length.out = 100)
y = dchisq(x,df = 1)
lines(x,y,col="red")

# I was curious how frequent highly variable bacteria were, some may be only present in a few samples, these may not be as biologically relevant to you, or they may be the focus of your investigation depending on your questions. For me, I was most curious about the lineage-specific bacteria and phylosymbiotic bacteria.
```


### Box and Whisker Plots - Individual. 
```{r}
library(ggplot2)
library(tidyverse)
library(ggthemes)  # for a mapping theme
library(ggalt)  # for custom map projections
library(ggrepel)  # for annotations
library(viridis)
bacteria_results <- read.csv("EVE_Results.csv")


#Visualize
p.piscicida_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (pseudoalteromonas.piscicida)), y = pseudoalteromonas.piscicida, fill = Species)) +
    geom_point(aes(y = pseudoalteromonas.piscicida, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
  #geom_text(aes(x=Species, y = pseudoalteromonas.piscicida, group=Status, label = paste(Genotype)),parse = TRUE ,
             # position =position_dodge(width = 1),size = 3, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Pseudoalteromonas piscicida", subtitle = "LRT = 102.9 | Beta = 0.055 | Alpha = 11.9 | pval = 3.53e-24") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))


endozoicomonas.spp_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (endozoicomonas.spp)), y = endozoicomonas.spp, fill = Species)) +
    geom_point(aes(y = endozoicomonas.spp, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Endozoicomonas spp.", subtitle = "LRT = 22.0 | Beta = 0.42 | Alpha = 7.05 | pval = 2.68e-6") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))
endozoicomonas.spp_species 


t.cechii_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (tetracoccus.cechii)), y = tetracoccus.cechii, fill = Species)) +
    geom_point(aes(y = tetracoccus.cechii, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Tetracoccus cechii", subtitle = "LRT = 4.92 | Beta = 0.86 | Alpha = 0.3 | pval = 0.026") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))
t.cechii_species 

roseicyclus.spp_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (roseicyclus.spp)), y = roseicyclus.spp, fill = Species)) +
    geom_point(aes(y = roseicyclus.spp, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Roseicyclus spp.", subtitle = "LRT = 9.02 | Beta = 0.59 | Alpha = 0.32 | pval = 0.0026") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))
roseicyclus.spp_species 

hyphomonas.spp_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (hyphomonas.spp)), y = hyphomonas.spp, fill = Species)) +
    geom_point(aes(y = hyphomonas.spp, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Hyphomonas spp.", subtitle = "LRT = 4.9 | Beta = 0.91 | Alpha = 0.36 | pval = 0.027") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))
hyphomonas.spp_species 

oleiphilus.spp_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (oleiphilus.spp)), y = oleiphilus.spp, fill = Species)) +
    geom_point(aes(y = oleiphilus.spp, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Oleiphilus spp.", subtitle = "LRT = 8.63 | Beta = 0.64 | Alpha = 0.39 | pval = 0.003") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))
oleiphilus.spp_species


p.veronii_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (pseudomonas.veronii)), y = pseudomonas.veronii, fill = Species)) +
    geom_point(aes(y = pseudomonas.veronii, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Pseudomonas veronii", subtitle = "LRT = 14.05 | Beta = 0.59 | Alpha = 3.74 | pval = 0.0001") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))
p.veronii_species

ruegeria.silicibacter.sp_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (ruegeria.silicibacter.sp)), y = ruegeria.silicibacter.sp, fill = Species)) +
    geom_point(aes(y = ruegeria.silicibacter.sp, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Ruegeria silicibacter sp.", subtitle = "LRT = 3.67 | Beta = 1.20 | Alpha = 0.67 | pval = 0.055") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))
ruegeria.silicibacter.sp_species

bellilinea.spp_species <- ggplot(data = bacteria_results, 
           aes(x = reorder(Species, (bellilinea.spp)), y = bellilinea.spp, fill = Species)) +
    geom_point(aes(y = bellilinea.spp, color = Species, shape=Status), 
               position =position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "bellilinea.spp", subtitle = "EVE: LRT = 44.96 | Beta = 0.199 | Alpha = 4.25 | pval = 2.01e-11 | Pearson: corr = 0.87 | pval = 0.009") +
scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                       values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
    #scale_color_brewer()+
    #scale_fill_brewer()+
    #theme_classic()+ 
  scale_y_continuous(breaks=c(-1,0,1,2,3,4,5))+
    theme(legend.position= "Null", legend.text=element_text(size=20), 
          legend.title = element_text(size=20), 
          axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
          axis.text.y = element_text(colour = "black", size=10), 
          axis.title.x = element_text(size=10), 
          axis.title.y = element_text(size=10), 
          strip.text.x = element_text(colour = "black", size = 20))
bellilinea.spp_species

```

### signifcantly correlated to Relative Risk
```{r}
#signifcantly correlated to Relative Risk
#positively
bellilinea.spp_RR_Corr #t = 4.0557, df = 5, p-value = 0.00977, cor = 0.876. More bellilinea is associated with susceptibility
pseudoalteromonas.piscicida_RR_Corr

```

### Correlation Test - Bacteria to Relative Risk
```{r}
library(car)
library(dplyr)
library(tidyverse)

bacteria_results <- read.csv("EVE_Results.csv")
Anova(model <- glm(nautella.italica~Status, data = bacteria_results))
TukeyHSD(aov(model))


avg.EVE_RR <- bacteria_results[,c(5:273)] %>% group_by(RelativeRisk) %>% summarise_all(list(mean))
#or
#avg.EV_RR <- read.csv("avg.EV_RR.csv")

cor.test(avg.EV_RR$endozoicomonas.spp,avg.EV_RR$RelativeRisk, method=c("pearson"))

cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson"))

cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$nautella.italica, method=c("pearson")) #cor = 0.8108507 t = 3.098, df = 5, p-value = 0.02691
cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$pseudoalteromonas.spp, method=c("pearson")) #cor = 0.9376406 t = 6.0316, df = 5, p-value = 0.001803
cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$thalassobius.spp, method=c("pearson")) #cor = 0.933 t = 5.8344, df = 5, p-value = 0.002092
cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$thalassobius.sp, method=c("pearson")) #cor = 0.9415741 t = 6.2511, df = 5, p-value = 0.001535
cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$bellilinea.spp, method=c("pearson")) #cor = 0.9415741 t = 6.2511, df = 5, p-value = 0.001535

cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$endozoicomonas.spp, method=c("pearson"))#ns

#summarize statuses and correlate Ruegeria to putative pathogens. 
bacteria_results <- read.csv("EVE_Results.csv")
avg.EVE_Status <- bacteria_results[,c(4,6:273)] %>% group_by(Status)  %>% mutate(across(4,6:273,list(mean))) 
avg.EVE_Status <- bacteria_results[,c(5:273)] %>% group_by(Status) %>%summarise(across(5:268, list(mean=mean))) 


avg.EVE_Status <- bacteria_results[,c(5:273)] %>% group_by(Status) %>% summarise_all(list(mean))
write.csv(avg.EVE_Status, file= "avg.EVE_Status_Bacteria.csv")
```

```{r}
bellilinea.spp_RR_Corr <- ggplot(avg.EV_RR, aes(x=RelativeRisk, y=bellilinea.spp)) +
  geom_point() +
  geom_text(aes(label = paste0("(", paste(round(bellilinea.spp, digits=2)), ")")), nudge_y = 0.25) +
  #label=sprintf("%0.2f", round(bellilinea.spp, digits = 2))+
  
  geom_smooth(method=glm) + geom_text(label=avg.EV_RR$Species, nudge_y = 0.5)+ 
  labs(y="Mean log10(abundance+1)", x = "Relative Risk of Disease Incidence", title = "Bellilinea spp.", subtitle = "EVE: LRT = 44.96 | Beta = 0.199 | Alpha = 4.25 | pval = 2.01e-11 | Pearson: corr = 0.87 | pval = 0.009") + 
  scale_y_continuous(breaks=c(-1,0,1,2,3,4,5)) #+  
  #annotate("text", label = "R= 0.87", x = 4.2, y = 4.8, size = 3, colour = "black") + annotate("text", label = "p= 0.009", x = 4.2, y = 4.3, size = 3, colour = "black") #+ scale_x_reverse()
bellilinea.spp_RR_Corr

pseudoalteromonas.piscicida_RR_Corr <- ggplot(avg.EV_RR, aes(x=RelativeRisk, y=pseudoalteromonas.piscicida)) +
  geom_point() +
  geom_text(aes(label = paste0("(", paste(round(pseudoalteromonas.piscicida, digits=2)), ")")), nudge_y = 0.25) +
  geom_smooth(method=glm) + geom_text(label=avg.EV_RR$Species, nudge_y = 0.5)+ 
  labs(y="Mean log10(abundance+1)", x = "Relative Risk of Disease Incidence", title = "Pseudoalteromonas piscicida", subtitle = "EVE: LRT = 102.9 | Beta =0.055 | Alpha = 11.9 | pval = 3.53e-24 | Pearson: corr = 0.67 | pval = 0.09") + 
  scale_y_continuous(breaks=c(0,1,2,3,4,5))

pseudoalteromonas.piscicida_RR_Corr

endozoicomonas.spp_RR_Corr <- ggplot(avg.EV_RR, aes(x=RelativeRisk, y=endozoicomonas.spp)) +
  geom_point() +
  geom_text(aes(label = paste0("(", paste(round(endozoicomonas.spp, digits=2)), ")")), nudge_y = 0.25) +
  geom_smooth(method=glm) + geom_text(label=avg.EV_RR$Species, nudge_y = 0.5)+ 
  labs(y="Mean log10(abundance+1)", x = "Relative Risk of Disease Incidence", title = "Endozoicomonas spp.", subtitle = "LRT = 22.0 | Beta = 0.42 | Alpha = 7.05 | pval = 2.68e-6 | Pearson: corr = -0.72 | pval = 0.066") + 
  scale_y_continuous(breaks=c(0,1,2,3,4,5))

endozoicomonas.spp_RR_Corr

ruegeria.silicibacter.sp.avg_RR_Corr <- ggplot(avg.EV_RR, aes(x=DysbiosisC2I, y=ruegeria.silicibacter.sp)) +
  geom_point() +
  geom_text(aes(label = paste0("(", paste(round(ruegeria.silicibacter.sp, digits=2)), ")")), nudge_y = 0.25) +
  geom_smooth(method=glm) + geom_text(label=avg.EV_RR$Species, nudge_y = 0.5)+ 
  labs(y="Mean log10(abundance+1)", x = "Microbial Dissimilarity Control to Disease-Infected", title = "Ruegeria silicibacter sp.", subtitle = "Pearson: corr = 0.89 | pval = 0.00654") + 
  scale_y_continuous(breaks=c(0,1,2,3,4,5))
ruegeria.silicibacter.sp.avg_RR_Corr


ruegeria.silicibacter.sp._pseudoalteromonas.spp_avg_RR_Corr <- ggplot(avg.EV_RR, aes(x=pseudoalteromonas.spp, y=ruegeria.silicibacter.sp)) +
  geom_point() +
  geom_text(aes(label = paste0("(", paste(round(ruegeria.silicibacter.sp, digits=2)), ")")), nudge_y = 0.25) +
  geom_smooth(method=glm) + geom_text(label=avg.EV_RR$Species, nudge_y = 0.5)+ 
  labs(y="ruegeria.silicibacter.sp abundance", x = "pseudoalteromonas.spp abundance", title = "Ruegeria silicibacter sp.", subtitle = "Pearson: corr = 0.94 | pval =  0.001803") + 
  scale_y_continuous(breaks=c(0,1,2,3,4,5))
ruegeria.silicibacter.sp._pseudoalteromonas.spp_avg_RR_Corr
```
```{r}
ggsave("ruegeria.silicibacter.sp.avg_RR_Corr.pdf")

```
```{r}
library(patchwork)

bellilinea.spp_species+bellilinea.spp_RR_Corr
```

### Basic line plot with points
```{r}


# Basic line plot with points
ruegeria.silicibacter.sp_statusLine <- ggplot(data=avg.EVE_Status, aes(x=Status, y=ruegeria.silicibacter.sp, group=1)) +
  geom_line(color="cornflowerblue")+
  geom_point(color="cornflowerblue")+
  ylim(0,3)

endozoicomonas.spp_statusLine <- ggplot(data=avg.EVE_Status, aes(x=Status, y=endozoicomonas.spp, group=1)) +
  geom_line(color="cornflowerblue")+
  geom_point(color="cornflowerblue")+
  ylim(0,3)

nautella.italica_statusLine <- ggplot(data=avg.EVE_Status, aes(x=Status, y=nautella.italica, group=1)) +
  geom_line(color="salmon")+
  geom_point(color="salmon")+
  ylim(0,3)
thalassobius.spp_statusLine <- ggplot(data=avg.EVE_Status, aes(x=Status, y=thalassobius.spp, group=1)) +
  geom_line(color="salmon")+
  geom_point(color="salmon")+
  ylim(0,3)
pseudoalteromonas.spp_statusLine <- ggplot(data=avg.EVE_Status, aes(x=Status, y=pseudoalteromonas.spp, group=1)) +
  geom_line(color="salmon")+
  geom_point(color="salmon")+
  ylim(0,3)

library(patchwork)
ruegeria.silicibacter.sp_statusLine+endozoicomonas.spp_statusLine
nautella.italica_statusLine+thalassobius.spp_statusLine+pseudoalteromonas.spp_statusLine
```
```{r}
ggsave("pseudoalteromonas.spp_statusLine.pdf")

```



### Box and Whisker Construction - For Loop
```{r}
# For Loop to create multiple bar graphs (of the same format) from a single csv file

#### Background info ####

# Desired boxplot figure output:
#   x-axis: coral species
#   y-axis: bacterium abundance
#   Colors/groups: based on species susceptibility
#   Individual figures for each: Bacteria species

# Load necessary libraries
library(ggplot2) 

# Read in data
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')

# Example figure
endozoicomonas.spp_species <- ggplot(data = bacteria_results, 
                                     aes(x = reorder(Species, (endozoicomonas.spp)), y = endozoicomonas.spp, fill = Species)) +
  geom_point(aes(y = endozoicomonas.spp, color = Species, shape=Status), 
             position =position_dodge(width = 0.47),size = 2, alpha = 2) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
  labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Endozoicomonas spp.", subtitle = "LRT = 22.0 | Beta = 0.42 | Alpha = 7.05 | pval = 2.68e-6") +
  scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                    values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                     values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))+
  theme(legend.position= "Null", legend.text=element_text(size=20), 
        legend.title = element_text(size=20), 
        axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
        axis.text.y = element_text(colour = "black", size=10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        strip.text.x = element_text(colour = "black", size = 20))
endozoicomonas.spp_species 

#### Original code  ####

# For loop prerequites (for this loop, not all loops):
#   1. All figures created should have the same format OR 
#      should have distinguishable features to be used in if/else statements 
#      to direct the loop through alternate chucks of code
#   2. All figures should come from the same csv file and columns should 
#      have a uniform naming convention over which the loop interates

# STEP 1: Prior to running any for loop create a folder on your machine into 
#         which you will put all of your completed figures
#         Ex) "C:/Users/Amy/Documents/Personal/For Nick/Bacterial_abundance_boxplots"

# STEP 2: Load any needed libraries to run code (see above)
#         (require() = If library is not loaded, load it. If it is loaded, do nothing.)
require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

# STEP 3: Set working directory to location of needed files and read in data files
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
sig_bacteria <- read.csv("Signficant EVE Bacteria.csv", header = T, fileEncoding = 'UTF-8-BOM')

library(dplyr)

# Bac_EVE_sig <- Bac_EVE_sig %>%
#   mutate(Pattern = ifelse(LRT > 5, "LS", "HV"))
# write.csv(Bac_EVE_sig, "/Users/nicholas.macknight/Desktop/Autumn16s/Autumn16s/Bac_EVE_sig.csv") # Moved Frequency and Pattern column closer to beginning of data for quick viewing within R. 

Bac_EVE_sig <- read.csv("/Users/nicholas.macknight/Desktop/Autumn16s/Autumn16s/Bac_EVE_sig.csv")

# STEP 4: Do some initial data formatting and object creation to make for loop smoother
# Create object for constant plot theme info 
# This chuck of code does not change regardless of the plot
my_fill <- scale_fill_manual(breaks = c("Ofav", "Past", "Mcav","Acer"), 
                                values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
my_col <- scale_color_manual(breaks = c("Ofav", "Past", "Mcav","Acer"), 
                     values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
my_theme <- theme(legend.position= "Null", legend.text=element_text(size=20), 
        legend.title = element_text(size=20), 
        axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
        axis.text.y = element_text(colour = "black", size=10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        strip.text.x = element_text(colour = "black", size = 20))

# Pull out all metadata columns from dataframe
metadata <- colnames(bacteria_results[,colnames(bacteria_results) %in% c("Bacteria","Species","Treatment","RelativeRisk","Status","Genotype")])
# Pull out all bacterial abundance columns from dataframe
bacteria <- colnames(bacteria_results[,!colnames(bacteria_results) %in% metadata])
# Pull out the relevant columns from the significant eve bacteria results
eve <- sig_bacteria[,colnames(sig_bacteria) %in% c("Bacteria","LRT_LS","alpha_phylosymbiosis","beta_HEV","pval","Expression")]
# rename columns for ease
colnames(eve)[2:5] <- c("LRT","Alpha","Beta","Pattern")

# Some of the Eve bacteria names end in a space
# Remove the trailing space and now the names should match up
eve$Bacteria <-  trimws(eve$Bacteria, which = "right") #cut out space after species name.

# Make sure we can match the bacteria names in the eve and bacteria vectors
bacteria_names1 <- eve$Bacteria # create vector of bacteria names
bacteria_names2 <- str_replace_all(bacteria, "[.]", " ") # create vector of names but replace '.' with spaces

bacteria_names1[!bacteria_names1 %in% bacteria_names2] # find mismatching cases. Ideally there would be 0. 
# THERE ARE 0

# ****************************************************************************************************
# For loop testing: Un-comment (remove hashtag for line immediately below) i = 1, this will allow you to test the code only on the first loop

#i = 1

# IMPORTANT NOTE: While testing loops skip any lines that start with "for (....){"
# Only run lines that are INSIDE the loop brackets "{}"
# This will allow you to test and fix code running only the first loop (i.e., creating 1 figure) 
# and will prevent R from crashing if the code is buggy.
# the reason you test a for loop is to make sure it works before doing many iterations. This is why the test loop is i=1, i.e do one loop.
#****************************************************************************************************

# STEP 7: The For Loop
# The "I loop"

for (i in 1:length(bacteria)){ #iterates loop over each column of bacteria
  trip <- bacteria[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  results <- eve[eve$Bacteria %in% trip_name,] # pulls out the relevant EVE result line
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(bacteria), " bacteria"))
  
  y_val <- paste0(trip) # Create variable for current bacteria name
  data <- bacteria_results[,names(bacteria_results) %in% c(y_val, metadata)] # create small df for current iteration
  colnames(data)[7] <- "y_val" #rename bacteria abundance column
  
  # PLOT CODE
  plot <- ggplot(data = data, aes(x = reorder(Species, (y_val)), y = y_val, fill = Species)) +
    geom_point(aes(y = y_val, color = Species, shape = Status), 
               position = position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", 
         title = paste0(str_to_sentence(trip_name)), 
         subtitle = paste0("LRT = ", round(results$LRT, digits = 1),
                           " | Beta = ", round(results$Beta, digits = 2),
                           " | Alpha = ", round(results$Alpha, digits = 2),
                           " | pval = ", formatC(results$pval, format = 'e', digits = 2))) +
        my_fill + my_col + my_theme 
  
  # setwd to folder for all the plots to be placed into
  if (results$Pattern == "HV"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/HV EVE Bacteria Box Plots")
    pdfname <- paste0(str_to_sentence(trip_name), "_HV_boxplot.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$Pattern == "LS") {
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Box Plots")
    pdfname <- paste0(str_to_sentence(trip_name), "_LS_boxplot.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  # DUNZO
}

# STEP 8: Check that pdfs are written to correct folder
# STEP 9: (Might required paid version of Adobe Acrobat? - If so, I can do this for you.)
# To combine into single pdf for easy scrolling:
#         In file explorer highlight all plots > right click > 
#         (IF available) "Combine files in Acrobat" > Combine > 
#         Rename "binder" > View > Page display > Two page scrolling
```

### Phylobacteriatic Correlation - Pegal's Lambda - For Loop - Figure 2
```{r}
# Example Phylobacteriatic correlation

trait<-sig_LS_Bac_avg$ruegeria.silicibacter.sp
names(trait)<-(sig_LS_Bac_avg$Species)
speciesTree <-read.tree(text="((P.porites:0.0952845,P.astreoides:0.101606):0.790892,((C.natans:0.122831,(M.cavernosa:0.0936041,(O.faveolata:0.174805,O.annularis:0.174805):.3):0.307649):0.658835,S.siderea:0.154932):0.790892);")
vcv.phylo(speciesTree)
phylosig(speciesTree, trait, method="lambda", test=TRUE, nsim=999)

#Prepare Loop

#Identify the phylobacteriatic tree
speciesTree <-read.tree(text="((P.porites:0.0952845,P.astreoides:0.101606):0.790892,((C.natans:0.122831,(M.cavernosa:0.0936041,(O.faveolata:0.174805,O.annularis:0.174805):.3):0.307649):0.658835,S.siderea:0.154932):0.790892);")
vcv.phylo(speciesTree)

sig_LS_Bac_avg <- sig_LS_Bac %>% group_by(Species) %>% summarize_if(is.numeric, mean, na.rm = T)
sig_LS_Bac_avg <- sig_LS_Bac_avg[,!colnames(sig_LS_Bac_avg) %in% c("RelativeRisk", "Genotype")]

corr_bacteria <- names(sig_LS_Bac_avg[,-1])
Species <- colnames(sig_LS_Bac_avg[,1])

all_outputs <- data.frame() # Create empty dataframe, outputs will be added to this with each loop iteration

#The For Loop
#i = 1
for (i in 1:length(corr_bacteria)){ #iterates loop over each column of bacteria
  trip <- corr_bacteria[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip)
  data <- sig_LS_Bac_avg[,colnames(sig_LS_Bac_avg) %in% c(Species, bacteria)] # create small df for current iteration
  colnames(data)[2] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Species # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(corr_bacteria), " bacteria"))
  
  # CORRELATION CODE
  # We chose Lineage Specific Bacterial Abundance as the trait we are testing for phylobacteriatic signal. Then, we do the test with 999 randomizations:
  output <- phylosig(speciesTree, trait, method="lambda", test=TRUE, nsim=999)
  
  #Put the output into a dataframe and bind it with the previous outputs
  this_output <- data.frame(bacteria = bacteria, lambda = output$lambda, logL = output$logL, logL0 = output$logL0, P = output$P)
  all_outputs <- rbind(all_outputs, this_output)

  # DUNZO
}

View(all_outputs)

sig_outputs <- all_outputs[all_outputs$P <= 0.05,] #these are phylosymbiotic bacteria!
length(sig_outputs$bacteria)
write.csv(sig_outputs, file="Phylosymbiotic Bacteria.csv")#replaced "." with " " in bacterial species names and capitalized Bacteria column name
sig_bacteria$Bacteria <-  trimws(sig_bacteria$Bacteria, which = "right") #cut out space after species name.
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
Phylosymbiotic_Bacteria <- read.csv("Phylosymbiotic Bacteria.csv")
Phylosymbiotic_Bacteria <- merge(Phylosymbiotic_Bacteria, sig_bacteria, by="Bacteria")
write.csv(Phylosymbiotic_Bacteria, file="Phylosymbiotic Bacteria.csv")


```

### Phylosymbiotic Bacteria Correlated to Dysbiosis
```{r}
#A phylosymbiotic bacteria is hypothesized to contribute beneficial and stabilizing functions to the coral host/microbiome. 
# We can explore this by using DYSBIOSIS metrics as a proxy for STABILITY. 
# I am also curious if some coral species do not have phylosymbiotic bacteria. This could be represented by low abundance in any bacteria marked as phylosymbiotic.

Phylosymbiotic_Bacteria <- read.csv("Phylosymbiotic Bacteria.csv")

avg.EV_RR <- read.csv("avg.EV_RR.csv") #avg.EVE_RR <- bacteria_results[,c(5:273)] %>% group_by(RelativeRisk) %>% summarise_all(list(mean))#this wasnt working so i read in a previous version of the file. 
avg.EV_RR_NA <- read.csv("avg.EV_RR_NA.csv") # Removed Dybsiosis metric for Mcav in Dysbiosis C2I ( no mcav were infected) and for ofav in DysbiosisC2E (no ofav survived exposure)

#tetracoccus.cechii
cor.test(avg.EV_RR$tetracoccus.cechii,avg.EV_RR$RelativeRisk, method=c("pearson")) #p-value = 0.2954
cor.test(avg.EV_RR$tetracoccus.cechii,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.9074, df = 5, p-value = 0.0335, cor = 0.7926804
cor.test(avg.EV_RR_NA$tetracoccus.cechii,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.4845, df = 4, p-value = 0.06788, cor = 0.7789662
cor.test(avg.EV_RR$tetracoccus.cechii,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #p-value = 0.3274

#roseicyclus.spp
cor.test(avg.EV_RR$roseicyclus.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #p-value = 0.3222
cor.test(avg.EV_RR$roseicyclus.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.8117, df = 5, p-value = 0.03747, cor = 0.7826751
cor.test(avg.EV_RR_NA$roseicyclus.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.3778, df = 4, p-value = 0.07617, cor = 0.7652916
cor.test(avg.EV_RR$roseicyclus.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #p-value = 0.3176

#marinobacter.salsuginis
cor.test(avg.EV_RR$marinobacter.salsuginis,avg.EV_RR$RelativeRisk, method=c("pearson")) #p-value = 0.3551
cor.test(avg.EV_RR$marinobacter.salsuginis,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.9548, df = 5, p-value = 0.03171, cor = 0.7974093
cor.test(avg.EV_RR_NA$marinobacter.salsuginis,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) # t = 2.3022, df = 4, p-value = 0.08274, cor = 0.7549162
cor.test(avg.EV_RR$marinobacter.salsuginis,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#aestuariibacter.spp
cor.test(avg.EV_RR$aestuariibacter.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$aestuariibacter.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.2662, df = 5, p-value = 0.07278, cor = 0.7118302
cor.test(avg.EV_RR_NA$aestuariibacter.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #NS
cor.test(avg.EV_RR$aestuariibacter.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#alteromonas.macleodii
cor.test(avg.EV_RR$alteromonas.macleodii,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$alteromonas.macleodii,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.8325, df = 5, p-value = 0.03657, cor = 0.7848989
cor.test(avg.EV_RR_NA$alteromonas.macleodii,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.2006, df = 4, p-value = 0.09259, cor = 0.7400356
cor.test(avg.EV_RR$alteromonas.macleodii,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#hyphomonas.spp
cor.test(avg.EV_RR$hyphomonas.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$hyphomonas.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 3.1565, df = 5, p-value = 0.02519, cor = 0.8160022 
cor.test(avg.EV_RR_NA$hyphomonas.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.4501, df = 4, p-value = 0.07044, cor = 0.7746705
cor.test(avg.EV_RR$hyphomonas.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#methylococcus.capsulatus
cor.test(avg.EV_RR$methylococcus.capsulatus,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$methylococcus.capsulatus,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.6023, df = 5, p-value = 0.04811, cor = 0.7584637
cor.test(avg.EV_RR_NA$methylococcus.capsulatus,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #NS
cor.test(avg.EV_RR$methylococcus.capsulatus,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#erythrobacter.citreus
cor.test(avg.EV_RR$erythrobacter.citreus,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$erythrobacter.citreus,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.9469, df = 5, p-value = 0.032, cor = 0.7966314
cor.test(avg.EV_RR_NA$erythrobacter.citreus,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.208, df = 4, p-value = 0.09183, cor = 0.7411541 
cor.test(avg.EV_RR$erythrobacter.citreus,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#novispirillum.spp
cor.test(avg.EV_RR$novispirillum.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$novispirillum.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 3.1719, df = 5, p-value = 0.02476, cor = 0.8173204
cor.test(avg.EV_RR_NA$novispirillum.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.4903, df = 4, p-value = 0.06746, cor = 0.7796865
cor.test(avg.EV_RR$novispirillum.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#paracoccus.sp*
cor.test(avg.EV_RR$paracoccus.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$paracoccus.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 3.7966, df = 5, p-value = 0.01267, cor = 0.8616579 
cor.test(avg.EV_RR_NA$paracoccus.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.9247, df = 4, p-value = 0.04304, cor = 0.8254526 
cor.test(avg.EV_RR$paracoccus.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#oleiphilus.spp
cor.test(avg.EV_RR$oleiphilus.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$oleiphilus.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 3.2234, df = 5, p-value = 0.02338, cor = 0.8216563
cor.test(avg.EV_RR_NA$oleiphilus.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.6025, df = 4, p-value = 0.05989, cor = 0.7929117
cor.test(avg.EV_RR$oleiphilus.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#micavibrio.sp
cor.test(avg.EV_RR$micavibrio.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$micavibrio.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 3.2965, df = 5, p-value = 0.02156, cor = 0.8275704
cor.test(avg.EV_RR_NA$micavibrio.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.7528, df = 4, p-value = 0.05123, cor = 0.809023 
cor.test(avg.EV_RR$micavibrio.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#pseudomonas.marincola
cor.test(avg.EV_RR$pseudomonas.marincola,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$pseudomonas.marincola,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 3.0894, df = 5, p-value = 0.02718, cor = 0.8100808
cor.test(avg.EV_RR_NA$pseudomonas.marincola,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.446, df = 4, p-value = 0.07075, cor = 0.7741492 
cor.test(avg.EV_RR$pseudomonas.marincola,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#lewinella.spp
cor.test(avg.EV_RR$lewinella.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$lewinella.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.8976, df = 5, p-value = 0.03388, cor = 0.7916843
cor.test(avg.EV_RR_NA$lewinella.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.4904, df = 4, p-value = 0.06746, cor =0.7796887
cor.test(avg.EV_RR$lewinella.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#cyanobacterium.spp*
cor.test(avg.EV_RR$cyanobacterium.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$cyanobacterium.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.6029, df = 5, p-value = 0.04808, cor = 0.7585396
cor.test(avg.EV_RR_NA$cyanobacterium.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.6029, df = 5, p-value = 0.04808, cor = 0.7585396 
cor.test(avg.EV_RR$cyanobacterium.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#granulosicoccus.coccoides
cor.test(avg.EV_RR$granulosicoccus.coccoides,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$granulosicoccus.coccoides,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.3287, df = 5, p-value = 0.06732, cor = 0.7213018 
cor.test(avg.EV_RR_NA$granulosicoccus.coccoides,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #NS
cor.test(avg.EV_RR$granulosicoccus.coccoides,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#alcanivorax.spp
cor.test(avg.EV_RR$alcanivorax.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$alcanivorax.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.7119, df = 5, p-value = 0.04218, cor = 0.7715509 
cor.test(avg.EV_RR_NA$alcanivorax.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #NS
cor.test(avg.EV_RR$alcanivorax.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#algimonas.porphyrae
cor.test(avg.EV_RR$algimonas.porphyrae,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$algimonas.porphyrae,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.9824, df = 5, p-value = 0.03071, cor = 0.8000977
cor.test(avg.EV_RR_NA$algimonas.porphyrae,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.209, df = 4, p-value = 0.09172, cor = 0.7413123 
cor.test(avg.EV_RR$algimonas.porphyrae,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#oleibacter.marinus
cor.test(avg.EV_RR$oleibacter.marinus,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$oleibacter.marinus,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.9746, df = 5, p-value = 0.03099, cor = 0.7993426 
cor.test(avg.EV_RR_NA$oleibacter.marinus,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.2408, df = 4, p-value = 0.08854, cor = 0.746063 
cor.test(avg.EV_RR$oleibacter.marinus,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#micavibrio.spp
cor.test(avg.EV_RR$micavibrio.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$micavibrio.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.773, df = 5, p-value = 0.03923, cor = 0.7784432
cor.test(avg.EV_RR_NA$micavibrio.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.3064, df = 4, p-value = 0.08236, cor = 0.7555015
cor.test(avg.EV_RR$micavibrio.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#reinekea.sp
cor.test(avg.EV_RR$reinekea.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$reinekea.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 3.3591, df = 5, p-value = 0.02013, cor = 0.8324311 
cor.test(avg.EV_RR_NA$reinekea.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.5159, df = 4, p-value = 0.06564, cor = 0.782799
cor.test(avg.EV_RR$reinekea.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#ruegeria.silicibacter.sp*
cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 4.4764, df = 5, p-value = 0.00654, cor = 0.8945984
cor.test(avg.EV_RR_NA$ruegeria.silicibacter.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 3.3375, df = 4, p-value = 0.0289, cor = 0.8577735 
cor.test(avg.EV_RR$ruegeria.silicibacter.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#kordiimonas.sp*
cor.test(avg.EV_RR$kordiimonas.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$kordiimonas.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 4.7466, df = 5, p-value = 0.00512, cor = 0.9046441 
cor.test(avg.EV_RR_NA$kordiimonas.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 3.6531, df = 4, p-value = 0.02171, cor = 0.87715
cor.test(avg.EV_RR$kordiimonas.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#alteromonas.marina*
cor.test(avg.EV_RR$alteromonas.marina,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$alteromonas.marina,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 4.402, df = 5, p-value = 0.007009, cor = 0.8915683
cor.test(avg.EV_RR_NA$alteromonas.marina,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 3.3147, df = 4, p-value = 0.02953, cor = 0.8562145 
cor.test(avg.EV_RR$alteromonas.marina,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#kordia.jejudonensis
cor.test(avg.EV_RR$kordia.jejudonensis,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$kordia.jejudonensis,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.7429, df = 5, p-value = 0.04065, cor = 0.7750858
cor.test(avg.EV_RR_NA$kordia.jejudonensis,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #NS
cor.test(avg.EV_RR$kordia.jejudonensis,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#roseivirga.sp
cor.test(avg.EV_RR$roseivirga.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$roseivirga.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 3.1261, df = 5, p-value = 0.02607, cor = 0.8133509 
cor.test(avg.EV_RR_NA$roseivirga.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 2.3612, df = 4, p-value = 0.07756, cor = 0.7630612
cor.test(avg.EV_RR$roseivirga.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#maricaulis.maris*
cor.test(avg.EV_RR$maricaulis.maris,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$maricaulis.maris,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 4.2509, df = 5, p-value = 0.008085, cor = 0.8850261 
cor.test(avg.EV_RR_NA$maricaulis.maris,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 3.1927, df = 4, p-value = 0.03313, cor = 0.8474505 
cor.test(avg.EV_RR$maricaulis.maris,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#thalassomonas.sp*
cor.test(avg.EV_RR$thalassomonas.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$thalassomonas.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 4.1678, df = 5, p-value = 0.008758, cor = 0.8811867
cor.test(avg.EV_RR_NA$thalassomonas.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 3.0786, df = 4, p-value = 0.03698, cor = 0.8385809
cor.test(avg.EV_RR$thalassomonas.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#reichenbachiella.sp
cor.test(avg.EV_RR$reichenbachiella.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$reichenbachiella.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 2.178, df = 5, p-value = 0.08132, cor = 0.6977459
cor.test(avg.EV_RR_NA$reichenbachiella.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #NS
cor.test(avg.EV_RR$reichenbachiella.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#neptuniibacter.sp*
cor.test(avg.EV_RR$neptuniibacter.sp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$neptuniibacter.sp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 4.783, df = 5, p-value = 0.004958, cor = 0.9058914
cor.test(avg.EV_RR_NA$neptuniibacter.sp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 3.5776, df = 4, p-value = 0.02322, cor = 0.8728644 
cor.test(avg.EV_RR$neptuniibacter.sp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#aliiglaciecola.glaciecola.lipolytica
cor.test(avg.EV_RR$aliiglaciecola.glaciecola.lipolytica,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$aliiglaciecola.glaciecola.lipolytica,avg.EV_RR$DysbiosisC2I, method=c("pearson")) # NS
cor.test(avg.EV_RR_NA$aliiglaciecola.glaciecola.lipolytica,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #NS
cor.test(avg.EV_RR$aliiglaciecola.glaciecola.lipolytica,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#celeribacter.neptunius*
cor.test(avg.EV_RR$celeribacter.neptunius,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$celeribacter.neptunius,avg.EV_RR$DysbiosisC2I, method=c("pearson")) #t = 4.6016, df = 5, p-value = 0.005832, cor = 0.8994316
cor.test(avg.EV_RR_NA$celeribacter.neptunius,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #t = 3.4402, df = 4, p-value = 0.02629, cor = 0.8645194
cor.test(avg.EV_RR$celeribacter.neptunius,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

#spongiibacter.spp
cor.test(avg.EV_RR$spongiibacter.spp,avg.EV_RR$RelativeRisk, method=c("pearson")) #
cor.test(avg.EV_RR$spongiibacter.spp,avg.EV_RR$DysbiosisC2I, method=c("pearson")) # NS
cor.test(avg.EV_RR_NA$spongiibacter.spp,avg.EV_RR_NA$DysbiosisC2I, method=c("pearson")) #NS
cor.test(avg.EV_RR$spongiibacter.spp,avg.EV_RR$DysbiosisC2E, method=c("pearson")) #

# The Nine Microbes significantly correlated to DysbiosisC2I in avg.EV_RR_NA
#Paracoccus sp. https://link.springer.com/article/10.1007/s00294-016-0658-3
#cyanobacterium.spp
#ruegeria.silicibacter.sp
#kordiimonas.sp
#alteromonas.marina
#maricaulis.maris
#thalassomonas.sp
#neptuniibacter.sp
#celeribacter.neptunius

```


### Correlation to RR and Dybsiosis - For Loop - Figure 3, Supp Data 2
```{r}
#EXAMPLE
# Example Pearson Correlation

cor.test(avg.EV_RR_NA$tetracoccus.cechii,avg.EV_RR_NA$RelativeRisk, method=c("pearson"))

#My code

#Identify the phylobacteriatic tree
results <- read.csv("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/EVE_results.csv")

library(tidyverse)
results_avg <- results %>% group_by(Species) %>% summarize_if(is.numeric, mean, na.rm = T)
results_avg <- results_avg[,!colnames(results_avg) %in% c("Bacteria","Treatment", "Status","Genotype")]

corr_bacteria <- names(results_avg[,-1:-4])
Species <- colnames(results_avg[,1])

all_outputs <- data.frame() # Create empty dataframe, outputs will be added to this with each loop iteration

#i = 1
for (i in 1:length(corr_bacteria)){ #iterates loop over each column of bacteria
  trip <- corr_bacteria[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip)
  data <- results_avg[,colnames(results_avg) %in% c(Species, "RelativeRisk", "DysbiosisC2I", "DysbiosisC2E", bacteria)] # create small df for current iteration
  colnames(data)[5] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Species # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(corr_bacteria), " bacteria"))
  
  # CORRELATION CODE
  # We chose Lineage Specific Bacterial Abundance as the trait we are testing for phylobacteriatic signal. Then, we do the test with 999 randomizations:
  RR_output <- cor.test(data$bacteria,data$RelativeRisk, method=c("pearson"))
  
  #Put the output into a dataframe and bind it with the previous outputs
  RR_output_sub <- data.frame(phenotype = RR_output$data.name, bacteria = bacteria, pvalue_RR = RR_output$p.value, cor_RR = RR_output$estimate, df_RR = RR_output$parameter, t_RR = RR_output$statistic)
  
  
  C2I_output <- cor.test(data$bacteria,data$DysbiosisC2I, method=c("pearson"))
 C2I_output_sub <- data.frame(phenotype = C2I_output$data.name, bacteria = bacteria, pvalue_C2I = C2I_output$p.value, cor_C2I = C2I_output$estimate, df_C2I = C2I_output$parameter, t_C2I = C2I_output$statistic)
  

 C2E_output <- cor.test(data$bacteria,data$DysbiosisC2E, method=c("pearson"))
  C2E_output_sub <- data.frame(phenotype = C2E_output$data.name, bacteria = bacteria, pvalue_C2E = C2E_output$p.value, cor_C2E = C2E_output$estimate, df_C2E = C2E_output$parameter, t_C2E = C2E_output$statistic)
  all_outputs <- rbind(all_outputs, RR_output_sub,C2I_output_sub,C2E_output_sub)
  # DUNZO
}

View(all_outputs)
write.csv(all_outputs, file="EVE Bacteria Correlated to RR and Dysbiosis.csv")
```
### Summary of Correlation Results
```{r}

Corr_results <- read.csv("EVE Bacteria Correlated to RR and Dysbiosis.csv") #Cleaned it up. 
sig_bacteria <- read.csv("Signficant EVE Bacteria.csv")
sig_bacteria$Bacteria <-  trimws(sig_bacteria$Bacteria, which = "right")
EVE_Bacteria_Correlates <- merge(Corr_results,sig_bacteria, by ="Bacteria")
write.csv(EVE_Bacteria_Correlates, file="EVE_Bacteria_Correlates.csv")

Phylosymbiotic_Bacteria <- read.csv("Phylosymbiotic Bacteria.csv")
Phylosymbiotic_correlates <- merge(Phylosymbiotic_Bacteria, Corr_results, by="Bacteria")
write.csv(Phylosymbiotic_correlates, file="Phylosymbiotic Correlates.csv")

#The Nine Microbes significantly correlated to DysbiosisC2I in avg.EV_RR_NA
# Paracoccus sp.
# cyanobacterium.spp
# ruegeria.silicibacter.sp
# kordiimonas.sp
# alteromonas.marina
# maricaulis.maris
# thalassomonas.sp
# neptuniibacter.sp
# celeribacter.neptunius

sig_LS_RR <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$phenotype == "data$bacteria and data$RelativeRisk" &
                                       EVE_Bacteria_Correlates$Expression == "LS" & EVE_Bacteria_Correlates$pvalue <0.05,]
# 7 LS bacteria are correlated to RR. All are highest abundance in Ofav. 
# bellilinea spp
# treponema brennaborense
# pandoraea spp
# uncultured candidatus planktophila sp
# marinitalea sucinacia
# treponema sp
# rubrobacter spp

sig_HV_RR <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$phenotype == "data$bacteria and data$RelativeRisk" &
                                       EVE_Bacteria_Correlates$Expression == "HV" & EVE_Bacteria_Correlates$pvalue <0.05,]
# 2 HV bacteria are correlated to RR. Both are in just two samples it seems. Not very biologically relevant. 
# bacillus decolorationis
# spirulina subsalsa

sig_LS_C2I <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$phenotype == "data$bacteria and data$DysbiosisC2I" &
                                       EVE_Bacteria_Correlates$Expression == "LS" & EVE_Bacteria_Correlates$pvalue <0.05,]
sig_LS_C2I <- na.omit(sig_LS_C2I)
# 25 LS bacteria are correlated to DysbiosisC2I
# ponticoccus sp
# leisingera aquimarina
# donghicola eburneus
# rhodovulum spp
# roseobacter spp
# kordiimonas sp
# thalassobius sp
# kordiimonas spp
# neptuniibacter sp
# celeribacter neptunius
# thalassobius aestuarii
# ruegeria silicibacter sp
# alteromonas marina
# pirellula spp
# maricaulis maris
# arthrobacter globiformis
# maritimibacter sp
# thalassomonas sp
# arthrobacter sp
# thalassobius gelatinovorus
# sneathiella spp
# paracoccus sp
# oceanicaulis spp
# alteromonas spp
# nautella italica

sig_HV_C2I <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$phenotype == "data$bacteria and data$DysbiosisC2I" &
                                       EVE_Bacteria_Correlates$Expression == "HV" & EVE_Bacteria_Correlates$pvalue <0.05,]
sig_HV_C2I <- na.omit(sig_HV_C2I)
# 2 HV bacteria are correlated to DysbiosisC2I
# acidaminobacter sp
# microbacterium aquimaris

sig_LS_C2E <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$phenotype == "data$bacteria and data$DysbiosisC2E" &
                                       EVE_Bacteria_Correlates$Expression == "LS" & EVE_Bacteria_Correlates$pvalue <0.05,]
sig_LS_C2E <- na.omit(sig_LS_C2E)
# 0 LS bacteria correlated to DysbiosisC2E

sig_HV_C2E <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$phenotype == "data$bacteria and data$DysbiosisC2E" &
                                       EVE_Bacteria_Correlates$Expression == "HV" & EVE_Bacteria_Correlates$pvalue <0.05,]
sig_HV_C2E <- na.omit(sig_HV_C2E)
# 2 HV bacteria are correlated to DysbiosisC2E
# pseudomonas fulva
# tepidiphilus petrobacter sp

# The 8 Phylosymbiotic bacteria correlated to DysbisoisC2I
# kordiimonas sp
# neptuniibacter sp
# celeribacter neptunius
# ruegeria silicibacter sp
# alteromonas marina
# maricaulis maris
# thalassomonas sp
# paracoccus sp

# The 3 Phylosymbiotic bacteria correlated to DysbisoisC2I AND in Joe's Bacterial Networks
# ruegeria silicibacter sp
# maricaulis maris
# thalassomonas sp

#Discussion: It remains to be seen if LS bacteria associated with these disease phenotypes are beneficial or dormant opportunists. This will be interesting to disentagle. Start by visualizing them side by side based on category (LS corr RR). I wonder if suscetpible species have LS dormant pathogens that make them vulnerable. 
```
### Bacterial Networks
```{r}
#These networks values are from Joe Mruzeks analysis where he assessed the network connectivity of the microbiome between treatments among species and identified the microbes central to that network connectivity. We can look for these microbes in our phylosymbiotic dataset. 

network <- read.csv("Bacterial Networks Unique.csv") #this is joe mruzeks data output from his network analysis.
Phylosymbiotic_Bacteria <- read.csv("Phylosymbiotic Bacteria.csv")
network$Bacteria <-  trimws(network$Bacteria, which = "right")
Phylosymbiotic_Bacteria$Bacteria <-  trimws(Phylosymbiotic_Bacteria$Bacteria, which = "right")
phylosymbiotic_network <- merge(network, Phylosymbiotic_Bacteria, by="Bacteria")
write.csv(phylosymbiotic_network, file="phylosymbiotic_network.csv")


```
### Correlation to RR and Dybsiosis - For Loop - Figure 3, Supp Data 2
```{r}
#EXAMPLE
# Example Pearson Correlation

cor.test(avg.EV_RR_NA$tetracoccus.cechii,avg.EV_RR_NA$RelativeRisk, method=c("pearson"))

#My code

#Identify the phylobacteriatic tree
results <- read.csv("EVE_results.csv")

results_avg <- results %>% group_by(Species) %>% summarize_if(is.numeric, mean, na.rm = T)
results_avg <- results_avg[,!colnames(results_avg) %in% c("Bacteria","Treatment", "Status","Genotype")]

corr_bacteria <- names(results_avg[,-1:-4])
Species <- colnames(results_avg[,1])

all_outputs <- data.frame() # Create empty dataframe, outputs will be added to this with each loop iteration

#i = 1
for (i in 1:length(corr_bacteria)){ #iterates loop over each column of bacteria
  trip <- corr_bacteria[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip)
  data <- results_avg[,colnames(results_avg) %in% c(Species, "RelativeRisk", "DysbiosisC2I", "DysbiosisC2E", bacteria)] # create small df for current iteration
  colnames(data)[5] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Species # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(corr_bacteria), " bacteria"))
  
  # CORRELATION CODE
  # We chose Lineage Specific Bacterial Abundance as the trait we are testing for phylobacteriatic signal. Then, we do the test with 999 randomizations:
  RR_output <- cor.test(data$bacteria,data$RelativeRisk, method=c("pearson"))
  
  #Put the output into a dataframe and bind it with the previous outputs
  RR_output_sub <- data.frame(phenotype = RR_output$data.name, bacteria = bacteria, pvalue_RR = RR_output$p.value, cor_RR = RR_output$estimate, df_RR = RR_output$parameter, t_RR = RR_output$statistic)
  
  
  C2I_output <- cor.test(data$bacteria,data$DysbiosisC2I, method=c("pearson"))
 C2I_output_sub <- data.frame(phenotype = C2I_output$data.name, bacteria = bacteria, pvalue_C2I = C2I_output$p.value, cor_C2I = C2I_output$estimate, df_C2I = C2I_output$parameter, t_C2I = C2I_output$statistic)
  

 C2E_output <- cor.test(data$bacteria,data$DysbiosisC2E, method=c("pearson"))
  C2E_output_sub <- data.frame(phenotype = C2E_output$data.name, bacteria = bacteria, pvalue_C2E = C2E_output$p.value, cor_C2E = C2E_output$estimate, df_C2E = C2E_output$parameter, t_C2E = C2E_output$statistic)
  all_outputs <- rbind(all_outputs, RR_output_sub,C2I_output_sub,C2E_output_sub)
  # DUNZO
}

View(all_outputs)
write.csv(all_outputs, file="EVE Bacteria Correlated to RR and Dysbiosis.csv")
```
### Box and Whisker Construction Individual - For Loop
```{r}
# For Loop to create multiple bar graphs (of the same format) from a single csv file

#### Background info ####

# Desired boxplot figure output:
#   x-axis: coral species
#   y-axis: bacterium abundance
#   Colors/groups: based on species susceptibility
#   Individual figures for each: Bacteria species

# Load necessary libraries
library(ggplot2) 

# Read in data
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')

# Example figure
endozoicomonas.spp_species <- ggplot(data = bacteria_results, 
                                     aes(x = reorder(Species, (endozoicomonas.spp)), y = endozoicomonas.spp, fill = Type))+
  geom_point(aes(y = endozoicomonas.spp, color = Type), 
             position =position_dodge(width = 0.47),size = 2, alpha = 2) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
  labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", title = "Endozoicomonas spp.", subtitle = "LRT = 22.0 | Beta = 0.42 | Alpha = 7.05 | pval = 2.68e-6") +
  scale_fill_manual(breaks = c("O.faveolata_c", "O.faveolata_i","O.annularis_c","O.annularis_e","O.annularis_i","C.natans_c","C.natans_e","C.natans_i", "S.siderea_c","S.siderea_e","S.siderea_i", "P.porites_c","P.porites_e","P.porites_i", "P.astreoides_c","P.astreoides_e","P.astreoides_i", "M.cavernosa_c","M.cavernosa_e"), 
                    values=c("mediumseagreen", "salmon","mediumseagreen","cornflowerblue","salmon", "mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","salmon"))+
  scale_color_manual(breaks = c("O.faveolata_c", "O.faveolata_i","O.annularis_c","O.annularis_e","O.annularis_i","C.natans_c","C.natans_e","C.natans_i", "S.siderea_c","S.siderea_e","S.siderea_i", "P.porites_c","P.porites_e","P.porites_i", "P.astreoides_c","P.astreoides_e","P.astreoides_i", "M.cavernosa_c","M.cavernosa_e"), 
                     values=c("mediumseagreen", "salmon","mediumseagreen","cornflowerblue","salmon", "mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","salmon"))+
  theme(legend.position= "Null", legend.text=element_text(size=20), 
        legend.title = element_text(size=20), 
        axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
        axis.text.y = element_text(colour = "black", size=10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        strip.text.x = element_text(colour = "black", size = 20))
endozoicomonas.spp_species 

#### Original code  ####

# For loop prerequites (for this loop, not all loops):
#   1. All figures created should have the same format OR 
#      should have distinguishable features to be used in if/else statements 
#      to direct the loop through alternate chucks of code
#   2. All figures should come from the same csv file and columns should 
#      have a uniform naming convention over which the loop interates

# STEP 1: Prior to running any for loop create a folder on your machine into 
#         which you will put all of your completed figures
#         Ex) "C:/Users/Amy/Documents/Personal/For Nick/Bacterial_abundance_boxplots"

# STEP 2: Load any needed libraries to run code (see above)
#         (require() = If library is not loaded, load it. If it is loaded, do nothing.)
require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

# STEP 3: Set working directory to location of needed files and read in data files
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
sig_bacteria <- read.csv("Signficant EVE Bacteria.csv", header = T, fileEncoding = 'UTF-8-BOM')

# STEP 4: Do some initial data formatting and object creation to make for loop smoother
# Create object for constant plot theme info 
# This chuck of code does not change regardless of the plot
my_fill <- scale_fill_manual(breaks =c("O.faveolata_c","O.faveolata_i","O.annularis_c","O.annularis_e","O.annularis_i","C.natans_c","C.natans_e","C.natans_i", "S.siderea_c","S.siderea_e","S.siderea_i", "P.porites_c","P.porites_e","P.porites_i", "P.astreoides_c","P.astreoides_e","P.astreoides_i", "M.cavernosa_c","M.cavernosa_e"), 
                    values=c("mediumseagreen", "salmon","mediumseagreen","cornflowerblue","salmon", "mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","salmon"))
my_col <- scale_color_manual(breaks = c("O.faveolata_c", "O.faveolata_i","O.annularis_c","O.annularis_e","O.annularis_i","C.natans_c","C.natans_e","C.natans_i", "S.siderea_c","S.siderea_e","S.siderea_i", "P.porites_c","P.porites_e","P.porites_i", "P.astreoides_c","P.astreoides_e","P.astreoides_i", "M.cavernosa_c","M.cavernosa_e"), 
                     values=c("mediumseagreen", "salmon","mediumseagreen","cornflowerblue","salmon", "mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","salmon"))
my_theme <- theme(legend.position= "Null", legend.text=element_text(size=20), 
        legend.title = element_text(size=20), 
        axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
        axis.text.y = element_text(colour = "black", size=10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        strip.text.x = element_text(colour = "black", size = 20))

# Pull out all metadata columns from dataframe
metadata <- colnames(bacteria_results[,colnames(bacteria_results) %in% c("Bacteria","Species","Treatment","RelativeRisk","Status","Genotype", "Type", "DysbiosisC2I", "DysbiosisC2E")])
# Pull out all bacterial abundance columns from dataframe
bacteria <- colnames(bacteria_results[,!colnames(bacteria_results) %in% metadata])
# Pull out the relevant columns from the significant eve bacteria results
eve <- sig_bacteria[,colnames(sig_bacteria) %in% c("Bacteria","LRT_LS","alpha_phylosymbiosis","beta_HEV","pval","Expression")]
# rename columns for ease
colnames(eve)[2:5] <- c("LRT","Alpha","Beta","Pattern")

# Some of the Eve bacteria names end in a space
# Remove the trailing space and now the names should match up
eve$Bacteria <-  trimws(eve$Bacteria, which = "right") #cut out space after species name.

# Make sure we can match the bacteria names in the eve and bacteria vectors
bacteria_names1 <- eve$Bacteria # create vector of bacteria names
bacteria_names2 <- str_replace_all(bacteria, "[.]", " ") # create vector of names but replace '.' with spaces

bacteria_names1[!bacteria_names1 %in% bacteria_names2] # find mismatching cases. Ideally there would be 0. 
# THERE ARE 0

# ****************************************************************************************************
# For loop testing: Un-comment (remove hashtag for line immediately below) i = 1, this will allow you to test the code only on the first loop

#i = 1

# IMPORTANT NOTE: While testing loops skip any lines that start with "for (....){"
# Only run lines that are INSIDE the loop brackets "{}"
# This will allow you to test and fix code running only the first loop (i.e., creating 1 figure) 
# and will prevent R from crashing if the code is buggy.
# the reason you test a for loop is to make sure it works before doing many iterations. This is why the test loop is i=1, i.e do one loop.
#****************************************************************************************************

# STEP 7: The For Loop
# The "I loop"

for (i in 1:length(bacteria)){ #iterates loop over each column of bacteria
  trip <- bacteria[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  results <- eve[eve$Bacteria %in% trip_name,] # pulls out the relevant EVE result line
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(bacteria), " bacteria"))
  
  y_val <- paste0(trip) # Create variable for current bacteria name
  data <- bacteria_results[,names(bacteria_results) %in% c(y_val, metadata)] # create small df for current iteration
  colnames(data)[10] <- "y_val" #rename bacteria abundance column
  
  # PLOT CODE
  plot <- ggplot(data = data, aes(x = reorder(Species, (y_val)), y = y_val, fill = Type)) +
    geom_point(aes(y = y_val, color = Type), 
               position = position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", 
         title = paste0(str_to_sentence(trip_name)), 
         subtitle = paste0("LRT = ", round(results$LRT, digits = 1),
                           " | Beta = ", round(results$Beta, digits = 2),
                           " | Alpha = ", round(results$Alpha, digits = 2),
                           " | pval = ", formatC(results$pval, format = 'e', digits = 2))) +
        my_fill + my_col + my_theme 

  # setwd to folder for all the plots to be placed into
  if (results$Pattern == "HV"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/HV EVE Bacteria Box Plots/Individual/")
    pdfname <- paste0(str_to_sentence(trip_name), "_HV_boxplot_Individual.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$Pattern == "LS") {
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Box Plots/Individual/")
    pdfname <- paste0(str_to_sentence(trip_name), "_LS_boxplot_Individual.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$Pattern == "PS") {
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/PS EVE Bacteria Box Plots/Individual/")
    pdfname <- paste0(str_to_sentence(trip_name), "_PS_boxplot_Individual.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  # DUNZO
}

# STEP 8: Check that pdfs are written to correct folder
# STEP 9: (Might required paid version of Adobe Acrobat? - If so, I can do this for you.)
# To combine into single pdf for easy scrolling:
#         In file explorer highlight all plots > right click > 
#         (IF available) "Combine files in Acrobat" > Combine > 
#         Rename "binder" > View > Page display > Two page scrolling
```
### Correlation Graphs
```{r}

setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
EVE_Bacteria_Correlates <- read.csv("EVE_Bacteria_Correlates.csv", header = T, fileEncoding = 'UTF-8-BOM')

bacteria_results_avg<- bacteria_results[,c(2,4,10:277)] %>% group_by(Species)  %>% summarize_all(list(mean)) 

plot <-  ggscatter(data, x = RelativeRisk, y = bellili, 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Relative Risk", ylab = "log10(abundance+1)",
          title = paste0(str_to_sentence(trip_name)), 
          subtitle = paste0("Cor = ", round(results$cor, digits = 1),
                           " | Pval = ", round(results$pvalue, digits = 2)) +
          my_fill + my_col + my_theme )
```

### Correlation Graphs - For Loop
```{r}
# For Loop to create multiple bar graphs (of the same format) from a single csv file

#### Background info ####

# Desired boxplot figure output:
#   x-axis: coral species
#   y-axis: bacterium abundance
#   Colors/groups: based on species susceptibility
#   Individual figures for each: Bacteria species

# Load necessary libraries
library(ggplot2) 
library("ggpubr")
require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

# Read in data
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
EVE_Bacteria_Correlates <- read.csv("EVE_Bacteria_Correlates.csv", header = T, fileEncoding = 'UTF-8-BOM')

# Example figure
library("ggpubr")
ggscatter(my_data, x = "mpg", y = "wt", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")


#### Original code ####

# For loop prerequites (for this loop, not all loops):
#   1. All figures created should have the same format OR 
#      should have distinguishable features to be used in if/else statements 
#      to direct the loop through alternate chucks of code
#   2. All figures should come from the same csv file and columns should 
#      have a uniform naming convention over which the loop interates

# STEP 1: Prior to running any for loop create a folder on your machine into 
#         which you will put all of your completed figures
#         Ex) "C:/Users/Amy/Documents/Personal/For Nick/Bacterial_abundance_boxplots"

# STEP 2: Load any needed libraries to run code (see above)
#         (require() = If library is not loaded, load it. If it is loaded, do nothing.)
require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

# STEP 3: Set working directory to location of needed files and read in data files
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
EVE_Bacteria_Correlates <- read.csv("EVE_Bacteria_Correlates.csv", header = T, fileEncoding = 'UTF-8-BOM')

# STEP 4: Do some initial data formatting and object creation to make for loop smoother
# Create object for constant plot theme info 
# This chuck of code does not change regardless of the plot
my_fill <- scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                                values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
my_col <- scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                     values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
my_theme <- theme(legend.position= "Null", legend.text=element_text(size=20), 
        legend.title = element_text(size=20), 
        axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
        axis.text.y = element_text(colour = "black", size=10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        strip.text.x = element_text(colour = "black", size = 20))

# Pull out correlation results for Relative Risk
RelativeRiskResults <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$phenotype == "data$bacteria and data$RelativeRisk", ]
# Keep only those that have a pvalue<0.05
SigRelativeRiskResults <- RelativeRiskResults[RelativeRiskResults$pvalue <= 0.05, ] 
# Calculate bacterial abundance averages based on species
bacteria_results_avg<- bacteria_results[,c(2,4,10:277)] %>% group_by(Species)  %>% summarize_all(list(mean)) 
# Pull out all metadata columns from dataframe
metadata <- colnames(bacteria_results_avg[,colnames(bacteria_results_avg) %in% c("Species","RelativeRisk")])
# Pull out all bacterial abundance columns from dataframe
bacteria <- colnames(bacteria_results_avg[,!colnames(bacteria_results_avg) %in% metadata])
# Pull out the relevant columns from the significant eve bacteria results
eve <- SigRelativeRiskResults[,colnames(SigRelativeRiskResults) %in% c("Bacteria","cor","pvalue")]
# rename columns for ease
#colnames(eve)[2:5] <- c("LRT","Alpha","Beta","Pattern")

# Some of the Eve bacteria names end in a space
# Remove the trailing space and now the names should match up
#eve$Bacteria <-  trimws(eve$Bacteria, which = "right") #cut out space after species name.

# Make sure we can match the bacteria names in the eve and bacteria vectors
bacteria_names1 <- eve$Bacteria # create vector of bacteria names
bacteria_names2 <- str_replace_all(bacteria, "[.]", " ") # create vector of names but replace '.' with spaces

bacteria_names_keep <- bacteria_names2[bacteria_names2 %in% bacteria_names1]

bacteria_names1[!bacteria_names1 %in% bacteria_names_keep] # find mismatching cases. Ideally there would be 0. 
# THERE ARE 0

# ****************************************************************************************************
# For loop testing: Un-comment (remove hashtag for line immediately below) i = 1, this will allow you to test the code only on the first loop

#i = 1

# IMPORTANT NOTE: While testing loops skip any lines that start with "for (....){"
# Only run lines that are INSIDE the loop brackets "{}"
# This will allow you to test and fix code running only the first loop (i.e., creating 1 figure) 
# and will prevent R from crashing if the code is buggy.
# the reason you test a for loop is to make sure it works before doing many iterations. This is why the test loop is i=1, i.e do one loop.
#****************************************************************************************************

# STEP 7: The For Loop
# The "I loop"

for (i in 1:length(bacteria_names_keep)){ #iterates loop over each column of bacteria
  trip <- bacteria_names_keep[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  trip_period_name <- str_replace_all(trip," ",".")
  results <- eve[eve$Bacteria %in% trip_name,] # pulls out the relevant EVE result line
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(bacteria_names_keep), " bacteria"))
  
  y_val <- paste0(trip_period_name) # Create variable for current bacteria name
  data <- bacteria_results_avg[,names(bacteria_results_avg) %in% c(y_val, metadata)] # create small df for current iteration
  colnames(data)[3] <- "y_val" #rename bacteria abundance column
  
  # PLOT CODE
  
  plot <- ggplot(data = data) +
    geom_smooth(method = "lm", aes(x = RelativeRisk, y = y_val), color="black")+
    geom_point(aes(x = RelativeRisk, y = y_val, col = Species), 
               #position = position_dodge(width = 0.47),
               size = 4) +
    labs(y = "log10(abundance+1)", x = "Relative Risk", col= "Species", 
         title = paste0(str_to_sentence(trip_name)), 
         subtitle = paste0("Cor = ", round(results$cor, digits = 2),
                           " | Pval = ", round(results$pvalue, digits = 3))) +
        my_col + my_theme +theme_classic() +theme(legend.position="none")

  # setwd to folder for all the plots to be placed into
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/Relative Risk/")
    pdfname <- paste0(str_to_sentence(trip_name), "_LS_CorrGraph_RR.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
     
  # DUNZO
}

# STEP 8: Check that pdfs are written to correct folder
# STEP 9: (Might required paid version of Adobe Acrobat? - If so, I can do this for you.)
# To combine into single pdf for easy scrolling:
#         In file explorer highlight all plots > right click > 
#         (IF available) "Combine files in Acrobat" > Combine > 
#         Rename "binder" > View > Page display > Two page scrolling
```
### Combining Box and Whisker and Correlation Graphs - For Loop - Incomplete
```{r}
# For Loop to create multiple bar graphs (of the same format) from a single csv file

#### Background info ####

# Desired boxplot figure output:
#   x-axis: coral species
#   y-axis: bacterium abundance
#   Colors/groups: based on species susceptibility
#   Individual figures for each: Bacteria species

# Load necessary libraries
library(ggplot2) 
library("ggpubr")
require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

# Read in data
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
EVE_Bacteria_Correlates <- read.csv("EVE_Bacteria_Correlates.csv", header = T, fileEncoding = 'UTF-8-BOM')

# Example figure
library("ggpubr")
ggscatter(my_data, x = "mpg", y = "wt", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")


#### Original code ####

# For loop prerequites (for this loop, not all loops):
#   1. All figures created should have the same format OR 
#      should have distinguishable features to be used in if/else statements 
#      to direct the loop through alternate chucks of code
#   2. All figures should come from the same csv file and columns should 
#      have a uniform naming convention over which the loop interates

# STEP 1: Prior to running any for loop create a folder on your machine into 
#         which you will put all of your completed figures
#         Ex) "C:/Users/Amy/Documents/Personal/For Nick/Bacterial_abundance_boxplots"

# STEP 2: Load any needed libraries to run code (see above)
#         (require() = If library is not loaded, load it. If it is loaded, do nothing.)
require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

# STEP 3: Set working directory to location of needed files and read in data files
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
EVE_Bacteria_Correlates <- read.csv("EVE_Bacteria_Correlates.csv", header = T, fileEncoding = 'UTF-8-BOM')

# STEP 4: Do some initial data formatting and object creation to make for loop smoother
# Create object for constant plot theme info 
# This chuck of code does not change regardless of the plot
my_fill <- scale_fill_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                                values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
my_col <- scale_color_manual(breaks = c("O.faveolata", "O.annularis", "C.natans", "S.siderea", "P.porites", "P.astreoides", "M.cavernosa"), 
                     values=c("orangered3", "indianred2", "cornflowerblue", "skyblue2", "mediumseagreen", "springgreen4", "green4"))
my_theme <- theme(legend.position= "Null", legend.text=element_text(size=20), 
        legend.title = element_text(size=20), 
        axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
        axis.text.y = element_text(colour = "black", size=10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        strip.text.x = element_text(colour = "black", size = 20))

# Pull out correlation results for Relative Risk
#RelativeRiskResults <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$phenotype == "data$bacteria and data$RelativeRisk", ]
# Keep only those that have a pvalue<0.05
SigCorrelates <- EVE_Bacteria_Correlates[EVE_Bacteria_Correlates$pvalue <= 0.05, ]
SigCorrelatesNA <- na.omit(SigCorrelates)
# Calculate bacterial abundance averages based on species
bacteria_results_avg<- bacteria_results[,c(2,4,8:277)] %>% group_by(Species)  %>% summarize_all(list(mean)) 
# Pull out all metadata columns from dataframe
metadata_corr <- colnames(bacteria_results_avg[,colnames(bacteria_results_avg) %in% c("Species","RelativeRisk","DysbiosisC2I","DysbiosisC2E")])
# Pull out all bacterial abundance columns from dataframe
bacteria <- colnames(bacteria_results_avg[,!colnames(bacteria_results_avg) %in% metadata_corr])
# Pull out the relevant columns from the significant eve bacteria results
eve <- SigCorrelatesNA[,colnames(SigCorrelatesNA) %in% c("Bacteria","Expression","phenotype","pvalue","cor")]
# rename columns for ease
colnames(eve)[2:5] <- c("Pattern","Phenotype", "Pval","R")

# Some of the Eve bacteria names end in a space
# Remove the trailing space and now the names should match up
#eve$Bacteria <-  trimws(eve$Bacteria, which = "right") #cut out space after species name.

# Make sure we can match the bacteria names in the eve and bacteria vectors
bacteria_names1 <- eve$Bacteria # create vector of bacteria names
bacteria_names2 <- str_replace_all(bacteria, "[.]", " ") # create vector of names but replace '.' with spaces

bacteria_names_keep <- bacteria_names2[bacteria_names2 %in% bacteria_names1]

bacteria_names1[!bacteria_names1 %in% bacteria_names_keep] # find mismatching cases. Ideally there would be 0. 
# THERE ARE 0

# ****************************************************************************************************
# For loop testing: Un-comment (remove hashtag for line immediately below) i = 1, this will allow you to test the code only on the first loop

#i = 1

# IMPORTANT NOTE: While testing loops skip any lines that start with "for (....){"
# Only run lines that are INSIDE the loop brackets "{}"
# This will allow you to test and fix code running only the first loop (i.e., creating 1 figure) 
# and will prevent R from crashing if the code is buggy.
# the reason you test a for loop is to make sure it works before doing many iterations. This is why the test loop is i=1, i.e do one loop.
#****************************************************************************************************

# STEP 7: The For Loop
# The "I loop"

for (i in 1:length(bacteria_names_keep)){ #iterates loop over each column of bacteria
  trip <- bacteria_names_keep[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  trip_period_name <- str_replace_all(trip," ",".")
  results <- eve[eve$Bacteria %in% trip_name,] # pulls out the relevant EVE result line
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(bacteria_names_keep), " bacteria"))
  
  y_val <- paste0(trip_period_name) # Create variable for current bacteria name
  data <- bacteria_results_avg[,names(bacteria_results_avg) %in% c(y_val, metadata_corr)] # create small df for current iteration
  colnames(data)[5] <- "y_val" #rename bacteria abundance column
  
  # PLOT CODE
  
  box <- ggplot(data = data, aes(x = reorder(Species, (y_val)), y = y_val, fill = Type)) +
    geom_point(aes(y = y_val, color = Type), 
               position = position_dodge(width = 0.47),size = 2, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", 
         title = paste0(str_to_sentence(trip_name)), 
         subtitle = paste0("LRT = ", round(results$LRT, digits = 1),
                           " | Beta = ", round(results$Beta, digits = 2),
                           " | Alpha = ", round(results$Alpha, digits = 2),
                           " | pval = ", formatC(results$pval, format = 'e', digits = 2))) +
        my_fill + my_col + my_theme 
  
  corr <- ggplot(data = data) +
    geom_smooth(method = "lm", aes(x = RelativeRisk, y = y_val), color="black")+
    geom_point(aes(x = RelativeRisk, y = y_val, col = Species), 
               #position = position_dodge(width = 0.47),
               size = 4) +
    labs(y = "log10(abundance+1)", x = "Relative Risk", col= "Species", 
         title = paste0(str_to_sentence(trip_name)), 
         subtitle = paste0("Cor = ", round(results$cor, digits = 2),
                           " | Pval = ", round(results$pvalue, digits = 3))) +
        my_col + my_theme +theme_classic() +theme(legend.position="none")
  
  plot <- grid.arrange(box, corr, nrow = 1)
    
   # setwd to folder for all the plots to be placed into
  if (results$phenotype == "data$bacteria and data$RelativeRisk") & (results$pattern == "LS"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/Relative Risk/")
    pdfname <- paste0(str_to_sentence(trip_name), "_LS_CorrGraph_RR.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$phenotype == "data$bacteria and data$DysbiosisC2I")& (results$pattern == "LS"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/DysbiosisC2I/")
    pdfname <- paste0(str_to_sentence(trip_name), "_LS_CorrGraph_C2I.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$phenotype == "data$bacteria and data$DysbiosisC2E")& (results$pattern == "LS"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/DysbiosisC2E/")
    pdfname <- paste0(str_to_sentence(trip_name), "_LS_CorrGraph_C2E.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
   
  if (results$phenotype == "data$bacteria and data$RelativeRisk") & (results$pattern == "HV"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/HV EVE Bacteria Correlation Plots/Relative Risk/")
    pdfname <- paste0(str_to_sentence(trip_name), "_HV_CorrGraph_RR.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$phenotype == "data$bacteria and data$DysbiosisC2I")& (results$pattern == "HV"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/DysbiosisC2I/")
    pdfname <- paste0(str_to_sentence(trip_name), "_HV_CorrGraph_C2I.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$phenotype == "data$bacteria and data$DysbiosisC2E")& (results$pattern == "HV"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/DysbiosisC2E/")
    pdfname <- paste0(str_to_sentence(trip_name), "_HV_CorrGraph_C2E.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
    
 if (results$phenotype == "data$bacteria and data$RelativeRisk") & (results$pattern == "PS"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/Relative Risk/")
    pdfname <- paste0(str_to_sentence(trip_name), "_PS_CorrGraph_RR.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$phenotype == "data$bacteria and data$DysbiosisC2I")& (results$pattern == "PS"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/DysbiosisC2I/")
    pdfname <- paste0(str_to_sentence(trip_name), "_PS_CorrGraph_C2I.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results$phenotype == "data$bacteria and data$DysbiosisC2E")& (results$pattern == "PS"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/LS EVE Bacteria Correlation Plots/DysbiosisC2E/")
    pdfname <- paste0(str_to_sentence(trip_name), "_PS_CorrGraph_C2E.pdf") # name pdf
    pdf(pdfname, width = 6, height = 6) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
     
  # DUNZO
}

# STEP 8: Check that pdfs are written to correct folder
# STEP 9: (Might required paid version of Adobe Acrobat? - If so, I can do this for you.)
# To combine into single pdf for easy scrolling:
#         In file explorer highlight all plots > right click > 
#         (IF available) "Combine files in Acrobat" > Combine > 
#         Rename "binder" > View > Page display > Two page scrolling
```
### Create hierarchical trees of HV, LS, and PS microbes.

Hierarchical Clustering
```{r, echo=FALSE}
library(minerva)
library(cluster)
library(dendextend)
library(factoextra)
```
``` {r}

Shared_Bacteria <- read.csv("Bacteria_EVE_input_transposed.csv")


rownames(Shared_Bacteria) <- Shared_Bacteria$Specified
gower.dist <- daisy(Shared_Bacteria[ ,11:68], metric = c("gower"))
```
```{r}
#DIVISIVE CLUSTERING
divisive.clust <- diana(as.matrix(gower.dist), 
                        diss = TRUE, keep.diss = TRUE)
plot(divisive.clust, main = "Divisive", xlab=)
as.matrix(gower.dist)
```
```{r}
#AGGLOMERATIVE CLUSTERING
aggl.clust.c <- hclust(gower.dist, method = "complete")
plot(aggl.clust.c,
     main = "Agglomerative, complete linkages", cex=0.5)
#names=c("annularis","faveolata","natans","siderastrea","astreoides","porites","cavernosa"))
#xlab=c("cavernosa", "annnularis", "faveolata", "astreoides", "porites", "siderea", "natans"))
rect.hclust(aggl.clust.c,k=2,border = 2:5)
```
```{r, echo=FALSE}
#Dendrogram Clustering
library("fpc")
library("ggplot2")
library("reshape2")
library("purrr")
library("dplyr")
library("dendextend")
```
```{r}
#Dendrogram of agglomerative hierarchy, k=3
dendro <- as.dendrogram(aggl.clust.c)
dendro.col <- dendro %>%
  set("branches_k_color", k = 2, value =   c("darkslategray", "darkslategray4", "darkslategray3")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", 
      value = c("darkslategray")) %>% 
  set("labels_cex", 0.5,
      xlab=c("cavernosa", "annnularis", "faveolata", "astreoides", "porites", "siderea", "natans"))
ggd1 <- as.ggdend(dendro.col)
ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram, k = 3")
ggplot(ggd1, labels = T) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")
```
### Hierarchical Trees for Phylosymbiotic Bacteria - Incomplete
```{r}
#Hierarchical Trees for Phylosymbiotic Bacteria
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
Phylosymbiotic_Bacteria_Transposed <- read.csv("Phylosymbiotic Bacteria_Transposed.csv")

rownames(Phylosymbiotic_Bacteria_Transposed) <- Phylosymbiotic_Bacteria_Transposed$Specified
gower.dist <- daisy(Phylosymbiotic_Bacteria_Transposed[ ,11:43], metric = c("gower"))
```
```{r}
#DIVISIVE CLUSTERING
divisive.clust <- diana(as.matrix(gower.dist), 
                        diss = TRUE, keep.diss = TRUE)
plot(divisive.clust, main = "Divisive", xlab=)
as.matrix(gower.dist)
```
```{r}
#AGGLOMERATIVE CLUSTERING
aggl.clust.c <- hclust(gower.dist, method = "complete")
plot(aggl.clust.c,
     main = "Agglomerative, complete linkages", cex=0.5)
#names=c("annularis","faveolata","natans","siderastrea","astreoides","porites","cavernosa"))
#xlab=c("cavernosa", "annnularis", "faveolata", "astreoides", "porites", "siderea", "natans"))
rect.hclust(aggl.clust.c,k=2,border = 2:5)
```
```{r, echo=FALSE}
#Dendrogram Clustering
library("fpc")
library("ggplot2")
library("reshape2")
library("purrr")
library("dplyr")
library("dendextend")
```
```{r}
#Dendrogram of agglomerative hierarchy, k=3
dendro <- as.dendrogram(aggl.clust.c)
dendro.col <- dendro %>%
  set("branches_k_color", k = 2, value =   c("darkslategray", "darkslategray4", "darkslategray3")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", 
      value = c("darkslategray")) %>% 
  set("labels_cex", 0.5,
      xlab=c("cavernosa", "annnularis", "faveolata", "astreoides", "porites", "siderea", "natans"))
ggd1 <- as.ggdend(dendro.col)
ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram, k = 3")
ggplot(ggd1, labels = T) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")
```

### Stacked Bar Plot - Phylosymbiosis
```{r}
# library
library(ggplot2)
library(dplyr)

setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
Phylosymbiotic_bacteria_transposed <- read.csv("Phylosymbiotic Bacteria_Transposed.csv", header = T, fileEncoding = 'UTF-8-BOM')

PS_Stacked_input <- read.csv("Phylosymbiotic Stacked Column input.csv")

plot <- ggplot(data = PS_Stacked_input, aes(x = reorder(Specified, (Order)), y = Abundance, fill=Bacteria)) +
  geom_bar(position="stack", stat="identity")+  #Position="stack" for real abundance, Position="fill" for relative abundance.
  scale_fill_manual(values=c33) +theme(axis.text.x = element_text(colour = "black", size=5,angle=90, hjust =1),
        legend.position= "Null")

```
### Stacked Bar Plot Loop - Lineage Specific 
```{r}
require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')

LS_EVE_transposed <- read.csv("LS EVE Bacteria transposed.csv")

metadata <- LS_EVE_transposed[,c("Species","Specified","RelativeRisk")]
LS_data <- LS_EVE_transposed[,c(11:118)]
LS_bacteria <- names(LS_data)

#v=rep(colnames(data), each=68)

all_outputs <- data.frame() # Create empty dataframe, outputs will be added to this with each loop iteration

#i = 1
for (i in 1:length(LS_bacteria)){ #iterates loop over each column of bacteria
  trip <- LS_bacteria[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip)
  data <- LS_EVE_transposed[,colnames(LS_EVE_transposed) %in% c("Species", "Specified", "RelativeRisk", bacteria)] # create small df for current iteration
  colnames(data)[4] <- "bacteria" #rename bacteria abundance column
  
  abundance <- pull(data, bacteria) # create a vector from data$bacteria
  names(abundance)<-data$Specified # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(LS_bacteria), " bacteria"))
  
  df1 <- cbind(metadata,rep(trip_name, each=68),data$bacteria)
  
  all_outputs <- rbind(all_outputs, df1)
  colnames(all_outputs)[4] <- "bacteria" #rename bacteria abundance column
  colnames(all_outputs)[5] <- "abundance" #rename bacteria abundance column
  # DUNZO
}

View(all_outputs)
write.csv(all_outputs, file="LS Bacteria Stacked Input.csv")

plot <- ggplot(data = all_outputs, aes(x = reorder(Specified, (RelativeRisk)), y = abundance, fill=bacteria)) +
  geom_bar(position="stack", stat="identity")+  #Position="stack" for real abundance, Position="fill" for relative abundance.
  theme(axis.text.x = element_text(colour = "black", size=5,angle=90, hjust =1),
        legend.position= "Null")

```
### Stacked Bar Plot - HV,LS,PS Bacteria
```{r}
#Color bacteria by EVE category (HV,LS,PS) in stacked column plot
#use species.percentages.csv for percentage values. 

require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
sig_bacteria<- read.csv("Signficant EVE Bacteria.csv", header = T, fileEncoding = 'UTF-8-BOM')
sp <- read.csv("Species.percentage.transposed.csv", header = T, fileEncoding = 'UTF-8-BOM')

sig_bacteria$Bacteria <-  trimws(sig_bacteria$Bacteria, which = "right") #cut out space after species name.
EVE_percentages <- merge(sig_bacteria, sp, by="Bacteria")
#write.csv(EVE_percentages, file= "EVE percentages.csv")
EVE_percentages_barinput <- read.csv("EVE percentages Bar input.csv")

metadata <- EVE_percentages_barinput[,c("Species","Specified","RelativeRisk")]
EVE_data <- EVE_percentages_barinput[,c(11:275)]
EVE_bacteria <- names(EVE_data)

EVE_bacteria <- str_replace_all(EVE_bacteria,"[.]"," ") # replaces '.' with space in bacteria name

#v=rep(colnames(data), each=68)

all_outputs <- data.frame() # Create empty dataframe, outputs will be added to this with each loop iteration

i = 1
for (i in 1:length(EVE_bacteria)){ #iterates loop over each column of bacteria
  trip <- EVE_bacteria[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip," ",".") # replaces '.' with space in bacteria name
 
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip_name)
  data <- EVE_percentages_barinput[,colnames(EVE_percentages_barinput) %in% c("Species", "Specified", "RelativeRisk", bacteria)] # create small df for current iteration
  colnames(data)[4] <- "bacteria" #rename bacteria abundance column
  
  abundance <- pull(data, bacteria) # create a vector from data$bacteria
  names(abundance)<-data$Specified # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(EVE_bacteria), " bacteria"))
  
  df1 <- cbind(metadata,rep(trip, each=68),data$bacteria)
  
  all_outputs <- rbind(all_outputs, df1)
  # DUNZO
}

colnames(all_outputs)[4] <- "Bacteria" #rename bacteria abundance column
colnames(all_outputs)[5] <- "Abundance" #rename bacteria abundance column
#View(all_outputs)
write.csv(all_outputs, file="EVE Bacteria Stacked Input.csv")
all_outputs <- read.csv("EVE Bacteria Stacked Input.csv")
ExpressionPattern <- sig_bacteria[,c("Bacteria","Expression")]
output <- merge(all_outputs,ExpressionPattern, by="Bacteria")
View(output)


plot <- ggplot(data = output, aes(x = reorder(Specified, (RelativeRisk)), y = Abundance, fill=Expression)) +
  geom_bar(position="stack", stat="identity")+  #Position="stack" for real abundance, Position="fill" for relative abundance.
  scale_fill_manual("legend", values = c("PS" = "steelblue1", "LS" = "mediumpurple1", "HV" = "palevioletred2"))+
  theme(axis.text.x = element_text(colour = "black", size=5,angle=90, hjust =1))
#saved as HV LS PS Stacked column"

```

### Signficant Abundance Differences of LS/PS Bacteria in Controls associated with Treatment Outcome
```{r}
QUESTION: Within a species - such as Porites astreoides - which PS/LS bacteria had the most difference in abundance between the controls whos genotypic pair got infected in disease treament vs controls whose genotypic pair in disease treatment remained healthy. Basically, we can see what is unique in disease naive samples whose genotypic fate is either resistant or susceptible.
-This will show which bacteria are present in resistant fragments and which bacteria are present in susceptible fragments. 

Steps:
  -Parse genotypes that got infected by those that remained resistant. 
  -Organize PS/LS bacteria abundance dataset and exclude HV bacteria and have a descriptive column indicating infected v exposed treatment outcome.
  -Prepare SIMPER analysis to compare infected v exposed treatment outcome. 
  -Box and whisker plot and ANOVA/tukey test most dissimilar bacterial abundances. 
  -Group bacteria by apparent abundance in resistant fragments as putatively "good" PS/LS bacteria and then another group for susceptible fragments as putatively co-evolved pathogenic bacteria.
  -Put these "good" and "bad" categorical groupings into WGCNA to see if a bacteria module is in CONSISTENT correlation with a majority of the bacteria in each categorical group
  - GO that module. This is a putative host-microbe interaction that is associated with whether that coral gets infected or not. 

#Porites astreoides
Genoytpe 1 was the only genotype to become infected. 

setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')
Past_bacteria_results <- bacteria_results[11:20,]
###SIMPER analysis
library(vegan)
coral.sim.Past=simper(Past_bacteria_results[,10:277],Past_bacteria_results$Type) #Past
summary(coral.sim.Past)

#If the simper is too large to view at once: I found this by going to the environment on the right and clicking "coral.sim" and then the magnfying glass near that
##coral.sim[["Control_Disease Band"]]
#or
##coral.sim[["Exposed astereoides_Exposed cavernosa"]][["overall"]]


#Send this simper data to an excel file 

results.test.coral.sim <- summary(coral.sim.Past)
df <- do.call(rbind.data.frame, results.test.coral.sim)
write.csv(df, "Past.EVE.simper.csv")

#####export file to excel using package xlsx
install.packages("xlsx")
library(xlsx)
install.packages("rJava")
install.packages("openxlsx")
library(openxlsx)
install.packages("writexl")

library(writexl)
library(rJava)
library(xlsx)
require(xlsx)
write.xlsx(df,"/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/Past.EVE.simper.xlsx",sheetName="Sheet1",col.names=TRUE,row.names=TRUE,append=FALSE,showNA=TRUE)

##Organize bacteria by expresion pattern. 
##Organize by those that are higher or lower in exposed v infected treatments. 
##Organize by those that are SIGNFICANTLY higher or lower in exposed v infected treatments. 
##Visualize by box and whisker by these categorical groups. So here are all the PS bacteria that seem to be pathogenic (higher in infected). Here are all the PS that seem to be probiotic (lower in infected relative to exposed). 

#I want to keep any row that contains "P.astreoides_e_P.astreoides_i" in the column "Comparison"
tibble::rownames_to_column(df, "Comparison") #turns the row names "P.astreoides_e_P.astreoides_i" into a column. 
df$Comparison <- rownames(df) 
Past_ei <- dplyr::filter(df, grepl('P.astreoides_e_P.astreoides_i', Comparison)) #filters for only a particular comparison.
Past_ei$Bacteria = substr(Past_ei$Comparison,31,50) #extracted only bacteria names. 
Past_ei$Bacteria <-  trimws(Past_ei$Bacteria, which = "right") #cut out space after species name.#remove blank spaces after bacteria name. 
require(stringr) # for string manipulation
Past_ei$Bacteria <- str_replace_all(Past_ei$Bacteria, "[.]", " ")

sig_bacteria <- read.csv("Signficant EVE Bacteria.csv")
sig_bacteria$Bacteria <-  trimws(sig_bacteria$Bacteria, which = "right")
Past_ei_exp <- merge(Past_ei, sig_bacteria, by="Bacteria")
colnames(Past_ei_exp)
Past_metadata <- Past_ei_exp[,c("Bacteria", "cumsum","LRT_LS", "alpha_phylosymbiosis","beta_HEV","Expression", "pval")]
Past_data <- Past_ei_exp[,c(38:47)]
Past_HighDissimilarity <- Past_ei_exp[Past_ei_exp$cumsum <= 0.5,]
Query <- Past_HighDissimilarity$Bacteria
one.way <- aov(fangia.hongkongensis ~ Status, data = Past_bacteria_results)
summary(one.way)
TukeyHSD(one.way)

all_outputs <- data.frame()

#i=1
for (i in 1:length(Query)){ #iterates loop over each column of bacteria
  trip <- Query[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip," ",".") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip_name)
  data <- Past_bacteria_results[,colnames(Past_bacteria_results) %in% c("Type", "Genotype", bacteria)] # create small df for current iteration
  colnames(data)[3] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Type # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(Query), " bacteria"))
  
  # ONE WAY ANOVA
  one.way <- aov(data$bacteria ~ data$Type, data = data)
  summary(one.way)
  
  
  #Put the output into a dataframe and bind it with the previous outputs
  one.way_output <- data.frame(bacteria =trip_name, AOVpvalue = summary(one.way)[[1]][["Pr(>F)"]][[1]], Test= "ONE WAY")
  
  #Tukey Test
   tukey <- TukeyHSD(one.way)
  tukey_results <- tukey$`data$Type`[,'p adj']
  tukey_output <- data.frame(bacteria=trip_name,TUKEYPval = tukey_results, Test= "TUKEY")

Test_output <- cbind(one.way_output, tukey_output)

all_outputs <- rbind(all_outputs,Test_output)
  # DUNZO
}

View(all_outputs)
write.csv(all_outputs, file="Past Dissimilar AOV TUKEY Results.csv")

#I think i just need to consider the controls and compare genotypes.
#I added a column that specifies the outcome of the paired fragment so that we can compare the controls based on the paired outcome. This was specified as Ce for control fragment whose pair was disease-exposed, and then Ci for a control fragment whose pair was disease-infected. So we will want to know the bacteria that are most dissimilar and then test for signficance between Ce v Ci within each species and then note the trend if its increasing or decreasing to suggest pathogen or probiotic. 

# Ce v Ci
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')

Ppor_bacteria_results <- bacteria_results[1:10,]
Past_bacteria_results <- bacteria_results[11:20,]
Cnat_bacteria_results <- bacteria_results[21:30,]
Mcav_bacteria_results <- bacteria_results[31:40,]
Ofav_bacteria_results <- bacteria_results[41:50,]
Oann_bacteria_results <- bacteria_results[51:58,]
Ssid_bacteria_results <- bacteria_results[59:68,]

###SIMPER analysis
library(vegan)
coral.sim=simper(bacteria_results[,12:279],bacteria_results$ControlOutcome) #Past
summary(coral.sim)

#If the simper is too large to view at once: I found this by going to the environment on the right and clicking "coral.sim" and then the magnfying glass near that
##coral.sim[["Control_Disease Band"]]
#or
##coral.sim[["Exposed astereoides_Exposed cavernosa"]][["overall"]]


#Send this simper data to an excel file 

results.test.coral.sim <- summary(coral.sim)
df <- do.call(rbind.data.frame, results.test.coral.sim)
write.csv(df, "EVE.byControlOutcome.simper.csv")


##Organize bacteria by expresion pattern. 
##Organize by those that are higher or lower in exposed v infected treatments. 
##Organize by those that are SIGNFICANTLY higher or lower in exposed v infected treatments. 
##Visualize by box and whisker by these categorical groups. So here are all the PS bacteria that seem to be pathogenic (higher in infected). Here are all the PS that seem to be probiotic (lower in infected relative to exposed). 


tibble::rownames_to_column(df, "Comparison") 
df$Comparison <- rownames(df) 

#keep only these rows. 
ControlOutcome <- dplyr::filter(df, grepl('O.annularis_Ce_O.annularis_Ci|C.natans_Ce_C.natans_Ci|S.siderea_Ci_S.siderea_Ce|P.porites_Ce_P.porites_Ci|P.astreoides_Ce_P.astreoides_Ci', Comparison)) #filters for only a particular comparison.

ControlOutcome$Bacteria <-gsub("O.annularis_Ce_O.annularis_Ci.|C.natans_Ce_C.natans_Ci.|S.siderea_Ci_S.siderea_Ce.|P.porites_Ce_P.porites_Ci.|P.astreoides_Ce_P.astreoides_Ci.","",as.character(ControlOutcome$Comparison))

require(stringr) # for string manipulation
ControlOutcome$Bacteria <- str_replace_all(ControlOutcome$Bacteria, "[.]", " ")

ControlOutcome$Bacteria <-  trimws(ControlOutcome$Bacteria, which = "right") #cut out space after species name.#remove blank spaces after bacteria name. 

sig_bacteria <- read.csv("Signficant EVE Bacteria.csv")
sig_bacteria$Bacteria <-  trimws(sig_bacteria$Bacteria, which = "right")
ControlOutcome_exp <- merge(ControlOutcome, sig_bacteria, by="Bacteria")
write.csv(ControlOutcome_exp, file="ControlOutcome_counts.csv") #Cleaned up comparison column to just say coral species name

ControlOutcome_counts <- read.csv("ControlOutcome_counts.csv")
Oann_ControlOutcome <- dplyr::filter(ControlOutcome_counts, grepl('O.annularis', Comparison))
Cnat_ControlOutcome <- dplyr::filter(ControlOutcome_counts, grepl('C.natans', Comparison))
Ssid_ControlOutcome <- dplyr::filter(ControlOutcome_counts, grepl('S.siderea', Comparison))
Past_ControlOutcome <- dplyr::filter(ControlOutcome_counts, grepl('P.astreoides', Comparison))
Ppor_ControlOutcome <- dplyr::filter(ControlOutcome_counts, grepl('P.porites', Comparison))

```
### Pastreoides Control Outcome AOV TUKEY
```{r}
colnames(ControlOutcome_exp)
Past_ControlOutcome_metadata <- Past_ControlOutcome[,c("Bacteria", "cumsum","LRT_LS", "alpha_phylosymbiosis","beta_HEV","Expression", "pval")]
Past_ControlOutcome_counts <- Past_ControlOutcome[,c(38:47)]
Past_HighDissimilarity <- Past_ControlOutcome[Past_ControlOutcome$cumsum <= 0.5,]
Query <- Past_HighDissimilarity$Bacteria
#one.way <- aov(endozoicomonas spp ~ Status, data = Past_ControlOutcome[,c(38:47)])
#summary(one.way)
#TukeyHSD(one.way)

all_outputs <- data.frame()

#i=1
for (i in 1:length(Query)){ #iterates loop over each column of bacteria
  trip <- Query[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip," ",".") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip_name)
  data <- Past_bacteria_results[,colnames(Past_bacteria_results) %in% c("ControlOutcome", "Genotype", bacteria)] # create small df for current iteration
  colnames(data)[3] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Status # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(Query), " bacteria"))
  
  # ONE WAY ANOVA
  one.way <- aov(data$bacteria ~ data$ControlOutcome, data = data)
  summary(one.way)
  
  
  #Put the output into a dataframe and bind it with the previous outputs
  one.way_output <- data.frame(bacteria =trip_name, AOVpvalue = summary(one.way)[[1]][["Pr(>F)"]][[1]], Test= "ONE WAY")
  
  #Tukey Test
  tukey <- TukeyHSD(one.way)
  tukey_results <- tukey$`data$ControlOutcome`[,'p adj']
  tukey_output <- data.frame(bacteria=trip_name,TUKEYPval = tukey_results, Test= "TUKEY")

Test_output <- cbind(one.way_output, tukey_output)

all_outputs <- rbind(all_outputs,Test_output)
  # DUNZO
}

View(all_outputs)
write.csv(all_outputs, file="Past Dissimilar AOV TUKEY Results.csv")

```
### P.porites - Control Outcome AOV TUKEY
```{r}
colnames(ControlOutcome_exp)
Ppor_ControlOutcome_metadata <- Ppor_ControlOutcome[,c("Bacteria", "cumsum","LRT_LS", "alpha_phylosymbiosis","beta_HEV","Expression", "pval")]
Ppor_ControlOutcome_counts <- Ppor_ControlOutcome[,c(28:37)]
Ppor_HighDissimilarity <- Ppor_ControlOutcome[Ppor_ControlOutcome$cumsum <= 0.5,]
Query <- Ppor_HighDissimilarity$Bacteria
#one.way <- aov(endozoicomonas spp ~ Status, data = Past_ControlOutcome[,c(38:47)])
#summary(one.way)
#TukeyHSD(one.way)

all_outputs <- data.frame()

#i=1
for (i in 1:length(Query)){ #iterates loop over each column of bacteria
  trip <- Query[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip," ",".") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip_name)
  data <- Ppor_bacteria_results[,colnames(Ppor_bacteria_results) %in% c("ControlOutcome", "Genotype", bacteria)] # create small df for current iteration
  colnames(data)[3] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Status # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(Query), " bacteria"))
  
  # ONE WAY ANOVA
  one.way <- aov(data$bacteria ~ data$ControlOutcome, data = data)
  summary(one.way)
  
  
  #Put the output into a dataframe and bind it with the previous outputs
  one.way_output <- data.frame(bacteria =trip_name, AOVpvalue = summary(one.way)[[1]][["Pr(>F)"]][[1]], Test= "ONE WAY")
  
  #Tukey Test
  tukey <- TukeyHSD(one.way)
  tukey_results <- tukey$`data$ControlOutcome`[,'p adj']
  tukey_output <- data.frame(bacteria=trip_name,TUKEYPval = tukey_results, Test= "TUKEY")

Test_output <- cbind(one.way_output, tukey_output)

all_outputs <- rbind(all_outputs,Test_output)
  # DUNZO
}

View(all_outputs)
write.csv(all_outputs, file="Ppor Dissimilar AOV TUKEY Results.csv")

```
### S.siderea - Control Outcome AOV TUKEY
```{r}
colnames(ControlOutcome_exp)
Ssid_ControlOutcome_metadata <- Ssid_ControlOutcome[,c("Bacteria", "cumsum","LRT_LS", "alpha_phylosymbiosis","beta_HEV","Expression", "pval")]
Ssid_ControlOutcome_counts <- Ssid_ControlOutcome[,c(86:95)]
Ssid_HighDissimilarity <- Ssid_ControlOutcome[Ssid_ControlOutcome$cumsum <= 0.5,]
Query <- Ssid_HighDissimilarity$Bacteria
#one.way <- aov(endozoicomonas spp ~ Status, data = Past_ControlOutcome[,c(38:47)])
#summary(one.way)
#TukeyHSD(one.way)

all_outputs <- data.frame()

#i=1
for (i in 1:length(Query)){ #iterates loop over each column of bacteria
  trip <- Query[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip," ",".") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip_name)
  data <- Ssid_bacteria_results[,colnames(Ssid_bacteria_results) %in% c("ControlOutcome", "Genotype", bacteria)] # create small df for current iteration
  colnames(data)[3] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Status # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(Query), " bacteria"))
  
  # ONE WAY ANOVA
  one.way <- aov(data$bacteria ~ data$ControlOutcome, data = data)
  summary(one.way)
  
  
  #Put the output into a dataframe and bind it with the previous outputs
  one.way_output <- data.frame(bacteria =trip_name, AOVpvalue = summary(one.way)[[1]][["Pr(>F)"]][[1]], Test= "ONE WAY")
  
  #Tukey Test
  tukey <- TukeyHSD(one.way)
  tukey_results <- tukey$`data$ControlOutcome`[,'p adj']
  tukey_output <- data.frame(bacteria=trip_name,TUKEYPval = tukey_results, Test= "TUKEY")

Test_output <- cbind(one.way_output, tukey_output)

all_outputs <- rbind(all_outputs,Test_output)
  # DUNZO
}

View(all_outputs)
write.csv(all_outputs, file="Ssid Dissimilar AOV TUKEY Results.csv")

```
### C.natans - Control Outcome AOV TUKEY
```{r}
colnames(ControlOutcome_exp)
Cnat_ControlOutcome_metadata <- Cnat_ControlOutcome[,c("Bacteria", "cumsum","LRT_LS", "alpha_phylosymbiosis","beta_HEV","Expression", "pval")]
Cnat_ControlOutcome_counts <- Cnat_ControlOutcome[,c(86:95)]
Cnat_HighDissimilarity <-Cnat_ControlOutcome[Cnat_ControlOutcome$cumsum <= 0.5,]
Query <- Cnat_HighDissimilarity$Bacteria
#one.way <- aov(endozoicomonas spp ~ Status, data = Past_ControlOutcome[,c(38:47)])
#summary(one.way)
#TukeyHSD(one.way)

all_outputs <- data.frame()

#i=1
for (i in 1:length(Query)){ #iterates loop over each column of bacteria
  trip <- Query[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip," ",".") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip_name)
  data <- Cnat_bacteria_results[,colnames(Cnat_bacteria_results) %in% c("ControlOutcome", "Genotype", bacteria)] # create small df for current iteration
  colnames(data)[3] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Status # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(Query), " bacteria"))
  
  # ONE WAY ANOVA
  one.way <- aov(data$bacteria ~ data$ControlOutcome, data = data)
  summary(one.way)
  
  
  #Put the output into a dataframe and bind it with the previous outputs
  one.way_output <- data.frame(bacteria =trip_name, AOVpvalue = summary(one.way)[[1]][["Pr(>F)"]][[1]], Test= "ONE WAY")
  
  #Tukey Test
  tukey <- TukeyHSD(one.way)
  tukey_results <- tukey$`data$ControlOutcome`[,'p adj']
  tukey_output <- data.frame(bacteria=trip_name,TUKEYPval = tukey_results, Test= "TUKEY")

Test_output <- cbind(one.way_output, tukey_output)

all_outputs <- rbind(all_outputs,Test_output)
  # DUNZO
}

#View(all_outputs)
write.csv(all_outputs, file="Cnat Dissimilar AOV TUKEY Results.csv")

```
### O.annularis - Control Outcome AOV TUKEY
```{r}
colnames(ControlOutcome_exp)
Oann_ControlOutcome_metadata <-Oann_ControlOutcome[,c("Bacteria", "cumsum","LRT_LS", "alpha_phylosymbiosis","beta_HEV","Expression", "pval")]
Oann_ControlOutcome_counts <- Oann_ControlOutcome[,c(86:95)]
Oann_HighDissimilarity <-Oann_ControlOutcome[Oann_ControlOutcome$cumsum <= 0.5,]
Query <- Oann_HighDissimilarity$Bacteria
#one.way <- aov(endozoicomonas spp ~ Status, data = Past_ControlOutcome[,c(38:47)])
#summary(one.way)
#TukeyHSD(one.way)

all_outputs <- data.frame()

#i=1
for (i in 1:length(Query)){ #iterates loop over each column of bacteria
  trip <- Query[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip," ",".") # replaces '.' with space in bacteria name
  
  # First, you need to define which trait you want to test and give names to each value according to species
  bacteria <- paste0(trip_name)
  data <- Oann_bacteria_results[,colnames(Oann_bacteria_results) %in% c("ControlOutcome", "Genotype", bacteria)] # create small df for current iteration
  colnames(data)[3] <- "bacteria" #rename bacteria abundance column
  
  trait <- pull(data, bacteria) # create a vector from data$bacteria
  names(trait)<-data$Status # Name vector elements by data$Species
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(Query), " bacteria"))
  
  # ONE WAY ANOVA
  one.way <- aov(data$bacteria ~ data$ControlOutcome, data = data)
  summary(one.way)
  
  
  #Put the output into a dataframe and bind it with the previous outputs
  one.way_output <- data.frame(bacteria =trip_name, AOVpvalue = summary(one.way)[[1]][["Pr(>F)"]][[1]], Test= "ONE WAY")
  
  #Tukey Test
  tukey <- TukeyHSD(one.way)
  tukey_results <- tukey$`data$ControlOutcome`[,'p adj']
  tukey_output <- data.frame(bacteria=trip_name,TUKEYPval = tukey_results, Test= "TUKEY")

Test_output <- cbind(one.way_output, tukey_output)

all_outputs <- rbind(all_outputs,Test_output)
  # DUNZO
}

#View(all_outputs)
write.csv(all_outputs, file="Oann Dissimilar AOV TUKEY Results.csv")

```
### Keeping only Sigdiff bacteria between Ce Ci
```{r}
Oann_tukeyresults <- read.csv("Oann Dissimilar AOV TUKEY Results.csv")
Cnat_tukeyresults <- read.csv("Cnat Dissimilar AOV TUKEY Results.csv")
Ssid_tukeyresults <- read.csv("Ssid Dissimilar AOV TUKEY Results.csv")
Past_tukeyresults <- read.csv("Past Dissimilar AOV TUKEY Results.csv")
Ppor_tukeyresults <- read.csv("Ppor Dissimilar AOV TUKEY Results.csv")

#keep only these rows. 
colnames(Oann_tukeyresults)[colnames(Oann_tukeyresults)=="X"] <- "Comparison"
Oann_tukeyresults <- dplyr::filter(Oann_tukeyresults, grepl('O.annularis_Ci-O.annularis_Ce', Comparison))
Oann_tukeysig <- Oann_tukeyresults %>% filter(TUKEYPval <= 0.1) #keep only sig diff bacteria.
Oann_tukeysig <- data.frame(Oann_tukeysig, Species= "O.annularis") #Adding a column that specifies the coral species name. 

colnames(Cnat_tukeyresults)[colnames(Cnat_tukeyresults)=="X"] <- "Comparison"
Cnat_tukeyresults <- dplyr::filter(Cnat_tukeyresults, grepl('C.natans_Ci-C.natans_Ce', Comparison))
Cnat_tukeysig <- Cnat_tukeyresults %>% filter(TUKEYPval <= 0.1) #keep only sig diff bacteria.
Cnat_tukeysig <- data.frame(Cnat_tukeysig, Species= "C.natans") #Adding a column that specifies the coral species name. 

colnames(Ssid_tukeyresults)[colnames(Ssid_tukeyresults)=="X"] <- "Comparison"
Ssid_tukeyresults <- dplyr::filter(Ssid_tukeyresults, grepl('S.siderea_Ci-S.siderea_Ce', Comparison))
Ssid_tukeysig <- Ssid_tukeyresults %>% filter(TUKEYPval <= 0.1) #keep only sig diff bacteria.
Ssid_tukeysig <- data.frame(Ssid_tukeysig, Species= "S.siderea") #Adding a column that specifies the coral species name. 

colnames(Past_tukeyresults)[colnames(Past_tukeyresults)=="X"] <- "Comparison"
Past_tukeyresults <- dplyr::filter(Past_tukeyresults, grepl('P.astreoides_Ci-P.astreoides_Ce', Comparison))
Past_tukeysig <- Past_tukeyresults %>% filter(TUKEYPval <= 0.1) #keep only sig diff bacteria.
Past_tukeysig <- data.frame(Past_tukeysig, Species= "P.astreoides") #Adding a column that specifies the coral species name.

colnames(Ppor_tukeyresults)[colnames(Ppor_tukeyresults)=="X"] <- "Comparison"
Ppor_tukeyresults <- dplyr::filter(Ppor_tukeyresults, grepl('P.porites_Ci-P.porites_Ce', Comparison))
Ppor_tukeysig <- Ppor_tukeyresults %>% filter(TUKEYPval <= 0.1) #keep only sig diff bacteria.
Ppor_tukeysig <- data.frame(Ppor_tukeysig, Species= "P.porites") #Adding a column that specifies the coral species name. 

#pull out unique bacteria. just to organize for individual graphs or to throw into WGCNA. 
Oann_CEVE <- select(Oann_tukeysig,c("bacteria","AOVpvalue","TUKEYPval","Species"))
Cnat_CEVE <- select(Cnat_tukeysig,c("bacteria","AOVpvalue","TUKEYPval","Species"))
Ssid_CEVE <- select(Ssid_tukeysig,c("bacteria","AOVpvalue","TUKEYPval","Species"))
Past_CEVE <- select(Past_tukeysig,c("bacteria","AOVpvalue","TUKEYPval","Species"))
Ppor_CEVE <- select(Ppor_tukeysig,c("bacteria","AOVpvalue","TUKEYPval","Species"))

write.csv(Oann_CEVE,file="Oann_CEVE.csv")
write.csv(Cnat_CEVE,file="Cnat_CEVE.csv")
write.csv(Ssid_CEVE,file="Ssid_CEVE.csv")
write.csv(Past_CEVE,file="Past_CEVE.csv")
write.csv(Ppor_CEVE,file="Ppor_CEVE.csv")

#I think I should make individual graphs of these bacteria where the pairs are indicated by symbols. 
#We can take a previous for loop, select for only bacteria in these species_CEVE datasets and the change the symbol of the controls to match their end phenotype and output them in a separate folder by species. Then work on WGCNA. 
```
### Individual Box and Whisker For Loop - AOV/TUKEY Significant Bacteria
```{r}
#I think I should make individual graphs of these bacteria where the pairs are indicated by symbols. 
#We can take a previous for loop, select for only bacteria in these species_CEVE datasets and the change the symbol of the controls to match their end phenotype and output them in a separate folder by species. Then work on WGCNA.

# Read in data
setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria")
#Signficant AOV/TUKEY Bacteria
Oann_CEVE <- read.csv("Oann_CEVE.csv")
Cnat_CEVE <- read.csv("Cnat_CEVE.csv")
Ssid_CEVE <- read.csv("Ssid_CEVE.csv")
Past_CEVE <- read.csv("Past_CEVE.csv")
Ppor_CEVE <- read.csv("Ppor_CEVE.csv")

bacteria_results <- read.csv("EVE_Results.csv", header = T, fileEncoding = 'UTF-8-BOM')

# Load necessary libraries
library(ggplot2)
require(stringr) # for string manipulation
require(dplyr) # for dataframe manipulation
require(ggplot2) # for plotting

# For Loop to create multiple bar graphs (of the same format) from a single csv fil

#### Background info ####

# Desired boxplot figure output:
#   x-axis: coral species
#   y-axis: bacterium abundance
#   Colors/groups: based on species susceptibility
#   Individual figures for each: Bacteria species

# STEP : Do some initial data formatting and object creation to make for loop smoother
# Create object for constant plot theme info 
# This chuck of code does not change regardless of the plot
my_fill <- scale_fill_manual(breaks =c("O.faveolata_c","O.faveolata_i","O.annularis_c","O.annularis_e","O.annularis_i","C.natans_c","C.natans_e","C.natans_i", "S.siderea_c","S.siderea_e","S.siderea_i", "P.porites_c","P.porites_e","P.porites_i", "P.astreoides_c","P.astreoides_e","P.astreoides_i", "M.cavernosa_c","M.cavernosa_e"), 
                    values=c("mediumseagreen", "salmon","mediumseagreen","cornflowerblue","salmon", "mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","salmon"))
my_col <- scale_color_manual(breaks = c("O.faveolata_c", "O.faveolata_i","O.annularis_c","O.annularis_e","O.annularis_i","C.natans_c","C.natans_e","C.natans_i", "S.siderea_c","S.siderea_e","S.siderea_i", "P.porites_c","P.porites_e","P.porites_i", "P.astreoides_c","P.astreoides_e","P.astreoides_i", "M.cavernosa_c","M.cavernosa_e"), 
                     values=c("mediumseagreen", "salmon","mediumseagreen","cornflowerblue","salmon", "mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","cornflowerblue","salmon","mediumseagreen","salmon"))
my_theme <- theme(legend.position="NULL",legend.text=element_text(size=20), 
        legend.title = element_text(size=12), 
        axis.text.x = element_text(colour = "black", size=10,angle=30, hjust =1), 
        axis.text.y = element_text(colour = "black", size=10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        strip.text.x = element_text(colour = "black", size = 8),
        plot.subtitle=element_text(size=8))

# Pull out all metadata columns from dataframe
metadata <- colnames(bacteria_results[,colnames(bacteria_results) %in% c("Bacteria","Species","Treatment","RelativeRisk","Status","Genotype", "Type", "DysbiosisC2I", "DysbiosisC2E","ControlOutcome","OutcomeSymbol","Specified")])
# Pull out all bacterial abundance columns from dataframe
bacteria <- colnames(bacteria_results[,!colnames(bacteria_results) %in% metadata])
# Pull out the relevant columns from the significant eve bacteria results
eve <- sig_bacteria[,colnames(sig_bacteria) %in% c("Bacteria","LRT_LS","alpha_phylosymbiosis","beta_HEV","pval","Expression")]
# rename columns for ease
colnames(eve)[2:5] <- c("LRT","Alpha","Beta","Pattern")

#Combine Species_CEVE results.
#First specify what species comparison was made. 
CEVE <- rbind(Oann_CEVE, Cnat_CEVE)
CEVE <- rbind(CEVE, Ssid_CEVE)
CEVE <- rbind(CEVE, Past_CEVE)
CEVE <- rbind(CEVE, Ppor_CEVE)
write.csv(CEVE, file="CEVE.csv") #EVE bacteria (LS/PS) that have signficantly different abundances in controls whose pair either was resistant or susceptible to disease exposure.

CEVE <- read.csv("CEVE.csv")
CEVE_Bacteria <- CEVE$bacteria

# Pull out all bacterial abundance Rows from dataframe
bacteria <- colnames(bacteria_results[,colnames(bacteria_results) %in% CEVE_Bacteria])

# Some of the Eve bacteria names end in a space
# Remove the trailing space and now the names should match up
eve$Bacteria <-  trimws(eve$Bacteria, which = "right") #cut out space after species name.

# Make sure we can match the bacteria names in the eve and bacteria vectors
bacteria_names1 <- eve$Bacteria # create vector of bacteria names
bacteria_names2 <- str_replace_all(bacteria, "[.]", " ") # create vector of names but replace '.' with spaces


bacteria_names1[!bacteria_names1 %in% bacteria_names2] # find mismatching cases. Ideally there would be 0. 

# ****************************************************************************************************
# For loop testing: Un-comment (remove hashtag for line immediately below) i = 1, this will allow you to test the code only on the first loop

#i = 1

# IMPORTANT NOTE: While testing loops skip any lines that start with "for (....){"
# Only run lines that are INSIDE the loop brackets "{}"
# This will allow you to test and fix code running only the first loop (i.e., creating 1 figure) 
# and will prevent R from crashing if the code is buggy.
# the reason you test a for loop is to make sure it works before doing many iterations. This is why the test loop is i=1, i.e do one loop.
#****************************************************************************************************

# STEP 7: The For Loop
# The "I loop"
i=1
for (i in 1:length(bacteria)){ #iterates loop over each column of bacteria
  trip <- bacteria[i] # One bacteria per trip through the loop
  trip_name <- str_replace_all(trip,"[.]"," ") # replaces '.' with space in bacteria name
  results <- eve[eve$Bacteria %in% trip_name,] # pulls out the relevant EVE result line
  results2 <- CEVE[CEVE$bacteria %in% trip,]
  results3 <- cbind(results,results2)
  
  print("-------------------------------------------------") # helps separate iterations in consule
  # Print state to help you keep track of current iteration and any error that may arrise
  print(paste0("Now running: ", str_to_sentence(trip_name), ", #", i, " of ", length(bacteria), " bacteria"))
  
  y_val <- paste0(trip) # Create variable for current bacteria name
  data <- bacteria_results[,names(bacteria_results) %in% c(y_val, metadata)] # create small df for current iteration
  colnames(data)[13] <- "y_val" #rename bacteria abundance column
  
  # PLOT CODE
  plot <- ggplot(data = data, aes(x = reorder(Species, (y_val)), y = y_val, fill = Type)) +
    geom_point(aes(y = y_val, color = Type, shape =OutcomeSymbol), 
               position = position_dodge(width = 0.47),size = 1, alpha = 2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
    labs(y = "log10(abundance+1)", x = NULL, fill = "Status", col= "Status", 
         title = paste0(str_to_sentence(trip_name)), 
         subtitle = paste0("LRT = ", round(results3$LRT, digits = 1),
                           " | B = ", round(results3$Beta, digits = 2),
                           " | a = ", round(results3$Alpha, digits = 2),
                           " | p = ", formatC(results3$pval, digits = 2),
                           "\n", (results3$Species),
                           " TukeyP= ", round(results3$TUKEYPval, digits = 3))) +
    scale_shape_manual(values=c(17,15,17,15))+
        my_fill + my_col + my_theme

  # setwd to folder for all the plots to be placed into
  if (results3$Species == "O.annularis"){
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/ControlOutcome/")
    pdfname <- paste0("O.annularis_",str_to_sentence(trip_name), str_to_sentence(results3$Pattern),"_boxplot_individual.pdf") # name pdf
    pdf(pdfname, width = 2.7, height = 2.7) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results3$Species == "P.astreoides") {
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/ControlOutcome/")
    pdfname <- paste0("P.astreoides_",str_to_sentence(trip_name), str_to_sentence(results3$Pattern),"_boxplot_individual.pdf") # name pdf
    pdf(pdfname, width = 2.7, height = 2.7) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results3$Species == "P.porites") {
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/ControlOutcome/")
    pdfname <- paste0("P.porites_",str_to_sentence(trip_name), str_to_sentence(results3$Pattern),"_boxplot_individual.pdf") # name pdf
    pdf(pdfname, width = 2.7, height = 2.7) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results3$Species == "S.siderea") {
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/ControlOutcome/")
    pdfname <- paste0("S.siderea_",str_to_sentence(trip_name), str_to_sentence(results3$Pattern),"_boxplot_individual.pdf") # name pdf
    pdf(pdfname, width = 2.7, height = 2.7) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  
  if (results3$Species == "C.natans") {
    setwd("~/Desktop/White Plague /Mechanisms/Mechanisms/EVE Bacteria/ControlOutcome/")
    pdfname <- paste0("C.natans_",str_to_sentence(trip_name), str_to_sentence(results3$Pattern),"_boxplot_individual.pdf") # name pdf
    pdf(pdfname, width = 2.7, height = 2.7) # Create blank pdf with W/H in inches
    print(plot) # print plot to pdf
    dev.off() # Honestly dunno what this does but it's necessary
  }
  # DUNZO
}

# STEP 8: Check that pdfs are written to correct folder
# STEP 9: (Might required paid version of Adobe Acrobat? - If so, I can do this for you.)
# To combine into single pdf for easy scrolling:
#         In file explorer highlight all plots > right click > 
#         (IF available) "Combine files in Acrobat" > Combine > 
#         Rename "binder" > View > Page display > Two page scrolling
```
