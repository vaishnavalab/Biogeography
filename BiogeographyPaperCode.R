###Load Packages
library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(vegan)
library(RColorBrewer)
library(DESeq2)
library(reshape2)
library(dplyr)
library(rio)
library(data.table)
library(VennDiagram)
library(microbiome)
library(goeveg)
library(venneuler)

###Set working directory
setwd("~/Desktop/Biogeography_new/")
list.files()

##DADA2
#load(file="Biogeography_Paper_DADA2.RData")
#samdf <- read.table("Biogeography_Mappingfile_newest.txt")
#colnames(samdf)<- c("SampleID","MouseID","Location","Organ","Genotype","Time","Treatment","Gender","Description")
#mappingfile <- samdf[,-1]
#rownames(mappingfile) <- samdf[,1]
#physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
#                   sample_data(mappingfile), 
#                   tax_table(taxa.plus))
#library("ape")
#random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
#physeq1 = merge_phyloseq(physeq, mappingfile, random_tree)
#save.image(file="BiogeographyPaper_physeq1.RData")

###Load Data
load(file="BiogeographyPaper_physeq1.RData")
theme_set(theme_bw())

###Remove Singletons
physeq1_trim<-prune_taxa(taxa_sums(physeq1) > 1, physeq1)
physeq1_trim

###Subset Data
WTF <- physeq1_trim %>%
  subset_samples(Description == "WTF_Characterization")
WTF

OGWT <- physeq1_trim %>%
  subset_samples(Description == "OGWT_Characterization")
OGWT

Combined_Layers<-merge_phyloseq(WTF,OGWT)

Kellys_Kvon_ABX_Experiment <- physeq1_trim %>%
  subset_samples(Description == "Kellys_Kvon_ABX_Experiment")
Kellys_Kvon_ABX_Experiment

Kellys_1st50 <- Kellys_Kvon_ABX_Experiment %>%
  subset_samples(Location == "1st50")
Kellys_1st50

Kellys_1st50_day0 <- Kellys_1st50 %>%
  subset_samples(Time == "day0")
Kellys_1st50_day0

Kellys_1st50_day10 <- Kellys_1st50 %>%
  subset_samples(Time == "day10")
Kellys_1st50_day10

Kellys_Luminal <- Kellys_Kvon_ABX_Experiment %>%
  subset_samples(Location == "Luminal")
Kellys_Luminal

Kellys_Luminal_day0 <- Kellys_Luminal %>%
  subset_samples(Time == "day0")
Kellys_Luminal_day0

Kellys_Luminal_day10 <- Kellys_Luminal %>%
  subset_samples(Time == "day10")
Kellys_Luminal_day10

Kellys_Layers<-merge_phyloseq(Kellys_1st50,Kellys_Luminal)
Kellys_Day10<-merge_phyloseq(Kellys_1st50_day10,Kellys_Luminal_day10)

Taconic_Untreated <- physeq1_trim %>%
  subset_samples(Description == "Taconic_Untreated")
Taconic_Untreated

Taconic_Vancomycin <- physeq1_trim %>%
  subset_samples(Description == "Taconic_Vancomycin")
Taconic_Vancomycin

Taconic_Ciprofloxacin <- physeq1_trim %>%
  subset_samples(Description == "Taconic_Ciprofloxacin")
Taconic_Ciprofloxacin

Taconic_ABX<-merge_phyloseq(Taconic_Untreated,Taconic_Vancomycin,Taconic_Ciprofloxacin)
Taconic_Untreated_Vancomycin<-merge_phyloseq(Taconic_Untreated,Taconic_Vancomycin)
Taconic_Untreated_Ciprofloxacin<-merge_phyloseq(Taconic_Untreated,Taconic_Ciprofloxacin)

#Figure 1

##Rarefied Diversity Table 

###Reads per sample
Analysis<-Taconic_ABX
sdt = data.table(as(sample_data(Analysis), "data.frame"),
                 TotalReads = sample_sums(Analysis), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt

###Remove Samples with low sequencing depth
Taconic_ABX = subset_samples(Taconic_ABX, MouseID != "TacC2")
Taconic_ABX

###Remove Outliers
Taconic_ABX_mod = subset_samples(Taconic_ABX, MouseID != "TacV1")
Taconic_ABX_mod = subset_samples(Taconic_ABX_mod, MouseID != "TacV3")
Taconic_ABX_mod = subset_samples(Taconic_ABX_mod, MouseID != "TacC3")

###Rarefy Table
set.seed(42)
rarefied_analysis = rarefy_even_depth(Taconic_ABX_mod)
rarefied_analysis

###Plot Diversity
p = plot_richness(rarefied_analysis, x="Description", color="MouseID", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 

rarefied_values<-p$data
rarefied_values
write.csv(rarefied_values, file = "TaconicABX_rarefied_values.csv")

##Bar plots

###Change ASV IDs
Analysis<-Taconic_ABX
taxa_names(Analysis) <- paste0("Seq", seq(ntaxa(Analysis)))

###Normalize dataset
subset=transform_sample_counts(Analysis,function(x)x/sum(x))

###Condense by taxonomic level
ps.class<-tax_glom(subset,taxrank="Class")

###Percent Abundance
otu.table<-otu_table(ps.class)
write.csv(otu.table,"TaconicABX.class.otu.table.csv")
tax.table<-tax_table(ps.class)
write.csv(tax.table,"TaconicABX.class.tax.table.csv")

##Differential Abundance between Untreated and Vancomycin

###Change ASV IDs
Analysis<-Taconic_Untreated_Vancomycin
taxa_names(Analysis) <- paste0("Seq", seq(ntaxa(Analysis)))

###Condense by taxonomic level
ps.class<-tax_glom(Analysis,taxrank="Class")

###Differential Abundance
dds = phyloseq_to_deseq2(ps.class, ~ Description)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.class)[rownames(sigtab_WD), ], "matrix"))
diff.abund
write.csv(diff.abund,"Untreated.vs.Vancomycin.class.csv")

##Differential Abundance between Untreated and Ciprofloxacin

###Change ASV IDs

Taconic_Untreated_Ciprofloxacin = subset_samples(Taconic_Untreated_Ciprofloxacin, MouseID != "TacC2")
Analysis<-Taconic_Untreated_Ciprofloxacin
taxa_names(Analysis) <- paste0("Seq", seq(ntaxa(Analysis)))

###Condense by taxonomic level
ps.class<-tax_glom(Analysis,taxrank="Class")

###Differential Abundance
dds = phyloseq_to_deseq2(ps.class, ~ Description)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.class)[rownames(sigtab_WD), ], "matrix"))
diff.abund
write.csv(diff.abund,"Untreated.vs.Ciprofloxacin.class.csv")

#Figure 2

##Combined Weighted Unifrac PCoA
Analysis<-Combined_Layers
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduW = ordinate(GP1, "PCoA","unifrac",weighted=TRUE)
pW = plot_ordination(GP1, orduW, color = "Location", shape = "Description", title = "Weighted Unifrac") + geom_point(size=6) + scale_colour_manual(values=c("red","blue"))
pW + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) #+ylim(-0.1,0.1) +xlim(-0.2,0.0)

###Permanova by Location
set.seed(42)
analysis_unifrac_weighted<-phyloseq::distance(Analysis,method = "unifrac", weighted=TRUE)
sampledf<- data.frame(sample_data(Analysis))
adonis(analysis_unifrac_weighted ~ Location, data = sampledf)

###Permanova by Description
adonis(analysis_unifrac_weighted ~ Description, data = sampledf)

##Unweighted Unifrac PCoA
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduU = ordinate(GP1, "PCoA","unifrac",weighted=FALSE)
pU = plot_ordination(GP1, orduU, color = "Location", shape = "Description", title = "Unweighted Unifrac") + geom_point(size=6) + scale_colour_manual(values=c("Red","Blue"))
pU + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  #+ ylim(-0.2,0.1) +xlim(-0.30,0.30) 
pU1 = plot_ordination(GP1, orduU, color = "Location", shape = "Description", title = "Unweighted Unifrac") + geom_point(size=5) 
pU1 + stat_ellipse(type = "t", level = 0.85) + theme_bw()  + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  + scale_colour_manual(values=c("red","blue"))

###Permanova by Location
set.seed(42)
analysis_unifrac_unweighted<-phyloseq::distance(Analysis,method = "unifrac", weighted=FALSE)
sampledf<- data.frame(sample_data(Analysis))
adonis(analysis_unifrac_unweighted ~ Location, data = sampledf)

###Permanova by Description
adonis(analysis_unifrac_unweighted ~ Description, data = sampledf)

###Enter Dataset for bar plots
Analysis<-WTF

###Set barplot Colors
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
ClassList = unique(tax_table(Analysis)[,"Class"])
ClassPalette = getPalette(length(ClassList))
names(ClassPalette) = ClassList

###Normalize dataset
subset=transform_sample_counts(Analysis,function(x)x/sum(x))

###Bar plots phylum level
p = plot_bar(subset, x= "Location", fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack") + theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

###Condense by taxonomic level
ps.class<-tax_glom(Analysis,taxrank="Class")

###Percent Abundance
otu.table<-otu_table(ps.class)
write.csv(otu.table,"WTF.class.otu.table.csv")
tax.table<-tax_table(ps.class)
write.csv(tax.table,"WTF.class.tax.table.csv")

###Change ASV IDs
taxa_names(Analysis) <- paste0("Seq", seq(ntaxa(Analysis)))

###Condense by taxonomic level
ps.phylum<-tax_glom(Analysis,taxrank="Phylum")
ps.class<-tax_glom(Analysis,taxrank="Class")

##Differential Abundance based on Location: Phylum
dds = phyloseq_to_deseq2(ps.phylum, ~ Location)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.phylum)[rownames(sigtab_WD), ], "matrix"))
diff.abund

###Plot counts for OTU determined differentially abundant and group by MouseID
data <- plotCounts(dds, "Seq4922", intgroup=c("Location","MouseID"), returnData=TRUE)
plot<-ggplot(data, aes(x=Location, y=count, color=MouseID, group=MouseID)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Phylum Proteobacteria"))

##Differential Abundance based on Location: Class
dds = phyloseq_to_deseq2(ps.class, ~ Location)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.class)[rownames(sigtab_WD), ], "matrix"))
diff.abund

###Plot counts for OTU determined differentially abundant and group by MouseID
data <- plotCounts(dds, "Seq5805", intgroup=c("Location","MouseID"), returnData=TRUE)
plot<-ggplot(data, aes(x=Location, y=count, color=MouseID, group=MouseID)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Class Gammaproteobacteria"))

###Make Excel sheet with normalized counts
dds <- estimateSizeFactors(dds)
dds.counts<-counts(dds, normalized=TRUE)
write.csv(dds.counts,"WTF.norm.counts.class.csv")

#Figure 3

##Rarefaction Curves
set.seed(42)
psdata<-WTF

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(psdata, c('Observed'), rep(c(100,1000,2000,4000,6000,8000,10000,12000,14000,16000,18000,20000,22000,24000,26000,28000,30000,32000), each = 50))

summary(rarefaction_curve_data)
###Summarize Alpha Diversity
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
###Add Sample Data
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')
###Plot
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Location,
    group = Sample
  )
) + geom_line(
) + geom_pointrange(
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
) + geom_point(size=1) + scale_colour_manual(values=c("red","blue")) + ggtitle("Alpha Rarefaction") + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 

##Rarefied Diversity Table
Analysis<-WTF_Layers
set.seed(42)
rarefied_analysis = rarefy_even_depth(Analysis)
rarefied_analysis<-Analysis
p = plot_richness(rarefied_analysis, x="Location", color="MouseID", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 
rarefied_values<-p$data
rarefied_values
write.csv(rarefied_values, file = "rarefied_diversity_values_WTF.csv")

###Chao1
p = plot_richness(Analysis, x="Location", color="MouseID", measures=c("Chao1"))
p + geom_point(size=5, alpha=0.7) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 
chao1_values<-p$data
chao1_values
write.csv(chao1_values, file = "chao1_values_WTF.csv")

###Evenness
even<-evenness(rarefied_analysis, index = "all", zeroes = TRUE)
even
write.csv(even, file = "even_WTF.csv")

#ASV Distribution
Analysis<-WTF

###Mouse 1
Mouse1 = subset_samples(Analysis,MouseID == "WTF1")
Mouse1<-prune_taxa(taxa_sums(Mouse1) > 0, Mouse1)
Mouse1 #1021

###Reads per sample
sdt = data.table(as(sample_data(Mouse1), "data.frame"),
                 TotalReads = sample_sums(Mouse1), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt

###Regional list of ASVs
Inner <- Mouse1 %>%
  subset_samples(Location == "1st100")
Inner
Mouse1.Inner.trim<-prune_taxa(taxa_sums(Inner) > 0, Inner)
Mouse1.Inner.trim
Inner.all.otu<-otu_table(Mouse1.Inner.trim)
Inner.all.melt<-psmelt(Inner.all.otu)
head(Inner.all.melt)
Inner.all.wide <- reshape(Inner.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Inner.all.wide)
Inner.all.wide[1:10, 1:2]
Mouse1.Inner.list<-Inner.all.wide[,1]
Mouse1.Inner.list

###Number of ASVs found in inner in each class
ps.class<-Mouse1.Inner.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #24
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #255
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #8
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #607
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #9
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #6

Outer <- Mouse1 %>%
  subset_samples(Location == "3rd50")
Outer
Mouse1.Outer.trim<-prune_taxa(taxa_sums(Outer) > 0, Outer)
Mouse1.Outer.trim
Outer.all.otu<-otu_table(Mouse1.Outer.trim)
Outer.all.melt<-psmelt(Outer.all.otu)
head(Outer.all.melt)
Outer.all.wide <- reshape(Outer.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Outer.all.wide)
Outer.all.wide[1:10, 1:2]
Mouse1.Outer.list<-Outer.all.wide[,1]
Mouse1.Outer.list

###Number of ASVs found in outer in each class
ps.class<-Mouse1.Outer.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #6
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #152
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #8
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #261
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #9
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #7

###ASVs unique to inner
allTaxa = taxa_names(Mouse1.Inner.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse1.Outer.list)]
Mouse1.Inner.unique = prune_taxa(myTaxa, Mouse1.Inner.trim)
Mouse1.Inner.unique #561

###Abundances of families unique to inner
Mouse1.Inner.unique.family<-tax_glom(Mouse1.Inner.unique,taxrank="Family")
taxa_names(Mouse1.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse1.Inner.unique.family)))
otu.table<-otu_table(Mouse1.Inner.unique.family)
write.csv(otu.table,"Mouse1.Inner.unique.otu.table.csv")
tax.table<-tax_table(Mouse1.Inner.unique.family)
write.csv(tax.table,"Mouse1.Inner.unique.tax.table.csv")

###Number of ASVs unique to inner in each class
ps.class<-tax_glom(Mouse1.Inner.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse1.Inner.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #23
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #138
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #1
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #386
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #1
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #3

###ASVs unique to outer
allTaxa = taxa_names(Mouse1.Outer.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse1.Inner.list)]
Mouse1.Outer.unique = prune_taxa(myTaxa, Mouse1.Outer.trim)
Mouse1.Outer.unique #89

###Abundances of families unique to outer
Mouse1.Outer.unique.family<-tax_glom(Mouse1.Outer.unique,taxrank="Family")
taxa_names(Mouse1.Outer.unique.family) <- paste0("Seq", seq(ntaxa(Mouse1.Outer.unique.family)))
otu.table<-otu_table(Mouse1.Outer.unique.family)
write.csv(otu.table,"Mouse1.Outer.unique.otu.table.csv")
tax.table<-tax_table(Mouse1.Outer.unique.family)
write.csv(tax.table,"Mouse1.Outer.unique.tax.table.csv")

###Number of ASVs unique to outer in each class
ps.class<-tax_glom(Mouse1.Outer.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse1.Outer.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #5
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #35
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #1
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #40
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #1
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #4

###ASVs Shared
allTaxa = taxa_names(Mouse1.Outer.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse1.Inner.list)]
Mouse1.Shared.unique1 = prune_taxa(myTaxa, Mouse1.Outer.trim)
Mouse1.Shared.unique1
allTaxa = taxa_names(Mouse1.Inner.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse1.Outer.list)]
Mouse1.Shared.unique2 = prune_taxa(myTaxa, Mouse1.Inner.trim)
Mouse1.Shared.unique2
Mouse1.Shared.unique<-merge_phyloseq(Mouse1.Shared.unique1,Mouse1.Shared.unique2)
Mouse1.Shared.unique #371

###Abundances of families shared
Mouse1.Shared.unique.family<-tax_glom(Mouse1.Shared.unique,taxrank="Family")
taxa_names(Mouse1.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse1.Shared.unique.family)))
otu.table<-otu_table(Mouse1.Shared.unique.family)
write.csv(otu.table,"Mouse1.Shared.unique.otu.table.csv")
tax.table<-tax_table(Mouse1.Shared.unique.family)
write.csv(tax.table,"Mouse1.Shared.unique.tax.table.csv")

###Number of ASVs shared in each class
ps.class<-tax_glom(Mouse1.Shared.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse1.Shared.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #1
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #117
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #7
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #221
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #8
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #3

###Plot unique vs. shared ASVs Mouse 1
venn=venneuler(c('A' = 561, 'B' = 89, 'A&B' = 371))
venn$labels = c("", "")
plot(venn, col=c("red", "blue")) 

###Mouse 2
Mouse2 = subset_samples(Analysis,MouseID == "WTF2")
Mouse2<-prune_taxa(taxa_sums(Mouse2) > 0, Mouse2)
Mouse2 #882

###Reads per sample
sdt = data.table(as(sample_data(Mouse2), "data.frame"),
                 TotalReads = sample_sums(Mouse2), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt

###Regional list of ASVs
Inner <- Mouse2 %>%
  subset_samples(Location == "1st100")
Inner
Mouse2.Inner.trim<-prune_taxa(taxa_sums(Inner) > 0, Inner)
Mouse2.Inner.trim
Inner.all.otu<-otu_table(Mouse2.Inner.trim)
Inner.all.melt<-psmelt(Inner.all.otu)
head(Inner.all.melt)
Inner.all.wide <- reshape(Inner.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Inner.all.wide)
Inner.all.wide[1:10, 1:2]
Mouse2.Inner.list<-Inner.all.wide[,1]
Mouse2.Inner.list

###Number of ASVs found in inner in each class
ps.class<-Mouse2.Inner.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #39
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #269
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #13
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #426
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #10
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #11

Outer <- Mouse2 %>%
  subset_samples(Location == "3rd50")
Outer
Mouse2.Outer.trim<-prune_taxa(taxa_sums(Outer) > 0, Outer)
Mouse2.Outer.trim
Outer.all.otu<-otu_table(Mouse2.Outer.trim)
Outer.all.melt<-psmelt(Outer.all.otu)
head(Outer.all.melt)
Outer.all.wide <- reshape(Outer.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Outer.all.wide)
Outer.all.wide[1:10, 1:2]
Mouse2.Outer.list<-Outer.all.wide[,1]
Mouse2.Outer.list

###Number of ASVs found in outer in each class
ps.class<-Mouse2.Outer.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #8
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #96
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #9
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #76
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #1
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #12

###ASVs unique to inner
allTaxa = taxa_names(Mouse2.Inner.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse2.Outer.list)]
Mouse2.Inner.unique = prune_taxa(myTaxa, Mouse2.Inner.trim)
Mouse2.Inner.unique #659

###Abundances of families unique to inner
Mouse2.Inner.unique.family<-tax_glom(Mouse2.Inner.unique,taxrank="Family")
taxa_names(Mouse2.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse2.Inner.unique.family)))
otu.table<-otu_table(Mouse2.Inner.unique.family)
write.csv(otu.table,"Mouse2.Inner.unique.otu.table.csv")
tax.table<-tax_table(Mouse2.Inner.unique.family)
write.csv(tax.table,"Mouse2.Inner.unique.tax.table.csv")

###Number of ASVs unique to inner in each class
ps.class<-tax_glom(Mouse2.Inner.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse2.Inner.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #36
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #195
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #6
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #366
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #9
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #5

###ASVs unique to outer
allTaxa = taxa_names(Mouse2.Outer.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse2.Inner.list)]
Mouse2.Outer.unique = prune_taxa(myTaxa, Mouse2.Outer.trim)
Mouse2.Outer.unique #62

###Abundances of families unique to outer
Mouse2.Outer.unique.family<-tax_glom(Mouse2.Outer.unique,taxrank="Family")
taxa_names(Mouse2.Outer.unique.family) <- paste0("Seq", seq(ntaxa(Mouse2.Outer.unique.family)))
otu.table<-otu_table(Mouse2.Outer.unique.family)
write.csv(otu.table,"Mouse2.Outer.unique.otu.table.csv")
tax.table<-tax_table(Mouse2.Outer.unique.family)
write.csv(tax.table,"Mouse2.Outer.unique.tax.table.csv")

###Number of ASVs unique to outer in each class
ps.class<-tax_glom(Mouse2.Outer.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse2.Outer.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #5
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #22
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #2
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #16
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #9
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #6

###ASVs shared
allTaxa = taxa_names(Mouse2.Outer.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse2.Inner.list)]
Mouse2.Shared.unique1 = prune_taxa(myTaxa, Mouse2.Outer.trim)
Mouse2.Shared.unique1
allTaxa = taxa_names(Mouse2.Inner.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse2.Outer.list)]
Mouse2.Shared.unique2 = prune_taxa(myTaxa, Mouse2.Inner.trim)
Mouse2.Shared.unique2
Mouse2.Shared.unique<-merge_phyloseq(Mouse2.Shared.unique1,Mouse2.Shared.unique2)
Mouse2.Shared.unique #161

###Abundances of families shared
Mouse2.Shared.unique.family<-tax_glom(Mouse2.Shared.unique,taxrank="Family")
taxa_names(Mouse2.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse2.Shared.unique.family)))
otu.table<-otu_table(Mouse2.Shared.unique.family)
write.csv(otu.table,"Mouse2.Shared.unique.otu.table.csv")
tax.table<-tax_table(Mouse2.Shared.unique.family)
write.csv(tax.table,"Mouse2.Shared.unique.tax.table.csv")

###Number of ASVs shared in each class
ps.class<-tax_glom(Mouse2.Shared.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse2.Shared.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #3
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #74
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #7
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #60
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #1
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #6

###Plot unique vs. shared ASVs Mouse 2
venn=venneuler(c('A' = 659, 'B' = 62, 'A&B' = 161))
venn$labels = c("", "")
plot(venn, col=c("red", "blue")) 

###Mouse 3
Mouse3 = subset_samples(Analysis,MouseID == "WTF3")
Mouse3<-prune_taxa(taxa_sums(Mouse3) > 0, Mouse3)
Mouse3 #1042

###Reads per sample
sdt = data.table(as(sample_data(Mouse3), "data.frame"),
                 TotalReads = sample_sums(Mouse3), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt

###Regional list of ASVs
Inner <- Mouse3 %>%
  subset_samples(Location == "1st100")
Inner
Mouse3.Inner.trim<-prune_taxa(taxa_sums(Inner) > 0, Inner)
Mouse3.Inner.trim
Inner.all.otu<-otu_table(Mouse3.Inner.trim)
Inner.all.melt<-psmelt(Inner.all.otu)
head(Inner.all.melt)
Inner.all.wide <- reshape(Inner.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Inner.all.wide)
Inner.all.wide[1:10, 1:2]
Mouse3.Inner.list<-Inner.all.wide[,1]
Mouse3.Inner.list

###Number of ASVs found in inner in each class
ps.class<-Mouse3.Inner.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #48
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #362
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #21
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #375
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #62
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #8

Outer <- Mouse3 %>%
  subset_samples(Location == "3rd50")
Outer
Mouse3.Outer.trim<-prune_taxa(taxa_sums(Outer) > 0, Outer)
Mouse3.Outer.trim
Outer.all.otu<-otu_table(Mouse3.Outer.trim)
Outer.all.melt<-psmelt(Outer.all.otu)
head(Outer.all.melt)
Outer.all.wide <- reshape(Outer.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Outer.all.wide)
Outer.all.wide[1:10, 1:2]
Mouse3.Outer.list<-Outer.all.wide[,1]
Mouse3.Outer.list

###Number of ASVs found in outer in each class
ps.class<-Mouse3.Outer.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #13
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #142
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #14
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #92
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #12
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #12

###ASVs unique to inner
allTaxa = taxa_names(Mouse3.Inner.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse3.Outer.list)]
Mouse3.Inner.unique = prune_taxa(myTaxa, Mouse3.Inner.trim)
Mouse3.Inner.unique #697

###Abundances of families unique to inner
Mouse3.Inner.unique.family<-tax_glom(Mouse3.Inner.unique,taxrank="Family")
taxa_names(Mouse3.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse3.Inner.unique.family)))
otu.table<-otu_table(Mouse3.Inner.unique.family)
write.csv(otu.table,"Mouse3.Inner.unique.otu.table.csv")
tax.table<-tax_table(Mouse3.Inner.unique.family)
write.csv(tax.table,"Mouse3.Inner.unique.tax.table.csv")

###Number of ASVs unique to inner in each class
ps.class<-tax_glom(Mouse3.Inner.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse3.Inner.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #39
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #251
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #12
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #301
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #51
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #4

###ASVs unique to outer
allTaxa = taxa_names(Mouse3.Outer.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse3.Inner.list)]
Mouse3.Outer.unique = prune_taxa(myTaxa, Mouse3.Outer.trim)
Mouse3.Outer.unique #111

###Abundances of families unique to outer
Mouse3.Outer.unique.family<-tax_glom(Mouse3.Outer.unique,taxrank="Family")
taxa_names(Mouse3.Outer.unique.family) <- paste0("Seq", seq(ntaxa(Mouse3.Outer.unique.family)))
otu.table<-otu_table(Mouse3.Outer.unique.family)
write.csv(otu.table,"Mouse3.Outer.unique.otu.table.csv")
tax.table<-tax_table(Mouse3.Outer.unique.family)
write.csv(tax.table,"Mouse3.Outer.unique.tax.table.csv")

###Number of ASVs unique to outer in each class
ps.class<-tax_glom(Mouse3.Outer.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse3.Outer.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #4
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #31
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #5
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #18
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #1
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #8

###ASVs shared
allTaxa = taxa_names(Mouse3.Outer.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse3.Inner.list)]
Mouse3.Shared.unique1 = prune_taxa(myTaxa, Mouse3.Outer.trim)
Mouse3.Shared.unique1
allTaxa = taxa_names(Mouse3.Inner.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse3.Outer.list)]
Mouse3.Shared.unique2 = prune_taxa(myTaxa, Mouse3.Inner.trim)
Mouse3.Shared.unique2
Mouse3.Shared.unique<-merge_phyloseq(Mouse3.Shared.unique1,Mouse3.Shared.unique2)
Mouse3.Shared.unique #234

###Abundances of families shared
Mouse3.Shared.unique.family<-tax_glom(Mouse3.Shared.unique,taxrank="Family")
taxa_names(Mouse3.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse3.Shared.unique.family)))
otu.table<-otu_table(Mouse3.Shared.unique.family)
write.csv(otu.table,"Mouse3.Shared.unique.otu.table.csv")
tax.table<-tax_table(Mouse3.Shared.unique.family)
write.csv(tax.table,"Mouse3.Shared.unique.tax.table.csv")

###Number of ASVs shared in each class
ps.class<-tax_glom(Mouse3.Shared.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse3.Shared.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #9
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #111
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #9
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #74
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #11
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #4

###Plot unique vs. shared ASVs Mouse 2
venn=venneuler(c('A' = 697, 'B' = 111, 'A&B' = 234))
venn$labels = c("", "")
plot(venn, col=c("red", "blue")) 

###Mouse 4
Mouse4 = subset_samples(Analysis,MouseID == "WTF4")
Mouse4<-prune_taxa(taxa_sums(Mouse4) > 0, Mouse4)
Mouse4 #1226

###Reads per sample
sdt = data.table(as(sample_data(Mouse4), "data.frame"),
                 TotalReads = sample_sums(Mouse4), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt

###Regional list of ASVs
Inner <- Mouse4 %>%
  subset_samples(Location == "1st100")
Inner
Mouse4.Inner.trim<-prune_taxa(taxa_sums(Inner) > 0, Inner)
Mouse4.Inner.trim
Inner.all.otu<-otu_table(Mouse4.Inner.trim)
Inner.all.melt<-psmelt(Inner.all.otu)
head(Inner.all.melt)
Inner.all.wide <- reshape(Inner.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Inner.all.wide)
Inner.all.wide[1:10, 1:2]
Mouse4.Inner.list<-Inner.all.wide[,1]
Mouse4.Inner.list

###Number of ASVs found in inner in each class
ps.class<-Mouse4.Inner.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #46
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #403
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #18
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #400
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #33
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #11

Outer <- Mouse4 %>%
  subset_samples(Location == "3rd50")
Outer
Mouse4.Outer.trim<-prune_taxa(taxa_sums(Outer) > 0, Outer)
Mouse4.Outer.trim
Outer.all.otu<-otu_table(Mouse4.Outer.trim)
Outer.all.melt<-psmelt(Outer.all.otu)
head(Outer.all.melt)
Outer.all.wide <- reshape(Outer.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Outer.all.wide)
Outer.all.wide[1:10, 1:2]
Mouse4.Outer.list<-Outer.all.wide[,1]
Mouse4.Outer.list

###Number of ASVs found in outer in each class
ps.class<-Mouse4.Outer.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #39
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #248
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #25
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #243
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #22
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #16

###ASVs unique to inner
allTaxa = taxa_names(Mouse4.Inner.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse4.Outer.list)]
Mouse4.Inner.unique = prune_taxa(myTaxa, Mouse4.Inner.trim)
Mouse4.Inner.unique #590

###Abundances of families unique to outer
Mouse4.Inner.unique.family<-tax_glom(Mouse4.Inner.unique,taxrank="Family")
taxa_names(Mouse4.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse4.Inner.unique.family)))
otu.table<-otu_table(Mouse4.Inner.unique.family)
write.csv(otu.table,"Mouse4.Inner.unique.otu.table.csv")
tax.table<-tax_table(Mouse4.Inner.unique.family)
write.csv(tax.table,"Mouse4.Inner.unique.tax.table.csv")

###Number of ASVs unique to inner in each class
ps.class<-tax_glom(Mouse4.Inner.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse4.Inner.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #23
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #246
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #4
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #256
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #19
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #5

###ASVs unique to outer
allTaxa = taxa_names(Mouse4.Outer.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse4.Inner.list)]
Mouse4.Outer.unique = prune_taxa(myTaxa, Mouse4.Outer.trim)
Mouse4.Outer.unique #256

###Abundances of families unique to outer
Mouse4.Outer.unique.family<-tax_glom(Mouse4.Outer.unique,taxrank="Family")
taxa_names(Mouse4.Outer.unique.family) <- paste0("Seq", seq(ntaxa(Mouse4.Outer.unique.family)))
otu.table<-otu_table(Mouse4.Outer.unique.family)
write.csv(otu.table,"Mouse4.Outer.unique.otu.table.csv")
tax.table<-tax_table(Mouse4.Outer.unique.family)
write.csv(tax.table,"Mouse4.Outer.unique.tax.table.csv")

###Number of ASVs unique to outer in each class
ps.class<-tax_glom(Mouse4.Outer.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse4.Outer.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #16
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #91
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #11
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #99
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #8
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #10

###ASVs Shared
allTaxa = taxa_names(Mouse4.Outer.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse4.Inner.list)]
Mouse4.Shared.unique1 = prune_taxa(myTaxa, Mouse4.Outer.trim)
Mouse4.Shared.unique1
allTaxa = taxa_names(Mouse4.Inner.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse4.Outer.list)]
Mouse4.Shared.unique2 = prune_taxa(myTaxa, Mouse4.Inner.trim)
Mouse4.Shared.unique2
Mouse4.Shared.unique<-merge_phyloseq(Mouse4.Shared.unique1,Mouse4.Shared.unique2)
Mouse4.Shared.unique #380

###Abundances of families shared
Mouse4.Shared.unique.family<-tax_glom(Mouse4.Shared.unique,taxrank="Family")
taxa_names(Mouse4.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse4.Shared.unique.family)))
otu.table<-otu_table(Mouse4.Shared.unique.family)
write.csv(otu.table,"Mouse4.Shared.unique.otu.table.csv")
tax.table<-tax_table(Mouse4.Shared.unique.family)
write.csv(tax.table,"Mouse4.Shared.unique.tax.table.csv")

###Number of ASVs Shared in each class
ps.class<-tax_glom(Mouse4.Shared.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse4.Shared.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #23
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #157
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #14
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #144
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #14
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #6

###Plot unique vs. shared ASVs Mouse 2
venn=venneuler(c('A' = 590, 'B' = 256, 'A&B' = 380))
venn$labels = c("", "")
plot(venn, col=c("red", "blue")) 

###Mouse 5
Mouse5 =  subset_samples(Analysis,MouseID == "WTF6")
Mouse5<-prune_taxa(taxa_sums(Mouse5) > 0, Mouse5)
Mouse5 #899

###Reads per sample
sdt = data.table(as(sample_data(Mouse5), "data.frame"),
                 TotalReads = sample_sums(Mouse5), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt

###Regional list of ASVs
Inner <- Mouse5 %>%
  subset_samples(Location == "1st100")
Inner
Mouse5.Inner.trim<-prune_taxa(taxa_sums(Inner) > 0, Inner)
Mouse5.Inner.trim
Inner.all.otu<-otu_table(Mouse5.Inner.trim)
Inner.all.melt<-psmelt(Inner.all.otu)
head(Inner.all.melt)
Inner.all.wide <- reshape(Inner.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Inner.all.wide)
Inner.all.wide[1:10, 1:2]
Mouse5.Inner.list<-Inner.all.wide[,1]
Mouse5.Inner.list

###Number of ASVs found in inner in each class
ps.class<-Mouse5.Inner.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #44
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #362
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #22
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #253
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #57
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #11

Outer <- Mouse5 %>%
  subset_samples(Location == "3rd50")
Outer
Mouse5.Outer.trim<-prune_taxa(taxa_sums(Outer) > 0, Outer)
Mouse5.Outer.trim
Outer.all.otu<-otu_table(Mouse5.Outer.trim)
Outer.all.melt<-psmelt(Outer.all.otu)
head(Outer.all.melt)
Outer.all.wide <- reshape(Outer.all.melt, idvar = "OTU", timevar = "Sample", direction = "wide")
dim(Outer.all.wide)
Outer.all.wide[1:10, 1:2]
Mouse5.Outer.list<-Outer.all.wide[,1]
Mouse5.Outer.list

###Number of ASVs found in inner in each class
ps.class<-Mouse5.Outer.trim
Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #9
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #144
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #11
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #80
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #12
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #16

###ASVs unique to inner
allTaxa = taxa_names(Mouse5.Inner.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse5.Outer.list)]
Mouse5.Inner.unique = prune_taxa(myTaxa, Mouse5.Inner.trim)
Mouse5.Inner.unique #592

###Abundances of families unique to inner
Mouse5.Inner.unique.family<-tax_glom(Mouse5.Inner.unique,taxrank="Family")
taxa_names(Mouse5.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse5.Inner.unique.family)))
otu.table<-otu_table(Mouse5.Inner.unique.family)
write.csv(otu.table,"Mouse5.Inner.unique.otu.table.csv")
tax.table<-tax_table(Mouse5.Inner.unique.family)
write.csv(tax.table,"Mouse5.Inner.unique.tax.table.csv")

###Number of ASVs unique to inner in each class
ps.class<-tax_glom(Mouse5.Inner.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse5.Inner.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #38
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #247
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #15
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #203
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #46
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #6
 
###ASVs unique to outer
allTaxa = taxa_names(Mouse5.Outer.trim)
myTaxa <- allTaxa[!(allTaxa %in% Mouse5.Inner.list)]
Mouse5.Outer.unique = prune_taxa(myTaxa, Mouse5.Outer.trim)
Mouse5.Outer.unique #93

###Abundances of families unique to outer
Mouse5.Outer.unique.family<-tax_glom(Mouse5.Outer.unique,taxrank="Family")
taxa_names(Mouse5.Outer.unique.family) <- paste0("Seq", seq(ntaxa(Mouse5.Outer.unique.family)))
otu.table<-otu_table(Mouse5.Outer.unique.family)
write.csv(otu.table,"Mouse5.Outer.unique.otu.table.csv")
tax.table<-tax_table(Mouse5.Outer.unique.family)
write.csv(tax.table,"Mouse5.Outer.unique.tax.table.csv")

###Number of ASVs unique to outer in each class
ps.class<-tax_glom(Mouse5.Outer.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse5.Outer.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #3
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #29
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria  #4
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #30
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #1
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #11

###ASVs Shared
allTaxa = taxa_names(Mouse5.Outer.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse5.Inner.list)]
Mouse5.Shared.unique1 = prune_taxa(myTaxa, Mouse5.Outer.trim)
Mouse5.Shared.unique1
allTaxa = taxa_names(Mouse5.Inner.trim)
myTaxa <- allTaxa[(allTaxa %in% Mouse5.Outer.list)]
Mouse5.Shared.unique2 = prune_taxa(myTaxa, Mouse5.Inner.trim)
Mouse5.Shared.unique2
Mouse5.Shared.unique<-merge_phyloseq(Mouse5.Shared.unique1,Mouse5.Shared.unique2)
Mouse5.Shared.unique #214

###Abundances of families shared
Mouse5.Shared.unique.family<-tax_glom(Mouse5.Shared.unique,taxrank="Family")
taxa_names(Mouse5.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse5.Shared.unique.family)))
otu.table<-otu_table(Mouse5.Shared.unique.family)
write.csv(otu.table,"Mouse5.Shared.unique.otu.table.csv")
tax.table<-tax_table(Mouse5.Shared.unique.family)
write.csv(tax.table,"Mouse5.Shared.unique.tax.table.csv")

###Number of ASVs Shared in each class
ps.class<-tax_glom(Mouse5.Shared.unique,taxrank="Class")
p = plot_bar(ps.class,fill = "Class") + scale_color_manual(values= ClassPalette) + scale_fill_manual(values= ClassPalette)
p + geom_bar(aes(color= Class, fill=Class), stat="identity", position="stack")
ps.class<-Mouse5.Shared.unique

Bacilli = subset_taxa(ps.class, Class=="Bacilli")
Bacilli #6
Bacteroidia = subset_taxa(ps.class, Class=="Bacteroidia")
Bacteroidia #115
Betaproteobacteria = subset_taxa(ps.class, Class=="Betaproteobacteria")
Betaproteobacteria #7
Clostridia = subset_taxa(ps.class, Class=="Clostridia")
Clostridia #50
Erysipelotrichia = subset_taxa(ps.class, Class=="Erysipelotrichia")
Erysipelotrichia #11
Gammaproteobacteria = subset_taxa(ps.class, Class=="Gammaproteobacteria")
Gammaproteobacteria #5

###Plot unique vs. shared ASVs Mouse 2
venn=venneuler(c('A' = 592, 'B' = 93, 'A&B' = 214))
venn$labels = c("", "")
plot(venn, col=c("red", "blue")) 

###Total number of reads in each mouse Inner sample
sdt = data.table(as(sample_data(Mouse1.Inner.trim), "data.frame"),
                 TotalReads = sample_sums(Mouse1.Inner.trim), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt #41081

sdt = data.table(as(sample_data(Mouse2.Inner.trim), "data.frame"),
                 TotalReads = sample_sums(Mouse2.Inner.trim), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt #33628

sdt = data.table(as(sample_data(Mouse3.Inner.trim), "data.frame"),
                 TotalReads = sample_sums(Mouse3.Inner.trim), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt #39755

sdt = data.table(as(sample_data(Mouse4.Inner.trim), "data.frame"),
                 TotalReads = sample_sums(Mouse4.Inner.trim), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt #45474

sdt = data.table(as(sample_data(Mouse5.Inner.trim), "data.frame"),
                 TotalReads = sample_sums(Mouse5.Inner.trim), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt #37052


###Number of reads of Inner unique ASVs at family level by mouse
Mouse1.Inner.unique 
Mouse1.Inner.unique.family<-tax_glom(Mouse1.Inner.unique,taxrank="Family")
Mouse1.Inner.unique.family 

taxa_names(Mouse1.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse1.Inner.unique.family)))
otu.table<-otu_table(Mouse1.Inner.unique.family)
write.csv(otu.table,"Mouse1.Inner.unique.family_OTU.csv")
tax.table<-tax_table(Mouse1.Inner.unique.family)
write.csv(tax.table,"Mouse1.Inner.unique.family_TAX.csv")

Mouse2.Inner.unique
Mouse2.Inner.unique.family<-tax_glom(Mouse2.Inner.unique,taxrank="Family")
Mouse2.Inner.unique.family 

taxa_names(Mouse2.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse2.Inner.unique.family)))
otu.table<-otu_table(Mouse2.Inner.unique.family)
write.csv(otu.table,"Mouse2.Inner.unique.family_OTU.csv")
tax.table<-tax_table(Mouse2.Inner.unique.family)
write.csv(tax.table,"Mouse2.Inner.unique.family_TAX.csv")

Mouse3.Inner.unique
Mouse3.Inner.unique.family<-tax_glom(Mouse3.Inner.unique,taxrank="Family")
Mouse3.Inner.unique.family 

taxa_names(Mouse3.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse3.Inner.unique.family)))
otu.table<-otu_table(Mouse3.Inner.unique.family)
write.csv(otu.table,"Mouse3.Inner.unique.family_OTU.csv")
tax.table<-tax_table(Mouse3.Inner.unique.family)
write.csv(tax.table,"Mouse3.Inner.unique.family_TAX.csv")

Mouse4.Inner.unique
Mouse4.Inner.unique.family<-tax_glom(Mouse4.Inner.unique,taxrank="Family")
Mouse4.Inner.unique.family 

taxa_names(Mouse4.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse4.Inner.unique.family)))
otu.table<-otu_table(Mouse4.Inner.unique.family)
write.csv(otu.table,"Mouse4.Inner.unique.family_OTU.csv")
tax.table<-tax_table(Mouse4.Inner.unique.family)
write.csv(tax.table,"Mouse4.Inner.unique.family_TAX.csv")

Mouse5.Inner.unique
Mouse5.Inner.unique.family<-tax_glom(Mouse5.Inner.unique,taxrank="Family")
Mouse5.Inner.unique.family 

taxa_names(Mouse5.Inner.unique.family) <- paste0("Seq", seq(ntaxa(Mouse5.Inner.unique.family)))
otu.table<-otu_table(Mouse5.Inner.unique.family)
write.csv(otu.table,"Mouse5.Inner.unique.family_OTU.csv")
tax.table<-tax_table(Mouse5.Inner.unique.family)
write.csv(tax.table,"Mouse5.Inner.unique.family_TAX.csv")

###Number of reads of Shared ASVs at family level by mouse
Mouse1.Shared.unique 
Mouse1.Shared.unique.family<-tax_glom(Mouse1.Shared.unique,taxrank="Family")
Mouse1.Shared.unique.family 

taxa_names(Mouse1.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse1.Shared.unique.family)))
otu.table<-otu_table(Mouse1.Shared.unique.family)
write.csv(otu.table,"Mouse1.Shared.unique.family_OTU.csv")
tax.table<-tax_table(Mouse1.Shared.unique.family)
write.csv(tax.table,"Mouse1.Shared.unique.family_TAX.csv")

Mouse2.Shared.unique
Mouse2.Shared.unique.family<-tax_glom(Mouse2.Shared.unique,taxrank="Family")
Mouse2.Shared.unique.family 

taxa_names(Mouse2.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse2.Shared.unique.family)))
otu.table<-otu_table(Mouse2.Shared.unique.family)
write.csv(otu.table,"Mouse2.Shared.unique.family_OTU.csv")
tax.table<-tax_table(Mouse2.Shared.unique.family)
write.csv(tax.table,"Mouse2.Shared.unique.family_TAX.csv")

Mouse3.Shared.unique
Mouse3.Shared.unique.family<-tax_glom(Mouse3.Shared.unique,taxrank="Family")
Mouse3.Shared.unique.family 

taxa_names(Mouse3.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse3.Shared.unique.family)))
otu.table<-otu_table(Mouse3.Shared.unique.family)
write.csv(otu.table,"Mouse3.Shared.unique.family_OTU.csv")
tax.table<-tax_table(Mouse3.Shared.unique.family)
write.csv(tax.table,"Mouse3.Shared.unique.family_TAX.csv")

Mouse4.Shared.unique
Mouse4.Shared.unique.family<-tax_glom(Mouse4.Shared.unique,taxrank="Family")
Mouse4.Shared.unique.family 

taxa_names(Mouse4.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse4.Shared.unique.family)))
otu.table<-otu_table(Mouse4.Shared.unique.family)
write.csv(otu.table,"Mouse4.Shared.unique.family_OTU.csv")
tax.table<-tax_table(Mouse4.Shared.unique.family)
write.csv(tax.table,"Mouse4.Shared.unique.family_TAX.csv")

Mouse5.Shared.unique
Mouse5.Shared.unique.family<-tax_glom(Mouse5.Shared.unique,taxrank="Family")
Mouse5.Shared.unique.family 

taxa_names(Mouse5.Shared.unique.family) <- paste0("Seq", seq(ntaxa(Mouse5.Shared.unique.family)))
otu.table<-otu_table(Mouse5.Shared.unique.family)
write.csv(otu.table,"Mouse5.Shared.unique.family_OTU.csv")
tax.table<-tax_table(Mouse5.Shared.unique.family)
write.csv(tax.table,"Mouse5.Shared.unique.family_TAX.csv")


#Figure 4

##Weighted Unifrac PCoA
Analysis<-Kellys_Layers
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduW = ordinate(GP1, "PCoA","unifrac",weighted=TRUE)
pW = plot_ordination(GP1, orduW, color = "Location", shape = "Time", title = "Weighted Unifrac") + geom_point(size=6) + scale_colour_manual(values=c("black","purple"))
pW + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) #+ylim(-0.1,0.1) +xlim(-0.2,0.0)

##Unweighted Unifrac PCoA
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduU = ordinate(GP1, "PCoA","unifrac",weighted=FALSE)
pU = plot_ordination(GP1, orduU, color = "Location", shape = "Time", title = "Unweighted Unifrac") + geom_point(size=6) + scale_colour_manual(values=c("black","purple"))
pU + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  #+ ylim(-0.2,0.1) +xlim(-0.30,0.30) 

##Rarefied Diversity Table 
Analysis<-Kellys_Layers

###Reads per sample
sdt = data.table(as(sample_data(Kellys_Layers), "data.frame"),
                 TotalReads = sample_sums(Kellys_Layers), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt

###Remove Mouse with lowest number of reads
Kellys_Layers_mod = subset_samples(Kellys_Layers, MouseID != "Kvon19")
Kellys_Layers_mod

###Rarefy Table
set.seed(42)
rarefied_analysis = rarefy_even_depth(Kellys_Layers_mod)
rarefied_analysis

###Plot Diversity
p = plot_richness(rarefied_analysis, x="Time", color="Location", measures=c("Observed","Shannon"))
p + geom_point(size=5, alpha=0.7) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 

rarefied_values<-p$data
rarefied_values
write.csv(rarefied_values, file = "KellysABX_rarefied_values.csv")

###Enter Dataset for bar plots
Analysis<-Kellys_Layers

###Set barplot Colors
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
PhylumList = unique(tax_table(Analysis)[,"Phylum"])
PhylumPalette = getPalette(length(PhylumList))
names(PhylumPalette) = PhylumList

##Mucosal
###Normalize dataset
subset=transform_sample_counts(Kellys_1st50,function(x)x/sum(x))

###Bar plots phylum level
p = plot_bar(subset, x= "Time", fill = "Phylum") + scale_color_manual(values= PhylumPalette) + scale_fill_manual(values= PhylumPalette)
p + geom_bar(aes(color= Phylum, fill=Phylum), stat="identity", position="stack") + theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

##Luminal
###Normalize dataset
subset=transform_sample_counts(Kellys_Luminal,function(x)x/sum(x))

###Bar plots phylum level
p = plot_bar(subset, x= "Time", fill = "Phylum") + scale_color_manual(values= PhylumPalette) + scale_fill_manual(values= PhylumPalette)
p + geom_bar(aes(color= Phylum, fill=Phylum), stat="identity", position="stack") + theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

###Condense by taxonomic level
ps.phylum<-tax_glom(Analysis,taxrank="Phylum")

###Percent Abundance
otu.table<-otu_table(ps.phylum)
write.csv(otu.table,"KellysABX.phylum.otu.table.csv")
tax.table<-tax_table(ps.phylum)
write.csv(tax.table,"KellysABX.phylum.tax.table.csv")

###Change ASV IDs
taxa_names(Kellys_Day10) <- paste0("Seq", seq(ntaxa(Kellys_Day10)))

###Condense by taxonomic level
ps.phylum<-tax_glom(Kellys_Day10,taxrank="Phylum")
ps.class<-tax_glom(Kellys_Day10,taxrank="Class")

##Differential Abundance based on Location: Phylum
dds = phyloseq_to_deseq2(ps.phylum, ~ Location)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.phylum)[rownames(sigtab_WD), ], "matrix"))
diff.abund

###Plot counts for OTU determined differentially abundant and group by MouseID
data <- plotCounts(dds, "Seq2762", intgroup=c("Location","MouseID"), returnData=TRUE)
plot<-ggplot(data, aes(x=Location, y=count, color=MouseID, group=MouseID)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Phylum Firmicutes"))

data <- plotCounts(dds, "Seq4554", intgroup=c("Location","MouseID"), returnData=TRUE)
plot<-ggplot(data, aes(x=Location, y=count, color=MouseID, group=MouseID)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Phylum Proteobacteria"))

##Differential Abundance based on Location: Class
dds = phyloseq_to_deseq2(ps.class, ~ Location)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.class)[rownames(sigtab_WD), ], "matrix"))
diff.abund

###Plot counts for OTU determined differentially abundant and group by MouseID
data <- plotCounts(dds, "Seq2762", intgroup=c("Location","MouseID"), returnData=TRUE)
plot<-ggplot(data, aes(x=Location, y=count, color=MouseID, group=MouseID)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Class Clostridia"))

###Plot counts for OTU determined differentially abundant and group by MouseID
data <- plotCounts(dds, "Seq3413", intgroup=c("Location","MouseID"), returnData=TRUE)
plot<-ggplot(data, aes(x=Location, y=count, color=MouseID, group=MouseID)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Class Bacilli"))

###Plot counts for OTU determined differentially abundant and group by MouseID
data <- plotCounts(dds, "Seq4554", intgroup=c("Location","MouseID"), returnData=TRUE)
plot<-ggplot(data, aes(x=Location, y=count, color=MouseID, group=MouseID)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Class Gammaproteobacteria"))

###Make Excel sheet with normalized counts
dds <- estimateSizeFactors(dds)
dds.counts<-counts(dds, normalized=TRUE)
write.csv(dds.counts,"KellysABX.norm.counts.class.csv")


