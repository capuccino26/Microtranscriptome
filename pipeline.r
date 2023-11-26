#Prepare Env
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.16")
BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")
library(dada2)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
library(edgeR)
library(gt)
library(tidyverse)
library(gtExtras)
library(webshot)

#Path
path <- "./FOLDER"
pathf <- sprintf("%s/READS/",path)

#Setting Samples
fnFs <- sort(list.files(pathf, pattern="R1", full.names = TRUE))
fnRs <- sort(list.files(pathf, pattern="R2", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

png(sprintf("%s/QUALITYF.png",path),width=4000, height=2000,res = 600)
plotQualityProfile(fnFs[1:2])
dev.off()
png(sprintf("%s/QUALITYR.png",path),width=4000, height=2000,res = 600)
plotQualityProfile(fnRs[1:2])
dev.off()

#Trimming and Filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)

#Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

png(sprintf("%s/ERRORF.png",path),width=2000, height=2000,res = 100)
plotErrors(errF, nominalQ=TRUE)
dev.off()
png(sprintf("%s/ERRORR.png",path),width=2000, height=2000,res = 100)
plotErrors(errR, nominalQ=TRUE)
dev.off()

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
#dada-class: object describing DADA2 denoising results
#16 sequence variants were inferred from 5220 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

#Merge paired end
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Sequence Table
seqtab <- makeSequenceTable(mergers)
##Stats
dim(seqtab)
#[1] 26 85
##85 ASVs das 26 amostras
table(nchar(getSequences(seqtab)))
#315 320 322 324 328 329 330 331 332 333 462 
#  5   1  18   9   6   5  13   2  21   4   1 
#A sequencia não foi removida dos primers

#Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 28 bimeras out of 85 input sequences.
dim(seqtab.nochim)
#[1] 26 57
sum(seqtab.nochim)/sum(seqtab)
#[1] 0.9817287

#Summary
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(as.data.frame(track), file=sprintf("%s/SUMMARY.csv",path))

#Set TaxID
train <- readDNAStringSet(sprintf("%s/Nematode_ITS2_1.0.0_idtaxa.fasta",path))
tax <- read_tsv(sprintf("%s/Nematode_ITS2_1.0.0_idtaxa.tax",path))
trainingSet <- LearnTaxa(train, names(train), tax)
dna <- DNAStringSet(getSequences(seqtab.nochim))
idtaxa <- IdTaxa(dna,trainingSet,strand = "both",threshold = 60,bootstraps = 100,processors = NULL,verbose = TRUE,type = "extended")
#Final analysis
taxid <- t(sapply(idtaxa, function(x) setNames(x$taxon, x$rank)))[, -1]
write.table(as.data.frame(taxid), file=sprintf("%s/RESULTS.csv",path))

#Phyloseq
samp_data <-  data.frame(row.names = sample.names,sample = sample.names)
samp_data[c("G1","G2","G3","G4","G5","G6","G7","G8","G9"),2]<-"1"
samp_data[c("G10","G11","G12","G13","G14","G15","G16","G17","G18","G19"),2]<-"2"
samp_data[c("G20","G21","G22","G23","G24","G25","G26"),2]<-"3"
colnames(samp_data)[2]<-"Group"
asvs <- paste0("ASV_", 1:length(dna))
rownames(taxid) <- asvs
colnames(seqtab.nochim) <- asvs
names(dna) <- asvs
physeq = phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),tax_table(taxid),sample_data(samp_data),dna)

#Estimate RICHNESS
estimate <- estimate_richness(physeq, split = TRUE, measures = NULL)

#Visual
theme_set(theme_bw())
png(sprintf("%s/RICHNESS.png",path),width=4000, height=4000,res = 600)
plot_richness(physeq, x="Group", measures=c("Shannon", "Simpson"), color="Group")
dev.off()

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(physeq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

png(sprintf("%s/Bray-Curtis.png",path),width=4000, height=4000,res = 600)
plot_ordination(ps.prop, ord.nmds.bray, color="Group", title="Bray NMDS")
dev.off()

#Barplot
top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
#Family
png(sprintf("%s/ABFAMILY.png",path),width=4000, height=4000,res = 600)
plot_bar(ps.top20, x="Group", fill="family") + facet_wrap(~Group, scales="free_x")
dev.off()
#Genus
png(sprintf("%s/ABGENUS.png",path),width=4000, height=4000,res = 600)
plot_bar(ps.top20, x="Group", fill="genus") + facet_wrap(~Group, scales="free_x")
dev.off()
#Species
png(sprintf("%s/ABSPECIES.png",path),width=4000, height=4000,res = 600)
plot_bar(ps.top20, x="Group", fill="species") + facet_wrap(~Group, scales="free_x")
dev.off()

#Execute phylotoedger.r
dge = phyloseq_to_edgeR(physeq, group="Group")
# Perform binary test
et = exactTest(dge)
# Extract values from test results
tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.001
sigtab = res[(res$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
write.table(as.data.frame(sigtab), file=sprintf("%s/DEGTAB.csv",path))

#Plot family
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(genus))
x = tapply(sigtabgen$logFC, sigtabgen$family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$family = factor(as.character(sigtabgen$family), levels = names(x))
x = tapply(sigtabgen$logFC, sigtabgen$genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$genus = factor(as.character(sigtabgen$genus), levels = names(x))
png(sprintf("%s/DEGFAMILY.png",path),width=4000, height=4000,res = 600)
ggplot(sigtabgen, aes(x = genus, y = logFC, color = family)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
dev.off()

#Plot Genus
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(genus))
x = tapply(sigtabgen$logFC, sigtabgen$genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$genus = factor(as.character(sigtabgen$genus), levels = names(x))
x = tapply(sigtabgen$logFC, sigtabgen$genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$genus = factor(as.character(sigtabgen$genus), levels = names(x))
png(sprintf("%s/DEGGENUS.png",path),width=4000, height=4000,res = 600)
ggplot(sigtabgen, aes(x = genus, y = logFC, color = genus)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
dev.off()

#Plot Species
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(genus))
x = tapply(sigtabgen$logFC, sigtabgen$species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$species = factor(as.character(sigtabgen$species), levels = names(x))
x = tapply(sigtabgen$logFC, sigtabgen$genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$genus = factor(as.character(sigtabgen$genus), levels = names(x))
png(sprintf("%s/DEGspecies.png",path),width=4000, height=4000,res = 600)
ggplot(sigtabgen, aes(x = genus, y = logFC, color = species)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
dev.off()

#Generate Tables
##Differential Abbundance
sigtabb <- subset(sigtab, select = -c(13,14,15,16,17,18,19,20))
sigtabb <- cbind(ASV = rownames(sigtabb), sigtabb)
rownames(sigtabb) <- 1:nrow(sigtabb)
gt(sigtabb) %>%  
    gt_theme_guardian() %>% 
    cols_width(.,everything()~px(150)) %>%
    tab_header(title = "Differential Abundance") %>%
	tab_options(row_group.as_column = TRUE) %>%
	tab_options(heading.align = "center") %>%
    cols_align(align = "center") %>%
    gtsave("tab_16.html")
webshot("tab_16.html" , "sigtabb.pdf", delay = 0.2,vwidth=2000)

##Cutoffs
gt(track) %>%  
    gt_theme_guardian() %>% 
    cols_width(.,everything()~px(150)) %>%
    tab_header(title = "Differential Abundance") %>%
	tab_options(row_group.as_column = TRUE) %>%
	tab_options(heading.align = "center") %>%
    cols_align(align = "center") %>%
    gtsave("tab_16.html")
webshot("tab_16.html" , "track.pdf", delay = 0.2,vwidth=2000)