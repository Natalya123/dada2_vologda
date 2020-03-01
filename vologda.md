# dada2_vologda
library(dada2); packageVersion("dada2")
path <- "~/storage/vologda/vologda/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# Trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "~/storage/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/storage/tax/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")


##############################################################################
# Metadata block #
##############################################################################

map <- read.csv2("~/storage/vologda/vologda/map_vologda_full.csv", sep = ';', row.names = 1)
mdat <- map[order(rownames(map)), ]
seqtab2 <- seqtab.nochim[order(rownames(seqtab.nochim)), ]

#viewer
View(cbind(as.vector(rownames(seqtab2)), as.vector(rownames(mdat))))




ps <- phyloseq(otu_table(seqtab2, taxa_are_rows=FALSE),
               sample_data(mdat),
               tax_table(taxa))
ps

ps <- prune_samples(sample_sums(ps)>=4500, ps)

ps
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Culture", title="Bray NMDS")


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
top20
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)




#Barplots
physeq2 = filter_taxa(ps, function(x) mean(x) > 0.1, TRUE)
physeq2
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
ph_ps <- tax_glom(physeq3, taxrank = 'Phylum')
data <- psmelt(ph_ps) # create dataframe from phyloseq object
data$Phylum <- as.character(data$Phylum) #convert to character
data$Phylum[data$Abundance < 0.01] <- "< 1% abund."
medians <- ddply(data, ~Phylum, function(x) c(median=median(x$Abundance)))
remainder <- medians[medians$median <= 0.01,]$Phylum
p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Phylum))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                               "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

write.csv(cbind(t(seqtab2), taxa), "seqtab.csv", quote=FALSE)
################################################



permanova.al <- function(ps, dist = "bray"){ # не обязательно использовать именно Брей-Кёртиса
  require(phyloseq)
  require(vegan)
  dist <- phyloseq::distance(ps, dist)
  metadata <- as(sample_data(ps), "data.frame")
  ad <- adonis2(dist ~ Lime, data = metadata)  #вместо Al известкование
  return(ad)
}
permanova.al(ps)

#######################################################################3


Des.soil.w.simper <- function(ps){
  require(DESeq2)
  require(vegan)
  require(tibble)
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local")
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(as.data.frame(as.data.frame(t(diagdds@assays@.xData$data$mu))), by=list(samp$Repeats), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])
  return(nice)
}
Des.soil.w.simper(ps) -> des_non_simp

