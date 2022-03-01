library(DESeq2)
library(plotly)
library(VennDiagram)
library(UpSetR)
library(grid)
# library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(plyr)

################################################################################

#Functions

save_pheatmap_pdf <- function(x, filename) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = 210, height = 297, units = "mm")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

heat_ds <- function(data, filt, sample_names) {
  genes <- filter(ens.all.out, ensgene %in% filt)
  heat.df <- as.data.frame(data)
  heat.df <- rownames_to_column(heat.df, "ensgene")
  heat.df <- left_join(heat.df, genes)
  heat.df <- filter(heat.df, !duplicated(gene))
  heat.df <- filter(heat.df, !is.na(gene))
  rownames(heat.df) <- heat.df$gene
  return(as.matrix(select(heat.df,all_of(sample_names))))
}

annoDE <- function(data, p_thresh=0.05, l2fc_thresh=1){
  res.df <- as.data.frame(data)
  res.filt <- res.df[complete.cases(res.df),]
  res.filt <- rownames_to_column(res.filt,"ensgene")
  res.filt <- filter(res.filt, padj<p_thresh & abs(res.filt$log2FoldChange) >l2fc_thresh)
  return(left_join(res.filt, ens.all.out))
}

volcano_int <- function(data) {
  g <- ggplot(data, aes(x=log2FoldChange, y=-log10(padj),
                             name=gene)) +
    geom_point(size =2, alpha=0.4)+
    theme_bw()+
    theme(legend.position = "none",
          panel.border = element_blank(), axis.line = element_line(size=1),
          axis.text = element_text(size=12, face="bold", colour = "black"),
          axis.title = element_text(size=12, face="bold", colour = "black"),
          legend.title = element_text(size=12, face="bold", colour = "black"),
          legend.text = element_text(size=10, face="bold", colour = "black"),
          axis.ticks = element_blank())
  ggplotly(g)
}

volcano_ny <- function(data) {
  g <- ggplot(data, aes(x=log2FoldChange, y=avg.memoryDP, color=analysis,
                      name=gene)) +
    geom_point(size =2, alpha=0.4)+
    theme_bw()+
    theme(panel.border = element_blank(), axis.line = element_line(size=1),
          axis.text = element_text(size=12, face="bold", colour = "black"),
          axis.title = element_text(size=12, face="bold", colour = "black"),
          legend.title = element_text(size=12, face="bold", colour = "black"),
          legend.text = element_text(size=10, face="bold", colour = "black"),
          axis.ticks = element_blank())
  ggplotly(g)
}

colname_index <- function(name, data) {
  len <- length(name)
  index <- numeric(len)
  for (i in 1:len){
    ind <- which(colnames(data)==name[i])
    if (!is.null(ind) & length(ind)>0) {
      index[i] <- ind
    }
  }
  return(index[index>0])
}

designate <- function(to_check, list_A, name_A, list_B, name_B){
  name_AB <- paste0(name_A, " & ",name_B)
  list_AB <- intersect(list_A, list_B)
  list_A_only <- list_A[!list_A %in% list_AB]
  list_B_only <- list_B[!list_B %in% list_AB]
  ds.desig <- data.frame(ensgene=to_check, designation="other")
  ds.desig[ds.desig$ensgene %in% list_A_only,2] <- name_A
  ds.desig[ds.desig$ensgene %in% list_B_only,2] <- name_B
  ds.desig[ds.desig$ensgene %in% list_AB,2] <- name_AB
  ds.desig$designation <- factor(ds.desig$designation, levels=c(name_A, name_B, name_AB, "other"))
  return(ds.desig)
}

deg <- function(data, interest, control, group="group", min_counts = 250){
  avg.i <- paste0("avg.",interest)
  avg.c <- paste0("avg.",control)
  DEres <- results(data, contrast=c(group,interest,control))
  degs.all <- left_join(annoDE(DEres),avg.spleen)
  degs <- filter(degs.all, eval(as.name(avg.c)) >= min_counts | eval(as.name(avg.i)) >= min_counts)
  return(degs)
}

return_set <- function(data, cs) {
  n <- length(cs)
  d.filt <- data
  d.filt$sum <- apply(data,1,sum)
  for (i in 1:n) {
    d.filt <- filter(d.filt, d.filt[[cs[i]]] ==1)
  }
  d.filt <- filter(d.filt, sum ==n)
  return(rownames(d.filt))
}

top_l2fc <- function(data, ascend=F, top=20){
  if (ascend==T) {
    d.sorted <- arrange(data, log2FoldChange)
    d.top <- d.sorted[1:top,]
  }
  else{
    d.sorted <- arrange(data, desc(log2FoldChange))
    d.top <- d.sorted[1:top,]
  }
  d.top$gene <- factor(d.top$gene, levels=rev(c(d.top$gene)))
  return(d.top)
}

lollipop_plot <- function(data, point_color="#e1af00"){
  if(grepl("down",deparse(substitute(data)))){
    ggplot(data, aes(x=log2FoldChange, y=gene, label=round(log2FoldChange, 2)))+
      geom_segment(aes(x=0, y=gene, xend=log2FoldChange,yend=gene), size=2,color=point_color)+
      geom_point(stat = "identity", size = 15, color=point_color)+
      geom_text(color="white",size=4, fontface="bold")+
      scale_x_continuous(limits = c(-10,0), expand = c(0,0))+
      scale_y_discrete(position = "right")+
      theme_bw()+
      theme(panel.background = element_rect(fill = "transparent", colour=NA),
            plot.background = element_rect(fill = "transparent", colour=NA),
            legend.background = element_rect(fill = "transparent", colour=NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            legend.key = element_rect(fill = "transparent", colour = NA),
            panel.border = element_blank(),
            axis.text.y = element_text(size=12, face="bold.italic", colour = "#262f46"),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())
  }
  else{
    ggplot(data, aes(x=log2FoldChange, y=gene, label=round(log2FoldChange, 2)))+
      geom_segment(aes(x=0, y=gene, xend=log2FoldChange,yend=gene), size=2,color=point_color)+
      geom_point(stat = "identity", size = 15, color=point_color)+
      geom_text(color="white",size=4, fontface="bold")+
      scale_x_continuous(limits = c(0,10), expand = c(0,0))+
      theme_bw()+
      theme(panel.background = element_rect(fill = "transparent", colour=NA),
            plot.background = element_rect(fill = "transparent", colour=NA),
            legend.background = element_rect(fill = "transparent", colour=NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            legend.key = element_rect(fill = "transparent", colour = NA),
            panel.border = element_blank(),
            axis.text.y = element_text(size=12, face="bold.italic", colour = "#262f46"),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())
  }

}


DEGbar_plot <- function(data, point_color="#e1af00"){
  if(grepl("down",deparse(substitute(data)))){
    ggplot(data, aes(x=log2FoldChange, y=gene, label=round(log2FoldChange, 2)))+
      geom_segment(aes(x=0, y=gene, xend=log2FoldChange,yend=gene), size=8,color=point_color)+
      geom_text(color="white",size=4, hjust=0, nudge_x = 0.1)+
      scale_x_continuous(limits = c(-10,0), expand = c(0,0))+
      scale_y_discrete(position = "right")+
      theme_bw()+
      theme(panel.background = element_rect(fill = "transparent", colour=NA),
            plot.background = element_rect(fill = "transparent", colour=NA),
            legend.background = element_rect(fill = "transparent", colour=NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            legend.key = element_rect(fill = "transparent", colour = NA),
            panel.border = element_blank(),
            axis.text.y = element_text(size=12, colour = "black", face="italic"),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())
  }
  else{
    ggplot(data, aes(x=log2FoldChange, y=gene, label=round(log2FoldChange, 2)))+
      geom_segment(aes(x=0, y=gene, xend=log2FoldChange,yend=gene), size=8,color=point_color)+
      geom_text(color="white",size=4, hjust=1, nudge_x = -0.1)+
      scale_x_continuous(limits = c(0,10), expand = c(0,0))+
      theme_bw()+
      theme(panel.background = element_rect(fill = "transparent", colour=NA),
            plot.background = element_rect(fill = "transparent", colour=NA),
            legend.background = element_rect(fill = "transparent", colour=NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            legend.key = element_rect(fill = "transparent", colour = NA),
            panel.border = element_blank(),
            axis.text.y = element_text(size=12, colour = "black", face="italic"),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())
  }
  
}


balloon_prep <- function(data){
  bal.ds.prep <- left_join(data, ens.all.out)
  bal.ds.prep <- select(bal.ds.prep, ensgene, gene, designation)
  
  bal.ds.table <- filter(avg.spleen, ensgene %in% bal.ds.prep$ensgene)
  bal.ds.table <- left_join(bal.ds.prep, bal.ds.table)
  bal.ds <- gather(bal.ds.table, "population", "avg.count", 4:9)
  bal.ds$population <- substring(bal.ds$population, 5)
  bal.ds$population <- mapvalues(bal.ds$population, from=c("naiveDN","memoryDN","naiveSP","memorySP","naiveDP","memoryDP"), to=c("ctrl-DN","mem-DN","ctrl-SP","mem-SP","ctrl-DP","mem-DP"))
  bal.ds$population <- factor(bal.ds$population, levels=c("ctrl-DN","mem-DN","ctrl-SP","mem-SP","ctrl-DP","mem-DP"))
  bal.gene.order <- arrange(bal.ds[bal.ds$population=="mem-DP",], desc(designation), avg.count)
  bal.gene.order <- bal.gene.order$gene
  bal.ds$gene <- factor(bal.ds$gene, levels=bal.gene.order)
  return(bal.ds)
}

balloon_plot <- function(data, point_color="Grey60"){
  ggplot(data, aes(x=gene, y=population)) +
    geom_point(aes(size=avg.count, shape=designation), color=point_color)+
    scale_shape_manual(values=c(16,15))+
    scale_size_area( max_size = 10, breaks = c(0,1000,2000,3000,4000,5000,6000)) +
    scale_y_discrete(limits=rev)+
    ylab("gene")+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent", colour=NA),
          plot.background = element_rect(fill = "transparent", colour=NA),
          legend.background = element_rect(fill = "transparent", colour=NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key = element_rect(fill = "transparent", colour = NA),
          panel.border = element_blank(), axis.line = element_line(size=0.5),
          axis.text = element_text(size=10, colour = "black"),
          axis.title = element_blank(),
          legend.title = element_text(size=10, face="bold", colour = "black"),
          legend.text = element_text(size=10, face="bold", colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face="italic"))
}

curate_balloon_prep <- function(data){
  bal.ds.prep <- left_join(data, ens.all.out)
  bal.ds.prep <- select(bal.ds.prep, ensgene, gene, designation, process)
  
  bal.ds.table <- filter(avg.spleen, ensgene %in% bal.ds.prep$ensgene)
  bal.ds.table <- left_join(bal.ds.prep, bal.ds.table)
  bal.ds <- gather(bal.ds.table, "population", "avg.count", 5:10)
  bal.ds$population <- substring(bal.ds$population, 5)
  bal.ds$population <- mapvalues(bal.ds$population, from=c("naiveDN","memoryDN","naiveSP","memorySP","naiveDP","memoryDP"), to=c("ctrl-DN","mem-DN","ctrl-SP","mem-SP","ctrl-DP","mem-DP"))
  bal.ds$population <- factor(bal.ds$population, levels=c("ctrl-DN","mem-DN","ctrl-SP","mem-SP","ctrl-DP","mem-DP"))
  bal.gene.order <- arrange(bal.ds[bal.ds$population=="mem-DP",], avg.count)
  bal.gene.order <- bal.gene.order$gene
  bal.ds$gene <- factor(bal.ds$gene, levels=bal.gene.order)
  return(bal.ds)
}


curate_balloon <- function(data, coloring){
  ggplot(data, aes(x=gene, y=population)) +
    geom_point(aes(size=avg.count/500, shape=designation, color=process))+
    scale_shape_manual(values=c(16,15))+
    scale_color_manual(values = coloring)+
    scale_size_identity(guide="legend") +
    scale_y_discrete(limits=rev)+
    ylab("gene")+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent", colour=NA),
          plot.background = element_rect(fill = "transparent", colour=NA),
          legend.background = element_rect(fill = "transparent", colour=NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key = element_rect(fill = "transparent", colour = NA),
          panel.border = element_blank(), axis.line = element_line(size=0.5),
          axis.text = element_text(size=10, colour = "black"),
          axis.title = element_blank(),
          legend.title = element_text(size=10, face="bold", colour = "black"),
          legend.text = element_text(size=10, face="bold", colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face="italic"))
}

################################################################################

#Global Variables
#################



redblue = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)


#Preparing Datasets
###################

#Reading in setup & counts from pooled files
counts <- read.csv("count-matrix.csv", row.names = 1)
meta <- read.csv("setup.csv", row.names = 1)

#Selecting spleen samples without second infection
meta.spleen.noLpn <- filter(meta, organ=="spleen" & !grepl("Lpn", infection) & exclude =="no")
names.allsamples <- rownames(meta.spleen.noLpn)
countdata <- select(counts, all_of(colname_index(names.allsamples, counts)))


#Retrieving gene information from ensembl gene name and saving it as csv

# ens.all <- rownames_to_column(counts, "ensgene")
# ens.all <- select(ens.all, ensgene)
# ensembl101 = useEnsembl("ensembl", version = 101)
# ensembl101 = useDataset("mmusculus_gene_ensembl", mart = ensembl101)
# ens.all.out = getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position",
#                                   "end_position","strand", "gene_biotype",
#                                   "external_gene_name", "description", "entrezgene_id"),
#                    filters = c("ensembl_gene_id"), values = ens.all$ensgene,
#                    mart = ensembl101)
# colnames(ens.all.out)[1] <- "ensgene"
# colnames(ens.all.out)[7] <- "gene"
# write.csv(ens.all.out,"ensembl_dataframe.csv", row.names = F)

#Reading in the gene info
ens.all.out <- read.csv("ensembl_dataframe_noGO.csv")
ens.all.out <- distinct(ens.all.out, ensgene, .keep_all = T)

ens.all.go <- read.csv("ensembl_dataframe.csv")



#GO-based categorization
ens.go.mmb <- filter(ens.all.go, grepl("integral component of membrane|cell surface|anchored component of plasma membrane", GO))
ens.go.mmb <- distinct(ens.go.mmb, ensgene, .keep_all = T)

ens.go.act <- filter(ens.all.go, grepl("activation|immune response|cytokine|positive regulation of cell cycle", GO))
ens.go.act <- distinct(ens.go.act, ensgene, .keep_all = T)

ens.go.tax <- filter(ens.all.go, grepl("taxis|migration|chemokine|cell adhesion", GO))
ens.go.tax <- distinct(ens.go.tax, ensgene, .keep_all = T)

ens.go.tkf <- filter(ens.all.go, GO=="DNA-binding transcription factor activity")
ens.go.tkf <- distinct(ens.go.tkf, ensgene, .keep_all = T)

ens.go.tkreg <- filter(ens.all.go, grepl("regulation of DNA-binding transcription factor activity", GO))
ens.go.tkreg <- filter(ens.go.tkreg, ensgene %in% setdiff(ens.go.tkreg$ensgene, ens.go.tkf$ensgene))
ens.go.tkreg <- distinct(ens.go.tkreg, ensgene, .keep_all = T)

ens.go.acttax <- filter(ens.go.act, ensgene %in% ens.go.tax$ensgene)


#Shaping the experimental condition = "coldata"
coldata <- meta.spleen.noLpn
coldata$condition <- factor(coldata$infection, levels = c("naive","memory"))
coldata$population <- factor(coldata$population, levels = c("DN", "SP", "DP"))


#DESeq for all of spleen samples
spleen.dds <- DESeqDataSetFromMatrix(countData = countdata,
                                               colData = coldata,
                                               design = ~population + condition + population:condition)
spleen.dds <- DESeq(spleen.dds)

#Normalized counts
spleen.estim <- estimateSizeFactors(spleen.dds)
spleen.counts <- counts(spleen.estim, normalized=T)
spleen.counts <- rownames_to_column(as.data.frame(spleen.counts), "ensgene")
spleen.counts <- left_join(spleen.counts, ens.all.out)
write.csv(spleen.counts, "spleen_normalized-counts.csv", row.names = F)

#Variables from Spleen DESeq
############################
spleen.ens <- spleen.counts$ensgene
spleen.vst <- varianceStabilizingTransformation(spleen.dds)
spleen.vst.assay <- assay(spleen.vst)




#PCA plot
spleen.pca <- plotPCA(spleen.vst, intgroup=c("condition","population"), returnData=TRUE)
percentVar <- round(100 * attr(spleen.pca, "percentVar"))
ggplot(spleen.pca, aes(PC1, PC2, shape=condition, fill = population)) +
  geom_point(size=4, stroke = 0.5) +
  scale_fill_manual(values = c("black","#55a0fb", "#ff8080"))+
  scale_shape_manual(values = c(21,22))+
  # scale_shape_manual(values = c(16,18,17))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(panel.background = element_rect(fill = "transparent", colour=NA),
        plot.background = element_rect(fill = "transparent", colour=NA),
        legend.background = element_rect(fill = "transparent", colour=NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        axis.line.x.bottom = element_line(color = '#262f46', size=0.5),
        axis.line.y.left   = element_line(color = '#262f46',size=0.5),
        panel.border = element_blank(),
        axis.text = element_text(size=10, colour = "#262f46"),
        axis.title = element_text(size=10, colour = "#262f46"),
        legend.title = element_text(size=10, colour = "#262f46"),
        legend.text = element_text(size=8, colour = "#262f46"),
        aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "#262f46", size=0.5))
ggsave("0_pca.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)


#Calculating averages of each group and creating dataset with it
inf.elements <- unique(meta.spleen.noLpn$infection)
pop.elements <- unique(meta.spleen.noLpn$population)
group.num <- length(inf.elements)*length(pop.elements)
group.names <- character(group.num)
avg.spleen.mat <- matrix(ncol = group.num, nrow = length(spleen.ens))
iter <- 1
for (inf in inf.elements) {
  for (pop in pop.elements){
    samples <- colname_index(rownames(filter(meta.spleen.noLpn, infection==inf, population==pop),spleen.counts),spleen.counts)
    group.names[iter] <- paste0("avg.",inf,pop)
    avg.spleen.mat[,iter] <- apply(spleen.counts[,samples],1,mean)
    iter <- iter+1
  }
}
colnames(avg.spleen.mat) <- group.names
avg.spleen <- data.frame(spleen.ens, avg.spleen.mat)
colnames(avg.spleen)[1] <- "ensgene"
avg.spleen <- distinct(avg.spleen, ensgene, .keep_all = T)


#Removing low count genes for total heatmap
heat.order <- meta.spleen.noLpn
heat.order$infection <- factor(heat.order$infection, levels = c("naive","memory"))
heat.order$population <- factor(heat.order$population, levels = c("DN","SP","DP"))
heat.order <- arrange(heat.order, population, infection)

avg.spleen.10 <- filter(avg.spleen, avg.naiveDN >=10 | avg.memoryDN >=10 | avg.naiveSP >=10 | avg.memorySP >=10 | avg.naiveDP >=10 | avg.memoryDP>=10)

heat.tot10 <- heat_ds(spleen.vst.assay, avg.spleen.10$ensgene, names.allsamples)
heat.tot10 <- heat.tot10[,colname_index(rownames(heat.order),heat.tot10)]

pheatmap(heat.tot10, scale = "row", annotation_names_row = F,show_rownames = F, cluster_cols = F,
         treeheight_row = 0,treeheight_col = 0, color = redblue)

###############################################################################


#Preparing dataframe for two-group analyses
###########################################

coldata$group <- paste0(coldata$condition,coldata$population)
coldata$group <- as.factor(coldata$group)
coldata.single <- select(coldata, group)
spleen.single.dds <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = coldata.single,
                                     design = ~group)
spleen.single.dds <- DESeq(spleen.single.dds)


#Getting the DEGs
############################

naiveDPpop.deg <- deg(spleen.single.dds,"naiveDP","naiveSP")

memDPpop.deg <- deg(spleen.single.dds,"memoryDP","memorySP")

memDPinf.deg <- deg(spleen.single.dds,"memoryDP","naiveDP")


#Venn diagram of the 3 DEG lists
deg.genes <- list(meDP_inf = memDPinf.deg$ensgene, meDP_pop = memDPpop.deg$ensgene, naDP_pop = naiveDPpop.deg$ensgene)
deg.sets <- c("meDP_inf","meDP_pop","naDP_pop")
venn.diagram(deg.genes, "0_venn-test.tiff")

#UpSet Plot
source("0_upset_deg.R")
deg.upset.plot #Save it as PDF with 8.27 x 7.67 inch in landscape

#Splitting DEGs into up- and downregulated genes
naiveDPpop.deg.up <- filter(naiveDPpop.deg, log2FoldChange > 0)
naiveDPpop.deg.down <- filter(naiveDPpop.deg, log2FoldChange < 0)

memDPpop.deg.up <- filter(memDPpop.deg, log2FoldChange > 0)
memDPpop.deg.down <- filter(memDPpop.deg, log2FoldChange < 0)

memDPinf.deg.up <- filter(memDPinf.deg, log2FoldChange > 0)
memDPinf.deg.down <- filter(memDPinf.deg, log2FoldChange < 0)

deg.genes.updown <- list(naDP_pop_up = naiveDPpop.deg.up$ensgene, naDP_pop_down = naiveDPpop.deg.down$ensgene,
                         meDP_pop_up = memDPpop.deg.up$ensgene, meDP_pop_down = memDPpop.deg.down$ensgene,
                         meDP_inf_up = memDPinf.deg.up$ensgene, meDP_inf_down = memDPinf.deg.down$ensgene)

#Top 20 up & down genes in each DEG list
naive.deg.up.20 <- top_l2fc(naiveDPpop.deg)
naive.deg.down.20 <- top_l2fc(naiveDPpop.deg,T)

memory.deg.up.20 <- top_l2fc(memDPpop.deg)
memory.deg.down.20 <- top_l2fc(memDPpop.deg,T)

inf.deg.up.20 <- top_l2fc(memDPinf.deg)
inf.deg.down.20 <- top_l2fc(memDPinf.deg,T)

lollipop_plot(naive.deg.up.20)
ggsave("0_naiveDEG-up20.eps", bg="transparent", width = 542, height = 532, units = "px")
DEGbar_plot(naive.deg.up.20)
ggsave("0_naiveDEG-up20.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)

lollipop_plot(naive.deg.down.20)
ggsave("0_naiveDEG-down20.pdf", bg="transparent", width = 210, height = 297, units = "mm")
DEGbar_plot(naive.deg.down.20)
ggsave("0_naiveDEG-down20.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)

lollipop_plot(memory.deg.up.20,"#3c9ab2")
ggsave("0_memoryDEG-up20.pdf", bg="transparent", width = 210, height = 297, units = "mm")
DEGbar_plot(memory.deg.up.20,"#3c9ab2")
ggsave("0_memoryDEG-up20.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)

lollipop_plot(memory.deg.down.20,"#3c9ab2")
ggsave("0_memoryDEG-down20.pdf", bg="transparent", width = 210, height = 297, units = "mm")
DEGbar_plot(memory.deg.down.20,"#3c9ab2")
ggsave("0_memoryDEG-down20.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)

lollipop_plot(inf.deg.up.20,"#fd6467")
ggsave("0_infDEG-up20.pdf", bg="transparent", width = 210, height = 297, units = "mm")
DEGbar_plot(inf.deg.up.20,"#fd6467")
ggsave("0_infDEG-up20.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)

lollipop_plot(inf.deg.down.20,"#fd6467")
ggsave("0_infDEG-down20.pdf", bg="transparent", width = 210, height = 297, units = "mm")
DEGbar_plot(inf.deg.down.20,"#fd6467")
ggsave("0_infDEG-down20.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)


#UpSet plot with up- and downregulated genes
#order of colors from highest setsize to lowest
upset(fromList(deg.genes.updown),order.by = "freq", keep.order = T, point.size = 4, line.size = 1,
      sets.bar.color = c("#e1af00","#3c9ab2","#e1af00","#3c9ab2","#fd6467"))
deg.updown.sets <- names(deg.genes.updown)[c(6:1)]
source("0_upset_deg_updown.R")
deg.updown.upset.plot #Save it as PDF with A4 in landscape
set.summary <- All_Freqs

#Retrieving Genes from the different intersects
deg.genes.vector <- unique(c(naiveDPpop.deg.up$ensgene,naiveDPpop.deg.down$ensgene,
                             memDPpop.deg.up$ensgene, memDPpop.deg.down$ensgene,
                             memDPinf.deg.up$ensgene, memDPinf.deg.down$ensgene))
deg.genes.colnames <- c("nu", "nd", "mu", "md", "iu", "id")
deg.genes.df <- as.data.frame(matrix(0, nrow = length(deg.genes.vector),ncol=length(deg.genes.colnames)))
colnames(deg.genes.df) <- deg.genes.colnames
deg.genes.df$ensgene <- deg.genes.vector
iter <- 0
for (ls in deg.genes.updown) {
  iter <- iter+1
  deg.genes.df[deg.genes.df$ensgene %in% ls,iter] <- 1
}
rownames(deg.genes.df) <- deg.genes.df$ensgene
deg.genes.df <- select(deg.genes.df, -ensgene)

set.combination <- set.summary[,1:6]
colnames(set.combination) <- deg.genes.colnames[6:1]
set.nr <- nrow(set.combination)
intersect.names <- character(set.nr)
set.genes <- vector(mode="list",length=0)
for (r in 1:set.nr) {
  set.sub <- set.combination[r,]
  set.cols <- colnames(select_if(set.sub, function(col) col==1))
  intersect.names[r] <- paste(set.cols, collapse="_")
  set.genes[[intersect.names[r]]] <- return_set(deg.genes.df,set.cols)
}

#Coloring Sets according to activation / migration

freq_acttax <- All_Freqs
freq_acttax$activation <- 0
freq_acttax$migration <- 0
freq_acttax$activation_migration <- 0
freq_acttax$transcription_factor <- 0
freq_acttax$tk_reg <- 0
for (setgene in set.genes) {
  set.transcription_factor <- intersect(setgene, ens.go.tkf$ensgene)
  freq_acttax[freq_acttax$freq==length(setgene),colname_index("transcription_factor",freq_acttax)] <- length(set.transcription_factor)
  set.tk_reg <- intersect(setgene, ens.go.tkreg$ensgene)
  freq_acttax[freq_acttax$freq==length(setgene),colname_index("tk_reg",freq_acttax)] <- length(set.tk_reg)
  set.acctax <- intersect(setgene, setdiff(ens.go.acttax$ensgene, set.transcription_factor))
  freq_acttax[freq_acttax$freq==length(setgene),colname_index("activation_migration",freq_acttax)] <- length(set.acctax)
  set.act <- setdiff(setdiff(intersect(setgene, ens.go.act$ensgene),set.acctax), set.transcription_factor)
  freq_acttax[freq_acttax$freq==length(setgene),colname_index("activation",freq_acttax)] <- length(set.act)
  set.tax <- setdiff(setdiff(intersect(setgene,ens.go.tax$ensgene),set.acctax), set.transcription_factor)
  freq_acttax[freq_acttax$freq==length(setgene),colname_index("migration",freq_acttax)] <- length(set.tax)
}
freq_acttax <- mutate(freq_acttax, other=freq-activation-migration-activation_migration-transcription_factor-tk_reg)
freq_acttax_gather <- gather(freq_acttax, "designation", "count", 10:15)
freq_acttax_gather$designation <- factor(freq_acttax_gather$designation,
                                         levels=c("other","tk_reg","transcription_factor","migration","activation_migration","activation"))

source("0_upset_deg_updown_acttax.R")
cairo_ps("0_upset-acttax.eps")
deg.updown.upset.acttax.plot
dev.off()


#Combine the 3 DEG lists and make heatmap
deg.all <- unique(c(naiveDPpop.deg$ensgene, memDPpop.deg$ensgene, memDPinf.deg$ensgene))

heat.deg_all <- heat_ds(spleen.vst.assay, deg.all, names.allsamples)
heat.deg_all <- heat.deg_all[,colname_index(rownames(heat.order),heat.deg_all)]
pheatmap(heat.deg_all, fontsize_row = 8, scale = "row", cluster_cols = F, show_rownames = F,
         treeheight_row = 0,treeheight_col = 0, color = redblue, border_color = NA)


pheatmap(heat.deg_all, fontsize_row = 5, scale = "row", cluster_cols = F,
         treeheight_row = 0,treeheight_col = 0, color = redblue, border_color = NA) #Save PDF as 8.27x24 inches


#Balloon Plots for combined DEG list
deg.all.tkf <- filter(ens.go.tkf, ensgene %in% deg.all)
deg.all.tkf.desig <- designate(deg.all.tkf$ensgene, ens.go.mmb$ensgene,"membrane",setdiff(ens.all.go$ensgene,ens.go.mmb$ensgene),"non-membrane")
deg.all.tkf.prep <- balloon_prep(deg.all.tkf.desig)
ggplot(deg.all.tkf.prep, aes(x=gene, y=population)) +
  geom_point(aes(size=avg.count, shape=designation), color="#941414")+
  scale_shape_manual(values=c(15))+
  scale_y_discrete(limits=rev)+
  scale_size_area( max_size = 10, breaks = c(0,1000,2000,3000,4000,5000,6000)) +
  ylab("gene")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent", colour=NA),
        plot.background = element_rect(fill = "transparent", colour=NA),
        legend.background = element_rect(fill = "transparent", colour=NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(), axis.line = element_line(size=0.5),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_blank(),
        legend.title = element_text(size=10, face="bold", colour = "black"),
        legend.text = element_text(size=10, face="bold", colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, face="italic"))
ggsave("0_balloon_deg-all_tkf.eps", bg="transparent", width = 225, height = 80, unit="mm", device = cairo_ps)


deg.all.tkreg <- filter(ens.go.tkreg, ensgene %in% deg.all)
deg.all.tkreg.desig <- designate(deg.all.tkreg$ensgene, ens.go.mmb$ensgene,"membrane",setdiff(ens.all.go$ensgene,ens.go.mmb$ensgene),"non-membrane")
balloon_plot(balloon_prep(deg.all.tkreg.desig), "#5FB5E5")
ggsave("0_balloon_deg-all_tkreg.eps", bg="transparent", width = 125, height = 80, unit="mm", device = cairo_ps)

deg.all.no_tkf <- setdiff(deg.all, deg.all.tkf$ensgene)
deg.all.no_tkf <- setdiff(deg.all.no_tkf, deg.all.tkreg$ensgene)
deg.all.no_tkf.desig <- designate(deg.all.no_tkf, ens.go.act$ensgene, "activation",ens.go.tax$ensgene, "migration")

deg.all.act <- filter(deg.all.no_tkf.desig, designation=="activation")
deg.all.act.desig <- designate(deg.all.act$ensgene, ens.go.mmb$ensgene,"membrane",setdiff(ens.all.go$ensgene,ens.go.mmb$ensgene),"non-membrane")
balloon_plot(balloon_prep(deg.all.act.desig), "#61bc78")
ggsave("0_balloon_deg-all_act.eps", bg="transparent", width = 450, height = 80, unit="mm", device = cairo_ps)

deg.all.acttax <- filter(deg.all.no_tkf.desig, designation=="activation & migration")
deg.all.acttax.desig <- designate(deg.all.acttax$ensgene, ens.go.mmb$ensgene,"membrane",setdiff(ens.all.go$ensgene,ens.go.mmb$ensgene),"non-membrane")
balloon_plot(balloon_prep(deg.all.acttax.desig), "#6e4a99")
ggsave("0_balloon_deg-all_acttax.eps", bg="transparent", width = 250, height = 80, unit="mm", device = cairo_ps)

deg.all.tax <- filter(deg.all.no_tkf.desig, designation=="migration")
deg.all.tax.desig <- designate(deg.all.tax$ensgene, ens.go.mmb$ensgene,"membrane",setdiff(ens.all.go$ensgene,ens.go.mmb$ensgene),"non-membrane")
balloon_plot(balloon_prep(deg.all.tax.desig), "#f27c48")
ggsave("0_balloon_deg-all_tax.eps", bg="transparent", width = 225, height = 80, unit="mm", device = cairo_ps)

deg.all.none <- filter(deg.all.no_tkf.desig, designation=="other")
deg.all.none.desig <- designate(deg.all.none$ensgene, ens.go.mmb$ensgene,"membrane",setdiff(ens.all.go$ensgene,ens.go.mmb$ensgene),"non-membrane")
balloon_plot(balloon_prep(deg.all.none.desig))
ggsave("0_balloon_deg-all_none.eps", bg="transparent", width = 1500, height = 80, unit="mm", device = cairo_ps, limitsize = F)


deg.all.ds$process <- "other"
deg.all.ds[deg.all.ds$ensgene %in% deg.all.tkreg$ensgene,colname_index("process",deg.all.ds)] <- "regulator of TF"
deg.all.ds[deg.all.ds$ensgene %in% deg.all.tkf$ensgene,colname_index("process",deg.all.ds)] <- "transcription factor"
deg.all.ds[deg.all.ds$ensgene %in% deg.all.act$ensgene,colname_index("process",deg.all.ds)] <- "activation"
deg.all.ds[deg.all.ds$ensgene %in% deg.all.acttax$ensgene,colname_index("process",deg.all.ds)] <- "activation & migration"
deg.all.ds[deg.all.ds$ensgene %in% deg.all.tax$ensgene,colname_index("process",deg.all.ds)] <- "migration"
deg.all.ds.desig <- designate(deg.all.ds$ensgene, ens.go.mmb$ensgene,"membrane",setdiff(ens.all.go$ensgene,ens.go.mmb$ensgene),"non-membrane")
deg.all.ds <- left_join(deg.all.ds,deg.all.ds.desig)


#Curated list of genes
#Coinhibitory
cois <- c("ENSMUSG00000026011","ENSMUSG00000026285","ENSMUSG00000071552","ENSMUSG00000038304","ENSMUSG00000034028","ENSMUSG00000015314","ENSMUSG00000030149")
deg.cois <- filter(deg.all.ds, deg.all.ds$ensgene %in% cois)
deg.cois$color <- "#61bc78"
deg.cois[deg.cois$process=="activation & migration",colname_index("color",deg.cois)] <- "#6e4a99"
#Cytokine Receptor
cykr <- c('ENSMUSG00000026770', 'ENSMUSG00000000791', 'ENSMUSG00000018341', 'ENSMUSG00000032089', 'ENSMUSG00000026068', 'ENSMUSG00000026070', 'ENSMUSG00000020009', 'ENSMUSG00000068227', 'ENSMUSG00000021756', 'ENSMUSG00000027947')
deg.cykr <- filter(deg.all.ds, deg.all.ds$ensgene %in% cykr)
deg.cykr$color <- "#61bc78"
deg.cykr[deg.cykr$process=="activation & migration",colname_index("color",deg.cykr)] <- "#6e4a99"
#Migration
mig <- c('ENSMUSG00000035042', 'ENSMUSG00000050335', 'ENSMUSG00000047880', 'ENSMUSG00000079227', 'ENSMUSG00000029530', 'ENSMUSG00000042284', 'ENSMUSG00000037405', 'ENSMUSG00000037944', 'ENSMUSG00000020689', 'ENSMUSG00000027009', 'ENSMUSG00000049103', 'ENSMUSG00000048521', 'ENSMUSG00000026581', 'ENSMUSG00000045092', 'ENSMUSG00000001281', 'ENSMUSG00000025809')
deg.mig <- filter(deg.all.ds, deg.all.ds$ensgene %in% mig)
deg.mig$color <- "#f27c48"
deg.mig[deg.mig$process=="activation & migration",colname_index("color",deg.mig)] <- "#6e4a99"

#Balloon plot with curated genes
curate_balloon(curate_balloon_prep(deg.cois), c("#61bc78","#6e4a99"))
ggsave("0_curate_balloon_cois.eps", bg="transparent", width = 125, height = 80, unit="mm", device = cairo_ps, limitsize = F)

curate_balloon(curate_balloon_prep(deg.cykr), c("#61bc78","#6e4a99"))
ggsave("0_curate_balloon_cykr.eps", bg="transparent", width = 150, height = 80, unit="mm", device = cairo_ps, limitsize = F)

curate_balloon(curate_balloon_prep(deg.mig), c("#6e4a99","#f27c48"))
ggsave("0_curate_balloon_mig.eps", bg="transparent", width = 240, height = 80, unit="mm", device = cairo_ps, limitsize = F)

population <- c("ctrl-DN","mem-DN","ctrl-SP","mem-SP","ctrl-DP","mem-DP")
population <- factor(population, levels=c("ctrl-DN","mem-DN","ctrl-SP","mem-SP","ctrl-DP","mem-DP"))

balloon_legend <- data.frame(gene=c(rep("leg",6)),population,avg.count=c(500,1000,2000,3000,4000,5000))
ggplot(balloon_legend, aes(x=gene, y=population)) +
  geom_point(aes(size=avg.count/500))+
  scale_size_identity(guide="legend") +
  scale_y_discrete(limits=rev)+
  ylab("gene")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent", colour=NA),
        plot.background = element_rect(fill = "transparent", colour=NA),
        legend.background = element_rect(fill = "transparent", colour=NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(), axis.line = element_line(size=0.5),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_blank(),
        legend.title = element_text(size=10, face="bold", colour = "black"),
        legend.text = element_text(size=10, face="bold", colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, face="italic"))
ggsave("0_curate_balloon_legend.eps", bg="transparent", width = 100, height = 80, unit="mm", device = cairo_ps, limitsize = F)


#Genes to highlight in volcano plot
curate_list <- c(deg.cois$ensgene,deg.cykr$ensgene,deg.mig$ensgene)
highlight_list <- c("ENSMUSG00000027009","ENSMUSG00000045092","ENSMUSG00000001281",
                    "ENSMUSG00000025809","ENSMUSG00000018341","ENSMUSG00000026070",
                    "ENSMUSG00000034028","ENSMUSG00000071552","ENSMUSG00000026285",
                    "ENSMUSG00000026011","ENSMUSG00000048521","ENSMUSG00000049103",
                    "ENSMUSG00000030149","ENSMUSG00000000791")


#Volcano plot of DP vs DP
ggplot(memDPinf.deg, aes(x=log2FoldChange, y=-log10(padj))) +
  xlim(-5,7.5)+
  scale_y_continuous(breaks = seq(0,160,20))+
  geom_point(size =2, color = "#FD6467")+
  geom_text(data=filter(filter(memDPinf.deg, ensgene %in% curate_list), !ensgene %in% highlight_list),
            aes(label=gene, y=-log10(padj), x=log2FoldChange), nudge_x = 0.1, nudge_y = 2, size=3, hjust="left")+
  geom_text(data=filter(filter(memDPinf.deg, ensgene %in% curate_list), ensgene %in% highlight_list),
            aes(label=gene, y=-log10(padj), x=log2FoldChange), nudge_x = 0.1, nudge_y = 2, size=3.5, hjust="left",fontface="bold")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent", colour=NA),
        plot.background = element_rect(fill = "transparent", colour=NA),
        legend.background = element_rect(fill = "transparent", colour=NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.position = "none",
        panel.border = element_blank(), axis.line = element_line(size=0.5),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"),
        legend.title = element_text(size=10, colour = "black"),
        legend.text = element_text(size=10, colour = "black"),
        aspect.ratio = 1,
        axis.ticks = element_blank())
ggsave("0_volcano_DP-DP.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)

#Volcano plot of mem DP vs SP
ggplot(memDPpop.deg, aes(x=log2FoldChange, y=-log10(padj))) +
  scale_y_continuous(breaks = seq(0,160,20))+
  scale_x_continuous(breaks=seq(-7.5,8.5,2.5), limits = c(-7.5,8.5))+
  geom_point(size =2, color = "#3C9AB2")+
  geom_text(data=filter(filter(memDPpop.deg, ensgene %in% curate_list), !ensgene %in% highlight_list),
            aes(label=gene, y=-log10(padj), x=log2FoldChange), nudge_x = 0.1, nudge_y = 2, size=3, hjust="left")+
  geom_text(data=filter(filter(memDPpop.deg, ensgene %in% curate_list), ensgene %in% highlight_list),
            aes(label=gene, y=-log10(padj), x=log2FoldChange), nudge_x = 0.1, nudge_y = 2, size=3.5, hjust="left",fontface="bold")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent", colour=NA),
        plot.background = element_rect(fill = "transparent", colour=NA),
        legend.background = element_rect(fill = "transparent", colour=NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.position = "none",
        panel.border = element_blank(), axis.line = element_line(size=0.5),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"),
        legend.title = element_text(size=10, colour = "black"),
        legend.text = element_text(size=10, colour = "black"),
        aspect.ratio = 1,
        axis.ticks = element_blank())
ggsave("0_volcano_memDP-SP.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)



#Volcano plot of naive DP vs SP
ggplot(naiveDPpop.deg, aes(x=log2FoldChange, y=-log10(padj))) +
  scale_y_continuous(breaks = seq(0,200,20))+
  scale_x_continuous(breaks=seq(-7.5,11,2.5), limits = c(-7.5,11))+
  geom_point(size =2, color = "#e1af00")+
  geom_text(data=filter(filter(naiveDPpop.deg, ensgene %in% curate_list), !ensgene %in% highlight_list),
            aes(label=gene, y=-log10(padj), x=log2FoldChange), nudge_x = 0.1, nudge_y = 2, size=3, hjust="left")+
  geom_text(data=filter(filter(naiveDPpop.deg, ensgene %in% curate_list), ensgene %in% highlight_list),
            aes(label=gene, y=-log10(padj), x=log2FoldChange), nudge_x = 0.1, nudge_y = 2, size=3.5, hjust="left",fontface="bold")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent", colour=NA),
        plot.background = element_rect(fill = "transparent", colour=NA),
        legend.background = element_rect(fill = "transparent", colour=NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.position = "none",
        panel.border = element_blank(), axis.line = element_line(size=0.5),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"),
        legend.title = element_text(size=10, colour = "black"),
        legend.text = element_text(size=10, colour = "black"),
        aspect.ratio = 1,
        axis.ticks = element_blank())
ggsave("0_volcano_ctrlDP-SP.eps", bg="transparent", width = 150, height = 150, units = "mm", device = cairo_ps)
