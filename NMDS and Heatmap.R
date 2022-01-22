#NMDS
attach(`Title of CSV`)
summary(`Title of CSV`)
str(`Title of CSV`)

install.packages("vegan")
library(vegan)

spe_log <- log1p (`Title of CSV`)  # log-transform raw species composition data
dist_bc <- vegdist (spe_log, method = 'bray') # calculate Bray-Curtis distance matrix
dist_bc_sqrt <- sqrt (dist_bc) # square-root BC distance matrix to make distances metric
clus_ward <- hclust (dist_bc_sqrt, method = 'ward.D2')  # calculate Ward's algorithm (using the correct ''method = 'ward.D2''')
plot(clus_ward)
barplot(clus_ward$height, names.arg=(nrow(`Title of CSV`)-1): 1) # show the number of cluster below each bars 

clus_ward_cut <- cutree (clus_ward, k = 5) # "cut the tree" - to which groups individual samples belong? Edit K based on barplot above; how many clusters?
write.csv(clus_ward_cut, row.names = FALSE ) #Copy to excel database

plot (clus_ward, cex = .5)  # argument cex reduced the size of the dendrogram leaf labels to make them readable
clus_in_dendro <- unique (clus_ward_cut[clus_ward$order]) # make sure to know which box is which cluster!
rect.hclust (clus_ward, k = 5, border = clus_in_dendro+1)
legend ('topright', legend = paste ('Cluster', 1:5), col = c(2,4,3,5,6), bg =c(2,4,3,5,6), pch = 0, bty = 'n')

NMDS <- metaMDS (dist_bc) #stress = around 0.1 is a good fit, errors may indicate you have bad samples, or samples were not enough for analysis
View(NMDS)
ordiplot (NMDS, type = 'n', main="")
points (NMDS, pch = 20, col = clus_ward_cut+1)
legend ('topright', legend = paste('Cluster', 1:5), pch = 20, col =c(2,4,3,5,6) , pt.bg = c(2,4,3,5,6))

  #high res plot
  tiff("NMDS.tiff", units = "in", width = 6, height = 5, res = 600)
  ordiplot (NMDS, type = 'n', main="")
  points (NMDS, pch = 20, col = clus_ward_cut+1)
  legend ('topright', legend = paste('Cluster', 1:3), pch = 20, col =c(2,4,3) , pt.bg = c(2,4,3))
  dev.off()

#Dendogram sample arrangement
write.csv(clus_ward_cut[clus_ward$order])
# Use this to arrange your samples based on cluster, so that you will be able to tell compositions of each cluster
# proceed to make pie graphs (composition), and heatmaps
  
####HEATMAP MAY ALSO BE DONE USING GRAPH PAD PRISM - much EASIER
#### Make 2 CSV files that are arranged based on cluster. One contains the counts, and the other one contains the parameters.
####CLUSTER HEATMAP
attach(`title`) #note that sample IDs should be replaced by order based on dendogram
install.packages("reshape")
install.packages ("ggplot2")
install.packages ("ggpubr")
install.packages ("ggdendro")
library(reshape)
library(ggplot2)
library(ggpubr)
library(ggdendro)

dendro.plot <- ggdendrogram(clus_ward) + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),plot.margin = unit(c(0,0.8,0,1.2),"cm"))
#heatmap for counts
data <- melt(`log`, id=c("Study","Cluster"))
View(data)
#highlight "data" to show 

all1 <- ggplot(data, aes(Study,variable, fill=value)) + geom_tile() +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(face="italic"),axis.ticks = element_blank(),plot.background = element_rect(fill="transparent",color=NA),panel.border = element_rect(color="black", fill=NA),plot.margin = unit(c(0.5,0.5,0.5,0),"cm")) +
  labs(y="",x="", title = "")  +
  scale_fill_gradient2(low="white",mid="cornsilk",high="firebrick1")
#guides(fill=guide_legend(title="Relative Abundance"))

#Sample Attributes
attach(`delta1`)
#heatmap component for parameter
data2 <- melt(`delta1`, id=c("Sample.ID","Number"))
View(data2)

levels <- ordered(c("Patient.Status","Tumor.Stage"))
data$value <- factor(data2$value, levels = c("Vagina","Cervix","Endometrium","Pregnant","Non-pregnant","Fertile","Infertile","follicular","luteal")) #change based on order of parameters

all2 <- ggplot(data2, aes(Sample.ID,ordered(variable, levels =rev(levels)), fill=value)) + geom_tile() +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill="transparent",color=NA),panel.border = element_rect(color="black", fill=NA), plot.margin = unit(c(0,0.5,0,-0.2),"cm")) + 
  labs(y="",x="", title = "")+
  scale_fill_manual(values = c("brown1","darkgoldenrod1","deeppink","mediumorchid","cornflowerblue","gold","firebrick1","darkseagreen1","darkslateblue"))

#Combine two components together
ggarrange(dendro.plot + theme(plot.margin = unit(c(0,-0.48,-1,3.53),"cm")),
          all2 + theme(plot.margin = unit(c(-0.3,0.5,-1.5,1.65),"cm")),
          all1,
          ncol=1,nrow=3,heights = c(2,0.5,5), legend ="none")
ggsave("Heatmap_Final_Sample.tiff", dpi = 600)

#Legend - print separate page
leg1 <-as_ggplot(get_legend(all1))
leg2 <- as_ggplot(get_legend(all2))
ggarrange(leg1+ theme(plot.margin = unit(c(0,-2,0,0),"cm")),
          leg2+ theme(plot.margin = unit(c(0,0,0,-2),"cm")),
          ncol = 2, nrow = 1)
ggsave("Heatmap_Final_Legend_Sample.tiff", dpi=600)

### Adapted from Renne Margaret Alcazar's code, special thanks to her!
### Prepared and annotated by Felippe Steven Louis Delos Reyes