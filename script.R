## SETUP WORKING DIRECTORY, LOAD PACKAGES AND DATA

# Set working directory
setwd("/Users/au231308/Dropbox/Projects/Physiology_retina/Github/")

# Load packages
Packages                 <-c("phyloseq","splitstackshape","RColorBrewer","ggsignif","ape","cowplot","readxl","geiger","phytools","ggplot2","nlme","coda","gridExtra","grid","phangorn")
lapply(Packages, library, character.only = TRUE) 


## EVOLUTIONARY MODEL OF RETINAL CAPILLARIZATION

# Load data and tree
tree<-read.tree("./Phylogeny.nex") 
df<-as.data.frame(read_excel("./Data.xlsx", sheet = "mammals"))

# Named vector for intraretinal capillarization (irc)
irc<-setNames(nm = df$sp,object = as.character(df$irc))
irc<-irc[!is.na(irc)]

# Any dublicates in the data set?
irc[duplicated(names(irc))]

# Prune phylogeny
tree.irc<-keep.tip(tree,tip = names(irc))
tree.irc<-ladderize(tree.irc)


# Fit different evolutionary models. 
fit.ER<-fitDiscrete(phy = tree.irc,dat = irc,model = "ER") # Equal rates model
fit.ARD<-fitDiscrete(phy = tree.irc,dat = irc,model = "ARD") # All rates different model

fit.ER
fit.ARD

# Q: Which model is better, using Akaike's weight?
aic.w(c(fit.ER$opt$aic,fit.ARD$opt$aic))
# A: fairly similar weights.

# Q: Which model is better, using likelihood ratio test?
k0<-length(fit.ER$rates)
k1<-length(fit.ARD$rates)
LR<--2*(logLik(fit.ER)-logLik(fit.ARD))
P_chisq<-pchisq(LR,df=k1-k0,lower.tail=FALSE)
P_chisq
# A: Model fits are not different. Use the simplest model (ER-model)


# Stochastic character mapping
simmap<-make.simmap(tree = tree.irc, x = irc,model = "ER",nsim = 10000)


# Summarize stochastic character maps
d.simmap<-describe.simmap(simmap)


## FIND TRANSITIONS

# Extract the coordinates of the phylogeny
coor<-as.data.frame(tree_layout(tree.irc)$edgeDT)

# Generate a matrix for storing phylogeny coordinates and baysian posterior probabilities
extract <- data.frame(
  node1 = rep(NA,1000),
  node2 = NA,
  bpp1 = NA,
  bpp2 = NA,
  time1 = NA,
  time2 = NA,
  time50 = NA,
  irc = NA,
  edge = NA
)

# Generate a character matrix (cmat).
# Col 1 is the probability of IRC is absent.
# Col 2 is the probability of IRC is present.

cmat<-data.frame(
  "no"=rep(0,length(irc)),
  "yes"=0)
colnames(cmat)<-c("0","1")

for (i in 1:length(tree.irc$tip.label)){
  cmat[i,which(colnames(cmat)==unname(irc[tree.irc$tip.label[i]]))]<-1
}

# Combined character matrix with inferred ancestral states 
ace<-rbind(cmat,d.simmap$ace)


# For each edge within the phylogeny, do the following
for (i in 1:dim(tree.irc$edge)[1]){
  extract$node1[i]<-node1<-tree.irc$edge[i,1]     # extract node 1 of the edge
  extract$node2[i]<-node2<-tree.irc$edge[i,2]     # extract node 2 of the edge
  node1n<-which(row.names(ace)==node1)
  node2n<-which(row.names(ace)==node2)
  
  trans<-ace[node1n,]<0.5&ace[node2n,]>0.5        # is there a transition between IRC being absent and IRC being present. 
  
  # if there is a transition between IRC being absent and present, do the following
  if(any(trans)){
    extract$irc[i]<-colnames(trans)[which(trans == T)]
    extract$bpp1[i]<-ace[node1,which(trans == T)]
    extract$bpp2[i]<-ace[node2,which(trans == T)]
    extract$time1[i]<-max(coor$xright)-coor[which(coor$V1==node1),]$xleft[1]
    extract$time2[i]<-max(coor$xright)-coor[which(coor$V2==node2),]$xright[1]
    extract$time50[i]<-(0.5-lm(c(extract$bpp1[i],extract$bpp2[i])~c(extract$time1[i],extract$time2[i]))$coefficients[1])/lm(c(extract$bpp1[i],extract$bpp2[i])~c(extract$time1[i],extract$time2[i]))$coefficients[2]
    extract$edge[i] <- which(tree.irc$edge[,1]==node1&tree.irc$edge[,2]==node2)
  }
}

extract.na<-na.omit(extract)
extract.na$num<-1:length(extract.na$node1)
extract.na


## PLOT SUMMARY OF STOCHASTIC CHARACTER MAPPING

# Plot as density map without tip-labels for in-text figure
map<-densityMap(trees = simmap,plot = F,legend = F)
map$cols[]<-colorRampPalette(c("black", "#D81B60"))(length(map$cols))
pdf(file = "./Mammal_densitymap.pdf",width = 2.5,height = 5,useDingbats = F)
par(mar = c(2,0,0,0))
plot.densityMap(map,legend = F,fsize =0.0001,lwd = 1)
edgelabels(edge=extract.na$edge[extract.na$irc=="0"],cex=2,bg = "black",pch = 21)
edgelabels(edge=extract.na$edge[extract.na$irc=="1"],cex=2,bg = "#D81B60",pch = 21)
dev.off()

# Plot as density map with tip-labels for supplementary figure 
pdf(file = "./Mammal_densitymap_labels.pdf",width = 2.5,height = 5,useDingbats = F)
par(mar = c(2,0,0,0))
plot.densityMap(map,legend = F,fsize=0.1,lwd = 1)
axisPhylo(lwd = 1,cex.axis=0.5,padj=-2.5)
title(xlab="Time (MYA)",font.lab=1,cex.lab=0.5,line = 1,adj = 0.35)
dev.off()



# QUESTIONS 

# Q: How many species in the data set
length(irc)


# Q: Phylogenetic signal in irc?
fitDiscrete(phy = tree.irc,dat = irc,model = "ER",transform = "lambda")$opt$lambda
# A: Yes!


# When did the MRCA of mammals live?
max(phyloseq::tree_layout(tree.irc)$edgeDT$xright)


# What is the probability that the most recent common ancestor of mammals was vascular?
d.simmap$ace[1,]






## RETINAL THICKNESS MAMMALS

# is retinal thickness higher in species with vascular retinas in mammals?
df.mammals<- data.frame(
  sp = df$sp,
  irc = df$irc, 
  thick = df$thick)

df.mammals<-na.omit(df.mammals)

df.mammals<-df.mammals[df.mammals$sp!="Pteropus_poliocephalus",] # Remove flying fox. See figure legend in publication. 

phylANOVA(
  tree = keep.tip(tree = tree,tip = df.mammals$sp),
  x = setNames(df.mammals$irc,df.mammals$sp), 
  y = setNames(df.mammals$thick,df.mammals$sp),
  nsim = 10000)




# Mammals Rug
ggplot(
  data = df.mammals,
  mapping = aes(x = thick,y=0,colour = as.character(irc))
)+
  geom_point(pch = 21,size = 2)+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text =element_blank(),
        axis.title =element_blank(),
        panel.spacing = unit(c(0, 0, 0, 0), "cm"))+
  
  scale_color_manual(
    values = c("black","#D81B60"))+
  scale_x_continuous(expand = c(0, 0),breaks = c(0,350,700),limits = c(0,700))
ggsave("./2_rug_mammals.pdf",height = .3,width = 3)

# Mammal distribution
ggplot(data = df.mammals,
       mapping = aes(x = thick,fill = as.character(irc)))+
  scale_fill_manual(
    values = c("black","#D81B60"))+
  theme_minimal()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(0,0,0,0), "mm"),
    axis.text = element_blank(),
    axis.title =element_blank()
  )+
  geom_density(alpha=.7)+
  scale_x_continuous(expand = c(0, 0),breaks = c(0,350,700),limits = c(0,700))
ggsave("./2_range_mammals.pdf",height = .3,width = 3)





## RETINAL THICKNESS FISHES

# is retinal thickness higher in species with choroid rete mirabile?
df.fish<-as.data.frame(read_excel("./Data.xlsx", sheet = "fishes"))
df.fish

phylANOVA(
  tree = keep.tip(tree = tree,tip = df.fish$sp),
  x = setNames(df.fish$crm,df.fish$sp), 
  y = setNames(df.fish$thick,df.fish$sp),
  nsim = 10000)




# Fish Rug
ggplot(
  data = df.fish,
  mapping = aes(x = thick,y=0,colour = as.character(crm))
)+
  geom_point(pch = 21,size = 2)+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text =element_blank(),
        axis.title =element_blank(),
        panel.spacing = unit(c(0, 0, 0, 0), "cm"))+
  
  scale_color_manual(
    values = c("black","#1E89E5"))+
  scale_x_continuous(expand = c(0, 0),breaks = c(0,350,700),limits = c(0,700))
ggsave("./1_rug_rff.pdf",height = .3,width = 3)

# fish distribution
ggplot(data = df.fish,
       mapping = aes(x = thick,fill = as.character(crm)))+
  scale_fill_manual(
    values = c("black","#1E89E5"))+
  theme_minimal()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(0,0,0,0), "mm"),
    axis.text = element_blank(),
    axis.title =element_blank()
  )+
  geom_density(alpha=.7)+
  scale_x_continuous(expand = c(0, 0),breaks = c(0,350,700),limits = c(0,700))
ggsave("./1_range_rff.pdf",height = .3,width = 3)







## RETINAL THICKNESS BIRDS
df.birds<-as.data.frame(read_excel("./Data.xlsx", sheet = "birds"))

# Birds Rug
ggplot(
  data = df.birds,
  mapping = aes(x = thick,y=0)
)+
  geom_point(pch = 21,size = 2,colour = "#004D40")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text =element_blank(),
        axis.title =element_blank(),
        panel.spacing = unit(c(0, 0, 0, 0), "cm"))+
  
  #scale_color_manual(
  #  values = c("black","#004D40"))+
  scale_x_continuous(expand = c(0, 0),breaks = c(0,350,700),limits = c(0,700))
ggsave("./3_rug_bird.pdf",height = .3,width = 3)

# Bird distribution
ggplot(data = df.birds,
       mapping = aes(x = thick))+
  geom_density(alpha=.7,fill="#004D40")+
  theme_minimal()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(0,0,0,0), "mm"),
    axis.text = element_blank(),
    axis.title =element_blank()
  )+
  
  scale_x_continuous(expand = c(0, 0),breaks = c(0,350,700),limits = c(0,700))
ggsave("./3_range_bird.pdf",height = .3,width = 3)





