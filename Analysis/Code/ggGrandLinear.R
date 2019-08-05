library(ggplot2)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)

getGrandLinearLayout = function(genome, chrorder=1:24){
  genome = as.data.frame(seqinfo(genome))
  genome$Chr = row.names(genome)
  row.names(genome) = 1:nrow(genome)
  genome = genome[chrorder,]
  genome = rbind(data.frame(seqlengths=0, isCircular=FALSE, genome=NA, Chr=NA), genome)
  genome[2:nrow(genome),"start"] = cumsum(genome[1:(nrow(genome)-1),"seqlengths"])+1
  genome$end = genome$seqlengths + genome$start - 1
  genome$mid = floor(genome$start + genome$end)/2
  genome$color = as.factor(1:nrow(genome)%%2)
  genome = genome[-c(1),-c(1:3)]
  return(genome)
}
getGrandLinearCoord = function(chr, pos, genome, chrorder=1:24){
  genome = getGrandLinearLayout(genome, chrorder)
  return(genome[match(chr, genome$Chr), "start"] + pos)
}
getGrandLinearColor = function(chr, pos, genome, chrorder=1:24){
  genome = getGrandLinearLayout(genome, chrorder)
  return(genome[match(chr, genome$Chr), "color"])
}

# ggGrandLinear = function(chr, pos, y, genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
#   data = data.frame(chr=chr, pos=getGrandLinearCoord(chr, pos, genome, chrorder), color=getGrandLinearColor(chr, pos, genome, chrorder), y=y)
#   GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
#   f=ggplot(data=data, aes(x=pos, y=y, color=color, group=chr)) + 
#    # scale_y_continuous(labels=scales::comma) +
#     scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
#     theme_minimal() +
#     theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
#     scale_color_manual(values=c("#6baed6", "#08519c")) +
#     geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end)) +
#     ylab("Circle Read Count") +
#     xlab("") +
#     ylab("") +
#     guides(colour=FALSE)
#   return(f)
# }

ggGrandLinear = function(chr, pos, y, genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
  data = data.frame(chr=chr, pos=getGrandLinearCoord(chr, pos, genome, chrorder), color=getGrandLinearColor(chr, pos, genome, chrorder), y=y)
  GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
  f=ggplot(data=data, aes(x=pos, y=y, color=color, group=chr)) + 
    # scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
    theme_classic() +
    theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
    #scale_color_manual(values=c("#6baed6", "#08519c")) +
    scale_color_manual(values=c("#000000", "#000000")) +
    geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end)) +
    xlab("") +
    ylab("") +
    guides(colour=FALSE) +
    theme(text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5))
  return(f)
}

ggGrandLinearNoChrSep = function(chr, pos, y, genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
  data = data.frame(chr=chr, pos=getGrandLinearCoord(chr, pos, genome, chrorder), color=getGrandLinearColor(chr, pos, genome, chrorder), y=y)
  GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
  f=ggplot(data=data, aes(x=pos, y=y)) + 
    scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
    theme_classic() +
    theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
    geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end), size=0.1) +
    xlab("") +
    ylab("") +
    theme(text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 12, face = "bold"), 
          plot.title = element_text(face ="bold", size = 12, hjust=0.5))
  return(f)
}

ggGrandLinearCond = function(chr, pos, y, condition, genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
  data = data.frame(chr=chr, pos=getGrandLinearCoord(chr, pos, genome, chrorder), color=getGrandLinearColor(chr, pos, genome, chrorder), y=y, condition=condition)
  GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
  f=ggplot(data=data, aes(x=pos, y=y, color=condition, group=interaction(chr, condition))) + 
    # scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
    theme_minimal() +
    theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
    scale_color_manual(values=c("steelblue", "firebrick")) +
    geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end)) +
    xlab("") +
    ylab("")
  return(f)
}

ggGrandLinearGenomeWideDensity = function(chr, pos, bw=2, genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
  data = data.frame(chr=chr, pos=getGrandLinearCoord(chr, pos, genome, chrorder), color=getGrandLinearColor(chr, pos, genome, chrorder))
  GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
  f=ggplot(data=data, aes(x=pos)) + 
    # scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
    theme_minimal() +
    theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
    scale_color_manual(values=c("#6baed6", "#08519c")) +
    geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end)) +
    xlab("") +
    ylab("") +
    geom_density(bw=bw) + 
    guides(color=FALSE)
  return(f)
}

ggGrandLinearGenomeWideHistogram = function(data, chr, pos, binwidthMb=25, genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
  
  data = as.data.frame(data)
  
  data$chr=data[,chr]
  data$pos=getGrandLinearCoord(data[,chr], data[,pos], genome, chrorder)
  data$color=getGrandLinearColor(data[,chr], data[,pos], genome, chrorder)
  
  GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
  
  breaks = unlist(mapply(function(a,b) seq(a,b,by=(binwidthMb*1000000)), GrandLinearLayout$start, GrandLinearLayout$end))
  
  f=ggplot(data=data, aes(x=pos, fill=color, group=chr)) + 
    # scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
    theme_minimal() +
    theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
    scale_fill_manual(values=c("#6baed6", "#08519c")) +
    geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end)) +
    xlab("") +
    ylab("") +
    geom_histogram(breaks = breaks) + 
    guides(fill=FALSE)
  return(f)
}

ggGrandLinearIntervals = function(data, chr, pos, binwidthMb=25, genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
  
  data = as.data.frame(data)
  
  data$chr=data[,chr]
  data$pos=getGrandLinearCoord(data[,chr], data[,pos], genome, chrorder)
  data$color=getGrandLinearColor(data[,chr], data[,pos], genome, chrorder)
  
  GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
  
  breaks = unlist(mapply(function(a,b) seq(a,b,by=(binwidthMb*1000000)), GrandLinearLayout$start, GrandLinearLayout$end))
  
  f=ggplot(data=data, aes(x=pos, fill=color, group=chr)) + 
    # scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
    theme_minimal() +
    theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
    scale_fill_manual(values=c("#6baed6", "#08519c")) +
    geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end)) +
    xlab("") +
    ylab("") +
    geom_histogram(breaks = breaks) + 
    guides(fill=FALSE)
  return(f)
}

ggGrandLinearBed = function(chr, start, end, y, genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
  data = data.frame(chr=chr, start=getGrandLinearCoord(chr, start, genome, chrorder), end=getGrandLinearCoord(chr, end, genome, chrorder), color=getGrandLinearColor(chr, start, genome, chrorder), y=y)
  GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
  f=ggplot(data=data, aes(x=start, xend=end, y=y, yend=y, color=color, group=chr)) + 
    geom_segment(size=2) + 
    # scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
    theme_classic() +
    theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
    #scale_color_manual(values=c("#6baed6", "#08519c")) +
    scale_color_manual(values=c("#000000", "#000000")) +
    geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end)) +
    xlab("") +
    ylab("") +
    guides(colour=FALSE) +
    theme(text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5))
  return(f)
}

ggGrandLinearColor = function(chr, pos, y, colors=c("#6baed6", "#08519c"), genome=BSgenome.Hsapiens.UCSC.hg19, chrorder=1:24){
  data = data.frame(chr=chr, pos=getGrandLinearCoord(chr, pos, genome, chrorder), color=getGrandLinearColor(chr, pos, genome, chrorder), y=y)
  GrandLinearLayout = getGrandLinearLayout(genome, chrorder)
  f=ggplot(data=data, aes(x=pos, y=y, color=color, group=chr)) + 
    # scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks = GrandLinearLayout$mid, labels=GrandLinearLayout$Chr) +
    theme_classic() +
    theme(axis.text=element_text(angle=90), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.title = element_text(size=10, hjust = 0.5, face="bold")) +
    scale_color_manual(values=colors) +
    geom_vline(color="lightgray", xintercept=c(GrandLinearLayout$start, GrandLinearLayout$end)) +
    xlab("") +
    ylab("") +
    guides(colour=FALSE) +
    theme(text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5))
  return(f)
}


