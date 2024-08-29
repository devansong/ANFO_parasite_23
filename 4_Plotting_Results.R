##################################################################
##################################################################
######ANFO parasite data analysis for Memphis Zoo ################
####Code by Anne Devan-Song. Bend, OR. 2023#######################
##################################################################
##################################################################


rm(list=ls())
graphics.off()

library(dplyr)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(viridis)
library(ggpubr)
setwd("~/Dropbox/Publications_Work/ANFO_Parasites/Data") 

df <- read.csv("data_MERGED.csv")


p2main <- ggplot(df, aes(x=SVL, y=Mass, color=Site))+
  geom_jitter(alpha=0.4, size=3)+
  scale_color_manual(values=c("#d45087", "#003f5c","#ffa600"))+ 
  xlab("Snout-vent Length (mm)")+ 
  ylab("Mass (g)")+ 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(color = guide_legend(title = "Site"))
# Marginal densities along x axis
x2dens <- axis_canvas(p2main, axis = "x")+
  geom_density(data = df, aes(x=SVL, fill = Site),
               alpha = 0.5, size = 0.2)+
  ggpubr::fill_palette(palette = c("#d45087", "#003f5c","#ffa600"))
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
y2dens <- axis_canvas(p2main, axis = "y", coord_flip = TRUE)+
  geom_density(data = df, aes(x = Mass, fill = Site),
               alpha = 0.5, size = 0.2)+
  coord_flip() +
  ggpubr::fill_palette(palette = c("#d45087", "#003f5c","#ffa600"))
p21 <- insert_xaxis_grob(p2main, x2dens, grid::unit(.2, "null"), position = "top")
p22<- insert_yaxis_grob(p21, y2dens, grid::unit(.2, "null"), position = "right")
ggdraw(p22)

png("SVLvsMASS.png", units="in", width=5, height=3.5, res=300)


df <- read.csv("data_MERGED.csv")
#df <- df[, -3]

library(reshape2)

parament <- melt(df, id.vars=c("Animal.ID", "Strongyle.LPG", "Mite.LPG", "Fecal.Mass.g.", "Site", "Sex", 
                                 "SVL", "Mass", "regline", "BCI"))

parament$variable <- factor(parament$variable, 
                            levels=c("Coccidia.EPG","Strongyle.EPG", "Ascarid.EPG","Cestode.EPG"))


plot <- ggplot(parament, aes(x=variable, y=log(value), fill=Site))+
  scale_fill_manual(values = c("#d45087", "#003f5c","#ffa600")) +
  scale_x_discrete(labels=c("Coccidia", "Strongyle", "Ascarid", "Cestode"))+
  ggbeeswarm::geom_quasirandom(shape = 21, size=1, 
                               color= "black", 
                               #fill= c("blue", "green", "red"), 
                               dodge.width = .75,alpha=0.5
                               ,show.legend = T)+
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color=NA, show.legend=F) +
  theme(panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        legend.position="right")+
  #ylim(0, 110)+
  ylab("Log(Parasite Load)")+
  xlab("")+
  ggtitle("")
plot

df <- read.csv("data_MERGED.csv")

df <- df[, c("Animal.ID",
             "Strongyle.EPG",
             "Coccidia.EPG",
             "Cestode.EPG",
             "Ascarid.EPG",
             "Site")]

colnames(df)[2] <- "Strongyle"
colnames(df)[3] <- "Coccidia"
colnames(df)[4] <- "Cestode"
colnames(df)[5] <- "Ascarid"


paronly <- df[, c(2:5)]
paronly[paronly > 0] <- 1 
newdf <- cbind(df[, 1], paronly, df[, 6])
colnames(newdf)[1] <- "Animal.ID"
colnames(newdf)[6] <- "Site"
newdf <- newdf[complete.cases(newdf),]

dat <- data.frame()

i="Zoo"
j="Strongyle"
for (i in unique(newdf$Site)){
  sub<- subset(newdf, Site == i)
  temp <- data.frame()
  for (j in (c("Coccidia", "Strongyle","Ascarid", "Cestode"))){
    subset <- sub[, j]
    prev <- sum(subset, na.rm=TRUE)/(length(subset))
    parasite <- colnames(subset)
    row <- as.data.frame(cbind(i, j, prev))
    temp <- rbind(temp, row)
  }
  dat <- rbind(dat, temp)
}
colnames(dat) <- c("Site", "Parasite", "Prevalence")
dat$Site <- as.factor(dat$Site)
dat$Parasite <- as.factor(dat$Parasite)
dat$Prevalence <- as.numeric(dat$Prevalence)

dat$Parasite <- factor(dat$Parasite, 
                            levels=c("Coccidia", "Strongyle", "Ascarid","Cestode"))
heatmap <- ggplot(dat, aes(Site, Parasite)) + 
  geom_tile(aes(fill = Prevalence), colour = "white") + 
  theme_classic()+
  scale_fill_viridis("Prevalence")+
  #scale_x_continuous(breaks = c(1:7))+
  ylab("")+
  xlab("")+
  theme(axis.line = element_line(colour = "white", 
                                 size = 0, linetype = "solid")
        #, 
        #axis.text.x = element_text(angle = 45, hjust=1)
        )

heatmap

