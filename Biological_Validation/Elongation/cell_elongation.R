
##########################################################
#          CELL ELONGATION AFTER 72 H AT HL              #
##########################################################

# Necessary package.

# install.packages("ggdist")
library(ggdist)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("MetBrewer")
library(MetBrewer)
# install.packages("ggrepel")
library(ggrepel)
# install.packages("cowplot")
library(cowplot)


# Cell size measured with the software ImageJ in confocal microscopy images for
# each condition.
len.ll <- read.table(file = "cell_length_ll.txt",as.is=T,header=F)[[1]]
len.hl <- read.table(file = "cell_length_hl.txt",as.is=T,header=F)[[1]]

# Preprocessing data.
ll <- data.frame(group="LL", value=len.ll)
hl <- data.frame(group="HL", value=len.hl)
elongation.data <- as.data.frame(rbind(ll,hl))
head(elongation.data)


# Normal boxplot.
boxplot(len.ll, len.hl,ylab="Cell length",names = c("LL","HL"),col=c("blue","red"),cex.axis=1.2,cex.lab=1.5)

# Statistical Studies.
shapiro.test(c(len.ll,len.hl)) # data doesn't follow a normal distribution
wilcox.test(x = len.ll,y = len.hl, alternative = "less")  # samples are statistically different


# Raincloud representation:

ggplot(as.data.frame(elongation.data), aes(x = reorder(group,value), y = value, fill = group)) +
  scale_fill_manual(values=c(met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[2]))+
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    ) + 
  stat_boxplot(aes(),geom = 'errorbar',linetype=1, width=0.1) +
  geom_boxplot(
    width = .25, 
    outlier.shape = NA,
    lwd=1
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  coord_cartesian(xlim = c(1.2, NA), clip = "off")+
  xlab(NULL) +
  ylab("Cell length")+
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20,vjust = +2),
        axis.text = element_text(size=17,face = "bold"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank())

ggsave2("Raincloud_plot_Elongation.png",width = 6,height = 6)

