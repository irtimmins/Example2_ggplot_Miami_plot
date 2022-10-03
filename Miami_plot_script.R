######################################################################

library(readr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(cowplot)

# read in summary statistic gwas data.

gwas.df <- read_csv("gwas_sum_stats_miami_high.csv",
                    col_types = cols(
                      SNP = col_character(),
                      CHR = col_double(),
                      BP = col_double(),
                      P2 = col_double(),
                      bp2 = col_double(),
                      chr2 = col_double(),
                      label = col_double()
                    ))

gwas.df2 <- read_csv("gwas_sum_stats_miami_low.csv",
                     col_types = cols(
                       SNP = col_character(),
                       CHR = col_double(),
                       BP = col_double(),
                       P2 = col_double(),
                       bp2 = col_double(),
                       chr2 = col_double(),
                       label = col_double()
                     ))

# reduce computational burden by sampling from dense regions of the plot,
# while keeping top hits intact.

sample1 <- sample(1:length(gwas.df$SNP), 100000)
sample1.logic <- is.element(1:length(gwas.df$SNP), sample1)
top.hits <- gwas.df$P2 < 10^(-2)

sample2 <- sample(1:length(gwas.df2$SNP), 100000)
sample2.logic <- is.element(1:length(gwas.df2$SNP), sample2)
top.hits2 <- abs(gwas.df2$P2) < 10^(-2)

gwas.df <- gwas.df[top.hits | sample1.logic,]
gwas.df2 <- gwas.df2[top.hits2 | sample2.logic,]



# identify position of chromosomes for x axis labels.

breaks.x <- rep(0,22)
names.x <- 1:22
  
for(i in 1:22){
    breaks.x[i] <- mean(gwas.df$bp2[gwas.df$CHR == i]) 
  } 

margin.set <- c(0,0.2,0,0.2)

plot.high <- ggplot(data = gwas.df, aes(bp2, y = P2, color = chr2)) +
    theme_classic()+
  theme(plot.title = element_text(color="black", size=7, hjust=0.5, vjust = -52), 
	axis.text=element_text(size=6),
       axis.title.y = element_text(size=6),
  	legend.position = "none",
	plot.margin = unit(margin.set, "cm"))+
  geom_point(size = 0.55, shape = 16)+
  scale_y_continuous(substitute(paste("-log10(", italic('p'), "value)")),  expand = c(0, 0), limits = c(-0.2, max(gwas.df$P2)*1.05),breaks = c(0,5,10,15,20,25), labels = c(0,5,10,15,20,25) )+
  scale_x_continuous(name =NULL,expand = c(0, 30000000), limits = c(0, max(gwas.df$bp2)*1.001),breaks=breaks.x, labels = names.x)+
  scale_color_gradient(low="grey80", high="grey88")+
  geom_hline(yintercept=-log10(5*10^(-8)), linetype="dashed", color = "red3", size=0.5)+
  geom_point(data = subset(gwas.df, label == 1), size = 0.55, shape = 21,fill = "blue", color = "grey40", stroke = 0.3 )
 

plot.low <- ggplot(data = gwas.df2, aes(bp2, y = P2, color = chr2)) +
  theme_classic() +
  theme(plot.title = element_text(color="black", size=7, hjust=0.5, vjust = -52),
        axis.text=element_text(size=6),
	 axis.title.y = element_text(size=6), 
	 legend.position = "none",
	 plot.margin = unit(margin.set, "cm"))+
  geom_point(size = 0.55, shape = 16)+
  scale_y_continuous(name = substitute(paste("-log10(", italic('p'), "value)")), expand = c(0, 0), limits = c(-max(gwas.df$P2)*1.05, 0.2), breaks = c(0,-5,-10,-15,-20,-25), labels = c(0,5,10,15,20,25))+
  scale_x_continuous(name =NULL, expand = c(0, 30000000), limits = c(0, max(gwas.df$bp2)*1.001),position = "top", breaks=breaks.x, labels = NULL) +
  scale_color_gradient(low="grey80", high="grey88")+
  geom_hline(yintercept=log10(5*10^(-8)), linetype="dashed", color = "red3", size=0.5)+
  geom_point(data = subset(gwas.df2, label == 1), size = 0.55, shape = 21,fill = "blue", color = "grey40", stroke = 0.3 )


# combine two plots using plot_grid().

pdf(file = "Miami_plot.pdf", width = 5, height = 3.6)
plot_grid(plot.high, plot.low, align = "v", ncol = 1, rel_widths=c(0.5,0.5))
dev.off()









