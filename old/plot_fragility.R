#plots fragility curves
library(ggplot2)
library(reshape)

setwd("/home/mhaas/PhD/Routines/fragility/")

data <- read.csv("dg4.csv")

Intensity <- data$iml
bt_1and4 <- data$bt1and4
bt_2 <- data$bt2
bt_3and6 <- data$bt3and6
bt_5 <- data$bt5

df <- data.frame(Intensity,bt_1and4,bt_2,bt_3and6,bt_5)
 
#plot data
df <- melt(df ,  id = 'Intensity', variable_name = 'Buildingtypes')
# plot on same grid, each series colored differently -- 
# good if the series have same scale
theme_set(theme_gray(base_size = 14))
p <- ggplot(df, aes(Intensity,value))
p + geom_line(aes(colour = Buildingtypes)) + ylab("Probability of exceedance")
 
ggsave("fragility_curves.png")