# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plotly)
data <- data.frame(NC_eval_3d[ , , 8])
rownames(data) <- as.character(seq(-0.9, 0.9, by = 0.1))
colnames(data) <- as.character(seq(0.1, 0.8, by = 0.1))

library(heatmaply)
MeanAE_NC_mlbd_org <- heatmaply(data, 
                              dendrogram = "none",
                              xlab = "TruncRate", ylab = "LatentR", 
                              main = "Mean Absolute Error for Ternary/Continuous (mlbd_org)",
                              #scale = "column",
                              margins = c(60,100,40,20),
                              grid_color = "white",
                              grid_width = 0.00001,
                              titleX = TRUE,
                              hide_colorbar = FALSE,
                              branches_lwd = 0.1,
                              label_names = c("LatentR:", "TrucRate:", "MeanAE:"),
                              fontsize_row = 5, fontsize_col = 5,
                              labCol = colnames(data),
                              labRow = rownames(data),
                              heatmap_layers = theme(axis.line=element_blank())
)
# save the widget
library(htmlwidgets)
saveWidget(MeanAE_NC_mlbd_org, file= "MeanAE_NC_mlbd_org.html")