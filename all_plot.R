# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plotly)
library(htmlwidgets)
library(heatmaply)
library(foreach)
library(doFuture)
obj = c("time_org", "time_ml", "time_approx",
        "MeanAE_org", "MeanAE_ml", "MeanAE_mlbd", "MeanAE_ml_org", "MeanAE_mlbd_org",
        "MaxAE_org", "MaxAE_ml", "MaxAE_mlbd", "MaxAE_ml_org", "MaxAE_mlbd_org")
for (i in 1:length(obj)) {
  for (type in c("BC", "BB", "TC", "TB", "TT", "NC", "NB", "NT", "NN")) {
    data <- data.frame(get(paste(type, "eval_3d", sep = "_"))[ , , i])
    rownames(data) <- as.character(seq(-0.9, 0.9, by = 0.1))
    if (type == "NC" | type == "NB" | type == "NT") {
      colnames(data) <- as.character(seq(0.1, 0.8, by = 0.1))
    } else if (type == "NN") {
      colnames(data) <- as.character(seq(0.2, 0.8, by = 0.1))
    } else {
      colnames(data) <- as.character(seq(0.1, 0.9, by = 0.1))
    }
    assign(paste(type, obj[i], sep = "_"), 
         heatmaply(data, dendrogram = "none",
                   xlab = "TruncRate", ylab = "LatentR", 
                   main = paste(type, obj[i], sep = " "),
                   margins = c(60,100,40,20),
                   grid_color = "white",
                   grid_width = 0.00001,
                   titleX = TRUE,
                   hide_colorbar = FALSE,
                   branches_lwd = 0.1,
                   label_names = c("LatentR:", "TrucRate:", paste(obj[i], ":", sep = "")),
                   fontsize_row = 10, fontsize_col = 10,
                   labCol = colnames(data),
                   labRow = rownames(data),
                   heatmap_layers = theme(axis.line=element_blank())
                   )
         )
    saveWidget(get(paste(type, obj[i], sep = "_")),
             file = paste(type, "_", obj[i], ".html", sep = ""))
  }
}