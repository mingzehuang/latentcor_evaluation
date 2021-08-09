library(plotly)
library(heatmaply)
R_pearson = cor(mtcars)
heatmaply(R_pearson, dendrogram = "none", main = "Pearson Correlation", margins = c(80,80,80,80),
                  grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Pearson correlation:"))
