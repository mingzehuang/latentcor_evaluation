library(png)
library(grid)
library(gridExtra)
library(plotly)
library(heatmaply)
library(latentcor)

mtcars_latentcor = estR(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin", "bin", "ter", "con"), showplot = TRUE)

heatmap_latentcor_mtcars = mtcars_latentcor$plotR

mtcars_pearson = cor(mtcars)
heatmap_pearson_mtcars = heatmaply(mtcars_pearson, dendrogram = "none", main = "Pearson correlation", margins = c(80,80,80,80),
                                   grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Pearson correlation:"), limits = c(-1, 1))

R_latentcor_mtcars = mtcars_latentcor$R

heatmap_diff = heatmaply(R_latentcor_mtcars - mtcars_pearson, dendrogram = "none", colors = heat.colors(100), main = "Difference (estimated latent correlation - pearson correlation)", margins = c(80,80,80,80),
                                    grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Difference:"))

img1 = rasterGrob(as.raster(readPNG("heatmap_latentcor_mtcars.png")), interpolate = FALSE)
img2 = rasterGrob(as.raster(readPNG("heatmap_pearson_mtcars.png")), interpolate = FALSE)
img3 = rasterGrob(as.raster(readPNG("heatmap_diff_mtcars.png")), interpolate = FALSE)
pdf(file = "all_heatmap.pdf", width = 16, height = 5)
grid.arrange(img1, img2, img3, ncol = 3)
dev.off()
