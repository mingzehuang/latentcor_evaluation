library(png)
library(grid)
library(gridExtra)
library(plotly)
library(heatmaply)
library(latentcor)
R_pearson = cor(mtcars)
heatmap_pearson_mtcars = heatmaply(R_pearson, dendrogram = "none", main = "Pearson correlation", margins = c(20,20,40,80),
                  grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Pearson correlation:"), limits = c(-1, 1), hide_colorbar = TRUE)
R_original = estR(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin", "bin", "ter", "con"), method = "original")$R
heatmap_original_mtcars = heatmaply(R_original, dendrogram = "none", main = "Estimated latent correlation (original)", margins = c(20,20,40,80),
                                    grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Latent correlation:"), limits = c(-1, 1), hide_colorbar = TRUE)
R_approx = estR(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin", "bin", "ter", "con"))$R
heatmap_approx_mtcars = heatmaply(R_approx, dendrogram = "none", main = "Estimated latent correlation (approx)", margins = c(70,20,40,0),
                                    grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Latent correlation:"), limits = c(-1, 1))
img1 = rasterGrob(as.raster(readPNG("heatmap_pearson_mtcars.png")), interpolate = FALSE)
img2 = rasterGrob(as.raster(readPNG("heatmap_original_mtcars.png")), interpolate = FALSE)
img3 = rasterGrob(as.raster(readPNG("heatmap_approx_mtcars.png")), interpolate = FALSE)
pdf(file = "all_heatmap.pdf", width = 16, height = 5)
grid.arrange(img1, img2, img3, ncol = 3)
dev.off()
