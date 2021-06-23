library(latentcor)

# Data generation
rhorep = rep(NA, 100); Rrep = matrix(NA, 100, 3)
for (rep in 1:100) {
  rho = runif(1, -1, 1)
  X = GenData(n = 1000, types = c("ter", "con"), rhos = rho, XP = list(c(.3, .5), NA))
  R_nc_org = estR(X, types = c("ter", "con"), method = "original")$R
  R_nc_approx = estR(X, types = c("ter", "con"), method = "approx")$R
  rhorep[rep] = rho; Rrep[rep, 1] = R_nc_org[2, 1]; Rrep[rep, 2] = R_nc_approx[2, 1]; Rrep[rep, 3] = cor(X)[2, 1]
}

# Plot ternary/continuous case estimation via original method.
R_nc_org = PlotPair(datapair = cbind(rhorep, Rrep[,1]), namepair = c("True latent correlation", "Estimated latent correlation"),
         title = "Ternary vs. continuous")
# Plot ternary/continuous case estimation via approximation method.
R_nc_approx = PlotPair(datapair = cbind(rhorep, Rrep[,2]), namepair = c("True latent correlation", "Estimated latent correlation (latentcor)"),
         title = "Ternary vs. continuous")
R_nc_pearson = PlotPair(datapair = cbind(rhorep, Rrep[,3]), namepair = c("True latent correlation", "Pearson correlation"),
         title = "Ternary vs. continuous")

cp1 = "exp"; cp2 = "cube"
for (tp1 in c("continuous", "binary", "ternary", "trunc")) {
  for (tp2 in c("continuous", "binary", "ternary", "trunc")) {
    for (md in c("original", "ml", "approx")) {
      c1 = c2 = NULL
      if (tp1 == "binary" | tp1 == "trunc") {
        c1 = 1
      } else if (tp1 == "ternary") {
        c1 = c(1, 2)
      }
      if (tp2 == "binary" | tp2 == "trunc") {
        c2 = 0
      } else if (tp2 == "ternary") {
        c2 = c(0, 1)
      }
      simdata = GenData(n=n, type1 = tp1, type2 = tp2, p1= p1, p2 = p2, rho = rho,
                        copula1 = cp1, copula2 = cp2, c1 = c1, c2 = c2)
      Sigma = simdata$Sigma; X1 = simdata$X1; X2 = simdata$X2
      assign(paste("R", cp1, cp2, tp1, tp2, md, sep = "_"),
             estR(X1 = X1, type1 = tp1, X2 = X2, type2 = tp2, method = md)$R)
      PlotPair(datapair = cbind(c(Sigma), c(get(paste("R", cp1, cp2, tp1, tp2, md, sep = "_")))),
               namepair = c("Sigma", paste("R", cp1, cp2, tp1, tp2, md, sep = "_")),
               title = "Latent correlation (True vs. Estimated)")
    }
  }
}


PlotPair = function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df = data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]))
}


pdf(file = "", width = , height = )
dev.off()