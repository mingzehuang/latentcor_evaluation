library(microbenchmark)
library(latentcor)
# The timing for simple simulation example
set.seed(1234)
X = GenData(n = 1000, types = c("ter", "con"), rhos = runif(1, -1, 1), XP = list(c(.3, .5), NA))$X
microbenchmark(latentcor(X, types = c("ter", "con"), method = "original"), latentcor(X, types = c("ter", "con")))
# Unit: milliseconds
# min     lq     mean    median     uq     max     neval
# 5.3444 5.8301 7.033555 6.06740 6.74975 20.8878   100
# 1.5049 1.6245 2.009371 1.73805 1.99820  5.0027   100

# The timing for mtcars example
microbenchmark(latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin", "bin", "ter", "con"), method = "original"),
               latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin", "bin", "ter", "con")))
# Unit: milliseconds
#  min       lq        mean      median        uq      max    neval
#  201.9872 215.6438   225.30385 221.5364 226.58330 411.4940   100
#   71.8457  75.1681   82.42531  80.1688  84.77845 238.3793    100
               