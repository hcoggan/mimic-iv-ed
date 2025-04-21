f = read.csv('packages.txt', header=FALSE, stringsAsFactors = FALSE)
z = install.packages(f[,1], repos='https://cran.rstudio.com') 