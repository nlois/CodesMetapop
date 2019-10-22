## Correlation plots for SIBER significance test

rm(list=ls())
graphics.off()

library(corrplot)

# Data and plot test
mydata = read.csv("https://wiki.q-researchsoftware.com/images/b/b9/Ownership.csv")
cordata <- cor(mydata)
corrplot(cordata,type = "upper", diag = TRUE)

# On SIBER data
siber <- read.csv("~/Documents/Proceedings Rockhopper/Data/SIBER.csv",
                  header = TRUE, row.names = 1)
data <- as.matrix(siber)
data

corrplot(data, method = "color", order = "hclust", is.corr = FALSE,na.label = "square",na.label.col = "black", tl.col = "black")


# On nicheROVER data
niche <- read.csv("~/Documents/Proceedings Rockhopper/Data/nicheROVER.csv",
                  header = TRUE, row.names = 1)
niche_lower <- read.csv("~/Documents/Proceedings Rockhopper/Data/nicheROVER_lower.csv",
                        header = TRUE, row.names = 1)
niche_upper <- read.csv("~/Documents/Proceedings Rockhopper/Data/nicheROVER_upper.csv",
                        header = TRUE, row.names = 1)


data_niche <- as.matrix(niche)
data_niche
data_niche_lower <- as.matrix(niche_lower)
data_niche_lower
data_niche_upper <- as.matrix(niche_upper)
data_niche_upper



corrplot(data_niche, method = "color", order = "hclust", tl.col = "black",is.corr = FALSE,na.label = "square",na.label.col = "black")

