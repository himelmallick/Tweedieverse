# Generate a matrix with feature and group data (from SDA package)
featureInfo <- matrix(runif(800, -2, 5), ncol = 40)
featureInfo[featureInfo<0] <- 0
feat.names <- paste("feature", 1:20, sep = '')
rownames(featureInfo) <- feat.names
sub.names <- paste('subject', 1:40, sep = '')
colnames(featureInfo) <- sub.names
groupInfo <- data.frame(grouping=matrix(sample(0:1, 40, replace = TRUE),
                                        ncol = 1))
rownames(groupInfo) <- colnames(featureInfo)

SEdata <- list(feature = featureInfo, group = groupInfo)

expect_error(Tweedieverse(SEdata))
