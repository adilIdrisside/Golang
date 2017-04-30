dist_mat = read.table ('Distances_R.txt',TRUE,sep=';')
dist_mat =dist_mat [, 1:16]
file_names = names (dist_mat)
file_names=file_names [2:16]
library (data.table)
setattr(dist_mat , "row.names", file_names)
dist_mat =dist_mat [,2:16]
distance <- dist(dist_mat , method = "euclidean")
fit <- hclust(distance, method="ward.D")
plot(fit)
groups <- cutree(fit, k=5)
rect.hclust(fit, k=5, border="red")
