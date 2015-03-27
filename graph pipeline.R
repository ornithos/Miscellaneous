dat <- read.csv("C:/Users/pc/Documents/revision.csv")
dat1 <- as.matrix(dat[,2:14])
n <- nrow(dat1)
cols <- ncol(dat1)
for(i in 1:cols) dat1[is.na(dat1[,i]),i] <- 0

dist <- cbind(expand.grid(1:n, 1:n), 0)
dist <- dist[dist[,1] < dist[ ,2],]
dist <- dist[order(dist[,1], dist[,2]),]
cd <- apply(dat1, 1, function(x) (dat1 %*% x) / (sqrt(sum(x^2))*apply(dat1, 1, function(a) sqrt(sum(a^2))) ) )
dist[,3] <- cd[lower.tri(cd)]

nodes <- data.frame(ID = 1:n, Label = dat[,1], Size = rowSums(dat1))
edges <- data.frame(Source = dist[,1], Target = dist[,2], Type = "Undirected", Weight = dist[,3])
write.csv(nodes, "nodes.csv"); write.csv(edges, "edges.csv")

# -> import these two files into Gephi
