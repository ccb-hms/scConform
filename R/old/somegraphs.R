###############################################################################
################### Test repository ###########################################
###############################################################################

#####################
# Build example graph
#####################
library(igraph)
t <- graph_from_literal(
    animale - +cane:gatto:topo, gatto - +british:persiano,
    gatto - +retriever,
    cane - +cocker:retriever, retriever - +golden:labrador
)
plot(t, layout = layout_as_tree(t, root = "animale"))
p <- c(0.2, 0.3, 0.2, 0.10, 0.15, 0.05)
names(p) <- c("cocker", "golden", "labrador", "british", "persiano", "topo")
p


######################
### Load old functions
######################
source("/Users/daniela/Documents/GitHub/ConfCell/R/old/conformal.R")
source("/Users/daniela/Documents/GitHub/ConfCell/R/old/utils_Hier.R")

# Check that scores are the same
nam <- V(t)$name
s <- sapply(nam, function(x) g(p, x, t))
s1 <- sapply(nam, function(x) .scores(p, x, t))
sum(s != s1) # ok

# Check ancestors
anc <- sapply(nam, function(x) ancestors(x, graph = t, include_self = T))
anc1 <- sapply(nam, function(x) .ancestors(x, onto = t, include_self = T))
tf <- rep(NA, length(anc))
for (i in 1:length(anc)) {
    tf[i] <- mean(anc[[i]] == anc1[[i]])
}
tf

anc <- sapply(nam, function(x) ancestors(x, graph = t, include_self = F))
anc1 <- sapply(nam, function(x) .ancestors(x, onto = t, include_self = F))
tf <- rep(NA, length(anc))
for (i in 1:length(anc)) {
    tf[i] <- mean(anc[[i]] == anc1[[i]])
}
tf

# Check children
child <- sapply(nam, function(x) children(x, graph = t, leaf = T))
child1 <- sapply(nam, function(x) .children(x, onto = t, leaf = T))
tf <- rep(NA, length(child))
for (i in 1:length(child)) {
    tf[i] <- mean(child[[i]] == child1[[i]])
}
tf

child <- sapply(nam, function(x) children(x, graph = t, leaf = F))
child1 <- sapply(nam, function(x) .children(x, onto = t, leaf = F))
tf <- rep(NA, length(child))
for (i in 1:length(child)) {
    tf[i] <- mean(child[[i]] == child1[[i]])
}
tf

# .predSets
pred_sets1(0.5, p, t)
.predSets(lambda = 0.5, pred = p, onto = t)

pred_sets1(0.7, p, t)
.predSets(lambda = 0.7, pred = p, onto = t)

pred_sets1(0.75, p, t)
.predSets(lambda = 0.75, pred = p, onto = t)

pred_sets1(0.8, p, t) # should be all
.predSets(lambda = 0.8, pred = p, onto = t)

pred_sets1(1, p, t)
.predSets(lambda = 1, pred = p, onto = t)

pred_sets1(0.6, p, t) # cane (golden, labrador, cocker)
.predSets(lambda = 0.6, pred = p, onto = t)

p1 <- c(0.3, 0.4, 0.3, 0, 0, 0)
names(p1) <- c("cocker", "golden", "labrador", "british", "persiano", "topo")
pred_sets1(0.7, p1, t)
.predSets(lambda = 0.7, pred = p1, onto = t)

#########################################
### Create a random dataset
#########################################
set.seed(1010)
n <- 2100

leaves <- c(4, 6, 7, 8, 9)
nclasses <- length(leaves)
lab <- sample(leaves, n, replace = T)
x <- (runif(n, 0, 3.1) + lab)
y <- rep(NA, length(lab))
y[lab == 4] <- "cocker"
y[lab == 6] <- "golden"
y[lab == 7] <- "labrador"
y[lab == 8] <- "british"
y[lab == 9] <- "persiano"
xydf <- data.frame(y = as.factor(y), x = x)
head(xydf)

# Divide in train, cal and test and fit model
library(VGAM)
train <- xydf[1:1000, ]
cal <- xydf[1001:1100, ]
test <- xydf[1101:2100, ]


fit <- vglm(y ~ x,
    family = multinomial(refLevel = "cocker"),
    data = train
)
lambdas <- seq(0.001, 0.999, length.out = 100)
# Compute predicted probabilities for calibration points
p.cal <- predict(fit, newdata = cal, type = "response")
library(doParallel)
library(foreach)
#
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

t3 <- graph_from_literal(
    animale - +cane:gatto, gatto - +british:persiano,
    gatto - +retriever,
    cane - +cocker:retriever, retriever - +golden:labrador
)
plot(t3, layout = layout_as_tree(t3, root = "animale"))
exportedFn <- c(".predSets", ".scores", ".children", ".ancestors")
system.time(sets1 <- foreach(lambda = lambdas, .packages = c("ConfCell", "igraph")) %dopar% {
    library(igraph)
    lapply(1:nrow(p.cal), function(i) pred_sets1(lambda, p.cal[i, ], t3))
})

#
system.time(l1 <- get_loss_table_mis(lambdas, sets1, as.character(cal$y), t3))
Rhat <- colMeans(l1)
plot(Rhat)
all(diff(Rhat) <= 0)

lhat <- get_lhat(l1, lambdas, 0.1)
p.test <- predict(fit, newdata = test, type = "response")
sets.test <- apply(p.test, 1, function(x) pred_sets1(lambda = lhat, x, t3))

t <- .getHierarchicalPredSets(p.cal, p.test, cal$y, onto = t3, alpha = 0.1, lambdas = lambdas)
t1 <- .getConformalPredSets(p.cal, p.test, as.character(cal$y), 0.1)
length(t1)
length(t$sets.test)

l.std <- sapply(t1, length)
l.crc <- sapply(t$sets.test, length)
barplot(table(l.std))
barplot(table(l.crc))
cvg <- cvg1 <- rep(NA, length(t))
for (i in 1:length(t1)) {
    cvg[i] <- test$y[i] %in% t1[[i]]
}

for (i in 1:length(t$sets.test)) {
    cvg1[i] <- test$y[i] %in% t$sets.test[[i]]
}
mean(cvg)
mean(cvg1)





sets <- foreach(j = lambdas, .export = exportedFn) %dopar% {
    library(igraph)
    lapply(
        1:nrow(p.cal),
        function(i) .predSets(lambda = j, pred = p.cal[i, ], onto = t3)
    )
}
