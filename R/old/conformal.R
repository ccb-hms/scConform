#### Function to get the conformal quantile
getConfQuant <- function(p, label, alpha) {
    # p prediction matrix for data in the calibration set (ncal x num labels)
    # label true label for data in the calibration set
    # alpha desired error rate
    true <- rep(NA, dim(p)[1])
    for (i in 1:dim(p)[1]) {
        true[i] <- p[i, label[i]]
    }

    cal.scores <- 1 - true

    # get adjusted quantile
    n <- dim(p)[1]
    q_level <- ceiling((n + 1) * (1 - alpha)) / n
    qhat <- quantile(cal.scores, q_level)
    return(list(scores = cal.scores, qhat = qhat))
}

#### Function to get the prediction sets for test data
getPredSets <- function(pred, label, qhat, getClass = T, acc.p = NULL, summary = F) {
    # pred prediction matrix for data in the test set (ntest x num labels)
    # label true label for data in the test set (needed to evaluate accuracy and coverage)
    # qhat conformal quantile

    # get predicted class from prediction matrix
    if (getClass) {
        acc.p <- apply(pred, 1, function(row) colnames(pred)[which.max(row)])
    }
    acc <- mean(acc.p == label)
    prediction_sets <- pred >= (1 - qhat)
    rs <- rowSums(prediction_sets) # to have summary on size of prediction sets

    # Get prediction set colnames
    pr.list <- lapply(1:nrow(prediction_sets), function(i) {
        colnames(prediction_sets)[prediction_sets[i, ]]
    })

    # Get coverage
    tf <- rep(NA, length(pr.list))
    for (i in 1:length(pr.list)) {
        tf[i] <- label[i] %in% pr.list[[i]]
    }
    coverage <- sum(tf) / length(pr.list)

    if (summary) {
        return(list(accuracy = acc, pred.sets = pr.list, coverage = coverage, summary = summary(rs)))
    } else {
        return(list(accuracy = acc, pred.sets = pr.list, coverage = coverage))
    }
}

#### Function to get class conditional coverage
class_sp_conf <- function(classes, pred.cal, pred.test, labels.cal, labels.test) {
    s <- rep(NA, length(classes))
    names(s) <- classes

    for (i in 1:length(classes)) {
        p <- as.matrix(pred.cal[labels.cal == classes[i], ])
        s[i] <- getConfQuant(p, rep(classes[i], nrow(p)), 0.05)$qhat
    }

    acc.p <- apply(pred.test, 1, function(row) colnames(pred.test)[which.max(row)])
    p.sets <- matrix(NA, nrow = nrow(pred.test), ncol = ncol(pred.test))
    colnames(p.sets) <- colnames(pred.test)
    for (i in 1:length(classes)) {
        p.sets[, classes[i]] <- pred.test[, classes[i]] >= (1 - s[classes[i]])
    }

    pr.list <- lapply(1:nrow(p.sets), function(i) {
        colnames(p.sets)[p.sets[i, ]]
    })

    # Get coverage
    tf <- rep(NA, length(pr.list))
    for (i in 1:length(pr.list)) {
        tf[i] <- labels.test[i] %in% pr.list[[i]]
    }
    return(list(pred.sets = pr.list, coverage = mean(tf)))
}
