########## Utils for resampling strategy

resample.two <- function(cal, test, cal.pred, test.pred, seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    s <- sample(1:nrow(test), round(nrow(test) / 2))
    test1 <- test[s, ]
    test2 <- test[-s, ]
    p.test1 <- test.pred[s, ]
    p.test2 <- test.pred[-s, ]
    cal_freq <- prop.table(table(cal$Y))
    pr.class1 <- apply(p.test1, 1, function(row) colnames(p.test1)[which.max(row)])
    pr.class2 <- apply(p.test2, 1, function(row) colnames(p.test2)[which.max(row)])
    test_freq1 <- prop.table(table(pr.class1))
    test_freq2 <- prop.table(table(pr.class2))
    des_freq1 <- round(test_freq1 * length(cal$Y))
    des_freq2 <- round(test_freq2 * length(cal$Y))

    idx1 <- idx2 <- NULL
    for (i in cl.types) {
        cat <- which(cal$Y == i)
        if (!is.na(des_freq1[i])) {
            idx_cat1 <- sample(cat, size = des_freq1[i], replace = TRUE)
            idx1 <- c(idx1, idx_cat1)
        }
        if (!is.na(des_freq2[i])) {
            idx_cat2 <- sample(cat, size = des_freq2[i], replace = TRUE)
            idx2 <- c(idx2, idx_cat2)
        }
    }

    return(list(
        cal1 = cal[idx1, ], cal2 = cal[idx2, ], p.cal1 = cal.pred[idx1, ],
        p.cal2 = cal.pred[idx2, ], test1 = test1, test2 = test2,
        p.test1 = p.test1, p.test2 = p.test2
    ))
}

# resample.oracle <- function(cal, test, cal.pred, test.pred, seed=NA){
#   if(!is.na(seed)) set.seed(seed)
#   s <- sample(1:nrow(test), round(nrow(test)/2))
#   test1 <- test[s,]
#   test2 <- test[-s,]
#   p.test1 <- test.pred[s,]
#   p.test2 <- test.pred[-s,]
#   cal_freq <- prop.table(table(cal$Y))
#   test_freq1 <- prop.table(table(test1$Y))
#   test_freq2 <- prop.table(table(test2$Y))
#   des_freq1 <- round(test_freq1 * length(cal$Y))
#   des_freq2 <- round(test_freq2 * length(cal$Y))
#
#   idx1 <- idx2 <- NULL
#   for (i in cl.types) {
#     cat <- which(cal$Y == i)
#     if(!is.na(des_freq1[i])){
#       idx_cat1 <- sample(cat, size = des_freq1[i], replace = TRUE)
#       idx1 <- c(idx1, idx_cat1)
#     }
#     if(!is.na(des_freq2[i])){
#       idx_cat2 <- sample(cat, size = des_freq2[i], replace = TRUE)
#       idx2 <- c(idx2, idx_cat2)
#     }
#   }
#
#   return(list(cal1=cal[idx1,], cal2=cal[idx2,], p.cal1=cal.pred[idx1,],
#               p.cal2=cal.pred[idx2,], test1=test1, test2=test2,
#               p.test1=p.test1, p.test2=p.test2))
# }

resample.oracle <- function(cal, test, cal.pred, test.pred, seed = NA) {
    if (!is.na(seed)) set.seed(seed)

    cal_freq <- prop.table(table(cal$Y))
    test_freq <- prop.table(table(test$Y))
    des_freq <- round(test_freq * length(cal$Y))

    idx <- NULL
    for (i in cl.types) {
        cat <- which(cal$Y == i)
        if (!is.na(des_freq[i])) {
            idx_cat <- sample(cat, size = des_freq[i], replace = TRUE)
            idx <- c(idx, idx_cat)
        }
    }

    return(list(cal = cal[idx, ], p.cal = cal.pred[idx, ]))
}

resample.freq <- function(cal, test, cal.pred, test.pred, freqs, seed = NA) {
    if (!is.na(seed)) set.seed(seed)

    cal_freq <- prop.table(table(cal$Y))
    test_freq <- freqs
    des_freq <- round(test_freq * length(cal$Y))

    idx <- NULL
    for (i in cl.types) {
        cat <- which(cal$Y == i)
        if (!is.na(des_freq[i])) {
            idx_cat <- sample(cat, size = des_freq[i], replace = TRUE)
            idx <- c(idx, idx_cat)
        }
    }

    return(list(cal = cal[idx, ], p.cal = cal.pred[idx, ]))
}





exportedFn <- c("pred_sets1", "children", "ancestors", "g")
conformal <- function(cal, test, p.cal, p.test, alpha = 0.1, graph = opi2,
    lambdas = seq(0.001, 0.999, length.out = 100),
    expFn = exportedFn) {
    conf.mln <- getConfQuant(p.cal, cal$Y, alpha)
    pred.mln <- getPredSets(p.test, test$Y, conf.mln$qhat,
        summary = T
    )
    tf.conf <- rep(NA, length(pred.mln$pred.sets))
    for (i in 1:length(pred.mln$pred.sets)) {
        tf.conf[i] <- test$Y[i] %in% pred.mln$pred.sets[[i]]
    }

    sets <- foreach(lambda = lambdas, .export = expFn) %dopar% {
        lapply(
            1:nrow(p.cal),
            function(i) pred_sets1(lambda, p.cal[i, ], graph)
        )
    }
    l <- get_loss_table_mis(lambdas, sets, as.character(cal$Y), graph)
    lhat <- get_lhat(l, lambdas, alpha)
    sets.test <- apply(p.test, 1, function(x) pred_sets1(lhat, x, graph))
    tf.crc <- rep(NA, nrow(test))
    for (j in 1:nrow(test)) {
        tf.crc[j] <- as.character(test$Y)[j] %in% sets.test[[j]]
    }
    return(list(tf.conf = tf.conf, tf.crc = tf.crc, sets.conf = pred.mln$pred.sets, sets.crc = sets.test))
}

conf.loo <- function(p.cal, p.test, Ycal, Ytest, alpha) {
    ntest <- nrow(p.test)
    emp.cov <- rep(NA, ntest)
    test.sets <- list()
    for (i in 1:ntest) {
        p.test.loo <- p.test[-i, ]
        cal_freq <- prop.table(table(Ycal))
        test_prfreq <- apply(p.test.loo, 1, function(row) colnames(p.test.loo)[which.max(row)])
        test_freq <- prop.table(table(test_prfreq))
        des_freq <- round(test_freq * length(Ycal))

        # Resample
        idx <- NULL
        for (j in cl.types) {
            cat <- which(Ycal == j)
            if (!is.na(des_freq[j])) {
                idx_cat <- sample(cat, size = des_freq[j], replace = TRUE)
                idx <- c(idx, idx_cat)
            }
        }
        p.cal.loo <- p.cal[idx, ]
        Ycal.loo <- Ycal[idx]

        # Conformal prediction for i obs
        q <- getConfQuant(p.cal.loo, Ycal.loo, alpha)
        p <- getPredSets(t(as.matrix(p.test[i, ])), Ytest[i], q$qhat)
        emp.cov[i] <- p$coverage
        test.sets[i] <- p$pred.sets
        # if(i%%100 == 0) cat(i)
    }
    return(list(emp.cov = mean(emp.cov), sets = test.sets))
}
