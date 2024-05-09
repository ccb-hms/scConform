
resample.two <- function(p.cal, p.test, y.cal, labels){
    s <- sample(1:nrow(test), round(nrow(test)/2))
    test1 <- p.test[s,]
    test2 <- p.test[-s,]
    # cal_freq <- prop.table(table(y.cal))
    # Compute predicted class
    pr.class1 <- apply(test1, 1, function(row) colnames(test1)[which.max(row)])
    pr.class2 <- apply(test2, 1, function(row) colnames(test2)[which.max(row)])
    test_freq1 <- prop.table(table(pr.class1))
    test_freq2 <- prop.table(table(pr.class2))
    # Transform to absolute frequencies
    des_freq1 <- round(test_freq1 * length(y.cal))
    des_freq2 <- round(test_freq2 * length(y.cal))

    idx1 <- idx2 <- NULL
    for (i in labels) {
      cat <- which(y.cal == i)
      if(!is.na(des_freq1[i])){
        idx_cat1 <- sample(cat, size = des_freq1[i], replace = TRUE)
        idx1 <- c(idx1, idx_cat1)
      }
      if(!is.na(des_freq2[i])){
        idx_cat2 <- sample(cat, size = des_freq2[i], replace = TRUE)
        idx2 <- c(idx2, idx_cat2)
      }
    }

    return(list(p.cal1=p.cal[idx1,],
                p.cal2=p.cal[idx2,],
                p.test1=p.test1,
                p.test2=p.test2,
                y.cal1=y.cal[idx1],
                y.cal2=y.cal[idx2])
           )
}

