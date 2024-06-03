# Function to implement the resampling strategy when calibration and test
# set supposedly have a different distribution of the cell types.
# Right now it implements a two-fold strategy, dividing randomly
# test data in two.

resample.two <- function(p.cal, p.test, y.cal, labels) {
    s <- sample(seq_len(nrow(p.test)), round(nrow(p.test) / 2))
    test1 <- p.test[s, ]
    test2 <- p.test[-s, ]

    # Compute predicted class
    pr.class1 <- apply(test1, 1, function(row) colnames(test1)[which.max(row)])
    pr.class2 <- apply(test2, 1, function(row) colnames(test2)[which.max(row)])
    test.freq1 <- prop.table(table(pr.class1))
    test.freq2 <- prop.table(table(pr.class2))
    # Transform to absolute frequencies
    des.freq1 <- round(test.freq1 * length(y.cal))
    des.freq2 <- round(test.freq2 * length(y.cal))

    idx1 <- idx2 <- NULL
    for (i in labels) {
        category <- which(y.cal == i)
        if (!is.na(des.freq1[i])) {
            idx.category1 <- sample(category, size = des.freq1[i], replace = TRUE)
            idx1 <- c(idx1, idx.category1)
        }
        if (!is.na(des.freq2[i])) {
            idx.category2 <- sample(category, size = des.freq2[i], replace = TRUE)
            idx2 <- c(idx2, idx.category2)
        }
    }

    return(
        list(
            p.cal1 = p.cal[idx1, ],
            p.cal2 = p.cal[idx2, ],
            p.test1 = test1,
            p.test2 = test2,
            y.cal1 = y.cal[idx1],
            y.cal2 = y.cal[idx2],
            idx = c(s, setdiff(seq_len(nrow(p.test)), s))
        ) # index in the original data
    )
}
