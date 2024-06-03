# Function to get conformal prediction sets with split conformal inference

.getConformalPredSets <- function(p.cal, p.test, y.cal, alpha) {
    # Get calibration scores (1-predicted probability for the true class)
    true <- rep(NA, dim(p.cal)[1])
    for (i in seq_len(dim(p.cal)[1])) {
        true[i] <- p.cal[i, y.cal[i]]
    }
    s <- 1 - true

    # Get adjusted quantile
    n <- nrow(p.cal)
    q_level <- ceiling((n + 1) * (1 - alpha)) / n
    qhat <- quantile(s, q_level)

    # Get prediction sets
    prediction_sets <- p.test >= 1 - qhat
    pr.list <- lapply(seq_len(nrow(prediction_sets)), function(i) {
        colnames(prediction_sets)[prediction_sets[i, ]]
    })
    return(pr.list)
}
