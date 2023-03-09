# Benjamini-Hochberg procedure

bh_apply_wide <- function(
        data,
        elec,
        alpha=0.05,
        tws=list(c(300, 500), c(600, 1000))) {
        preds <- unique(data$Spec)
        keep_ts <- c()
        for (tw in tws) {
                keep_ts <- c(keep_ts, seq(tw[1], tw[2]))
                for (p in preds) {
                        df <- data[Spec == p & Timestamp >= tw[1] &
                                               Timestamp <= tw[2], c(4:29)]
                        uncorrected <- unlist(df)
                        corrected <- p.adjust(uncorrected, method = "fdr")
                        corrected_matrix <- matrix(corrected,
                                                   nrow = nrow(df), ncol = 26)

                        data[Spec == p & Timestamp >= tw[1] &
                                         Timestamp <= tw[2],
                                         c(4:29)] <- data.table(
                                                        corrected_matrix)
                        data[Spec == p & Timestamp >= tw[1] &
                                         Timestamp <= tw[2],
                                         paste0(colnames(data)[4:29],
                                                "_sig")] <- data.table(
                                                    corrected_matrix < alpha)
                }
        }
        sigcols <- grep("_sig", colnames(data))
        data[!(Timestamp %in% keep_ts), c(c(4:29), sigcols)] <- 0
        data
}