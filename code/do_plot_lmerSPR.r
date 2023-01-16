source("plot_lmerSPR.R")

produce_spr_plots <- function(
    path
) {
    lmerSPR <- fread(paste0("../data/", path, ".csv"))

    # make dirs
    system(paste0("mkdir ../plots/", path))

    if (grepl("PrecritRT", path)) {
        leg_vals <- c("#000000", "red", "#E349F6", "#00FFFF")
        leg_labs <- c("Intercept", "PrecritRT", "Plausibility",
            "Distractor Cloze")
    } else {
        leg_vals <- c("#000000", "#E349F6", "#00FFFF")
        leg_labs <- c("Intercept", "Plausibility", "Distractor Cloze")
    }

    # Observed Data
    plot_lmerSPR(lmerSPR, "logRT", yunit = "logRT", title = "Observed RTs",
        ylims = c(5.50, 5.80), name = path)

    # Estimated RTs
    estimates <- colnames(lmerSPR)[which(grepl("est_", colnames(lmerSPR)))]
    for (x in estimates) {
        plot_lmerSPR(lmerSPR, x, yunit = "logRT", title = "Estimated RTs",
            ylims = c(5.50, 5.80), name = path)
    }

    # Residuals
    residuals <- colnames(lmerSPR)[which(grepl("res_", colnames(lmerSPR)))]
    for (x in residuals) {
        plot_lmerSPR(lmerSPR, x, yunit = "logRT",
            title = "Residuals: Plausibility + Cloze Distractor",
            ylims = c(0.12, -0.12), name = path)
    }

    # Coefficients
    coefs_cols <- colnames(lmerSPR)[grep("coef", colnames(lmerSPR))]
    SE_cols <- colnames(lmerSPR)[grep("SE_", colnames(lmerSPR))]
    cols <- c("Condition", "Item", "Subject", "Region", coefs_cols, SE_cols)
    data <- lmerSPR[,..cols]
    data1 <- melt(data, id.vars=c("Condition", "Region"),
        measure.vars = coefs_cols, variable.name = "Coefficient",
        value.name = "coefficients")
    data1$SE <- melt(data, id.vars = c("Condition", "Region"),
        measure.vars = SE_cols, variable.name = "Coefficient",
        value.name = "SE")$SE
    data1 <- data1[, lapply(.SD, mean), by = list(Region, Coefficient),
        .SDcols = c("coefficients", "SE")]
    plot_lmerSPR(data1, "coefficients", yunit = "SPR Coefficients",
        title = "Coefficients", grouping = "Coefficient", name = path,
        leg_vals = leg_vals, leg_labs = leg_labs)

    # Z-values
    zval_cols <- setdiff(colnames(lmerSPR)[grep("zval_",
        colnames(lmerSPR))], "zval_1")
    pval_cols <- setdiff(colnames(lmerSPR)[grep("pval_",
        colnames(lmerSPR))], "pval_1")
    cols <- c("Condition", "Item", "Subject", "Region", zval_cols, pval_cols)
    data <- lmerSPR[,..cols]
    data1 <- melt(data, id.vars=c("Condition", "Region"),
        measure.vars = zval_cols, variable.name = "Zvalue",
        value.name = "zvalue")
    data1$pvalue <- melt(data, id.vars=c("Condition", "Region"),
        measure.vars = pval_cols, variable.name = "Pvalue",
        value.name="pvalue")$pvalue
    data1 <- data1[, lapply(.SD, mean), by = list(Region, Zvalue),
        .SDcols = c("zvalue", "pvalue")]
    plot_lmerSPR(data1, "zvalue", yunit="Z-values",
        title="Inferential Statistics", grouping="Zvalue",
        ylims=c(-2, 7), name=path,
        leg_vals = leg_vals[2:length(leg_labs)], 
        leg_labs = leg_labs[2:length(leg_labs)])
}

produce_spr_plots("rRTs_Plaus_Clozedist")

produce_spr_plots("rRTs_PrecritRT_Plaus_Clozedist")