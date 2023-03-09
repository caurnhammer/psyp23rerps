##
# Christoph Aurnhammer <aurnhammer@coli.uni-saarland.de>
# EEG plotting options for (lme)rRTs
##

library(data.table)
library(ggplot2)

se <- function(
    x,
    na.rm = FALSE
) {
    if (na.rm == TRUE) {
        sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))
    } else {
        sd(x) / sqrt(length(x))
    }
}

plot_lmerSPR <- function(
    data,
    DV,
    yunit,
    title,
    ylims = NULL,
    grouping = "Condition",
    name,
    leg_labs,
    leg_vals,
    omit_legend = FALSE,
    save_legend = FALSE
) {
    if (!(DV %in% c("coefficients", "zvalue"))){
        pm <- aggregate(data[[DV]] ~ Region + data[[grouping]] + Subject,
            data, FUN = mean)
        colnames(pm)[c(2,4)] <- c("group", "DV")

        plusminus <- aggregate(DV ~ Region + group, pm, FUN = mean)
        plusminus$SE <- aggregate(DV ~ Region + group, pm, FUN = se)[, 3]
        plusminus <- data.table(plusminus)
    } else if (DV == "zvalue") {
        plusminus <- data
        colnames(plusminus)[c(2, 3)] <- c("group", "DV")
        plusminus$sig <- plusminus$pvalue < 0.05
        df <- plusminus[, c("Region", "group", "pvalue", "sig")]
        df$posit <- rep(seq(ylims[1] - 3.5, ylims[1] - 0.3,
            length = length(unique(plusminus$group))),
            each = length(unique(plusminus$Region)))
        df$sig <- factor(df$sig, levels = c("TRUE", "FALSE"),
            labels = c("Significant", "Insignificant"))
        df$Region <- plusminus$Region
    } else {
        plusminus <- data
        colnames(plusminus)[c(2, 3)] <- c("group", "DV")
    }

    plusminus$Region <- factor(plusminus$Region,
        levels = c("Pre-critical", "Critical", "Spillover", "Post-spillover"))
    if (DV == "zvalue") {
        df$Region <- plusminus$Region
    }

    p <- ggplot(plusminus, aes(x = Region, y = DV, color = group,
        group = group)) + geom_point(size = 2.5, shape = "cross") +
        geom_line(size = 0.5)
    p <- p + theme_minimal()
    if (!(DV %in% c("zvalue", "coefficients"))) {
        p <- p + geom_errorbar(aes(ymin = DV - SE, ymax = DV + SE),
            width = .1, size = 0.3)
    }
    if (DV == "coefficients") {
        p <- p + geom_errorbar(aes(ymin = DV - SE, ymax = DV + SE),
            width = .1, size = 0.3)
        p <- p + scale_color_manual(name = "Coefficients",
            values = leg_vals, labels = leg_labs)
        # p <- p + guides(color = guide_legend(nrow = length(leg_labs), byrow = TRUE))
    } else if (DV == "zvalue") {
        p <- p + geom_hline(yintercept = 0, linetype = "dashed")
        p <- p + scale_color_manual(name = "Z-value",
            values = leg_vals, labels = leg_labs)
        p <- p + geom_point(data = df, aes(x = Region, y = posit, shape = sig),
            size = 2.5)
        p <- p + scale_shape_manual(values = c(20, 32),
            name = "Corrected p-values",
            labels = c("Significant", "Nonsignificant"))
    }
    else { # RTs, Residuals
        p <- p + scale_color_manual(name = "Condition",
            labels = c("A: Plausible", "B: Less plausible, attraction",
                        "C: Implausible, no attraction"),
            values = c("black", "red", "blue"))
        # p <- p + guides(color = guide_legend(nrow = 3, byrow = TRUE))
    }

    if ((is.vector(ylims) == TRUE) & (DV != "zvalue")) {
        p <- p + ylim(ylims[1], ylims[2])
    }

    p <- p + theme(plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        legend.box = "vertical",
        legend.spacing.y = unit(-0.2, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10, -10, -10, -10))
    p <- p + labs(x = "Region", y = yunit, title = title)

    # Option to omit legend and save it separately.
    if (omit_legend) {
        if (save_legend) {
            lgnd <- get_legend(p)
            # file_trimmed <- strtrim(name, nchar(file) - 4)
            ggsave(paste0("../plots/", name, "/RT_", DV, "_wavelegend.pdf"),
                    lgnd, device = cairo_pdf,
                    width = 3.6, height = 0.5)
        }
        p <- p + theme(legend.position = "none")
    } 

    file <- paste0("../plots/", name, "/RT_", DV, ".pdf")
    ggsave(file, p, width = 3, height = 3)
}
