source("plot_rERP.r")
source("benjamini-hochberg.r")

make_plots <- function(
    file,
    elec = c("F3", "Fz", "F4", "C3", "Cz", "C4", "P3", "Pz", "P4"),
    predictor = "Intercept"
) {
    # make dirs
    system(paste0("mkdir ../plots/", file))
    system(paste0("mkdir ../plots/", file, "/Waveforms"))
    system(paste0("mkdir ../plots/", file, "/Topos"))

    # nine elecs
    elec_nine <-  c("F3", "Fz", "F4", "C3", "Cz", "C4", "P3", "Pz", "P4")

    # MODELS
    # Coefficents
    mod <- fread(paste0("../data/", file, "_models.csv"))
    mod$Spec <- factor(mod$Spec, levels = predictor)
    coef <- mod[Type == "Coefficient", ]
    coef$Condition <- coef$Spec
    if (predictor[2] != "ReadingTime") {
        model_labs <- c("Intercept", "Plausibility", "Distractor Cloze")
        model_vals <- c("black", "#E349F6", "#00FFFF")
        plot_single_elec(coef, "C3",
            file = paste0("../plots/", file, "/Waveforms/Coefficients_C3.pdf"),
            modus = "Coefficient", ylims = c(10, -5.5),
            leg_labs = model_labs, leg_vals = model_vals,
            title = "Coefficients",
            omit_legend = TRUE, save_legend = TRUE)
        plot_single_elec(coef, "Pz", 
            file = paste0("../plots/", file, "/Waveforms/Coefficients_Pz.pdf"),
            modus = "Coefficient", ylims = c(10, -5.5),
            leg_labs = model_labs, leg_vals = model_vals,
            title = "Coefficients",
            omit_legend = TRUE, save_legend = FALSE)
    } else {
        model_labs <- c("Intercept", "Reading Time")
        model_vals <- c("black", "#E349F6")
        plot_single_elec(coef, "Pz",
            file = paste0("../plots/", file, "/Waveforms/Coefficients_Pz.pdf"),
            modus = "Coefficient", ylims = c(11, -5.5),
            title = "Coefficients",
            leg_labs = model_labs, leg_vals = model_vals)
    }

    # Models: t-value
    time_windows <- list(c(300, 500), c(600, 1000))
    tval <- mod[Type == "t-value" & Spec != "Intercept", ]
    sig <- mod[Type == "p-value" & Spec != "Intercept", ]
    colnames(sig) <- gsub("_CI", "_sig", colnames(sig))

    sig_corr <- bh_apply_wide(sig, elec, alpha=0.05, tws = time_windows)
    sigcols <- grepl("_sig", colnames(sig_corr))
    tval <- cbind(tval, sig_corr[,..sigcols])
    mod$Spec <- factor(mod$Spec, levels = predictor)
    tval$Condition <- tval$Spec

    plot_nine_elec(tval, elec_nine,
        file = paste0("../plots/", file, "/Waveforms/t-values.pdf"),
        title = "Inferential statistics",
        modus = "t-value", ylims = c(7, -5), tws = time_windows,
        leg_labs = model_labs[2:length(model_vals)],
        leg_vals = model_vals[2:length(model_vals)])

    ## DATA
    eeg <- fread(paste0("../data/", file, "_data.csv"))
    eeg$Condition <- factor(plyr::mapvalues(eeg$Condition, c(2, 1, 3),
                             c("B", "A", "C")), levels = c("A", "B", "C"))

    data_labs <- c("A: Plausible", "B: Less plausible, attraction",
                   "C: Implausible, no attraction")
    data_vals <- c("black", "red", "blue")

    # Observed
    obs <- eeg[Type == "EEG", ]
    plot_full_elec(data = obs, e = elec_all, 
        file = paste0("../plots/", file, "/Waveforms/Observed_Full.pdf"),
        title = "Observed ERPs", modus = "Condition",
        ylims = c(9, -5), leg_labs = data_labs, leg_vals = data_vals)
    plot_single_elec(obs, "Pz",
            file = paste0("../plots/", file, "/Waveforms/Observed_Pz.pdf"),
            modus = "Condition", ylims = c(10.5, -5.5),
            leg_labs = data_labs, leg_vals = data_vals)

    plot_topo(obs, file = paste0("../plots/", file, "/Topos/Observed"),
                tw = c(250, 400), cond_man = "B", cond_base = "A",
                add_title = "\nObserved", omit_legend = TRUE,
                save_legend = TRUE)
    plot_topo(obs, file = paste0("../plots/", file, "/Topos/Observed"),
                tw = c(300, 500), cond_man = "B", cond_base = "A",
                add_title = "\nObserved", omit_legend = TRUE,
                save_legend = TRUE)
    plot_topo(obs, file = paste0("../plots/", file, "/Topos/Observed"),
                tw = c(600, 1000), cond_man = "B", cond_base = "A",
                add_title = "\nObserved", omit_legend = TRUE)
    plot_topo(obs, file = paste0("../plots/", file, "/Topos/Observed"),
                tw = c(250, 400), cond_man = "C", cond_base = "A",
                add_title = "\nObserved", omit_legend = TRUE)
    plot_topo(obs, file = paste0("../plots/", file, "/Topos/Observed"),
                tw = c(300, 500), cond_man = "C", cond_base = "A",
                add_title = "\nObserved", omit_legend = TRUE,
                save_legend = TRUE)
    plot_topo(obs, file = paste0("../plots/", file, "/Topos/Observed"),
                tw = c(600, 1000), cond_man = "C", cond_base = "A",
                add_title = "\nObserved", omit_legend = TRUE)

    # Estimated
    est <- eeg[Type == "est",]
    pred <- c("Intercept", "Tar. Plaus.", "Dist. Cloze",
              "Dist. Cloze + Tar. Plaus.")
    predictor <- c("Intercept", "target plausibility", "distractor cloze",
              "distractor cloze + target plausibility")
    for (i in seq(1, length(unique(est$Spec)))) {
        spec <- unique(est$Spec)[i]
        est_set <- est[Spec == spec, ]
        spec <- unique(est_set$Spec)
        name <- gsub("\\[|\\]|:|,| ", "", spec)
        if (i < length(unique(est$Spec))) {
            plot_single_elec(est_set, "Pz",
                file = paste0("../plots/", file, "/Waveforms/EstimatedPz_",
                name, ".pdf"), modus = "Condition", ylims = c(10.5, -5.5),
                leg_labs = data_labs, leg_vals = data_vals,
                title = paste("Isolated contribution of", predictor[i]),
                omit_legend = TRUE, save_legend = FALSE)
        }
        else {
            plot_single_elec(est_set, "Pz",
                file = paste0("../plots/", file, "/Waveforms/EstimatedPz_",
                name, ".pdf"), modus = "Condition", ylims = c(10.5, -5.5),
                leg_labs = data_labs, leg_vals = data_vals,
                title = paste("Estimated ERPs"),
                omit_legend = TRUE, save_legend = TRUE)
        }

        plot_topo(est_set,
            file = paste0("../plots/", file, "/Topos/Estimated_", name),
            tw = c(600, 1000), cond_man = "B", cond_base = "A",
            add_title = paste("\nEstimate", pred[i]), omit_legend = TRUE)
    }

    # Residual
    res <- eeg[Type == "res", ]
    spec <- unique(est$Spec)[max(length(unique(est$Spec)))]
    res_set <- res[Spec == spec, ]
    plot_single_elec(res_set, "Pz",
        file = paste0("../plots/", file, "/Waveforms/ResidualPz_",
        name, ".pdf"), modus = "Condition", ylims = c(4.7, -4),
        leg_labs = data_labs, leg_vals = data_vals,
        title = paste("Residuals (Observed - Estimated)"),
        omit_legend = TRUE, save_legend = FALSE)
}

elec_all <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "FC5",
                "FC1", "FC2", "FC6", "C3", "Cz", "C4", "CP5", "CP1",
                "CP2", "CP6", "P7", "P3", "Pz", "P4", "P8", "O1", "Oz", "O2")

make_plots("rERPs_Plaus_Clozedist", elec_all,
     predictor = c("Intercept", "Plaus", "Cloze_distractor"))

make_plots("rERPs_Plaus_Clozedist_across",
    predictor = c("Intercept", "Plaus", "Cloze_distractor"))

make_plots("rERPS_ReadingTime", elec_all,
   predictor = c("Intercept", "ReadingTime"))
