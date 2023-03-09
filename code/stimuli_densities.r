# Christoph Aurnhammer, 2022
# Dissertation, Chapter 4
# Create density plots of Design 2 stimulus properties

source("plot_rERP.r")

quad_density <- function(path, file, leg_labs, leg_vals) {
    dt <- fread(path)

    cloze   <- items_and_means(dt, "Cloze")
    cloze_d <- items_and_means(dt, "Cloze_distractor")
    plaus   <- items_and_means(dt, "Plaus")
    plaus_d <- items_and_means(dt, "Plaus_distractor")

    cl <- plot_density(cloze[[1]], cloze[[2]],
            "Density", "Probability", "Cloze",
            leg_labs, leg_vals, c(0, 21), c(0, 0.25, 0.5, 0.75, 1))
    cl_d <- plot_density(cloze_d[[1]], cloze_d[[2]],
            "Density", "Probability", "Cloze_distractor",
            leg_labs, leg_vals, c(0, 21), c(0, 0.25, 0.5, 0.75, 1))
    pl <- plot_density(plaus[[1]], plaus[[2]],
            "Density", "Rating",  "Plaus",
            leg_labs, leg_vals, c(0, 1.5), c(1:7))
    pl_d <- plot_density(plaus_d[[1]], plaus_d[[2]],
            "Density", "Rating", "Plaus_distractor",
            leg_labs, leg_vals, c(0, 1.5), c(1:7))

    leg <- get_legend(pl)
    nl <- theme(legend.position = "none")
    ttlpos <- theme(plot.title = element_text(hjust = 0.5))
    gg <- arrangeGrob(arrangeGrob(
            cl + nl + ggtitle("Target Cloze Probability") + ttlpos,
            cl_d + nl + ggtitle("Distractor Cloze Probability") + ttlpos,
            pl + nl + ggtitle("Target Plausibility") + ttlpos,
            pl_d + nl + ggtitle("Distractor Plausibility") + ttlpos,
            layout_matrix = matrix(1:4, ncol = 2)),
            leg, heights = c(10, 1))
    ggsave(file, gg, device = cairo_pdf, width = 7, height = 7)
}

system(paste0("mkdir ../plots/Stimuli"))

quad_density(
    "../data/adbc23_stimuli.csv",
    "../plots/Stimuli/adbc23_Densities.pdf",
    c("A: Plausible", "B: Less plausible, attraction", "C: Implausible, no attraction"),
    c("black", "red", "blue")
)
