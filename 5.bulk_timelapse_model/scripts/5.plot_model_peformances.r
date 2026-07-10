for (pkg in c("ggplot2", "dplyr", "patchwork", "ggplotify")) {
    suppressPackageStartupMessages(
        suppressWarnings(
            library(
                pkg,
                character.only = TRUE,
                quietly = TRUE,
                warn.conflicts = FALSE
            )
        )
    )
}

performances_file_path <- file.path("..", "results", "model_performances.parquet")
df <- arrow::read_parquet(performances_file_path)
df$feature <- gsub("all_terminal_features", "All features", df$feature)
df$feature <- gsub("Terminal_Cytoplasm_Intensity_UpperQuartileIntensity_AnnexinV", "AnnexinV\nsingle feature", df$feature)
df$shuffled <- gsub("shuffled", "Shuffled", df$shuffled)
df$shuffled <- gsub("not_Shuffled", "Not shuffled", df$shuffled)
df$shuffled <- factor(df$shuffled, levels = c("Not shuffled", "Shuffled"))


width <- 8
height <- 6
options(repr.plot.width = width, repr.plot.height = height)

df

mse_plot <- (
    ggplot(df, aes(x = feature, y = mse, fill = shuffled))
    + geom_bar(stat = "identity", position = position_dodge())
    + theme_bw()

    + labs(
        x = "Feature Set",
        y = "Mean Squared Error (MSE)"
    )
    + theme(text = element_text(size = 18))
    + guides(fill = guide_legend(title = "Model Type"))
    + facet_wrap(~data_split, nrow = 1, scales = "free_x")
)
mse_plot

r2_plot <- (
    ggplot(df, aes(x = feature, y = r2, fill = shuffled))
    + geom_bar(stat = "identity", position = position_dodge())
    + theme_bw()

    + labs(
        x = "Feature Set",
        y = expression(R^2)
    )
    + ylim(min(df$r2) - 0.1, 1)
    + theme(text = element_text(size = 18))
    + guides(fill = guide_legend(title = "Model Type"))
    + facet_wrap(~data_split, nrow = 1, scales = "free_x")


)
r2_plot

width <- 12
height <- 6
dpi <- 600
options(repr.plot.width = width, repr.plot.height = height, repr.plot.res = dpi)
final_plot <- (
    mse_plot
    + r2_plot
    + plot_layout(nrow = 1, guides = "collect")
    + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 24))
)
png(
    filename = file.path("..", "figures", "model_performances.png"),
    width = width,
    height = height,
    units = "in",
    res = dpi
)
final_plot
dev.off()
final_plot
