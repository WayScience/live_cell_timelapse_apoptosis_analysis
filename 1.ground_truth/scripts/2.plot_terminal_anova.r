suppressPackageStartupMessages(suppressWarnings({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
}))

tukey_results_file_path <- file.path(
    "../data/0.ground_truth/tukey_results.parquet"
)
intensity_features_file_path <- file.path(
    "../data/0.ground_truth/annexinv_intensity_features_df.parquet"
)

figures_dir <- file.path("../figures/")
if (!dir.exists(figures_dir)) {
    dir.create(figures_dir)
}

# Read the intensity features
intensity_features_df <- arrow::read_parquet(intensity_features_file_path)

tukey_df <- arrow::read_parquet(tukey_results_file_path)
head(tukey_df)

# get the intensity feature only
tukey_df <- tukey_df %>% filter(feature == "Intensity_MedianIntensity_AnnexinV")

# tidy long
intensity_features_df <- intensity_features_df %>%
    pivot_longer(
        cols = colnames(intensity_features_df)[-1],
        names_to = "feature",
        values_to = "value"
    )

# select only annexin features
intensity_features_df$channel <- gsub("Intensity_", "", intensity_features_df$feature)
intensity_features_df$channel <- sub(".*_(.*)", "\\1", intensity_features_df$feature)
intensity_features_df <- intensity_features_df %>% filter(
    channel == "AnnexinV"
)
head(intensity_features_df)

intensity_features_df$Metadata_dose <- as.character(intensity_features_df$Metadata_dose)
intensity_features_df$Metadata_dose <- factor(
    intensity_features_df$Metadata_dose,
    levels = c(
        '0',
        '0.61',
        '1.22',
        '2.44',
        '4.88',
        '9.77',
        '19.53',
        '39.06',
        '78.13',
        '156.25'
)
    )

# plot the intensity_features_df
width <- 10
height <- 10
options(repr.plot.width = width, repr.plot.height = height)
intensity_plot <- (
    ggplot(intensity_features_df, aes(x = Metadata_dose, y = value, fill = Metadata_dose))
    + geom_boxplot(aes(group=Metadata_dose), outlier.size = 0.5, outlier.colour = "gray")
    + labs(
        x = "Dose (µM)",
        y = "Intensity",
    )
    + facet_wrap(~ feature, scales = "free_y", ncol = 2)
        + theme(
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size = 14)
    )
)
ggsave(
    file.path(figures_dir, "intensity_features_boxplot.png"),
    plot = intensity_plot,
    width = width,
    height = height,
    dpi = 600
)
intensity_plot

# sort the data by group1 and group2
tukey_df <- tukey_df %>%
    arrange(group1, group2)
# make sure the group1 and group2 are factors
tukey_df$comparison <- paste0(
  tukey_df$group1, "_", tukey_df$group2
)
tukey_df$channel <- sub(".*_(.*)", "\\1", tukey_df$feature)

width <- 15
height <- 15
options(repr.plot.width = width, repr.plot.height = height)
tukey_plot <- (
    ggplot(tukey_df, aes(x = meandiff, y = comparison, col = channel))
    + geom_point(size = 3)
    + geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2)
    + geom_hline(yintercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "Group Comparison",
        y = "Mean Difference"
    )
    + theme_bw()
    + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    + facet_wrap(feature ~ ., scales = "free_y")

        + facet_wrap(~ feature, scales = "free_y", ncol = 2)
        + theme(
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size = 14)
    )
)
ggsave(
    file.path(figures_dir, "tukey_plot.png"),
    plot = tukey_plot,
    width = width,
    height = height,
    dpi = 600
)
tukey_plot

