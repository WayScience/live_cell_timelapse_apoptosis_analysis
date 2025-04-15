suppressPackageStartupMessages(suppressWarnings({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
}))

map_results_file_path <- file.path(
    "../data/0.ground_truth/map.parquet"
)

map_df <- arrow::read_parquet(map_results_file_path)
head(map_df)

map_df$Metadata_dose <- factor(
    map_df$Metadata_dose,
    levels = c(
        "0.61",
        "1.22",
        "2.44",
        "4.88",
        "9.77",
        "19.53",
        "39.06",
        "78.13",
        "156.25"
))


map_plot <- (
    ggplot(map_df, aes(x=Metadata_dose, y=mean_average_precision, fill=Metadata_dose))
    + geom_bar(stat="identity")
    + ylim(0,1)
    + theme_bw()
    + labs(
        title="Mean Average Precision by Dose",
        x="Dose (uM)",
        y="Mean Average Precision"
    )
    + theme(
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none"
    )
)
ggsave("../figures/mean_average_precision_by_dose.png", plot=map_plot, width=6, height=4)
map_plot

# write the ground truth map_df to a csv file
ground_truth_df <- data.frame(
    Metadata_dose = map_df$Metadata_dose,
    mean_average_precision = map_df$mean_average_precision
)
ground_truth_df$apoptosis <- "control"
# change to control negative or positive
ground_truth_df <- ground_truth_df %>% mutate(
    apoptosis = ifelse(mean_average_precision > 0.8, "positive", "negative")
)
# drop the mean_average_precision column
ground_truth_df <- ground_truth_df %>% select(-mean_average_precision)
# sort by Metadata_dose
# add the 0 dose row
ground_truth_df <- ground_truth_df %>% add_row(
    Metadata_dose = "0",
    apoptosis = "control"
)
ground_truth_df$Metadata_dose <- factor(
    ground_truth_df$Metadata_dose,
    levels = c(
        "0",
        "0.61",
        "1.22",
        "2.44",
        "4.88",
        "9.77",
        "19.53",
        "39.06",
        "78.13",
        "156.25"
    )
)

ground_truth_df
write.csv(
    ground_truth_df,
    file = "../data/0.ground_truth/ground_truth.csv",
    row.names = FALSE
)
