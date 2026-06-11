packages <- c("ggplot2", "dplyr", "patchwork", "viridis", "platetools")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    suppressPackageStartupMessages(
        suppressWarnings(
            library(pkg, character.only = TRUE)
        )
    )
  }
}
source("../../utils/r_themes.r")


train_splits_df <- arrow::read_parquet("../data_splits/train.parquet")
test_splits_df <- arrow::read_parquet("../data_splits/test.parquet")
train_test_wells_df <- arrow::read_parquet("../data_splits/train_test_wells.parquet")

columns_to_keep <- c("Metadata_Well", "Metadata_dose")
train_splits_df <- train_splits_df %>% select(all_of(columns_to_keep))
test_splits_df <- test_splits_df %>% select(all_of(columns_to_keep))
# combine train and test splits for plotting
df <- rbind(
    train_splits_df,
    test_splits_df
)

df <- merge(
    df,
    train_test_wells_df,
)
# replace the "-" with ""
df$Metadata_Well <- gsub("-", "", df$Metadata_Well)
df$data_split <- gsub("train", "Train", df$data_split)
df$data_split <- gsub("test", "Test", df$data_split)
head(df)

train_test_split_platemap <- (
    raw_map(data = df$data_split,
        well = df$Metadata_Well,
        plate = 96)
    + theme_dark()
    # change the legend title
    + labs(fill = "Data Split")
)


train_test_split_platemap_dose <- (
    raw_map(data = factor(df$Metadata_dose),
            well = df$Metadata_Well,
            plate = 96)
    + theme_dark()
    # use color palette for dose
    + scale_fill_manual(values = color_palette_dose)
    # set the legend title to "Dose"
    + labs(fill = "Staurosporine Dose (nM)")
    # make the legend two columns
    + guides(fill = guide_legend(ncol = 2))
)


final_plot <- (
    train_test_split_platemap / train_test_split_platemap_dose
    # + plot_layout(guides = "collect")
    # add annotations for A and B and make them big
    + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 20))
)
final_plot
ggsave("../figures/train_test_split_platemap.png", final_plot, width = 8, height = 10, dpi = 600)
