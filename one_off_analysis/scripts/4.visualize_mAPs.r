suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(arrow)))
suppressWarnings(suppressPackageStartupMessages(library(argparse)))

percent_cell_mAP_file_path <- file.path("../results/mAP_cell_percentages.parquet")
across_channels_mAP_file_path <- file.path("../results/mAP_across_channels.parquet")

percent_cell_mAP <- arrow::read_parquet(percent_cell_mAP_file_path)
across_channels_mAP <- arrow::read_parquet(across_channels_mAP_file_path)
dim(percent_cell_mAP)
dim(across_channels_mAP)

percent_cell_mAP <- percent_cell_mAP %>% filter(Metadata_treatment == "Thapsigargin 10 uM")
percent_cell_mAP <- percent_cell_mAP %>% group_by(shuffle, percentage_of_cells,Metadata_Time, Metadata_treatment) %>%
  summarise(mAP = mean(mean_average_precision))
head(percent_cell_mAP)

width <- 15
height <- 15
options(repr.plot.width = width, repr.plot.height = height)
percent_cell_line_plot <- (
    ggplot(data = percent_cell_mAP, aes(x = Metadata_Time, y = mAP))
    + geom_line(aes(
        group = Metadata_treatment,
        color = Metadata_treatment))
    + facet_grid(percentage_of_cells ~ shuffle)
    + ylim(0, 1)

)
percent_cell_line_plot



# get the first time point only
percent_cell_mAP <- percent_cell_mAP %>% filter(Metadata_Time == "09")

percent_cell_mAP <- percent_cell_mAP %>% group_by(Metadata_treatment, shuffle, percentage_of_cells) %>%
  summarise(mAP = mean(mean_average_precision))

width <- 30
height <- 30
options(repr.plot.width = width, repr.plot.height = height)
percent_cell_plot <- (
    ggplot(data = percent_cell_mAP, aes(x = Metadata_treatment, y = mAP, fill=Metadata_treatment))
    + geom_bar(stat = "identity", position = "dodge")
    + facet_wrap(shuffle~percentage_of_cells, scales = "free")
)
percent_cell_plot

across_channels_mAP <- arrow::read_parquet(across_channels_mAP_file_path)

across_channels_mAP <- across_channels_mAP %>%
  group_by(Metadata_treatment, Metadata_Time, shuffle, Channel) %>%
  summarise(mAP = mean(mean_average_precision))
unique(across_channels_mAP$Metadata_Time)
# selct only 'LPS 1 ug/ml + Nigericin 5uM',
across_channels_mAP <- across_channels_mAP %>%
  filter(Metadata_treatment == 'LPS 1 ug/ml + Nigericin 5uM')

width <- 15
height <- 15
options(repr.plot.width = width, repr.plot.height = height)
channels_mAP_plot <- (
    ggplot(data = across_channels_mAP, aes(x = Metadata_Time, y = mAP))
    + geom_line(aes(group = Metadata_treatment, color = Metadata_treatment))
    + facet_wrap(shuffle ~ Channel)
    + ylim(0, 1)

)
channels_mAP_plot


