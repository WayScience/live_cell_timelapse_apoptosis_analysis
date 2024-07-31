suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(argparse)))

# Create an ArgumentParser object
parser <- ArgumentParser(description = "UMAP Visualization Script")

# Define the arguments
parser$add_argument("--data_mode", type = "character", help = "data to plot", required = TRUE)
# Parse the arguments
args <- parser$parse_args()

data_mode <- args$data_mode

# set paths
umap_file_path <- file.path("../../data/umap/",paste0(data_mode,"_umap_transformed.parquet"))
umap_file_path <- normalizePath(umap_file_path)
figures_path <- file.path(paste0("../figures/",data_mode,"/"))
if (!dir.exists(figures_path)) {
  dir.create(figures_path)
}

umap_df <- arrow::read_parquet(umap_file_path)


head(umap_df,1)

# add nM to the dose column
umap_df$Metadata_dose <- paste0(umap_df$Metadata_dose, " nM")

# make the dose a factor with levels
umap_df$Metadata_dose <- factor(umap_df$Metadata_dose, levels = c(
    "0 nM",
    "0.61 nM",
    "1.22 nM",
    "2.44 nM",
    "4.88 nM",
    "9.77 nM",
    "19.53 nM",
    "39.06 nM",
    "78.13 nM",
    "156.25 nM"
    )
    )



# make time a factor with levels
# replace the "T000" with ""
umap_df$Metadata_Time <- gsub("T00", "", umap_df$Metadata_Time)
# make time an integer
umap_df$Metadata_Time <- as.integer(umap_df$Metadata_Time)
# change the Metadata Time columnd to minutes
umap_df$Metadata_Time <- ((umap_df$Metadata_Time)-1) * 30
# add "min" to the time column
umap_df$Metadata_Time <- paste0(umap_df$Metadata_Time, " min")
# make the metadata time column a factor with levels
umap_df$Metadata_Time <- factor(umap_df$Metadata_Time, levels = c(
    "0 min",
    "30 min",
    "60 min",
    "90 min",
    "120 min",
    "150 min",
    "180 min",
    "210 min",
    "240 min",
    "270 min",
    "300 min",
    "330 min",
    "360 min"
))


# make a ggplot of the umap
width <- 30
height <- 20
options(repr.plot.width = width, repr.plot.height = height)
umap_plot <- (
    ggplot(data = umap_df, aes(x = UMAP0, y = UMAP1, color = Metadata_dose))
    + geom_point(size = 0.2)
    + theme_bw()
    + facet_grid(Metadata_dose~Metadata_Time)

    + labs( x = "UMAP0", y = "UMAP1")
    + theme(
        legend.position = "none",
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        )
)
umap_plot
# save
ggsave(paste0("../figures/",data_mode,"/umap_plot_time.png"), plot = umap_plot, width = width, height = height, dpi = 600)

# set temporal colour palette of 13 hues of blue
temporal_palette <- c(
    "#008CF5", "#0079E7", "#0066D9", "#0053CB", "#0040BD", "#002D9F", "#001A91", "#000781", "#000570", "#000460", "#000350", "#000240", "#000130"
)
# calculate the centroid of each UMAP cluster dose and time wise
umap_df_centroids <- umap_df %>% group_by(Metadata_dose, Metadata_Time) %>% summarise(
    UMAP0_centroid = mean(UMAP0),
    UMAP1_centroid = mean(UMAP1)
)

width <- 15
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
# plot the centroids per dose over time
umap_centroid_plot <- (
    ggplot(data = umap_df_centroids, aes(x = UMAP0_centroid, y = UMAP1_centroid, color = Metadata_Time))
    + geom_point(size = 5)
    + theme_bw()
    + labs( x = "UMAP0", y = "UMAP1", title = "Centroids of UMAP space per dose of Staurosporine over time")
    # add custom colors
    + scale_color_manual(values = c(
        "0 min" = temporal_palette[1],
        "30 min" = temporal_palette[2],
        "60 min" = temporal_palette[3],
        "90 min" = temporal_palette[4],
        "120 min" = temporal_palette[5],
        "150 min" = temporal_palette[6],
        "180 min" = temporal_palette[7],
        "210 min" = temporal_palette[8],
        "240 min" = temporal_palette[9],
        "270 min" = temporal_palette[10],
        "300 min" = temporal_palette[11],
        "330 min" = temporal_palette[12],
        "360 min" = temporal_palette[13]
    ))


    # change legend title
    + guides(color = guide_legend(title = "Time (min)", ncol = 2), size = 5)

    + theme(
        strip.text.x = element_text(size = 24),
        strip.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24, hjust = 0.5),
        plot.title = element_text(size = 24, hjust = 0.5)
        )
    + facet_wrap(~Metadata_dose,nrow = 2)

)
umap_centroid_plot
# save
ggsave(paste0("../figures/",data_mode,"/umap_centroid_plot.png"), plot = umap_centroid_plot, width = width, height = height, dpi = 600)


# get the first two, middle, and last two doses
dose_levels <- c("0 nM", "0.61 nM", "9.77 nM", "78.13 nM", "156.25 nM")
umap_df <- umap_df %>% filter(Metadata_dose %in% dose_levels)

# make a ggplot of the umap
width <- 30
height <- 10
options(repr.plot.width = width, repr.plot.height = height)
umap_plot <- (
    ggplot(data = umap_df, aes(x = UMAP0, y = UMAP1, color = Metadata_dose))
    + geom_point(size = 0.2)
    + theme_bw()
    + facet_grid(Metadata_dose~Metadata_Time)

    + labs( x = "UMAP0", y = "UMAP1")
    + theme(
        legend.position = "none",
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),

        )


)
umap_plot
# save
ggsave(paste0("../figures/",data_mode,"/umap_plot_time_part_of_doses.png"), plot = umap_plot, width = width, height = height, dpi = 600)

# set temporal colour palette of 13 hues of blue
temporal_palette <- c(
    "#008CF5", "#0079E7", "#0066D9", "#0053CB", "#0040BD", "#002D9F", "#001A91", "#000781", "#000570", "#000460", "#000350", "#000240", "#000130"
)
# calculate the centroid of each UMAP cluster dose and time wise
umap_df_centroids <- umap_df %>% group_by(Metadata_dose, Metadata_Time) %>% summarise(
    UMAP0_centroid = mean(UMAP0),
    UMAP1_centroid = mean(UMAP1)
)

width <- 20
height <- 4
options(repr.plot.width = width, repr.plot.height = height)
# plot the centroids per dose over time
umap_centroid_plot <- (
    ggplot(data = umap_df_centroids, aes(x = UMAP0_centroid, y = UMAP1_centroid, color = Metadata_Time))
    + geom_point(size = 5)
    + theme_bw()
    + labs( x = "UMAP0", y = "UMAP1", title = "Centroids of UMAP space per dose of Staurosporine over time")
    # add custom colors
    + scale_color_manual(values = c(
        "0 min" = temporal_palette[1],
        "30 min" = temporal_palette[2],
        "60 min" = temporal_palette[3],
        "90 min" = temporal_palette[4],
        "120 min" = temporal_palette[5],
        "150 min" = temporal_palette[6],
        "180 min" = temporal_palette[7],
        "210 min" = temporal_palette[8],
        "240 min" = temporal_palette[9],
        "270 min" = temporal_palette[10],
        "300 min" = temporal_palette[11],
        "330 min" = temporal_palette[12],
        "360 min" = temporal_palette[13]
    ))


    # change legend title
    + guides(color = guide_legend(title = "Time (min)", ncol = 2), size = 5)

    + theme(
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, hjust = 0.5),
        plot.title = element_text(size = 24, hjust = 0.5)
        )
    + facet_grid(~Metadata_dose)

)
umap_centroid_plot
# save
ggsave(paste0("../figures/",data_mode,"/umap_centroid_plot_part of doses.png"), plot = umap_centroid_plot, width = width, height = height, dpi = 600)
