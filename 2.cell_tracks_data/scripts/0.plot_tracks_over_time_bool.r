suppressPackageStartupMessages(suppressWarnings({
    library("ggplot2")
    library(dplyr)
    library(tidyr)
    library(ComplexHeatmap)
    library(tibble)
    library(RColorBrewer)
    library(circlize)
}))


# show all columns
options(repr.matrix.max.cols=200, repr.matrix.max.rows=100)

profile_path <- file.path("../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs.parquet")
figure_path <- file.path("../figures/")
if (!dir.exists(figure_path)) {
    dir.create(figure_path, recursive = TRUE)
}
df <- arrow::read_parquet(profile_path)
# get the metadata columns only
metadata_df <- df %>%
  select(contains("Metadata"))

metadata_df$unique_cell <- paste0(
    metadata_df$Metadata_track_id, "_",

    metadata_df$Metadata_Well, "_",
    metadata_df$Metadata_FOV, "_",
    metadata_df$Metadata_dose
    )
# sort by Metadata_Well
metadata_df <- metadata_df %>%
  arrange(Metadata_Well, Metadata_FOV, Metadata_dose)
head(metadata_df,1)

plotting_df <- metadata_df %>%
  select(c(
    "Metadata_Time",
    "Metadata_dose",
    "unique_cell"
    ))
plotting_df$Metadata_dose <- as.numeric(plotting_df$Metadata_dose)
plotting_df$values <- 1
# sort by metadata dose
plotting_df <- plotting_df %>%
  arrange(Metadata_dose)
# drop dose column
plotting_df <- plotting_df %>%
  select(-Metadata_dose)


# pivot wide such that the Metadata_id in the columns
plotting_df <- plotting_df %>%
  pivot_wider(
    names_from = unique_cell,
    values_from = values
  )
head(plotting_df)

# replace NA with 0
plotting_df[is.na(plotting_df)] <- 0
# turn the values into characters
# remove the dose column
plotting_df$Metadata_dose <- NULL
# remove the Metadata_Time column


# remove the Metadata_Time column
mat <- t(as.matrix(plotting_df))
colnames(mat) <- plotting_df$Metadata_Time
# drop the first row
mat <- mat[-1, ]
# add another column that is the sum of the columns
mat <- as.data.frame(mat)
# make the values numeric
mat <- mat %>%
  mutate(across(where(is.character), as.numeric))
# add a column that is the sum of the columns
mat$sum <- rowSums(mat)
# remove the last column
# sort the matrix by the sum column
# remove the last column
mat <- mat %>%
  select(-c("sum"))
# make the values character
mat <- mat %>%
  mutate(across(where(is.numeric), as.character))
# convert to a matrix
mat <- as.matrix(mat)

rows <- rownames(mat)
# get only the well names
well_dose <- rows %>%
  strsplit("_") %>%
  sapply(function(x) x[4]) %>%
  as.data.frame() %>%
rename(Metadata_dose = ".") %>%
mutate(Metadata_dose = as.character(Metadata_dose))

# convert to a factor
well_dose$Metadata_dose <- factor(
    well_dose$Metadata_dose,
    levels = c(
        "0.0", "0.61", "1.22", "2.44", "4.88",
        "9.77", "19.53", "39.06", "78.13", "156.25"
    )
)


row_ha <- rowAnnotation(
    Dose = well_dose$Metadata_dose,

    annotation_legend_param = list(
        title_position = "topcenter",
        title_gp = gpar(fontsize = 16, angle = 0, fontface = "bold", hjust = 0.5),
        labels_gp = gpar(fontsize = 16,
            # make annotation bar text bigger
            legend = gpar(fontsize = 16),
            annotation_name = gpar(fontsize = 16),
            legend_height = unit(20, "cm"),
            legend_width = unit(1, "cm"),
            # make legend taller
            legend_height = unit(10, "cm"),
            legend_width = unit(1, "cm"),
            legend_key = gpar(fontsize = 16)
            )
    ),


    col = list(
        Dose = c(
            "0.0" = "#57F2F2",
            "0.61" = "#63D6D6",
            "1.22" = "#65BABA",
            "2.44" = "#68A3A3",
            "4.88" = "#668A8A",
            "9.77" = "#5E7070",
            "19.53" = "#4B5757",
            "39.06" = "#2F3D3D",
            "78.13" = "#182424",
            "156.25" = "#030A0A"
        )
    ),



    annotation_name_side = "bottom",
    # make font size bigger
    annotation_name_gp = gpar(fontsize = 16)

)


# make 0 in the matrix Cell Absent and 1 Cell Present but not in the row names or column names
mat[mat == 0.0] <- "Cell Absent"
mat[mat == 1.0] <- "Cell Present"

colors = structure(
    c("#c8c8c8","#2a2a2a"),
    names = c("Cell Absent", "Cell Present")
)

width <- 10
height <- 10
options(repr.plot.width=width, repr.plot.height=height)
ht_opt$message = FALSE
heatmap <- Heatmap(
    mat,
    cluster_rows = TRUE,    # Cluster rows
    # cluster_columns = FALSE, # Cluster columns
    show_row_names = FALSE,  # Show row names
    show_column_names = TRUE, # Show column names
    column_names_gp = gpar(fontsize = 16), # Column name label formatting
    row_names_gp = gpar(fontsize = 14),    # Row name label formatting
    right_annotation = row_ha,

    heatmap_legend_param = list(
                title = "Track\nBoolean",
                title_position = "topcenter",
                title_gp = gpar(fontsize = 16, angle = 0, fontface = "bold", hjust = 0.5),
                labels_gp = gpar(fontsize = 16),
                legend_height = unit(6.6, "cm")
                ),
        # set color for 0 and 1
        col = colors,

)

png(filename = paste0(figure_path, "cell_tracks_over_time_heatmap.png"), width = width, height = height, units = "in", res = 600)
heatmap
dev.off()
heatmap
